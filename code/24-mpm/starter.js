// Chapter 24 - The Material Point Method (starter)
// Explicit MPM with APIC transfers on a 64 x 64 grid, following the eight steps of the MPM course
// (Jiang et al. 2016, section 10.5). The APIC transfers of mass and momentum from chapter 23 are done
// for you, so the blocks already drift with their initial velocities. Your job, in the TODO blocks:
//   1. scatter the stress forces to the grid in P2G,
//   2. the grid velocity update (forces, gravity, walls),
//   3. the F update in G2P,
// then the neo-Hookean stress and the water equation of state.

const scene = Sim.select('scene', ['falling jelly', 'dam break', 'jelly into water'], 'falling jelly', { reset: true });
const model = Sim.select('jelly model', ['fixed corotated', 'neo-Hookean'], 'fixed corotated');
const youngs = Sim.logSlider("jelly Young's modulus E (Pa)", 3, 6.5, 2e4);
const poisson = Sim.slider('jelly Poisson ratio nu', 0, 0.45, 0.3, { step: 0.01 });
const bulk = Sim.logSlider('water bulk modulus K (Pa)', 3, 6.5, 3e4);
const stepMs = Sim.logSlider('time step h (ms)', -1, 1, 1);
const walls = Sim.select('walls', ['slip', 'sticky'], 'slip');

// ---------------------------------------------------------------------------------------------
// Grid: N x N cells on the unit square, nodes at (i dx, j dx) for i, j = 0..N.
// ---------------------------------------------------------------------------------------------
const N = 64;
const dx = 1 / N;                  // grid spacing (m)
const invDx = 1 / dx;
const NN = N + 1;                  // nodes per row
const BOUND = 3;                   // nodes closer than this to the border are wall nodes
const gravity = -9.81;
const rhoJelly = 1000, rhoWater = 1000;   // densities (kg/m^3; per metre of depth in 2D)

const gm = new Float64Array(NN * NN);     // node mass
const gvx = new Float64Array(NN * NN);    // node momentum, then velocity
const gvy = new Float64Array(NN * NN);
const gfx = new Float64Array(NN * NN);    // node force
const gfy = new Float64Array(NN * NN);

// ---------------------------------------------------------------------------------------------
// Particles, in flat arrays.
// ---------------------------------------------------------------------------------------------
const MAXP = 6000;
const JELLY = 0, WATER = 1;
const px = new Float64Array(MAXP), py = new Float64Array(MAXP);
const pvx = new Float64Array(MAXP), pvy = new Float64Array(MAXP);
const pC = new Float64Array(4 * MAXP);     // APIC affine velocity C_p = B_p D_p^-1 (row-major 2x2)
const pF = new Float64Array(4 * MAXP);     // deformation gradient F_p (row-major 2x2)
const pJ = new Float64Array(MAXP);         // volume ratio J for water (water does not need all of F)
const pm = new Float64Array(MAXP);         // mass m_p
const pV0 = new Float64Array(MAXP);        // initial volume V_p^0 (an area in 2D)
const pmat = new Uint8Array(MAXP);
const pcol = new Array(MAXP);
let count = 0;

function addBlock(x0, y0, x1, y1, vx, vy, mat, color) {
  const s = dx / 2;                        // 2 x 2 particles per cell
  for (let x = x0 + s / 2; x < x1; x += s) {
    for (let y = y0 + s / 2; y < y1; y += s) {
      if (count >= MAXP) return;
      const p = count++;
      px[p] = x; py[p] = y; pvx[p] = vx; pvy[p] = vy;
      pC[4 * p] = pC[4 * p + 1] = pC[4 * p + 2] = pC[4 * p + 3] = 0;
      pF[4 * p] = 1; pF[4 * p + 1] = 0; pF[4 * p + 2] = 0; pF[4 * p + 3] = 1;
      pJ[p] = 1;
      pV0[p] = s * s;
      pm[p] = (mat === WATER ? rhoWater : rhoJelly) * pV0[p];
      pmat[p] = mat;
      pcol[p] = color;
    }
  }
}

let exploded = false;

function reset() {
  count = 0;
  exploded = false;
  const c = Sim.colors;
  if (scene.value === 'falling jelly') {
    addBlock(0.12, 0.05, 0.32, 0.25, 0, 0, JELLY, c.red);
    addBlock(0.22, 0.45, 0.42, 0.65, 0, 0, JELLY, c.orange);
    addBlock(0.66, 0.30, 0.86, 0.50, -2.5, 1.5, JELLY, c.purple);
  } else if (scene.value === 'dam break') {
    addBlock(0.05, 0.05, 0.37, 0.62, 0, 0, WATER, c.blue);
  } else {
    addBlock(0.05, 0.05, 0.95, 0.22, 0, 0, WATER, c.blue);
    addBlock(0.40, 0.50, 0.60, 0.70, 0, -1, JELLY, c.red);
  }
}

// Lame parameters from Young's modulus and Poisson's ratio.
function lameMu() { return youngs.value / (2 * (1 + poisson.value)); }
function lameLambda() { const nu = poisson.value; return youngs.value * nu / ((1 + nu) * (1 - 2 * nu)); }

// ---------------------------------------------------------------------------------------------
// Quadratic B-spline weights and their gradients for the 3 x 3 nodes around a particle.
// fx = x/dx - base is in [0.5, 1.5); node k sits at (base + k) dx, so x_i - x_p = (k - fx) dx.
// ---------------------------------------------------------------------------------------------
const wx = new Float64Array(3), wy = new Float64Array(3);
const dwx = new Float64Array(3), dwy = new Float64Array(3);
let baseX = 0, baseY = 0, fracX = 0, fracY = 0;

function computeWeights(x, y) {
  baseX = Math.floor(x * invDx - 0.5);
  baseY = Math.floor(y * invDx - 0.5);
  fracX = x * invDx - baseX;
  fracY = y * invDx - baseY;
  quadratic(fracX, wx, dwx);
  quadratic(fracY, wy, dwy);
}

function quadratic(f, w, dw) {
  w[0] = 0.5 * (1.5 - f) * (1.5 - f);
  w[1] = 0.75 - (f - 1) * (f - 1);
  w[2] = 0.5 * (f - 0.5) * (f - 0.5);
  dw[0] = (f - 1.5) * invDx;               // derivatives with respect to x_p (1/m)
  dw[1] = -2 * (f - 1) * invDx;
  dw[2] = (f - 0.5) * invDx;
}

// ---------------------------------------------------------------------------------------------
// Constitutive models. Each writes the first Piola-Kirchhoff stress P(F) into P (row-major 2x2).
// ---------------------------------------------------------------------------------------------
const P = new Float64Array(4);

// Fixed corotated (course eq. 52): P = 2 mu (F - R) + lambda (J - 1) J F^-T.
function fixedCorotatedP(F00, F01, F10, F11, mu, la, P) {
  // Closest rotation R (polar decomposition). In 2D it is a rotation by this angle.
  const theta = Math.atan2(F10 - F01, F00 + F11);
  const c = Math.cos(theta), s = Math.sin(theta);
  const J = F00 * F11 - F01 * F10;
  // J F^-T is the cofactor matrix of F: [[F11, -F10], [-F01, F00]].
  P[0] = 2 * mu * (F00 - c) + la * (J - 1) * F11;
  P[1] = 2 * mu * (F01 + s) + la * (J - 1) * -F10;
  P[2] = 2 * mu * (F10 - s) + la * (J - 1) * -F01;
  P[3] = 2 * mu * (F11 - c) + la * (J - 1) * F00;
}

// Neo-Hookean (course eq. 48): P = mu (F - F^-T) + lambda log(J) F^-T.
function neoHookeanP(F00, F01, F10, F11, mu, la, P) {
  // TODO (exercise 3): write P = mu (F - F^-T) + lambda log(J) F^-T into P[0..3].
  // F^-T is the cofactor matrix [[F11, -F10], [-F01, F00]] divided by J. Clamp J away from zero.
  P[0] = 0; P[1] = 0; P[2] = 0; P[3] = 0;
}

// Water: energy psi(J) = K/2 (J - 1)^2, so P F^T = K (J - 1) J I = -J p I with pressure p = -K (J - 1).
// Returns the scalar s with P F^T = s I.
function waterPFt(J, K) {
  // TODO (exercise 4): return the scalar s. Until then water feels no pressure and collapses like dust.
  return 0;
}

// ---------------------------------------------------------------------------------------------
// One explicit MPM step.
// ---------------------------------------------------------------------------------------------
function step(h) {
  if (exploded) return;
  gm.fill(0); gvx.fill(0); gvy.fill(0); gfx.fill(0); gfy.fill(0);
  const mu = lameMu(), la = lameLambda(), K = bulk.value;
  const neo = model.value === 'neo-Hookean';

  // Steps 1 and 4: P2G of mass and APIC momentum, and the stress forces f_i = -sum_p V_p^0 P F^T grad w_ip.
  for (let p = 0; p < count; p++) {
    computeWeights(px[p], py[p]);
    const k = 4 * p;
    const C00 = pC[k], C01 = pC[k + 1], C10 = pC[k + 2], C11 = pC[k + 3];

    // A = V_p^0 P(F_p) F_p^T, the particle's contribution to the force (units: N m in 2D).
    let A00, A01, A10, A11;
    if (pmat[p] === WATER) {
      const s = pV0[p] * waterPFt(pJ[p], K);
      A00 = s; A01 = 0; A10 = 0; A11 = s;
    } else {
      const F00 = pF[k], F01 = pF[k + 1], F10 = pF[k + 2], F11 = pF[k + 3];
      if (neo) neoHookeanP(F00, F01, F10, F11, mu, la, P);
      else fixedCorotatedP(F00, F01, F10, F11, mu, la, P);
      const V = pV0[p];
      A00 = V * (P[0] * F00 + P[1] * F01);
      A01 = V * (P[0] * F10 + P[1] * F11);
      A10 = V * (P[2] * F00 + P[3] * F01);
      A11 = V * (P[2] * F10 + P[3] * F11);
    }

    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const node = (baseX + a) * NN + (baseY + b);
        const w = wx[a] * wy[b];
        const gx = dwx[a] * wy[b], gy = wx[a] * dwy[b];     // grad w_ip
        const rx = (a - fracX) * dx, ry = (b - fracY) * dx; // x_i - x_p
        const mw = w * pm[p];
        gm[node] += mw;
        gvx[node] += mw * (pvx[p] + C00 * rx + C01 * ry);
        gvy[node] += mw * (pvy[p] + C10 * rx + C11 * ry);
        // TODO (exercise 1): the stress force on this node, f_i -= A grad w_ip, into gfx and gfy.
      }
    }
  }

  // Steps 2, 3 and 5: grid velocities on active nodes (m_i > 0), forces and gravity, walls.
  const sticky = walls.value === 'sticky';
  for (let i = 0; i < NN; i++) {
    for (let j = 0; j < NN; j++) {
      const node = i * NN + j;
      const m = gm[node];
      if (m <= 0) { gvx[node] = 0; gvy[node] = 0; continue; }   // inactive node
      let vx = gvx[node] / m, vy = gvy[node] / m;
      // TODO (exercise 1): add the force and gravity, v_i += h (f_i / m_i + g).
      // TODO (exercise 1): wall nodes have i < BOUND, i > N - BOUND, j < BOUND or j > N - BOUND.
      //   Sticky: set the velocity to zero. Slip: zero only a component that points into the wall.
      gvx[node] = vx; gvy[node] = vy;
    }
  }
  applyMouse(h);

  // Steps 6, 7 and 8: G2P of velocity and C_p, the F update F <- (I + h grad v) F, advection.
  for (let p = 0; p < count; p++) {
    computeWeights(px[p], py[p]);
    let vx = 0, vy = 0, B00 = 0, B01 = 0, B10 = 0, B11 = 0;
    let G00 = 0, G01 = 0, G10 = 0, G11 = 0;     // velocity gradient grad v_p = sum_i v_i (grad w_ip)^T
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const node = (baseX + a) * NN + (baseY + b);
        const w = wx[a] * wy[b];
        const gx = dwx[a] * wy[b], gy = wx[a] * dwy[b];
        const rx = (a - fracX) * dx, ry = (b - fracY) * dx;
        const ux = gvx[node], uy = gvy[node];
        vx += w * ux; vy += w * uy;
        B00 += w * ux * rx; B01 += w * ux * ry;
        B10 += w * uy * rx; B11 += w * uy * ry;
        G00 += ux * gx; G01 += ux * gy;
        G10 += uy * gx; G11 += uy * gy;
      }
    }
    const k = 4 * p;
    const s = 4 * invDx * invDx;                // D_p^-1 = 4 / dx^2 for quadratic B-splines
    pC[k] = s * B00; pC[k + 1] = s * B01; pC[k + 2] = s * B10; pC[k + 3] = s * B11;
    pvx[p] = vx; pvy[p] = vy;

    // TODO (exercise 1): F <- (I + h grad v) F, with grad v = [[G00, G01], [G10, G11]].
    //   Write the new F into pF[k..k+3]. Water only needs J: J <- det(I + h grad v) J, into pJ[p].

    px[p] += h * vx;
    py[p] += h * vy;
    clampInside(p);
  }
}

// ---------------------------------------------------------------------------------------------
// Scaffolding: keep particles inside the grid, mouse push, explosion check, drawing.
// ---------------------------------------------------------------------------------------------
function clampInside(p) {
  const lo = 1.0 * dx, hi = 1 - 1.0 * dx;
  if (!(px[p] > lo)) px[p] = lo; else if (px[p] > hi) px[p] = hi;   // also catches NaN
  if (!(py[p] > lo)) py[p] = lo; else if (py[p] > hi) py[p] = hi;
}

// While the mouse is held down, grid nodes near it are pushed away from it.
function applyMouse(h) {
  if (!Sim.mouse.down || !Sim.mouse.inside) return;
  const r = 0.08, strength = 40;           // radius (m), acceleration (m/s^2)
  for (let i = 0; i < NN; i++) {
    for (let j = 0; j < NN; j++) {
      const node = i * NN + j;
      if (gm[node] <= 0) continue;
      const ex = i * dx - Sim.mouse.x, ey = j * dx - Sim.mouse.y;
      const d = Math.hypot(ex, ey);
      if (d < r && d > 1e-6) {
        gvx[node] += h * strength * ex / d;
        gvy[node] += h * strength * ey / d;
      }
    }
  }
}

const FRAME = 1 / 60, MAX_SUBSTEPS = 60;
function substeps() { return Math.min(MAX_SUBSTEPS, Math.max(1, Math.round(FRAME / (stepMs.value / 1000)))); }

function draw(g) {
  const wallColor = '#e3e0d6', b = BOUND * dx;
  g.rect(0, 0, 1, b, { fill: wallColor });
  g.rect(0, 1 - b, 1, b, { fill: wallColor });
  g.rect(0, 0, b, 1, { fill: wallColor });
  g.rect(1 - b, 0, b, 1, { fill: wallColor });
  g.rect(b, b, 1 - 2 * b, 1 - 2 * b, { stroke: Sim.colors.muted, width: 1 });

  let vmax = 0;
  for (let p = 0; p < count; p++) {
    const v = Math.hypot(pvx[p], pvy[p]);
    if (!(v < 1e3)) { vmax = Infinity; break; }
    if (v > vmax) vmax = v;
  }
  if (!exploded && vmax > 50) exploded = true;   // nothing in these scenes moves faster than ~5 m/s
  g.points(px, py, count, { radius: 1.6, colors: pcol });

  const h = stepMs.value / 1000;
  const mu = lameMu(), la = lameLambda();
  const cJelly = Math.sqrt((la + 2 * mu) / rhoJelly);
  const cWater = Math.sqrt(bulk.value / rhoWater);
  const c = scene.value === 'falling jelly' ? cJelly : scene.value === 'dam break' ? cWater : Math.max(cJelly, cWater);
  Sim.show('particles', count);
  Sim.show('wave speed c', c.toFixed(2) + ' m/s');
  Sim.show('dx / c', (dx / c * 1000).toFixed(2) + ' ms');
  Sim.show('h c / dx', (h * c / dx).toFixed(3));
  Sim.show('h vmax / dx', isFinite(vmax) ? (h * vmax / dx).toFixed(3) : 'inf');
  const n = substeps();
  Sim.show('substeps / frame', n + (n * h < FRAME * 0.99 ? ' (slow motion)' : ''));
  Sim.plot('max speed (m/s)', Math.min(vmax, 20), { color: Sim.colors.red, min: 0 });
  if (exploded) {
    g.text('Exploded! Lower h or E, then press R.', 0.5, 0.55, { color: Sim.colors.red, size: 18, align: 'center' });
  } else if (Sim.time < 3) {
    g.text('Click and hold to push the material', 0.5, 0.94, { color: Sim.colors.muted, size: 12, align: 'center' });
  }
}

Sim.run({
  plotCorner: 'top-right',
  title: 'Explicit MPM: jelly and water (starter)', view: [0, 0, 1, 1], reset, step, draw, grid: false,
  dt: () => stepMs.value / 1000, stepsPerFrame: substeps,
});
