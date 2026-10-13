// Chapter 24 - APIC (starter)
// The transfers of chapter 23, plus two richer ones: RPIC, where each particle also carries a spin,
// and APIC, where each particle carries a whole affine velocity field v_p + C_p (x - x_p).
// Scene "rotating disk": a rigidly spinning disk held together by a torque-free pull to the centre.
// Scene "colliding blobs": two blobs meet off-centre and splash; a grid pressure keeps them apart.
// Angular momentum is conserved only until the first particle reaches a wall (shown in the info box).
// The info box shows the total angular momentum, counting the particles' spin for RPIC and APIC.

const scene = Sim.select('scene', ['rotating disk', 'colliding blobs'], 'rotating disk', { reset: true });
const transfer = Sim.select('transfer', ['PIC', 'FLIP 0.95', 'FLIP', 'RPIC', 'APIC'], 'APIC', { reset: true });
const cellsAcross = Sim.select('grid cells across', ['16', '32', '64', '128'], '64', { reset: true });
const perAxis = Sim.slider('particles per cell, per axis', 1, 3, 2, { step: 1, reset: true });
const colourBy = Sim.select('colour by', ['speed', 'velocity the grid cannot see'], 'speed');
Sim.button('add particle noise', addNoise);

const density = 1000;          // kg/m^3 (the 2D slab is 1 m thick)
const spinRate = 1.5;          // rad/s, rotating disk
const blobSpeed = 0.5;         // m/s, each colliding blob
const soundSpeed = 3.0;        // m/s, sets the stiffness of the grid pressure
const bound = 2;               // wall boundary: nodes this close to the edge get wall conditions

// ---------------------------------------------------------------------------------------------
// State: flat Float64Arrays. Particle p: px[p], py[p], pvx[p], pvy[p], pm[p];
// APIC matrix C_p = [[pC[4p], pC[4p+1]], [pC[4p+2], pC[4p+3]]]; RPIC spin angular momentum pS[p].
// Grid node (i, j) at position (i dx, j dx), stored at index i + j * nodesX.
// ---------------------------------------------------------------------------------------------
let n, dx, invDx, nodesX, nodeCount;
let gm, gvx, gvy, gvxOld, gvyOld;
let count, px, py, pvx, pvy, pm, pC, pS, hidden, marker;
let L0, E0, stepSize = 1 / 60, dpReport = '', wallTime = -1, simTime = 0;

const wx = new Float64Array(3), wy = new Float64Array(3);

// Quadratic B-spline weights of nodes base, base + 1, base + 2 along one axis (chapter 23).
// fx = x / dx - base is in [0.5, 1.5); the distances to the three nodes are fx, fx - 1, fx - 2 cells.
function quadraticWeights(fx, w) {
  w[0] = 0.5 * (1.5 - fx) * (1.5 - fx);
  w[1] = 0.75 - (fx - 1) * (fx - 1);
  w[2] = 0.5 * (fx - 0.5) * (fx - 0.5);
}

// D_p = sum_i w_ip (x_i - x_p)(x_i - x_p)^T, averaged over many random particle positions,
// plus the largest deviation of any single D_p / dx^2 from the prediction (1/4) I.
function checkDp() {
  const rand = Sim.random(42), trials = 1000;
  let d00 = 0, d01 = 0, d11 = 0, worst = 0;
  for (let t = 0; t < trials; t++) {
    const x = (8 + 10 * rand()) * dx, y = (8 + 10 * rand()) * dx;
    const gx = x * invDx, gy = y * invDx;
    const bx = Math.floor(gx - 0.5), by = Math.floor(gy - 0.5);
    quadraticWeights(gx - bx, wx);
    quadraticWeights(gy - by, wy);
    let a00 = 0, a01 = 0, a11 = 0;
    // TODO (exercise 2): loop over the 3x3 stencil as in particlesToGrid and accumulate
    // a00, a01, a11 = the entries of sum_i w_ip (x_i - x_p)(x_i - x_p)^T for the point (x, y).
    const s = 1 / (dx * dx);
    d00 += a00 * s / trials; d01 += a01 * s / trials; d11 += a11 * s / trials;
    worst = Math.max(worst, Math.abs(a00 * s - 0.25), Math.abs(a01 * s), Math.abs(a11 * s - 0.25));
  }
  return `${d00.toFixed(4)}, ${d01.toFixed(4)}, ${d11.toFixed(4)} (worst ${worst.toExponential(0)})`;
}

function reset() {
  n = +cellsAcross.value; dx = 1 / n; invDx = n;
  nodesX = n + 1; nodeCount = nodesX * nodesX;
  gm = new Float64Array(nodeCount); gvx = new Float64Array(nodeCount); gvy = new Float64Array(nodeCount);
  gvxOld = new Float64Array(nodeCount); gvyOld = new Float64Array(nodeCount);

  // Each blob: centre x, centre y, radius, velocity x, velocity y
  const blobs = scene.value === 'rotating disk'
    ? [[0.5, 0.5, 0.25, 0, 0]]
    : [[0.3, 0.45, 0.1, blobSpeed, 0], [0.7, 0.55, 0.1, -blobSpeed, 0]];
  const k = perAxis.value, xs = [], ys = [], us = [], vs = [];
  for (let j = 0; j < n * k; j++) {
    for (let i = 0; i < n * k; i++) {
      const x = (i + 0.5) * dx / k, y = (j + 0.5) * dx / k;
      for (const [cx, cy, r, u, v] of blobs) {
        if ((x - cx) * (x - cx) + (y - cy) * (y - cy) >= r * r) continue;
        xs.push(x); ys.push(y);
        if (scene.value === 'rotating disk') { us.push(-spinRate * (y - cy)); vs.push(spinRate * (x - cx)); }
        else { us.push(u); vs.push(v); }
      }
    }
  }
  count = xs.length;
  px = new Float64Array(xs); py = new Float64Array(ys);
  pvx = new Float64Array(us); pvy = new Float64Array(vs);
  pm = new Float64Array(count).fill(density * dx * dx / (k * k));
  pC = new Float64Array(4 * count);
  pS = new Float64Array(count);
  hidden = new Float64Array(count);
  marker = new Uint8Array(count);
  for (let p = 0; p < count; p++) {
    if (scene.value !== 'rotating disk') continue;
    // Start each method with the exact local field of the rotation: APIC with its velocity
    // gradient, RPIC with the matching spin K_p omega.
    if (transfer.value === 'APIC') { pC[4 * p + 1] = -spinRate; pC[4 * p + 2] = spinRate; }
    if (transfer.value === 'RPIC') pS[p] = pm[p] * 0.5 * dx * dx * spinRate;
    if (Math.abs(py[p] - 0.5) < 0.6 * dx && px[p] > 0.5) marker[p] = 1;   // a spoke, to see the rotation
  }
  L0 = angularMomentum(); E0 = kineticEnergy();
  dpReport = checkDp();
  wallTime = -1; simTime = 0;
}

// ---------------------------------------------------------------------------------------------
// Transfers
// ---------------------------------------------------------------------------------------------

// P2G: m_i = sum_p w_ip m_p and (m v)_i = sum_p w_ip m_p (v_p + C_p (x_i - x_p))   (APIC eq. 8).
// For quadratic weights D_p = dx^2/4 I, so C_p = B_p D_p^-1 = (4 / dx^2) B_p, and we store C_p.
// RPIC uses the rigid local field v_p + omega_p x (x_i - x_p) with omega_p = spin / K_p (eqs. 6-7).
function particlesToGrid() {
  gm.fill(0); gvx.fill(0); gvy.fill(0);
  const mode = transfer.value;
  for (let p = 0; p < count; p++) {
    const gx = px[p] * invDx, gy = py[p] * invDx;
    const bx = Math.floor(gx - 0.5), by = Math.floor(gy - 0.5);
    quadraticWeights(gx - bx, wx);
    quadraticWeights(gy - by, wy);
    let c00 = 0, c01 = 0, c10 = 0, c11 = 0;
    if (mode === 'APIC') {
      // TODO (exercise 1): read this particle's C_p from pC into c00, c01, c10, c11.
    } else if (mode === 'RPIC') {
      // TODO (exercise 4, challenge): turn the spin pS[p] into an angular velocity omega with the
      // inertia K_p of the stencil, and set c01, c10 so that C r = omega x r.
    }
    for (let b = 0; b < 3; b++) {
      for (let a = 0; a < 3; a++) {
        const w = wx[a] * wy[b], node = (bx + a) + (by + b) * nodesX;
        const rx = (bx + a) * dx - px[p], ry = (by + b) * dx - py[p];   // x_i - x_p
        const mw = w * pm[p];
        gm[node] += mw;
        gvx[node] += mw * (pvx[p] + c00 * rx + c01 * ry);
        gvy[node] += mw * (pvy[p] + c10 * rx + c11 * ry);
      }
    }
  }
  for (let i = 0; i < nodeCount; i++) {
    if (gm[i] > 0) { gvx[i] /= gm[i]; gvy[i] /= gm[i]; }
    gvxOld[i] = gvx[i]; gvyOld[i] = gvy[i];
  }
}

// Grid update: forces and walls, as in chapter 23.
function updateGrid(h) {
  if (scene.value === 'rotating disk') {
    const w2 = spinRate * spinRate;
    for (let j = 0; j < nodesX; j++) {
      for (let i = 0; i < nodesX; i++) {
        const node = i + j * nodesX;
        if (gm[node] === 0) continue;
        gvx[node] -= h * w2 * (i * dx - 0.5);
        gvy[node] -= h * w2 * (j * dx - 0.5);
      }
    }
  } else {
    applyPressure(h);
  }
  for (let j = 0; j < nodesX; j++) {
    for (let i = 0; i < nodesX; i++) {
      const node = i + j * nodesX;
      if (i < bound && gvx[node] < 0) gvx[node] = 0;
      if (i > n - bound && gvx[node] > 0) gvx[node] = 0;
      if (j < bound && gvy[node] < 0) gvy[node] = 0;
      if (j > n - bound && gvy[node] > 0) gvy[node] = 0;
    }
  }
}

// G2P: v_p = sum_i w_ip v_i for PIC, RPIC and APIC; FLIP adds the grid change instead.
// APIC also gathers B_p = sum_i w_ip v_i (x_i - x_p)^T (eq. 10); RPIC gathers the spin (eq. 4).
function gridToParticles() {
  const mode = transfer.value;
  const alpha = mode === 'FLIP' ? 1 : mode === 'FLIP 0.95' ? 0.95 : 0;
  const DInv = 4 / (dx * dx);
  for (let p = 0; p < count; p++) {
    const gx = px[p] * invDx, gy = py[p] * invDx;
    const bx = Math.floor(gx - 0.5), by = Math.floor(gy - 0.5);
    quadraticWeights(gx - bx, wx);
    quadraticWeights(gy - by, wy);
    let picX = 0, picY = 0, dvX = 0, dvY = 0;
    let b00 = 0, b01 = 0, b10 = 0, b11 = 0;
    for (let b = 0; b < 3; b++) {
      for (let a = 0; a < 3; a++) {
        const w = wx[a] * wy[b], node = (bx + a) + (by + b) * nodesX;
        const vx = gvx[node], vy = gvy[node];
        const rx = (bx + a) * dx - px[p], ry = (by + b) * dx - py[p];
        picX += w * vx; picY += w * vy;
        dvX += w * (vx - gvxOld[node]); dvY += w * (vy - gvyOld[node]);
        // TODO (exercise 1): accumulate B_p = sum_i w_ip v_i (x_i - x_p)^T into b00, b01, b10, b11.
      }
    }
    pvx[p] = alpha * (pvx[p] + dvX) + (1 - alpha) * picX;
    pvy[p] = alpha * (pvy[p] + dvY) + (1 - alpha) * picY;
    hidden[p] = Math.hypot(pvx[p] - picX, pvy[p] - picY);
    if (mode === 'APIC') {
      // TODO (exercise 1): store C_p = B_p D_p^-1 in pC (DInv holds the inverse of dx^2 / 4).
    } else if (mode === 'RPIC') {
      // TODO (exercise 4, challenge): store the spin sum_i w_ip (x_i - x_p) x m_p v_i in pS[p].
      // In 2D, a x b = a.x b.y - a.y b.x; it can be written with the entries of B_p.
    }
  }
}

// Grid pressure (chapter 23): p_i = K max(0, rho_i / rho0 - 1), f_i = -dx^2 (grad p)_i.
// Central differences make the forces sum to zero and exert no net torque.
let pressure = null;
function applyPressure(h) {
  if (!pressure || pressure.length !== nodeCount) pressure = new Float64Array(nodeCount);
  const K = density * soundSpeed * soundSpeed;
  for (let i = 0; i < nodeCount; i++) pressure[i] = K * Math.max(0, gm[i] / (density * dx * dx) - 1);
  for (let j = 1; j < n; j++) {
    for (let i = 1; i < n; i++) {
      const node = i + j * nodesX;
      if (gm[node] === 0) continue;
      const gpx = (pressure[node + 1] - pressure[node - 1]) / (2 * dx);
      const gpy = (pressure[node + nodesX] - pressure[node - nodesX]) / (2 * dx);
      gvx[node] -= h * dx * dx * gpx / gm[node];
      gvy[node] -= h * dx * dx * gpy / gm[node];
    }
  }
}

function advect(h) {
  const lo = dx, hi = 1 - dx, near = (bound + 1) * dx;
  for (let p = 0; p < count; p++) {
    px[p] = Math.min(hi, Math.max(lo, px[p] + h * pvx[p]));
    py[p] = Math.min(hi, Math.max(lo, py[p] + h * pvy[p]));
    // Once a particle's stencil touches the wall nodes, the walls can exert a torque.
    if (wallTime < 0 && (px[p] < near || px[p] > 1 - near || py[p] < near || py[p] > 1 - near)) wallTime = simTime;
  }
  simTime += h;
}

function simulate(h) {
  particlesToGrid();
  updateGrid(h);
  gridToParticles();
  advect(h);
}

function step(dt) {
  const sub = scene.value === 'rotating disk' ? 1 : Math.ceil(dt / (0.4 * dx / soundSpeed));
  stepSize = dt / sub;
  for (let s = 0; s < sub; s++) simulate(stepSize);
}

function addNoise() {
  const rand = Sim.random(Sim.frame + 1);
  for (let p = 0; p < count; p++) {
    const a = 2 * Math.PI * rand();
    pvx[p] += 0.3 * Math.cos(a); pvy[p] += 0.3 * Math.sin(a);
  }
}

// ---------------------------------------------------------------------------------------------
// Measurements and drawing
// ---------------------------------------------------------------------------------------------
// Total angular momentum about (0.5, 0.5): orbital part plus the particles' own spin.
// APIC eq. 11: the spin of particle p is m_p times the skew part of B_p, in 2D m_p (B10 - B01),
// with B_p = C_p D_p = (dx^2 / 4) C_p. RPIC stores its spin directly.
function angularMomentum() {
  let L = 0;
  for (let p = 0; p < count; p++) {
    L += pm[p] * ((px[p] - 0.5) * pvy[p] - (py[p] - 0.5) * pvx[p]);
    L += pm[p] * 0.25 * dx * dx * (pC[4 * p + 2] - pC[4 * p + 1]);
    L += pS[p];
  }
  return L;
}
function kineticEnergy() {
  let E = 0;
  for (let p = 0; p < count; p++) E += 0.5 * pm[p] * (pvx[p] * pvx[p] + pvy[p] * pvy[p]);
  return E;
}

const ramp = [];
for (let i = 0; i <= 32; i++) {
  const t = i / 32, a = [205, 214, 230], b = [24, 64, 160];
  ramp.push(`rgb(${Math.round(a[0] + (b[0] - a[0]) * t)},${Math.round(a[1] + (b[1] - a[1]) * t)},${Math.round(a[2] + (b[2] - a[2]) * t)})`);
}
let colours = [];

function draw(g) {
  if (n <= 64) {
    for (let i = 0; i <= n; i++) {
      g.line(i * dx, 0, i * dx, 1, { color: '#e6e4dc', width: 0.5 });
      g.line(0, i * dx, 1, i * dx, { color: '#e6e4dc', width: 0.5 });
    }
  }
  g.rect(0, 0, 1, 1, { color: Sim.colors.muted, width: 1.5 });
  if (colours.length !== count) colours = new Array(count);
  const bySpeed = colourBy.value === 'speed';
  const scale = bySpeed ? (scene.value === 'rotating disk' ? spinRate * 0.25 : blobSpeed) : 0.3;
  for (let p = 0; p < count; p++) {
    const s = bySpeed ? Math.hypot(pvx[p], pvy[p]) : hidden[p];
    colours[p] = marker[p] && bySpeed ? Sim.colors.orange : ramp[Math.min(32, Math.round(32 * s / scale))];
  }
  g.points(px, py, count, { radius: n > 64 ? 1 : 1.5, colors: colours });

  const L = angularMomentum(), E = kineticEnergy();
  Sim.show('particles / grid nodes', `${count} / ${nodeCount}`);
  Sim.show('time step h (ms)', stepSize * 1000, 2);
  Sim.show('L / L0', L / L0, 6);
  Sim.show('L / L0 - 1', (L / L0 - 1).toExponential(2));
  Sim.show('kinetic energy / start', E / E0, 4);
  Sim.show('first wall contact', wallTime < 0 ? 'none yet' : wallTime.toFixed(2) + ' s');
  Sim.show('D_p / dx^2: xx, xy, yy', dpReport);
  Sim.plot('L / L0', L / L0, { color: Sim.colors.blue, min: 0, max: 1.05 });
  Sim.plot('KE / KE0', Math.min(E / E0, 2), { color: Sim.colors.red, min: 0 });
}

Sim.run({ title: 'APIC transfers (starter)', view: [-0.02, -0.32, 1.02, 1.18], reset, step, draw, grid: false });
