// Chapter 22 - Particles Meet Grids: PIC and FLIP (reference solution)
// Particles carry mass and velocity. Every step they are splatted onto a background grid (P2G),
// the grid velocities are updated, and the result is read back (G2P) with PIC, FLIP or a blend.
// Scene "rotating disk": a disk spinning rigidly, held together by a torque-free pull toward the centre,
// so any change in angular momentum is the transfer's fault.
// Scene "blob into wall": a blob flies into the right wall; a simple grid pressure makes it bounce.

const scene = Sim.select('scene', ['rotating disk', 'blob into wall'], 'rotating disk', { reset: true });
const transfer = Sim.select('transfer', ['PIC', 'FLIP', 'FLIP/PIC blend'], 'PIC', { reset: true });
const flipFraction = Sim.slider('FLIP fraction (blend only)', 0, 1, 0.95, { step: 0.01 });
const cellsAcross = Sim.select('grid cells across', ['16', '32', '64', '128'], '64', { reset: true });
const perAxis = Sim.slider('particles per cell, per axis', 1, 3, 2, { step: 1, reset: true });
const colourBy = Sim.select('colour by', ['speed', 'velocity the grid cannot see'], 'speed');
Sim.button('add particle noise', addNoise);

const density = 1000;          // kg/m^3 (the 2D slab is 1 m thick)
const spinRate = 1.5;          // rad/s, rotating disk
const blobSpeed = 2.0;         // m/s, blob into wall
const soundSpeed = 3.0;        // m/s, sets the stiffness of the grid pressure
const bound = 2;               // wall boundary: nodes this close to the edge get wall conditions

// ---------------------------------------------------------------------------------------------
// State: flat Float64Arrays. Particle p: px[p], py[p], pvx[p], pvy[p], pm[p].
// Grid node (i, j) at position (i dx, j dx), stored at index i + j * nodesX.
// ---------------------------------------------------------------------------------------------
let n, dx, invDx, nodesX, nodeCount;
let gm, gvx, gvy, gvxOld, gvyOld;
let count, px, py, pvx, pvy, pm, hidden, marker;
let L0, E0, P0, stepSize = 1 / 60;

const wx = new Float64Array(3), wy = new Float64Array(3);

// Quadratic B-spline weights of the three nodes base, base + 1, base + 2 along one axis.
// fx = x / dx - base is the particle's distance from node 'base' in cells, in [0.5, 1.5).
// Its distances to the three nodes are fx, fx - 1 and fx - 2.
function quadraticWeights(fx, w) {
  w[0] = 0.5 * (1.5 - fx) * (1.5 - fx);
  w[1] = 0.75 - (fx - 1) * (fx - 1);
  w[2] = 0.5 * (fx - 0.5) * (fx - 0.5);
}

function reset() {
  n = +cellsAcross.value; dx = 1 / n; invDx = n;
  nodesX = n + 1; nodeCount = nodesX * nodesX;
  gm = new Float64Array(nodeCount); gvx = new Float64Array(nodeCount); gvy = new Float64Array(nodeCount);
  gvxOld = new Float64Array(nodeCount); gvyOld = new Float64Array(nodeCount);

  const k = perAxis.value, cx = scene.value === 'rotating disk' ? 0.5 : 0.3, cy = 0.5;
  const radius = scene.value === 'rotating disk' ? 0.25 : 0.15;
  const xs = [], ys = [];
  for (let j = 0; j < n * k; j++) {
    for (let i = 0; i < n * k; i++) {
      const x = (i + 0.5) * dx / k, y = (j + 0.5) * dx / k;
      if ((x - cx) * (x - cx) + (y - cy) * (y - cy) < radius * radius) { xs.push(x); ys.push(y); }
    }
  }
  count = xs.length;
  px = new Float64Array(xs); py = new Float64Array(ys);
  pvx = new Float64Array(count); pvy = new Float64Array(count);
  pm = new Float64Array(count).fill(density * dx * dx / (k * k));   // each particle owns its share of a cell
  hidden = new Float64Array(count);
  marker = new Uint8Array(count);
  for (let p = 0; p < count; p++) {
    if (scene.value === 'rotating disk') {
      pvx[p] = -spinRate * (py[p] - cy); pvy[p] = spinRate * (px[p] - cx);
    } else {
      pvx[p] = blobSpeed; pvy[p] = 0;
    }
    if (scene.value === 'rotating disk' && Math.abs(py[p] - cy) < 0.6 * dx && px[p] > cx) marker[p] = 1;   // a spoke
  }
  L0 = angularMomentum(); E0 = kineticEnergy(); P0 = momentumX();
}

// ---------------------------------------------------------------------------------------------
// The three transfer stages
// ---------------------------------------------------------------------------------------------

// P2G: m_i = sum_p w_ip m_p and (m v)_i = sum_p w_ip m_p v_p, then v_i = (m v)_i / m_i.
function particlesToGrid() {
  gm.fill(0); gvx.fill(0); gvy.fill(0);
  for (let p = 0; p < count; p++) {
    const gx = px[p] * invDx, gy = py[p] * invDx;
    const bx = Math.floor(gx - 0.5), by = Math.floor(gy - 0.5);
    quadraticWeights(gx - bx, wx);
    quadraticWeights(gy - by, wy);
    for (let b = 0; b < 3; b++) {
      for (let a = 0; a < 3; a++) {
        const w = wx[a] * wy[b], node = (bx + a) + (by + b) * nodesX;
        const mw = w * pm[p];
        gm[node] += mw;
        gvx[node] += mw * pvx[p];
        gvy[node] += mw * pvy[p];
      }
    }
  }
  for (let i = 0; i < nodeCount; i++) {
    if (gm[i] > 0) { gvx[i] /= gm[i]; gvy[i] /= gm[i]; }    // momentum -> velocity; empty nodes stay 0
    gvxOld[i] = gvx[i]; gvyOld[i] = gvy[i];                  // kept for FLIP
  }
}

// Grid update: forces and wall boundary conditions.
// The disk is held together by a pull toward the centre, a = -omega^2 (x - c): exactly the centripetal
// acceleration of a rigid rotation at omega, and it exerts no torque. (A real disk gets this force
// from its internal stress; supplying it from outside means only the transfers are being tested.)
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
      if (i < bound && gvx[node] < 0) gvx[node] = 0;           // no flow out through the walls
      if (i > n - bound && gvx[node] > 0) gvx[node] = 0;
      if (j < bound && gvy[node] < 0) gvy[node] = 0;
      if (j > n - bound && gvy[node] > 0) gvy[node] = 0;
    }
  }
}

// G2P. PIC: v_p = sum_i w_ip v_i. FLIP: v_p += sum_i w_ip (v_i - v_i_old). Blend: mix the two.
function gridToParticles() {
  const alpha = transfer.value === 'PIC' ? 0 : transfer.value === 'FLIP' ? 1 : flipFraction.value;
  for (let p = 0; p < count; p++) {
    const gx = px[p] * invDx, gy = py[p] * invDx;
    const bx = Math.floor(gx - 0.5), by = Math.floor(gy - 0.5);
    quadraticWeights(gx - bx, wx);
    quadraticWeights(gy - by, wy);
    let picX = 0, picY = 0, dvX = 0, dvY = 0;
    for (let b = 0; b < 3; b++) {
      for (let a = 0; a < 3; a++) {
        const w = wx[a] * wy[b], node = (bx + a) + (by + b) * nodesX;
        picX += w * gvx[node]; picY += w * gvy[node];
        dvX += w * (gvx[node] - gvxOld[node]); dvY += w * (gvy[node] - gvyOld[node]);
      }
    }
    const flipX = pvx[p] + dvX, flipY = pvy[p] + dvY;
    pvx[p] = alpha * flipX + (1 - alpha) * picX;
    pvy[p] = alpha * flipY + (1 - alpha) * picY;
    hidden[p] = Math.hypot(pvx[p] - picX, pvy[p] - picY);    // the part of v_p the grid cannot represent
  }
}

// ---------------------------------------------------------------------------------------------
// A stand-in for real material forces (chapter 24): pressure from the grid density,
// p_i = K max(0, rho_i / rho0 - 1), and nodal force f_i = -dx^2 (grad p)_i by central differences.
// ---------------------------------------------------------------------------------------------
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
  const lo = 1 * dx, hi = 1 - 1 * dx;     // keep the 3x3 stencil inside the grid
  for (let p = 0; p < count; p++) {
    px[p] = Math.min(hi, Math.max(lo, px[p] + h * pvx[p]));
    py[p] = Math.min(hi, Math.max(lo, py[p] + h * pvy[p]));
  }
}

function simulate(h) {
  particlesToGrid();
  updateGrid(h);
  gridToParticles();
  advect(h);
}

function step(dt) {
  // The disk has no forces, so one transfer per frame. The blob needs h below dx / c.
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
function angularMomentum() {   // about the centre of the box (0.5, 0.5)
  let L = 0;
  for (let p = 0; p < count; p++) L += pm[p] * ((px[p] - 0.5) * pvy[p] - (py[p] - 0.5) * pvx[p]);
  return L;
}
function kineticEnergy() {
  let E = 0;
  for (let p = 0; p < count; p++) E += 0.5 * pm[p] * (pvx[p] * pvx[p] + pvy[p] * pvy[p]);
  return E;
}
function momentumX() {
  let P = 0;
  for (let p = 0; p < count; p++) P += pm[p] * pvx[p];
  return P;
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
  Sim.show('kinetic energy / start', E / E0, 4);
  if (scene.value === 'rotating disk') {
    Sim.show('angular momentum / start', L / L0, 4);
    Sim.plot('L / L0', L / L0, { color: Sim.colors.blue, min: 0, max: 1 });
  } else {
    Sim.show('x momentum / start', momentumX() / P0, 4);
    Sim.plot('x momentum', momentumX() / P0, { color: Sim.colors.blue, min: -1, max: 1 });
  }
  Sim.plot('KE / KE0', Math.min(E / E0, 2), { color: Sim.colors.red, min: 0 });
}

Sim.run({ title: 'PIC and FLIP transfers', view: [-0.02, -0.32, 1.02, 1.18], reset, step, draw, grid: false });
