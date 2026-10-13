// Chapter 07 - XPBD: Compliance (reference solution)
// An elastic chain hangs from a fixed point with a heavy weight at the bottom. Every link has the same
// compliance alpha (metres of stretch per newton). Compare the measured stretch with the exact static
// stretch, and see which solver reaches it, and with how many iterations.

const method = Sim.select('solver', ['PBD', 'XPBD'], 'XPBD', { reset: true });
const links = Sim.slider('links n', 1, 20, 10, { step: 1, reset: true });
const compliance = Sim.logSlider('compliance alpha (m/N)', -6, -2.5, 1e-3);
const iterations = Sim.slider('iterations per step', 1, 200, 20, { step: 1 });
const weight = Sim.slider('end weight M (kg)', 0.1, 10, 2, { step: 0.1 });
const stepMs = Sim.slider('time step h (ms)', 4, 33, 1000 / 60, { step: 0.1 });
const damping = Sim.slider('damping beta (N s/m)', 0, 200, 0, { step: 1 });
Sim.button('kick the weight sideways', () => { vx[n] += 2; });

const gravity = 9.81;       // m/s^2
const chainLength = 1.5;    // total rest length (m)
const beadMass = 0.1;       // mass of every particle except the top one and the weight (kg)

let n, restLength;
let x, y, px, py, vx, vy, w, lambda;
let viewBottom = -2.2;

function reset() {
  n = links.value;
  restLength = chainLength / n;
  x = new Float64Array(n + 1); y = new Float64Array(n + 1);
  px = new Float64Array(n + 1); py = new Float64Array(n + 1);
  vx = new Float64Array(n + 1); vy = new Float64Array(n + 1);
  w = new Float64Array(n + 1);
  lambda = new Float64Array(n);
  for (let i = 0; i <= n; i++) {
    x[i] = 0;
    y[i] = -i * restLength;            // hanging straight down, every link at rest length
    w[i] = 1 / beadMass;
  }
  w[0] = 0;                            // the top particle is fixed: infinite mass
  w[n] = 1 / weight.value;
  viewBottom = -2.2;
}

// Link j joins particle j (above) and particle j + 1 (below).
// C = |x_b - x_a| - L, with gradient n for x_b and -n for x_a, where n is the unit vector from a to b.
function solveLink(j, alphaTilde, gamma) {
  const a = j, b = j + 1;
  const dx = x[b] - x[a], dy = y[b] - y[a];
  const len = Math.hypot(dx, dy);
  if (len < 1e-12) return;
  const nx = dx / len, ny = dy / len;
  const C = len - restLength;
  const wSum = w[a] + w[b];            // sum of w_i |grad_i C|^2, since |n| = 1
  // Damping: the rate of change of C over this step, times h. grad C . (x - x_prev)
  const dC = nx * ((x[b] - px[b]) - (x[a] - px[a])) + ny * ((y[b] - py[b]) - (y[a] - py[a]));
  const dLambda = (-C - alphaTilde * lambda[j] - gamma * dC) / ((1 + gamma) * wSum + alphaTilde);
  lambda[j] += dLambda;
  x[a] -= w[a] * nx * dLambda; y[a] -= w[a] * ny * dLambda;
  x[b] += w[b] * nx * dLambda; y[b] += w[b] * ny * dLambda;
}

function step(h) {
  w[n] = 1 / weight.value;
  const xpbd = method.value === 'XPBD';
  // PBD is XPBD with zero compliance (and no damping term).
  const alphaTilde = xpbd ? compliance.value / (h * h) : 0;
  const gamma = xpbd ? compliance.value * damping.value / h : 0;

  // 1. predict
  for (let i = 0; i <= n; i++) {
    px[i] = x[i]; py[i] = y[i];
    if (w[i] === 0) continue;
    vy[i] -= h * gravity;
    x[i] += h * vx[i]; y[i] += h * vy[i];
  }
  // 2. project, with the multipliers reset at the start of every step
  lambda.fill(0);
  for (let it = 0; it < iterations.value; it++) {
    for (let j = 0; j < n; j++) solveLink(j, alphaTilde, gamma);
  }
  // 3. update velocities
  for (let i = 0; i <= n; i++) {
    vx[i] = (x[i] - px[i]) / h; vy[i] = (y[i] - py[i]) / h;
  }
}

// Link j carries the weight plus every bead below it: (M + (n - 1 - j) m) g.
function linkLoad(j) {
  return (weight.value + (n - 1 - j) * beadMass) * gravity;
}

function exactStretch() {
  let s = 0;
  for (let j = 0; j < n; j++) s += compliance.value * linkLoad(j);
  return s;
}

function measuredStretch() {
  let s = 0;
  for (let j = 0; j < n; j++) s += Math.hypot(x[j + 1] - x[j], y[j + 1] - y[j]) - restLength;
  return s;
}

function draw(g) {
  const h = stepMs.value / 1000;
  const exact = exactStretch(), measured = measuredStretch();

  // Keep the weight in view when the chain stretches a lot.
  const target = Math.min(-2.2, y[n] - 0.4);
  viewBottom += (target - viewBottom) * 0.1;
  Sim.setView([-1.3, viewBottom, 1.3, 0.3]);

  g.rect(-0.3, 0, 0.6, 0.06, { fill: '#d9d6cc' });
  // Where the bottom would be at rest length, and where the exact stretch puts it.
  g.line(-0.5, -chainLength, 0.5, -chainLength, { color: Sim.colors.muted, dash: [4, 4], width: 1 });
  g.text('rest length', 0.55, -chainLength, { color: Sim.colors.muted, size: 11 });
  g.line(-0.5, -chainLength - exact, 0.5, -chainLength - exact, { color: Sim.colors.green, dash: [4, 4], width: 1 });
  g.text('exact XPBD stretch', 0.55, -chainLength - exact, { color: Sim.colors.green, size: 11 });

  const pts = [];
  for (let i = 0; i <= n; i++) pts.push({ x: x[i], y: y[i] });
  g.polyline(pts, { color: Sim.colors.ink, width: 2 });
  for (let i = 1; i < n; i++) g.circle(x[i], y[i], 0.025, { fill: Sim.colors.blue });
  const s = 0.06 + 0.03 * Math.cbrt(weight.value);
  g.rect(x[n] - s, y[n] - s, 2 * s, 2 * s, { fill: Sim.colors.red });

  // Stiffness ratio S = (w_a + w_b) h^2 / alpha for the stiffest-looking link.
  let maxS = 0;
  for (let j = 0; j < n; j++) maxS = Math.max(maxS, (w[j] + w[j + 1]) * h * h / compliance.value);

  Sim.show('exact stretch (sum alpha * load)', (100 * exact).toFixed(2) + ' cm');
  Sim.show('measured stretch', (100 * measured).toFixed(2) + ' cm');
  Sim.show('measured / exact', (measured / exact).toFixed(3));
  Sim.show('top link force |lambda| / h^2', (Math.abs(lambda[0]) / (h * h)).toFixed(2) + ' N');
  Sim.show('top link load (exact)', linkLoad(0).toFixed(2) + ' N');
  Sim.show('largest S = w h^2 / alpha', maxS < 100 ? maxS.toFixed(2) : maxS.toExponential(1));
  Sim.plot('stretch (cm)', 100 * measured, { color: Sim.colors.blue, min: 0 });
  Sim.plot('exact (cm)', 100 * exact, { color: Sim.colors.green, min: 0 });
}

Sim.run({ title: 'A hanging chain: PBD vs XPBD', view: [-1.3, -2.2, 1.3, 0.3], reset, step, draw, dt: () => stepMs.value / 1000 });
