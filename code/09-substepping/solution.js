// Chapter 09 - Small Steps (reference solution)
// Two identical chains with a heavy weight. Both do the same number of constraint projections per frame:
// left:  one step of h, with N Gauss-Seidel iterations
// right: N substeps of h/N, with one iteration each
// The plot shows how far each chain is from its exact static stretch, as a percentage of its length.

const count = Sim.slider('n (iterations or substeps)', 1, 40, 10, { step: 1 });
const ratio = Sim.logSlider('mass ratio M / m', 0, 3, 100, { reset: true });
const compliance = Sim.logSlider('compliance alpha (m/N)', -8, -3, 1e-6);
Sim.button('swing both weights', () => { left.vx[links] += 1.5; right.vx[links] += 1.5; });
Sim.button('time 1000 frames of each', benchmark);

const gravity = 9.81;       // m/s^2
const links = 20;
const chainLength = 1.5;    // m
const restLength = chainLength / links;
const beadMass = 0.1;       // kg, every particle except the top one and the weight
const h = 1 / 60;           // frame time step (s)

let left, right, viewBottom = -2.2, timing = 'press the button';

function makeChain(x0) {
  const c = {
    x: new Float64Array(links + 1), y: new Float64Array(links + 1),
    px: new Float64Array(links + 1), py: new Float64Array(links + 1),
    vx: new Float64Array(links + 1), vy: new Float64Array(links + 1),
    w: new Float64Array(links + 1), lambda: new Float64Array(links),
  };
  for (let i = 0; i <= links; i++) {
    c.x[i] = x0; c.y[i] = -i * restLength; c.w[i] = 1 / beadMass;
  }
  c.w[0] = 0;                                   // fixed support
  c.w[links] = 1 / (ratio.value * beadMass);    // the heavy weight
  return c;
}

function reset() {
  left = makeChain(-0.7);
  right = makeChain(0.7);
  viewBottom = -2.2;
}

function predict(c, dt) {
  for (let i = 0; i <= links; i++) {
    c.px[i] = c.x[i]; c.py[i] = c.y[i];
    if (c.w[i] === 0) continue;
    c.vy[i] -= dt * gravity;
    c.x[i] += dt * c.vx[i]; c.y[i] += dt * c.vy[i];
  }
}

// One Gauss-Seidel sweep of XPBD distance constraints, top to bottom.
function solveSweep(c, alphaTilde) {
  const x = c.x, y = c.y, w = c.w, lambda = c.lambda;
  for (let j = 0; j < links; j++) {
    const a = j, b = j + 1;
    const dx = x[b] - x[a], dy = y[b] - y[a];
    const len = Math.hypot(dx, dy);
    const nx = dx / len, ny = dy / len;
    const C = len - restLength;
    const dLambda = (-C - alphaTilde * lambda[j]) / (w[a] + w[b] + alphaTilde);
    lambda[j] += dLambda;
    x[a] -= w[a] * nx * dLambda; y[a] -= w[a] * ny * dLambda;
    x[b] += w[b] * nx * dLambda; y[b] += w[b] * ny * dLambda;
  }
}

function updateVelocities(c, dt) {
  for (let i = 0; i <= links; i++) {
    c.vx[i] = (c.x[i] - c.px[i]) / dt; c.vy[i] = (c.y[i] - c.py[i]) / dt;
  }
}

// Left chain: one big step, N iterations.
function stepIterations(c, dt, N) {
  predict(c, dt);
  c.lambda.fill(0);
  const alphaTilde = compliance.value / (dt * dt);
  for (let it = 0; it < N; it++) solveSweep(c, alphaTilde);
  updateVelocities(c, dt);
}

// Right chain: N small steps of dt / N, one iteration each.
function stepSubsteps(c, dt, N) {
  const hs = dt / N;
  const alphaTilde = compliance.value / (hs * hs);
  for (let s = 0; s < N; s++) {
    predict(c, hs);
    c.lambda.fill(0);
    solveSweep(c, alphaTilde);
    updateVelocities(c, hs);
  }
}

function step(dt) {
  left.w[links] = right.w[links] = 1 / (ratio.value * beadMass);
  stepIterations(left, dt, count.value);
  stepSubsteps(right, dt, count.value);
}

function exactStretch() {
  let s = 0;
  const M = ratio.value * beadMass;
  for (let j = 0; j < links; j++) s += compliance.value * (M + (links - 1 - j) * beadMass) * gravity;
  return s;
}

function stretch(c) {
  let s = 0;
  for (let j = 0; j < links; j++) s += Math.hypot(c.x[j + 1] - c.x[j], c.y[j + 1] - c.y[j]) - restLength;
  return s;
}

// Run each method on a copy of the current state for 1000 frames and report the time per frame.
function benchmark() {
  const frames = 1000, N = count.value;
  const results = [];
  for (const run of [stepIterations, stepSubsteps]) {
    const c = makeChain(0);
    const t0 = performance.now();
    for (let f = 0; f < frames; f++) run(c, h, N);
    results.push((performance.now() - t0) / frames * 1000);
  }
  timing = `${results[0].toFixed(1)} µs vs ${results[1].toFixed(1)} µs (n = ${N})`;
}

function drawChain(g, c, color) {
  const pts = [];
  for (let i = 0; i <= links; i++) pts.push({ x: c.x[i], y: c.y[i] });
  g.polyline(pts, { color: Sim.colors.ink, width: 1.5 });
  for (let i = 1; i < links; i++) g.dot(c.x[i], c.y[i], 2.5, color);
  const s = 0.05 + 0.025 * Math.cbrt(ratio.value * beadMass);
  g.rect(c.x[links] - s, c.y[links] - s, 2 * s, 2 * s, { fill: color });
}

function draw(g) {
  const target = Math.min(-2.2, Math.min(left.y[links], right.y[links]) - 0.4);
  viewBottom += (target - viewBottom) * 0.1;
  Sim.setView([-1.6, viewBottom, 1.6, 0.35]);

  const exact = exactStretch();
  g.rect(-1.0, 0, 2.0, 0.05, { fill: '#d9d6cc' });
  g.line(-1.3, -chainLength - exact, 1.3, -chainLength - exact, { color: Sim.colors.green, dash: [4, 4], width: 1 });
  drawChain(g, left, Sim.colors.orange);
  drawChain(g, right, Sim.colors.blue);
  g.text('1 step x N iterations', -0.7, 0.2, { color: Sim.colors.orange, align: 'center' });
  g.text('N substeps x 1 iteration', 0.7, 0.2, { color: Sim.colors.blue, align: 'center' });

  const errL = 100 * (stretch(left) - exact) / chainLength;
  const errR = 100 * (stretch(right) - exact) / chainLength;
  Sim.show('weight M', (ratio.value * beadMass).toFixed(2) + ' kg');
  Sim.show('exact stretch', (100 * exact).toFixed(3) + ' cm');
  Sim.show('iterations: stretch error', errL.toFixed(3) + ' %');
  Sim.show('substeps: stretch error', errR.toFixed(3) + ' %');
  Sim.show('projections per frame, each', count.value * links);
  Sim.show('time per frame, iter vs sub', timing);
  Sim.plot('iterations %', errL, { color: Sim.colors.orange, min: 0 });
  Sim.plot('substeps %', errR, { color: Sim.colors.blue, min: 0 });
}

Sim.run({ title: 'Iterations or substeps?', view: [-1.6, -2.2, 1.6, 0.35], reset, step, draw, dt: h });
