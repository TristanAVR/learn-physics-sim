// Chapter 33 - Verifying a Simulator (reference solution)
// A convergence-study tool. It runs a test problem to a fixed time T at eleven step sizes with the chosen
// integrator, measures the error against the exact (spring) or a very fine reference (pendulum) solution,
// and draws error against step size on a log-log plot. The slope of that line is the order of accuracy.
// Three regression tests run alongside; plant a bug and see which of them catches it.

const integ = Sim.select('integrator', ['explicit Euler', 'symplectic Euler', 'midpoint (RK2)', 'velocity Verlet', 'RK4'],
  'symplectic Euler', { reset: true });
const problem = Sim.select('test problem', ['pendulum (fine reference)', 'spring (exact)'], 'pendulum (fine reference)', { reset: true });
const bug = Sim.select('planted bug', ['none', 'damping sign flipped', 'Verlet uses old acceleration',
  'off-by-one final time', 'g = 9.8 in the integrator'], 'none', { reset: true });
const hAnim = Sim.logSlider('animation step h (s)', -3, -0.7, 0.05);

const G = 9.81, Lp = 1.0, theta0 = 1.0;        // pendulum: 1 m, released from 1 rad
const W = 2 * Math.PI;                         // spring: period 1 s
const T = 1.7;                                 // compare here: mid-swing, not at a turning point
const expectedOrder = { 'explicit Euler': 1, 'symplectic Euler': 1, 'midpoint (RK2)': 2, 'velocity Verlet': 2, 'RK4': 4 };

// ---------------------------------------------------------------------------------------------
// Accelerations. The planted bugs live here and in the integrators, as real bugs would.
// ---------------------------------------------------------------------------------------------
function pendulumAcc(q, v) {
  const g = bug.value === 'g = 9.8 in the integrator' ? 9.8 : G;
  return -(g / Lp) * Math.sin(q);
}
function springAcc(zeta) {
  const sign = bug.value === 'damping sign flipped' ? +1 : -1;
  return (q, v) => -W * W * q + sign * 2 * zeta * W * v;
}

// One step of each integrator for q'' = acc(q, v). Returns [q, v].
function stepWith(name, q, v, h, acc) {
  if (name === 'explicit Euler') {
    const a = acc(q, v);
    return [q + h * v, v + h * a];
  }
  if (name === 'symplectic Euler') {
    const v1 = v + h * acc(q, v);
    return [q + h * v1, v1];
  }
  if (name === 'midpoint (RK2)') {
    const a = acc(q, v), qm = q + 0.5 * h * v, vm = v + 0.5 * h * a;
    return [q + h * vm, v + h * acc(qm, vm)];
  }
  if (name === 'velocity Verlet') {
    const a0 = acc(q, v), q1 = q + h * v + 0.5 * h * h * a0;
    if (bug.value === 'Verlet uses old acceleration') return [q1, v + h * a0];
    const a1 = acc(q1, v + h * a0);                  // velocity predicted for the damping term
    return [q1, v + 0.5 * h * (a0 + a1)];
  }
  // RK4
  const k1q = v, k1v = acc(q, v);
  const k2q = v + 0.5 * h * k1v, k2v = acc(q + 0.5 * h * k1q, v + 0.5 * h * k1v);
  const k3q = v + 0.5 * h * k2v, k3v = acc(q + 0.5 * h * k2q, v + 0.5 * h * k2v);
  const k4q = v + h * k3v, k4v = acc(q + h * k3q, v + h * k3v);
  return [q + h / 6 * (k1q + 2 * k2q + 2 * k3q + k4q), v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)];
}

// Run from (q0, v0) to time Tend in N steps. (The off-by-one bug takes one step too few.)
function runTo(name, acc, q0, v0, Tend, N) {
  const h = Tend / N, steps = bug.value === 'off-by-one final time' ? N - 1 : N;
  let q = q0, v = v0;
  for (let i = 0; i < steps; i++) [q, v] = stepWith(name, q, v, h, acc);
  return [q, v];
}

// ---------------------------------------------------------------------------------------------
// The convergence study
// ---------------------------------------------------------------------------------------------
let study, tests, anim;

function exactOrReference() {
  if (problem.value === 'spring (exact)') return Math.cos(W * T);
  // No closed form for the pendulum at 1 rad: RK4 with 40,000 steps (error far below 1e-12).
  let q = theta0, v = 0;
  const N = 40000, h = T / N, acc = (q, v) => -(G / Lp) * Math.sin(q);
  for (let i = 0; i < N; i++) [q, v] = stepWith('RK4', q, v, h, acc);
  return q;
}

// Errors at T for N = 10, 20, 40, ..., 10240 steps.
function convergenceStudy(name) {
  const acc = problem.value === 'spring (exact)' ? springAcc(0) : pendulumAcc;
  const q0 = problem.value === 'spring (exact)' ? 1 : theta0;
  const ref = exactOrReference(), points = [];
  for (let k = 0; k <= 10; k++) {
    const N = 10 * Math.pow(2, k);
    const [q] = runTo(name, acc, q0, 0, T, N);
    points.push({ h: T / N, err: Math.abs(q - ref) });
  }
  return points;
}

// Least-squares slope of log10(err) against log10(h) over the four finest points that are still above
// round-off (1e-11): the asymptotic regime.
function fitSlope(points) {
  const use = points.filter((p) => p.err > 1e-11 && Number.isFinite(p.err)).slice(-4);
  if (use.length < 2) return NaN;
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (const p of use) {
    const x = Math.log10(p.h), y = Math.log10(p.err);
    sx += x; sy += y; sxx += x * x; sxy += x * y;
  }
  const n = use.length;
  return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

// ---------------------------------------------------------------------------------------------
// Regression tests, each with a tolerance and a closed form to compare with.
// ---------------------------------------------------------------------------------------------
function runTests(name) {
  const out = [];
  // T1: the measured order matches the integrator's order to within 0.2.
  const slope = study.slope, want = expectedOrder[name];
  out.push({ name: 'order of accuracy', pass: Math.abs(slope - want) < 0.2,
    detail: 'slope ' + slope.toFixed(2) + ', expected ' + want });

  // T2: energy of the undamped spring after 100 s at h = 0.01 s, against the theory for this integrator.
  // Energies are averaged over the first and last 100 steps (one period) to remove the bounded wobble.
  const h = 0.01, n = 10000, x = h * W, acc = springAcc(0);
  let q = 1, v = 0, e0 = 0, e1 = 0;
  const E = (q, v) => 0.5 * v * v + 0.5 * W * W * q * q;
  for (let i = 0; i < n; i++) {
    if (i < 100) e0 += E(q, v);
    [q, v] = stepWith(name, q, v, h, acc);
    if (i >= n - 100) e1 += E(q, v);
  }
  const drift = e1 / e0 - 1, m = n - 100;
  const predicted = {
    'explicit Euler': Math.pow(1 + x * x, m) - 1,                              // |lambda|^2 = 1 + x^2
    'symplectic Euler': 0, 'velocity Verlet': 0,                               // bounded, no drift
    'midpoint (RK2)': Math.pow(1 + Math.pow(x, 4) / 4, m) - 1,
    'RK4': Math.pow(1 - Math.pow(x, 6) / 72 + Math.pow(x, 8) / 576, m) - 1,
  }[name];
  out.push({ name: 'energy drift', pass: Math.abs(drift - predicted) <= 1e-4 + 0.01 * Math.abs(predicted),
    detail: 'drift ' + Sim.format(drift) + ', theory ' + Sim.format(predicted) });

  // T3: a damped spring (zeta = 0.1) at T with h = 0.001 s against the closed form, tolerance 0.01.
  const z = 0.1, wd = W * Math.sqrt(1 - z * z);
  const exact = Math.exp(-z * W * T) * (Math.cos(wd * T) + z * W / wd * Math.sin(wd * T));
  const [qd] = runTo(name, springAcc(z), 1, 0, T, 1700);
  out.push({ name: 'damped decay', pass: Math.abs(qd - exact) < 0.01,
    detail: 'x(T) ' + qd.toFixed(4) + ', exact ' + exact.toFixed(4) });
  return out;
}

function reset() {
  const points = convergenceStudy(integ.value);
  study = { points, slope: fitSlope(points) };
  tests = runTests(integ.value);
  anim = { t: 0, acc: 0, q: problem.value === 'spring (exact)' ? 1 : theta0, v: 0,
    rq: problem.value === 'spring (exact)' ? 1 : theta0, rv: 0 };
}

// The animation: the chosen integrator at a coarse step next to a fine RK4 run.
function step(dt) {
  const isSpring = problem.value === 'spring (exact)';
  const acc = isSpring ? springAcc(0) : pendulumAcc;
  const refAcc = isSpring ? ((q) => -W * W * q) : ((q) => -(G / Lp) * Math.sin(q));
  anim.acc += dt;
  while (anim.acc >= hAnim.value) {
    [anim.q, anim.v] = stepWith(integ.value, anim.q, anim.v, hAnim.value, acc);
    anim.acc -= hAnim.value;
  }
  for (let i = 0; i < 20; i++) [anim.rq, anim.rv] = stepWith('RK4', anim.rq, anim.rv, dt / 20, refAcc);
  anim.t += dt;
}

// ---------------------------------------------------------------------------------------------
// Drawing: the animation on the left, the log-log plot on the right (drawn with polylines).
// ---------------------------------------------------------------------------------------------
const plotBox = { x0: 0.35, x1: 3.35, y0: -1.1, y1: 0.8, lhMin: -4, lhMax: -0.5, leMin: -16, leMax: 1 };
function px(lh) { return plotBox.x0 + (lh - plotBox.lhMin) / (plotBox.lhMax - plotBox.lhMin) * (plotBox.x1 - plotBox.x0); }
function py(le) {
  const c = Math.max(plotBox.leMin, Math.min(plotBox.leMax, le));
  return plotBox.y0 + (c - plotBox.leMin) / (plotBox.leMax - plotBox.leMin) * (plotBox.y1 - plotBox.y0);
}

function draw(g) {
  const C = Sim.colors, isSpring = problem.value === 'spring (exact)';
  // Animation, top right: the chosen integrator (red) at a coarse step and a fine reference (grey).
  const ox = 2.85, oy = 1.68, len = 0.55;
  if (isSpring) {
    g.line(2.2, 1.0, 2.2, 1.7, { color: C.gray, width: 3 });
    [[anim.rq, C.gray, 1.2], [anim.q, C.red, 1.5]].forEach(([q, col, y]) => {
      const x = 2.85 + 0.25 * Math.max(-2.2, Math.min(2.2, Number.isFinite(q) ? q : 0));
      g.line(2.2, y, x, y, { color: col, width: 1, dash: [3, 3] });
      g.rect(x - 0.06, y - 0.06, 0.12, 0.12, { fill: col });
    });
  } else {
    g.rect(ox - 0.12, oy, 0.24, 0.03, { fill: '#d9d6cc' });
    for (const [q, col] of [[anim.rq, C.gray], [anim.q, C.red]]) {
      if (!Number.isFinite(q)) continue;
      const bx = ox + len * Math.sin(q), by = oy - len * Math.cos(q);
      g.line(ox, oy, bx, by, { color: col, width: 1.5 });
      g.circle(bx, by, 0.05, { fill: col });
    }
  }
  g.text('red: h = ' + hAnim.value.toFixed(3) + ' s,  grey: reference', 3.35, 0.93, { color: C.muted, size: 11, align: 'right' });

  // Plot frame, decades, labels.
  const B = plotBox;
  g.rect(B.x0, B.y0, B.x1 - B.x0, B.y1 - B.y0, { fill: '#ffffff' });
  for (let lh = B.lhMin; lh <= B.lhMax; lh += 1) {
    g.line(px(lh), B.y0, px(lh), B.y1, { color: C.grid, width: 1 });
    g.text('1e' + lh, px(lh), B.y0 - 0.06, { color: C.muted, size: 10, align: 'center' });
  }
  for (let le = B.leMin; le <= B.leMax; le += 2) {
    g.line(B.x0, py(le), B.x1, py(le), { color: C.grid, width: 1 });
    g.text('1e' + le, B.x0 - 0.04, py(le), { color: C.muted, size: 10, align: 'right' });
  }
  g.rect(B.x0, B.y0, B.x1 - B.x0, B.y1 - B.y0, { color: C.axis, width: 1 });
  g.text('step size h (s)', 0.5 * (B.x0 + B.x1), B.y0 - 0.16, { color: C.muted, size: 11, align: 'center' });
  g.text('error at T = ' + T + ' s', B.x0, B.y1 + 0.06, { color: C.muted, size: 11 });

  // Reference slopes 1, 2 and 4 through the finest valid point.
  const pts = study.points.filter((p) => p.err > 0 && Number.isFinite(p.err));
  const anchor = pts.filter((p) => p.err > 1e-11).slice(-1)[0];
  if (anchor) {
    const ax = Math.log10(anchor.h), ay = Math.log10(anchor.err);
    for (const [s, col] of [[1, C.blue], [2, C.green], [4, C.purple]]) {
      const lh1 = Math.min(B.lhMax, ax + (B.leMax - ay) / s);      // stop where the line leaves the box
      const ey = ay + s * (lh1 - ax);
      g.line(px(ax), py(ay), px(lh1), py(ey), { color: col, width: 1, dash: [5, 4] });
      g.text('slope ' + s, px(lh1) + 0.03, py(ey) - 0.06, { color: col, size: 10, align: 'right' });
    }
  }
  // Measured errors.
  const line = pts.map((p) => ({ x: px(Math.log10(p.h)), y: py(Math.log10(p.err)) }));
  if (line.length > 1) g.polyline(line, { color: C.red, width: 2 });
  for (const p of line) g.dot(p.x, p.y, 3.5, C.red);

  Sim.show('integrator', integ.value);
  Sim.show('measured order (fit)', Number.isFinite(study.slope) ? study.slope.toFixed(3) : '-');
  Sim.show('expected order', expectedOrder[integ.value]);
  for (const t of tests) Sim.show((t.pass ? 'PASS ' : 'FAIL ') + t.name, t.detail);
}

Sim.run({ title: 'Convergence study', view: [-0.25, -1.45, 3.45, 1.75], reset, step, draw, dt: 1 / 60, grid: false });
