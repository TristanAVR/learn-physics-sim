// Chapter 09 - Iterative and Direct Solvers (starter)
// Two identical long, stiff chains with a heavy weight, advanced with XPBD.
// left:  Gauss-Seidel, the usual XPBD loop, with a chosen number of iterations per step
// right: one direct solve per step of the linearised system, which is tridiagonal for a chain
//        (optionally with substeps: a linearised solve needs them when the chain moves fast)
// The button at the bottom measures how many sweeps Jacobi and Gauss-Seidel need on a uniform rigid chain.
// The direct solver needs a tridiagonal solve: implement thomasSolve below.

const links = Sim.slider('links n', 2, 100, 50, { step: 1, reset: true });
const iterations = Sim.slider('Gauss-Seidel iterations', 1, 500, 50, { step: 1 });
const ratio = Sim.logSlider('mass ratio M / m', 0, 3, 100, { reset: true });
const compliance = Sim.logSlider('compliance alpha (m/N)', -9, -4, 1e-7);
const directSubsteps = Sim.slider('direct: substeps per frame', 1, 16, 1, { step: 1 });
Sim.button('swing both weights', () => { left.vx[left.n] += 1.5; right.vx[right.n] += 1.5; });
Sim.button('count sweeps to 1 % (uniform rigid chain)', () => { countText = countSweeps(links.value, 1); });

const gravity = 9.81;       // m/s^2
const chainLength = 2.0;    // m
const beadMass = 0.02;      // kg, every particle except the top one and the weight
const h = 1 / 60;           // s

let left, right, viewBottom = -2.8, countText = 'press the button';

function makeChain(x0, n) {
  const c = {
    n, rest: chainLength / n,
    x: new Float64Array(n + 1), y: new Float64Array(n + 1),
    px: new Float64Array(n + 1), py: new Float64Array(n + 1),
    vx: new Float64Array(n + 1), vy: new Float64Array(n + 1),
    w: new Float64Array(n + 1), lambda: new Float64Array(n),
    // scratch for the direct solver: link directions, constraint values, the three diagonals, rhs, solution
    nx: new Float64Array(n), ny: new Float64Array(n), C: new Float64Array(n),
    a: new Float64Array(n), b: new Float64Array(n), c: new Float64Array(n),
    d: new Float64Array(n), dl: new Float64Array(n),
  };
  for (let i = 0; i <= n; i++) { c.x[i] = x0; c.y[i] = -i * c.rest; c.w[i] = 1 / beadMass; }
  c.w[0] = 0;
  c.w[n] = 1 / (ratio.value * beadMass);
  return c;
}

function reset() {
  left = makeChain(-0.8, links.value);
  right = makeChain(0.8, links.value);
  viewBottom = -2.8;
}

function predict(c, dt) {
  for (let i = 0; i <= c.n; i++) {
    c.px[i] = c.x[i]; c.py[i] = c.y[i];
    if (c.w[i] === 0) continue;
    c.vy[i] -= dt * gravity;
    c.x[i] += dt * c.vx[i]; c.y[i] += dt * c.vy[i];
  }
}

function updateVelocities(c, dt) {
  for (let i = 0; i <= c.n; i++) { c.vx[i] = (c.x[i] - c.px[i]) / dt; c.vy[i] = (c.y[i] - c.py[i]) / dt; }
}

// ---------------------------------------------------------------------------------------------
// Left chain: XPBD with Gauss-Seidel sweeps over the links.
// ---------------------------------------------------------------------------------------------
function stepGaussSeidel(c, dt, iters) {
  predict(c, dt);
  c.lambda.fill(0);
  const at = compliance.value / (dt * dt);
  const x = c.x, y = c.y, w = c.w;
  for (let it = 0; it < iters; it++) {
    for (let j = 0; j < c.n; j++) {
      const dx = x[j + 1] - x[j], dy = y[j + 1] - y[j];
      const len = Math.hypot(dx, dy);
      const nx = dx / len, ny = dy / len;
      const dLambda = (-(len - c.rest) - at * c.lambda[j]) / (w[j] + w[j + 1] + at);
      c.lambda[j] += dLambda;
      x[j] -= w[j] * nx * dLambda; y[j] -= w[j] * ny * dLambda;
      x[j + 1] += w[j + 1] * nx * dLambda; y[j + 1] += w[j + 1] * ny * dLambda;
    }
  }
  updateVelocities(c, dt);
}

// ---------------------------------------------------------------------------------------------
// Right chain: one Newton step per (sub)step. Linearise every link at the predicted positions and
// solve (J W J^T + alphaTilde) dLambda = -C - alphaTilde lambda for all links at once.
// For a chain, link j only shares a particle with links j - 1 and j + 1, so the matrix is tridiagonal.
// ---------------------------------------------------------------------------------------------
function stepDirect(c, dt, substeps) {
  const hs = dt / substeps;
  for (let s = 0; s < substeps; s++) {
    predict(c, hs);
    c.lambda.fill(0);
    newtonStep(c, hs);
    updateVelocities(c, hs);
  }
}

function newtonStep(c, dt) {
  const n = c.n, at = compliance.value / (dt * dt);
  const x = c.x, y = c.y, w = c.w;
  for (let j = 0; j < n; j++) {
    const dx = x[j + 1] - x[j], dy = y[j + 1] - y[j];
    const len = Math.hypot(dx, dy);
    c.nx[j] = dx / len; c.ny[j] = dy / len;
    c.C[j] = len - c.rest;
  }
  for (let j = 0; j < n; j++) {
    // grad C_j is -n_j on particle j and +n_j on particle j + 1.
    c.b[j] = w[j] + w[j + 1] + at;
    c.a[j] = j > 0 ? -w[j] * (c.nx[j - 1] * c.nx[j] + c.ny[j - 1] * c.ny[j]) : 0;          // shares particle j with link j - 1
    c.c[j] = j < n - 1 ? -w[j + 1] * (c.nx[j] * c.nx[j + 1] + c.ny[j] * c.ny[j + 1]) : 0;  // shares particle j + 1 with link j + 1
    c.d[j] = -c.C[j] - at * c.lambda[j];
  }
  thomasSolve(c.a, c.b, c.c, c.d, c.dl, n);
  // dx_i = w_i sum_j grad_i C_j dLambda_j: particle i ends link i - 1 and starts link i.
  for (let i = 1; i <= n; i++) {
    let gx = c.nx[i - 1] * c.dl[i - 1], gy = c.ny[i - 1] * c.dl[i - 1];
    if (i < n) { gx -= c.nx[i] * c.dl[i]; gy -= c.ny[i] * c.dl[i]; }
    x[i] += w[i] * gx; y[i] += w[i] * gy;
  }
  for (let j = 0; j < n; j++) c.lambda[j] += c.dl[j];
}

// Solve a tridiagonal system: a[i] x[i-1] + b[i] x[i] + c[i] x[i+1] = d[i], for i = 0 .. n-1.
// (a[0] and c[n-1] are ignored.) Forward elimination, then back substitution: O(n).
const cp = new Float64Array(101), dp = new Float64Array(101);
function thomasSolve(a, b, c, d, x, n) {
  // TODO: implement the Thomas algorithm (forward elimination, then back substitution).
  // Use the scratch arrays cp and dp for the modified super-diagonal and right-hand side.
  // This placeholder ignores the off-diagonals: it is a single Jacobi sweep from zero, under-relaxed
  // by a factor 0.5, because a full Jacobi step overshoots and blows the chain up.
  for (let i = 0; i < n; i++) x[i] = 0.5 * d[i] / b[i];
}

// ---------------------------------------------------------------------------------------------
// Convergence experiment: a rigid chain of n links hanging straight down, bead masses equal except
// the last, which is massRatio times heavier. Right after the prediction only the top link is stretched,
// by h^2 g. Solve J W J^T dLambda = -C exactly (Thomas), then count the Jacobi and Gauss-Seidel sweeps
// needed to get within 1 % of it, and the observed error reduction per sweep.
// ---------------------------------------------------------------------------------------------
function countSweeps(n, massRatio) {
  const w = new Float64Array(n + 1).fill(1 / beadMass);
  w[0] = 0; w[n] = 1 / (massRatio * beadMass);
  const a = new Float64Array(n), b = new Float64Array(n), c = new Float64Array(n), d = new Float64Array(n);
  for (let j = 0; j < n; j++) { b[j] = w[j] + w[j + 1]; a[j] = -w[j]; c[j] = -w[j + 1]; }
  d[0] = -h * h * gravity;
  const exact = new Float64Array(n);
  thomasSolve(a, b, c, d, exact, n);
  let norm = 0;
  for (let j = 0; j < n; j++) norm += exact[j] * exact[j];
  norm = Math.sqrt(norm);

  function run(jacobi) {
    let z = new Float64Array(n), zNew = new Float64Array(n);
    const errs = [];
    for (let k = 1; k <= 200000; k++) {
      for (let j = 0; j < n; j++) {
        const src = jacobi ? z : zNew;       // Gauss-Seidel uses the values updated in this sweep
        let s = d[j];
        if (j > 0) s -= a[j] * src[j - 1];
        if (j < n - 1) s -= c[j] * z[j + 1];
        zNew[j] = s / b[j];
      }
      const t = z; z = zNew; zNew = t;
      if (!jacobi) zNew.set(z);
      let e = 0;
      for (let j = 0; j < n; j++) e += (z[j] - exact[j]) * (z[j] - exact[j]);
      errs.push(Math.sqrt(e) / norm);
      if (errs[k - 1] < 0.01) {
        const rate = k > 2 ? Math.sqrt(errs[k - 1] / errs[k - 3]) : NaN;
        return { k, rate };
      }
    }
    return null;
  }
  const gs = run(false), jac = run(true);
  if (!gs || !jac) return 'did not converge: is thomasSolve right?';
  return `n = ${n}: GS ${gs.k} sweeps (x${gs.rate.toFixed(5)} per sweep), Jacobi ${jac.k} (x${jac.rate.toFixed(5)})`;
}

function step(dt) {
  left.w[left.n] = right.w[right.n] = 1 / (ratio.value * beadMass);
  stepGaussSeidel(left, dt, iterations.value);
  stepDirect(right, dt, directSubsteps.value);
}

function exactStretch(n) {
  const M = ratio.value * beadMass;
  let s = 0;
  for (let j = 0; j < n; j++) s += compliance.value * (M + (n - 1 - j) * beadMass) * gravity;
  return s;
}

function stretch(c) {
  let s = 0;
  for (let j = 0; j < c.n; j++) s += Math.hypot(c.x[j + 1] - c.x[j], c.y[j + 1] - c.y[j]) - c.rest;
  return s;
}

function drawChain(g, c, color) {
  const pts = [];
  for (let i = 0; i <= c.n; i++) pts.push({ x: c.x[i], y: c.y[i] });
  g.polyline(pts, { color: Sim.colors.ink, width: 1.5 });
  for (let i = 1; i < c.n; i++) g.dot(c.x[i], c.y[i], 2, color);
  const s = 0.05 + 0.04 * Math.cbrt(ratio.value * beadMass);
  g.rect(c.x[c.n] - s, c.y[c.n] - s, 2 * s, 2 * s, { fill: color });
}

function draw(g) {
  const target = Math.max(-40, Math.min(-2.8, Math.min(left.y[left.n], right.y[right.n]) - 0.4));
  if (isFinite(target)) viewBottom += (target - viewBottom) * 0.1;
  Sim.setView([-1.8, viewBottom, 1.8, 0.4]);

  const exact = exactStretch(left.n);
  g.rect(-1.2, 0, 2.4, 0.05, { fill: '#d9d6cc' });
  g.line(-1.5, -chainLength - exact, 1.5, -chainLength - exact, { color: Sim.colors.green, dash: [4, 4], width: 1 });
  drawChain(g, left, Sim.colors.orange);
  drawChain(g, right, Sim.colors.blue);
  g.text('Gauss-Seidel', -0.8, 0.22, { color: Sim.colors.orange, align: 'center' });
  g.text('direct (tridiagonal)', 0.8, 0.22, { color: Sim.colors.blue, align: 'center' });

  const eL = 1000 * (stretch(left) - exact), eR = 1000 * (stretch(right) - exact);
  Sim.show('exact stretch', (1000 * exact).toFixed(3) + ' mm');
  Sim.show('Gauss-Seidel stretch error', eL.toFixed(3) + ' mm');
  Sim.show('direct stretch error', Math.abs(eR) < 1e5 ? eR.toFixed(4) + ' mm' : 'exploded: press R');
  Sim.show('sweeps', countText);
  Sim.plot('Gauss-Seidel error (mm)', eL, { color: Sim.colors.orange, min: 0 });
  Sim.plot('direct error (mm)', eR, { color: Sim.colors.blue, min: 0 });
}

Sim.run({ title: 'Gauss-Seidel vs a direct solve (starter)', view: [-1.8, -2.8, 1.8, 0.4], reset, step, draw, dt: h });
