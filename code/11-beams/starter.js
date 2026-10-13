// Chapter 11 - Bending, Twisting, Buckling (starter)
// A rod made of particles. Distance constraints keep the segments at length l, so right now it hangs like a chain.
// Your job: write the bending constraint on every triple of particles (solveBend) and the series-rule compliance
// of each joint (jointCompliance). Then the rod should stand out from the wall like a real cantilever.
// Cantilever: compare the tip with P L^3/(3 EI) + w L^4/(8 EI). Buckling: push down on the top and watch the rod
// stay straight below P_cr = pi^2 EI / (2L)^2 and bow out above it.

const scene = Sim.select('scene', ['cantilever', 'buckling'], 'cantilever', { reset: true });
const rodType = Sim.select('rod', ['uniform', 'stiff base, soft tip'], 'uniform', { reset: true });
const rule = Sim.select('joint compliance', ['series rule', 'naive (average EI)'], 'series rule', { reset: true });
const segments = Sim.slider('segments N', 2, 30, 10, { step: 1, reset: true });
const load = Sim.slider('load P (N)', 0, 6, 0, { step: 0.05 });
const gravityOn = Sim.checkbox('self-weight (cantilever)', true);
const substepMode = Sim.select('substeps per frame', ['auto', '10', '30', '100'], 'auto');

const L = 1.0;          // rod length (m)
const EI = 1.25;        // bending stiffness of the soft part (N m^2)
const stiffRatio = 4;   // the stiff base is this many times stiffer
const mu = 0.1;         // mass per length (kg/m)
const gravity = 9.81;
const damping = 3.0;    // velocity damping rate (1/s), only so the rod settles
const iterations = 4;   // solver sweeps per substep
const imperfection = 0.002; // tiny sideways force at the tip in buckling mode (N)

// Nodes: 0 is a fixed "ghost" behind the clamp, 1 is the clamped root, 2..N+1 are free. Node k (k >= 1) sits at
// arc length s = (k - 1) l. Segment k joins nodes k and k+1; segment 0 (ghost to root) is the rigid clamp.
let N, l, n, x, y, px, py, vx, vy, w, lambdaBend, lambdaDist, alphaJoint;

function segmentEI(k) {
  if (k === 0) return Infinity; // the clamp does not bend
  const s = (k - 0.5) * l;      // arc length of the segment's midpoint
  if (rodType.value === 'stiff base, soft tip' && s < 0.5 * L) return stiffRatio * EI;
  return EI;
}

// Compliance of the joint at node k, which sits between segment k-1 and segment k.
function jointCompliance(k) {
  const EIa = segmentEI(k - 1), EIb = segmentEI(k);
  if (rule.value === 'naive (average EI)') return l / (0.5 * (EIa + EIb));
  // TODO: the series rule. The joint stands for half of segment k-1 and half of segment k; compliances in series
  // add, and a rigid half (EI = Infinity, like the clamp) contributes nothing.
  return 0;
}

function reset() {
  N = segments.value;
  l = L / N;
  n = N + 2;
  x = new Float64Array(n); y = new Float64Array(n);
  px = new Float64Array(n); py = new Float64Array(n);
  vx = new Float64Array(n); vy = new Float64Array(n);
  w = new Float64Array(n);
  lambdaBend = new Float64Array(n);
  lambdaDist = new Float64Array(n);
  alphaJoint = new Float64Array(n);
  const vertical = scene.value === 'buckling';
  for (let k = 0; k < n; k++) {
    const s = (k - 1) * l;
    x[k] = vertical ? 0 : s;
    y[k] = vertical ? s : 0;
    const m = k === n - 1 ? 0.5 * mu * l : mu * l; // lumped masses, half at the free end
    w[k] = k <= 1 ? 0 : 1 / m;
  }
  for (let k = 1; k <= N; k++) alphaJoint[k] = jointCompliance(k);
  Sim.setView(vertical ? [-1.0, -0.55, 1.0, 1.2] : [-0.2, -0.9, 1.3, 0.55]);
}

// Keep nodes k and k+1 at distance l (rigid: EA is so large that we treat it as infinite).
function solveDistance(k, h) {
  const a = k, b = k + 1;
  const wSum = w[a] + w[b];
  if (wSum === 0) return;
  const dx = x[b] - x[a], dy = y[b] - y[a];
  const len = Math.hypot(dx, dy);
  const C = len - l;
  const nx = dx / len, ny = dy / len;
  const dlambda = -C / wSum;
  lambdaDist[k] += dlambda;
  x[a] -= w[a] * nx * dlambda; y[a] -= w[a] * ny * dlambda;
  x[b] += w[b] * nx * dlambda; y[b] += w[b] * ny * dlambda;
}

// Bending constraint on the triple (k-1, k, k+1): C = turning angle between the two segments (rest angle 0).
function solveBend(k, h) {
  const a = k - 1, b = k, c = k + 1;
  // TODO:
  // 1. e1 = x_b - x_a, e2 = x_c - x_b. C = signed angle from e1 to e2 = atan2(cross(e1, e2), dot(e1, e2)).
  // 2. Gradients: dC/dx_a = perp(e1)/|e1|^2, dC/dx_c = perp(e2)/|e2|^2, dC/dx_b = -(dC/dx_a + dC/dx_c),
  //    with perp(e) = (-e.y, e.x).
  // 3. The XPBD update with alphaTilde = alphaJoint[k] / h^2 and the accumulated multiplier lambdaBend[k]:
  //    dlambda = (-C - alphaTilde * lambda) / (sum of w |grad|^2 + alphaTilde); move each particle by w * grad * dlambda.
}

// The bending joints get stiffer as 1/l and the particles lighter as l, so a joint's own oscillation gets faster
// as N^2. A Gauss-Seidel sweep only finds the right equilibrium if the substep resolves it: we keep h omega <= 0.5,
// with omega^2 = (w |grad C|^2) / alpha = 6 / (mu l^3 alpha) for the stiffest joint.
function substepCount(dt) {
  if (substepMode.value !== 'auto') return +substepMode.value;
  let omegaMax = 0;
  for (let k = 1; k <= N; k++) {
    if (alphaJoint[k] > 0) omegaMax = Math.max(omegaMax, Math.sqrt(6 / (mu * l * l * l * alphaJoint[k])));
  }
  return Math.min(800, Math.max(20, Math.ceil(omegaMax * dt / 0.5)));
}

function step(dt) {
  const ns = substepCount(dt), h = dt / ns;
  const tip = n - 1;
  const buckling = scene.value === 'buckling';
  for (let s = 0; s < ns; s++) {
    for (let k = 0; k < n; k++) {
      if (w[k] === 0) continue;
      let fx = 0, fy = 0;
      if (!buckling && gravityOn.value) fy -= gravity / w[k];      // weight m g
      if (k === tip) {
        fy -= load.value;                                         // load P, always straight down
        if (buckling) fx -= imperfection;   // a tiny push to the left
      }
      vx[k] += h * w[k] * fx; vy[k] += h * w[k] * fy;
      const d = 1 / (1 + damping * h);
      vx[k] *= d; vy[k] *= d;
      px[k] = x[k]; py[k] = y[k];
      x[k] += h * vx[k]; y[k] += h * vy[k];
    }
    lambdaBend.fill(0); lambdaDist.fill(0);
    for (let it = 0; it < iterations; it++) {
      for (let k = 1; k <= N; k++) solveDistance(k, h);
      for (let k = 1; k <= N; k++) solveBend(k, h);
    }
    for (let k = 0; k < n; k++) {
      if (w[k] === 0) continue;
      vx[k] = (x[k] - px[k]) / h; vy[k] = (y[k] - py[k]) / h;
    }
  }
}

// Small-deflection beam theory for the cantilever: kappa = M / EI(s), theta = integral of kappa,
// deflection = integral of theta. For a uniform rod this gives P L^3/(3 EI) + w L^4/(8 EI) at the tip.
function analyticCantilever() {
  const P = load.value, q = gravityOn.value ? mu * gravity : 0;
  const steps = 400, ds = L / steps, pts = [{ x: 0, y: 0 }];
  let theta = 0, defl = 0;
  for (let i = 0; i < steps; i++) {
    const s = (i + 0.5) * ds, r = L - s;
    const M = P * r + 0.5 * q * r * r;
    const stiff = rodType.value === 'stiff base, soft tip' && s < 0.5 * L ? stiffRatio * EI : EI;
    const thetaMid = theta + 0.5 * ds * M / stiff;
    theta += ds * M / stiff;
    defl += ds * thetaMid;
    pts.push({ x: s + 0.5 * ds, y: -defl });
  }
  return { pts, tip: defl };
}

// Euler's critical load for a uniform fixed-free column (effective length K L with K = 2).
function eulerLoad() {
  if (rodType.value !== 'uniform') return NaN;
  return Math.PI * Math.PI * EI / (4 * L * L);
}

function draw(g) {
  const buckling = scene.value === 'buckling';
  const tip = n - 1;
  // clamp
  if (buckling) g.rect(-0.25, -0.05, 0.5, 0.05, { fill: '#d9d6cc' });
  else g.rect(-0.12, -0.25, 0.12, 0.5, { fill: '#d9d6cc' });

  const pts = [];
  for (let k = 1; k < n; k++) pts.push({ x: x[k], y: y[k] });
  g.polyline(pts, { color: Sim.colors.blue, width: 5, alpha: 0.8 });
  for (let k = 1; k < n; k++) {
    const s = (k - 1) * l;
    const stiffPart = rodType.value === 'stiff base, soft tip' && s < 0.5 * L - 1e-9;
    g.dot(x[k], y[k], 3.5, stiffPart ? Sim.colors.purple : Sim.colors.ink);
  }

  if (!buckling) {
    const a = analyticCantilever();
    g.polyline(a.pts, { color: Sim.colors.red, width: 1.5, dash: [6, 4] });
    const tipDefl = -y[tip];
    Sim.show('tip deflection (sim)', (1000 * tipDefl).toFixed(1) + ' mm');
    Sim.show('beam theory', (1000 * a.tip).toFixed(1) + ' mm');
    Sim.show('sim / theory', a.tip > 0 ? (tipDefl / a.tip).toFixed(4) : '-');
    Sim.show('substeps per frame', substepCount(1 / 60));
    Sim.plot('tip deflection (mm)', 1000 * tipDefl, { color: Sim.colors.blue, min: 0 });
  } else {
    const Pcr = eulerLoad();
    Sim.show('P_cr = pi^2 EI/(2L)^2', isFinite(Pcr) ? Pcr.toFixed(3) + ' N' : 'non-uniform rod');
    if (isFinite(Pcr)) Sim.show('P / P_cr', (load.value / Pcr).toFixed(3));
    Sim.show('tip sideways', (-1000 * x[tip]).toFixed(1) + ' mm');
    Sim.plot('tip sideways (mm)', -1000 * x[tip], { color: Sim.colors.blue, min: 0 });
  }

  const arrowLen = Math.min(0.35, 0.06 * load.value + (load.value > 0 ? 0.05 : 0));
  if (arrowLen > 0) g.arrow(x[tip], y[tip] + arrowLen, 0, -arrowLen + 0.01, { color: Sim.colors.red, width: 2 });
  g.text(buckling ? 'load P stays vertical' : 'dashed: beam theory',
         buckling ? 0.15 : -0.15, buckling ? 1.15 : -0.8, { color: Sim.colors.muted });
}

Sim.run({ title: 'Cantilever and column (starter)', view: [-0.2, -0.9, 1.3, 0.55], reset, step, draw, grid: 0.1, plotCorner: 'top-right' });
