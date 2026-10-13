// Chapter 13 - Cosserat Rods (reference solution)
// A rod made of rigid segments. Each segment has a centre (x, y) AND an angle theta, so the rod carries its
// material frame with it. Two constraints per joint:
//   stretch-shear: the end of segment i must meet the start of segment i+1 (couples positions and angles)
//   bend:          theta[i+1] - theta[i] must equal the rest angle (rest curvature times length), compliance l / EI
// Drag any segment with the mouse (grab it off-centre to twist it), curl it with the rest curvature, or crack the whip.

const restCurvature = Sim.slider('rest curvature (1/m)', -12, 12, 0, { step: 0.1 });
const bendStiffness = Sim.logSlider('bending stiffness EI (N m^2)', -2, 1, 0.5);
const segments = Sim.slider('segments N', 4, 24, 16, { step: 1, reset: true });
const gravityOn = Sim.checkbox('gravity', true);
const whip = Sim.checkbox('whip: shake the base', false, { reset: true });

const L = 1.0;           // length of the free rod (m)
const mu = 0.1;          // mass per length (kg/m)
const gravity = 9.81;
const damping = 1.0;     // velocity damping rate (1/s)
const iterations = 2;    // solver sweeps per substep
const dragCompliance = 1e-4;
const whipAmplitude = 0.12, whipFrequency = 2.0; // base motion (m, Hz)

// Segment 0 is the clamped base (kinematic: infinite mass and inertia), segments 1..N are free.
// Segment i has centre (x, y), angle th (its direction d = (cos th, sin th)), length l.
let N, l, x, y, th, px, py, pth, vx, vy, om, w, wI, lambdaJoint, lambdaBend, drag, time;

function reset() {
  N = segments.value;
  l = L / N;
  const n = N + 1;
  x = new Float64Array(n); y = new Float64Array(n); th = new Float64Array(n);
  px = new Float64Array(n); py = new Float64Array(n); pth = new Float64Array(n);
  vx = new Float64Array(n); vy = new Float64Array(n); om = new Float64Array(n);
  w = new Float64Array(n); wI = new Float64Array(n);
  lambdaJoint = new Float64Array(n); lambdaBend = new Float64Array(n);
  const m = mu * l, inertia = m * l * l / 12; // thin rod about its centre
  // Lay the segments out along the rest shape, starting at the wall (x = 0) heading right.
  let ex = 0, ey = 0;
  x[0] = -0.5 * l; y[0] = 0; th[0] = 0;
  for (let i = 1; i < n; i++) {
    th[i] = restCurvature.value * (i - 0.5) * l;
    x[i] = ex + 0.5 * l * Math.cos(th[i]); y[i] = ey + 0.5 * l * Math.sin(th[i]);
    ex += l * Math.cos(th[i]); ey += l * Math.sin(th[i]);
    w[i] = 1 / m; wI[i] = 1 / inertia;
  }
  drag = null;
  time = 0;
}

// The length of rod that joint i (between segment i and i+1) stands for. The base is rigid, so the first joint
// only stands for half of segment 1 (the series rule from chapter 12).
function jointLength(i) { return i === 0 ? 0.5 * l : l; }

// Stretch-shear: the end point of segment i, x_i + (l/2) d_i, must coincide with the start point of segment
// i+1, x_j - (l/2) d_j. Thin rods are so stiff in stretch and shear (EA, GA) that we make it rigid (alpha = 0).
function solveStretchShear(i, h) {
  const j = i + 1;
  const rax = 0.5 * l * Math.cos(th[i]), ray = 0.5 * l * Math.sin(th[i]);    // offset of the end of i
  const rbx = -0.5 * l * Math.cos(th[j]), rby = -0.5 * l * Math.sin(th[j]);  // offset of the start of j
  const dx = (x[j] + rbx) - (x[i] + rax), dy = (y[j] + rby) - (y[i] + ray);
  const C = Math.hypot(dx, dy);
  if (C < 1e-12) return;
  const nx = dx / C, ny = dy / C;
  // Generalised inverse masses: 1/m + (r x n)^2 / I for each body (chapter 11).
  const rna = rax * ny - ray * nx, rnb = rbx * ny - rby * nx;
  const wa = w[i] + wI[i] * rna * rna, wb = w[j] + wI[j] * rnb * rnb;
  if (wa + wb === 0) return;
  const dlambda = -C / (wa + wb);
  lambdaJoint[i] += dlambda;
  const pxImp = dlambda * nx, pyImp = dlambda * ny;   // positional impulse on j; minus it on i
  x[j] += w[j] * pxImp; y[j] += w[j] * pyImp; th[j] += wI[j] * (rbx * pyImp - rby * pxImp);
  x[i] -= w[i] * pxImp; y[i] -= w[i] * pyImp; th[i] -= wI[i] * (rax * pyImp - ray * pxImp);
}

// Bend: C = (th_j - th_i) - restAngle. In 2D the "Darboux vector" is a single number: the relative rotation of
// the two frames. The rest angle is rest curvature times the length the joint stands for.
function solveBend(i, h) {
  const j = i + 1;
  const len = jointLength(i);
  const C = th[j] - th[i] - restCurvature.value * len;
  const alphaTilde = (len / bendStiffness.value) / (h * h);
  const wsum = wI[i] + wI[j];
  const dlambda = (-C - alphaTilde * lambdaBend[i]) / (wsum + alphaTilde);
  lambdaBend[i] += dlambda;
  th[i] -= wI[i] * dlambda;
  th[j] += wI[j] * dlambda;
}

// Mouse: a soft constraint pulling a point fixed in the grabbed segment's frame towards the mouse.
// It is the one-body version of solveStretchShear (the world has infinite mass).
function solveDrag(h) {
  const i = drag.seg, c = Math.cos(th[i]), s = Math.sin(th[i]);
  const rx = c * drag.lx - s * drag.ly, ry = s * drag.lx + c * drag.ly;
  const dx = Sim.mouse.x - (x[i] + rx), dy = Sim.mouse.y - (y[i] + ry);
  const C = Math.hypot(dx, dy);
  if (C < 1e-9) return;
  const nx = dx / C, ny = dy / C;
  const rn = rx * ny - ry * nx;
  const wa = w[i] + wI[i] * rn * rn;
  const alphaTilde = dragCompliance / (h * h);
  const dlambda = (C - alphaTilde * drag.lambda) / (wa + alphaTilde);
  drag.lambda += dlambda;
  const pX = dlambda * nx, pY = dlambda * ny;
  x[i] += w[i] * pX; y[i] += w[i] * pY; th[i] += wI[i] * (rx * pY - ry * pX);
}

// Keep h * omega_joint below 0.5 for the bend joints: omega^2 = (wI_i + wI_j) / alpha = 24 EI / (mu l^4).
function substepCount(dt) {
  const omega = Math.sqrt(24 * bendStiffness.value / (mu * l * l * l * l)) * Math.SQRT2; // the root joint is 2x stiffer
  return Math.min(1000, Math.max(10, Math.ceil(omega * dt / 0.5)));
}

function step(dt) {
  if (Sim.mouse.pressed) {
    let best = -1, bestD = 0.1;
    for (let i = 1; i <= N; i++) {
      const d = Math.hypot(Sim.mouse.x - x[i], Sim.mouse.y - y[i]);
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best > 0) {
      const c = Math.cos(th[best]), s = Math.sin(th[best]);
      const dx = Sim.mouse.x - x[best], dy = Sim.mouse.y - y[best];
      drag = { seg: best, lx: c * dx + s * dy, ly: -s * dx + c * dy, lambda: 0 };
    }
  }
  if (!Sim.mouse.down) drag = null;

  const ns = substepCount(dt), h = dt / ns;
  const decay = 1 / (1 + damping * h);
  for (let s = 0; s < ns; s++) {
    time += h;
    if (whip.value) y[0] = whipAmplitude * Math.sin(2 * Math.PI * whipFrequency * time);
    for (let i = 1; i <= N; i++) {
      if (gravityOn.value) vy[i] -= h * gravity;
      vx[i] *= decay; vy[i] *= decay; om[i] *= decay;
      px[i] = x[i]; py[i] = y[i]; pth[i] = th[i];
      x[i] += h * vx[i]; y[i] += h * vy[i]; th[i] += h * om[i];
    }
    lambdaJoint.fill(0); lambdaBend.fill(0);
    if (drag) drag.lambda = 0;
    for (let it = 0; it < iterations; it++) {
      for (let i = 0; i < N; i++) solveStretchShear(i, h);
      for (let i = 0; i < N; i++) solveBend(i, h);
      if (drag) solveDrag(h);
    }
    for (let i = 1; i <= N; i++) {
      vx[i] = (x[i] - px[i]) / h; vy[i] = (y[i] - py[i]) / h; om[i] = (th[i] - pth[i]) / h;
    }
  }
}

function draw(g) {
  // the rest shape, for reference (a circular arc of radius 1 / kappa)
  const rest = [{ x: 0, y: y[0] }];
  let ex = 0, ey = y[0];
  for (let i = 1; i <= N; i++) {
    const a = restCurvature.value * (i - 0.5) * l;
    ex += l * Math.cos(a); ey += l * Math.sin(a);
    rest.push({ x: ex, y: ey });
  }
  g.polyline(rest, { color: Sim.colors.gray, width: 1, dash: [4, 4] });

  g.rect(-0.12, y[0] - 0.12, 0.12, 0.24, { fill: '#d9d6cc' });
  for (let i = 1; i <= N; i++) {
    const c = Math.cos(th[i]), s = Math.sin(th[i]);
    const hx = 0.5 * l * c, hy = 0.5 * l * s;
    g.line(x[i] - hx, y[i] - hy, x[i] + hx, y[i] + hy,
           { color: i % 2 ? Sim.colors.blue : Sim.colors.teal, width: 6 });
    // the material frame: d1 (normal) drawn as a short tick, so you can see each segment's orientation
    g.line(x[i], y[i], x[i] - 0.35 * l * s, y[i] + 0.35 * l * c, { color: Sim.colors.orange, width: 1.5 });
  }
  if (drag) {
    const i = drag.seg, c = Math.cos(th[i]), s = Math.sin(th[i]);
    const gx = x[i] + c * drag.lx - s * drag.ly, gy = y[i] + s * drag.lx + c * drag.ly;
    g.line(gx, gy, Sim.mouse.x, Sim.mouse.y, { color: Sim.colors.red, width: 1, dash: [3, 3] });
    g.dot(gx, gy, 3, Sim.colors.red);
  }

  const tipX = x[N] + 0.5 * l * Math.cos(th[N]), tipY = y[N] + 0.5 * l * Math.sin(th[N]);
  g.dot(tipX, tipY, 3, Sim.colors.ink);
  Sim.show('substeps per frame', substepCount(1 / 60));
  Sim.show('tip angle', (th[N] * 180 / Math.PI).toFixed(1) + ' deg');
  if (restCurvature.value === 0 && !whip.value) {
    const theory = gravityOn.value ? mu * gravity * Math.pow(L, 4) / (8 * bendStiffness.value) : 0;
    Sim.show('tip deflection', (-1000 * (tipY - y[0])).toFixed(1) + ' mm');
    Sim.show('w L^4 / (8 EI)', (1000 * theory).toFixed(1) + ' mm');
  }
  let gap = 0;
  for (let i = 0; i < N; i++) {
    const j = i + 1;
    gap = Math.max(gap, Math.hypot(x[j] - 0.5 * l * Math.cos(th[j]) - x[i] - 0.5 * l * Math.cos(th[i]),
                                   y[j] - 0.5 * l * Math.sin(th[j]) - y[i] - 0.5 * l * Math.sin(th[i])));
  }
  Sim.show('largest joint gap', (1000 * gap).toFixed(2) + ' mm');
}

Sim.run({ title: 'A rod of rigid segments', view: [-0.25, -0.85, 1.25, 0.55], reset, step, draw, grid: 0.1 });
