// Chapter 21 - The Deformation Gradient (starter)
// A square of material, sampled by points and cut into triangles.
// Mode "set F by hand": every point moves by x = F X, with F built from the sliders.
// Mode "flow field": every point is carried by a velocity field and integrates its own F.
// Colours show J = det F: blue where the material is compressed, red where it is stretched,
// purple where it is inverted (J < 0).

const mode = Sim.select('mode', ['set F by hand', 'flow field'], 'set F by hand', { reset: true });
const stretchX = Sim.slider('F: stretch x', 0.2, 2.5, 1.5, { step: 0.01 });
const stretchY = Sim.slider('F: stretch y', -1.5, 2.5, 1, { step: 0.01 });
const shear = Sim.slider('F: shear', -1.5, 1.5, 0, { step: 0.01 });
const angleDeg = Sim.slider('F: rotation (degrees)', -180, 180, 0, { step: 1 });
const field = Sim.select('flow: velocity field', ['rigid rotation', 'shear flow', 'expansion', 'vortex'],
  'rigid rotation', { reset: true });
const rate = Sim.slider('flow: rate (1/s)', 0.1, 3, 1, { step: 0.05, reset: true });
const stepMs = Sim.slider('flow: time step h (ms)', 1, 300, 16, { step: 1 });
const colourBy = Sim.select('colour by', ['J = det F', 'strain |E|'], 'J = det F');

// ---------------------------------------------------------------------------------------------
// 2x2 matrices are plain arrays [a, b, c, d] meaning [[a, b], [c, d]] (row by row).
// ---------------------------------------------------------------------------------------------
function mul(A, B) {
  return [A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3],
          A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3]];
}
function transpose(A) { return [A[0], A[2], A[1], A[3]]; }
function det(A) { return A[0] * A[3] - A[1] * A[2]; }
function inverse(A) {
  const d = det(A);
  return [A[3] / d, -A[1] / d, -A[2] / d, A[0] / d];
}
function rotation(theta) {
  const c = Math.cos(theta), s = Math.sin(theta);
  return [c, -s, s, c];
}

// Green strain E = (F^T F - I) / 2. Zero for any rotation.
function greenStrain(F) {
  const C = mul(transpose(F), F);
  return [0.5 * (C[0] - 1), 0.5 * C[1], 0.5 * C[2], 0.5 * (C[3] - 1)];
}
function frobenius(A) { return Math.hypot(A[0], A[1], A[2], A[3]); }

// Eigen-decomposition of a symmetric 2x2 matrix S = [[p, q], [q, r]] = V diag(l1, l2) V^T,
// with V a rotation by phi and l1 >= l2.
function symmetricEigen(S) {
  const p = S[0], q = 0.5 * (S[1] + S[2]), r = S[3];
  const mean = 0.5 * (p + r), radius = Math.hypot(0.5 * (p - r), q);
  const phi = 0.5 * Math.atan2(2 * q, p - r);
  return { l1: mean + radius, l2: mean - radius, phi: phi };
}

// ---------------------------------------------------------------------------------------------
// The three exercises
// ---------------------------------------------------------------------------------------------

// Deformation gradient of a triangle from its rest corners (X0, X1, X2) and current corners
// (x0, x1, x2), each an {x, y}. The edge matrices hold edge vectors as columns:
// Dm = [X1 - X0, X2 - X0], Ds = [x1 - x0, x2 - x0], and F maps one to the other: Ds = F Dm.
function triangleF(X0, X1, X2, x0, x1, x2) {
  // TODO (exercise 1): build Dm and Ds as 2x2 arrays [a, b, c, d] whose COLUMNS are the edge
  // vectors, then return F = Ds Dm^-1 (use mul and inverse above).
  return [1, 0, 0, 1];
}

// Polar decomposition F = R S, with R a rotation and S symmetric.
// R^T F must be symmetric, which fixes the rotation angle: theta = atan2(c - b, a + d).
function polar(F) {
  // TODO (exercise 2): find the angle theta of the rotation R, then S = R^T F.
  // Hint: write R^T F for R = rotation(theta) and ask when it is symmetric.
  const theta = 0;
  const R = rotation(theta);
  const S = mul(transpose(R), F);
  return { R: R, S: S, theta: theta };
}

// One explicit step of dF/dt = (grad v) F for the point whose F is stored at Fs[4k .. 4k+3]:
// F <- (I + h L) F, where L = grad v is the velocity gradient at the point.
function integrateF(Fs, k, L, h) {
  const o = 4 * k;
  // TODO (exercise 3): read F = [Fs[o], Fs[o+1], Fs[o+2], Fs[o+3]], compute (I + h L) F and
  // write the result back into the same four slots.
}

// SVD F = U diag(s1, s2) V^T from the polar decomposition: S = V diag(s1, s2) V^T, U = R V.
// U and V are rotations; s2 is negative when the material is inverted (det F < 0).
function svd(F) {
  const pd = polar(F);
  const e = symmetricEigen(pd.S);
  return { s1: e.l1, s2: e.l2, U: mul(pd.R, rotation(e.phi)), V: rotation(e.phi), R: pd.R, theta: pd.theta };
}

// ---------------------------------------------------------------------------------------------
// Velocity fields: velocity v(x, y) and its gradient L = grad v = [[dvx/dx, dvx/dy], [dvy/dx, dvy/dy]]
// ---------------------------------------------------------------------------------------------
const vortexRadius = 0.6; // m
function velocity(x, y, out) {
  const w = rate.value;
  if (field.value === 'rigid rotation') { out.x = -w * y; out.y = w * x; }
  else if (field.value === 'shear flow') { out.x = w * y; out.y = 0; }
  else if (field.value === 'expansion') { out.x = w * x; out.y = w * y; }
  else {
    const f = w / (1 + (x * x + y * y) / (vortexRadius * vortexRadius));
    out.x = -f * y; out.y = f * x;
  }
}
function velocityGradient(x, y) {
  const w = rate.value;
  if (field.value === 'rigid rotation') return [0, -w, w, 0];
  if (field.value === 'shear flow') return [0, w, 0, 0];
  if (field.value === 'expansion') return [w, 0, 0, w];
  // vortex: v = f(r^2) (-y, x) with f = w / (1 + r^2 / r0^2); fp = df/d(r^2)
  const a2 = vortexRadius * vortexRadius, q = 1 + (x * x + y * y) / a2;
  const f = w / q, fp = -w / (a2 * q * q);
  return [-2 * x * y * fp, -f - 2 * y * y * fp, f + 2 * x * x * fp, 2 * x * y * fp];
}

// ---------------------------------------------------------------------------------------------
// The material: an n x n lattice of points over a unit square, two triangles per lattice cell.
// ---------------------------------------------------------------------------------------------
const n = 16, side = 1.0;
const count = (n + 1) * (n + 1);
const restX = new Float64Array(count), restY = new Float64Array(count);
const posX = new Float64Array(count), posY = new Float64Array(count);
const pointF = new Float64Array(4 * count);
const tris = [];
for (let j = 0; j <= n; j++) {
  for (let i = 0; i <= n; i++) {
    restX[i + j * (n + 1)] = -0.5 * side + side * i / n;
    restY[i + j * (n + 1)] = -0.5 * side + side * j / n;
  }
}
for (let j = 0; j < n; j++) {
  for (let i = 0; i < n; i++) {
    const a = i + j * (n + 1), b = a + 1, c = a + n + 1, d = c + 1;
    tris.push([a, b, d], [a, d, c]);
  }
}
const centre = (n / 2) + (n / 2) * (n + 1);
const pointColors = new Array(count);
let steps = 0;

function reset() {
  for (let k = 0; k < count; k++) {
    posX[k] = restX[k]; posY[k] = restY[k];
    pointF[4 * k] = 1; pointF[4 * k + 1] = 0; pointF[4 * k + 2] = 0; pointF[4 * k + 3] = 1;
  }
  steps = 0;
}

function handF() {
  // Stretch and shear first (upper-triangular part), then rotate.
  const K = [stretchX.value, shear.value, 0, stretchY.value];
  return mul(rotation(angleDeg.value * Math.PI / 180), K);
}

const vel = { x: 0, y: 0 };
function step() {
  if (mode.value !== 'flow field') return;
  const h = stepMs.value / 1000;
  for (let k = 0; k < count; k++) {
    const L = velocityGradient(posX[k], posY[k]);
    velocity(posX[k], posY[k], vel);
    integrateF(pointF, k, L, h);      // F and x use the velocity field at the OLD position
    posX[k] += h * vel.x;
    posY[k] += h * vel.y;
  }
  steps++;
}

// ---------------------------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------------------------
// Diverging colour map: log2(J) from -1 (blue) through 0 (pale) to +1 (red); J <= 0 is purple.
const ramp = [];
for (let i = 0; i <= 64; i++) {
  const t = i / 64 * 2 - 1, u = Math.abs(t);
  const pale = [238, 236, 228], end = t < 0 ? [47, 111, 223] : [217, 71, 58];
  const c = pale.map((p, m) => Math.round(p + (end[m] - p) * u));
  ramp.push(`rgb(${c[0]},${c[1]},${c[2]})`);
}
function colourJ(J) {
  if (!(J > 0)) return '#5b2a86';
  const t = Math.max(-1, Math.min(1, Math.log2(J)));
  return ramp[Math.round((t + 1) * 32)];
}
function colourStrain(e) {
  if (!isFinite(e)) return '#5b2a86';
  return ramp[32 + Math.round(Math.min(1, e / 0.5) * 32)];
}
function colourFor(F) {
  return colourBy.value === 'J = det F' ? colourJ(det(F)) : colourStrain(frobenius(greenStrain(F)));
}

const P = (k) => ({ x: restX[k], y: restY[k] });
const Q = (k) => ({ x: posX[k], y: posY[k] });

function drawMesh(g) {
  let firstF = null;
  for (const t of tris) {
    const F = triangleF(P(t[0]), P(t[1]), P(t[2]), Q(t[0]), Q(t[1]), Q(t[2]));
    if (!firstF) firstF = F;
    g.polyline([posX[t[0]], posY[t[0]], posX[t[1]], posY[t[1]], posX[t[2]], posY[t[2]]],
      { closed: true, fill: colourFor(F), color: 'rgba(60,60,60,0.25)', width: 0.6 });
  }
  return firstF;
}

function drawEllipse(g, F, cx, cy, r) {
  const circle = [], ellipse = [];
  for (let i = 0; i <= 64; i++) {
    const a = i / 64 * 2 * Math.PI, X = r * Math.cos(a), Y = r * Math.sin(a);
    circle.push(cx + X, cy + Y);
    ellipse.push(cx + F[0] * X + F[1] * Y, cy + F[2] * X + F[3] * Y);
  }
  g.polyline(circle, { color: Sim.colors.muted, dash: [4, 4], width: 1 });
  g.polyline(ellipse, { color: Sim.colors.ink, width: 2 });
  // Principal stretches: the columns of U, scaled by the singular values.
  const s = svd(F);
  g.arrow(cx, cy, r * s.s1 * s.U[0], r * s.s1 * s.U[2], { color: Sim.colors.orange, width: 2 });
  g.arrow(cx, cy, r * s.s2 * s.U[1], r * s.s2 * s.U[3], { color: Sim.colors.green, width: 2 });
  return s;
}

function showMatrix(label, A) {
  Sim.show(label, `[${A[0].toFixed(3)}, ${A[1].toFixed(3)}; ${A[2].toFixed(3)}, ${A[3].toFixed(3)}]`);
}

function draw(g) {
  if (mode.value === 'set F by hand') {
    const F = handF();
    for (let k = 0; k < count; k++) {
      posX[k] = F[0] * restX[k] + F[1] * restY[k];
      posY[k] = F[2] * restX[k] + F[3] * restY[k];
    }
    g.polyline([-0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5], { closed: true, color: Sim.colors.muted, dash: [5, 4], width: 1 });
    const Ftri = drawMesh(g);
    g.points(posX, posY, count, { radius: 1.6, color: Sim.colors.ink });
    const s = drawEllipse(g, F, 0, 0, 0.3);
    const pd = polar(F);
    const RS = mul(pd.R, pd.S);
    showMatrix('F', F);
    Sim.show('J = det F', det(F), 3);
    Sim.show('R: angle (deg)', (pd.theta * 180 / Math.PI).toFixed(1));
    Sim.show('principal stretches', `${s.s1.toFixed(3)}, ${s.s2.toFixed(3)}`);
    showMatrix('Green strain E', greenStrain(F));
    Sim.show('|R S - F|', frobenius([RS[0] - F[0], RS[1] - F[1], RS[2] - F[2], RS[3] - F[3]]));
    showMatrix('F of triangle 0', Ftri);
  } else {
    drawMesh(g);
    for (let k = 0; k < count; k++) pointColors[k] = colourFor(pointF.subarray(4 * k, 4 * k + 4));
    g.points(posX, posY, count, { radius: 2.2, colors: pointColors });
    // Images of the material axes at the centre point: the columns of its F.
    const o = 4 * centre, cx = posX[centre], cy = posY[centre], r = 0.3;
    g.arrow(cx, cy, r * pointF[o], r * pointF[o + 2], { color: Sim.colors.orange, width: 2 });
    g.arrow(cx, cy, r * pointF[o + 1], r * pointF[o + 3], { color: Sim.colors.green, width: 2 });
    let sumJ = 0, minJ = Infinity, maxJ = -Infinity, sumE = 0;
    for (let k = 0; k < count; k++) {
      const F = pointF.subarray(4 * k, 4 * k + 4), J = det(F);
      sumJ += J; minJ = Math.min(minJ, J); maxJ = Math.max(maxJ, J);
      sumE += frobenius(greenStrain(F));
    }
    const Fc = Array.from(pointF.subarray(o, o + 4));
    showMatrix('F (centre point)', Fc);
    Sim.show('J (centre point)', det(Fc), 5);
    Sim.show('J min / max', `${minJ.toFixed(4)} / ${maxJ.toFixed(4)}`);
    Sim.show('mean strain |E|', sumE / count, 5);
    Sim.show('steps', steps);
    Sim.plot('J centre', Math.min(det(Fc), 3), { color: Sim.colors.red });
    Sim.plot('mean |E|', Math.min(sumE / count, 3), { color: Sim.colors.purple, min: 0 });
  }
}

Sim.run({
  title: 'The deformation gradient (starter)', view: [-2.2, -1.6, 2.2, 1.6], reset, step, draw,
  dt: () => stepMs.value / 1000, grid: 0.25,
});
