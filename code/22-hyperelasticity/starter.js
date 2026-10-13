// Chapter 22 - Stress and Hyperelasticity (starter)
// A block of jelly meshed with linear triangles. Each triangle computes its deformation gradient
// F = Ds Dm^-1, a first Piola-Kirchhoff stress P(F) from the chosen material model, and the
// nodal forces f = -V0 P Dm^-T. Symplectic Euler with as many substeps as stability needs.
// Drag the block with the mouse. "spin" sets it rotating; "slam" throws it at the floor.

const model = Sim.select('material model', ['linear', 'St Venant-Kirchhoff', 'fixed corotated', 'neo-Hookean'],
  'fixed corotated');
const youngs = Sim.logSlider("Young's modulus E (Pa)", 4.5, 6, 2e5);
const poisson = Sim.slider("Poisson's ratio nu", 0, 0.45, 0.3, { step: 0.01 });
const damping = Sim.slider('damping (1/s)', 0, 5, 1, { step: 0.1 });
const gravityOn = Sim.checkbox('gravity', true);
Sim.button('spin (8 rad/s)', () => spin(8));
Sim.button('slam into the floor', () => { for (let k = 0; k < count; k++) vy[k] = -slamSpeed; });

const density = 1000;     // kg/m^3; the 2D block is a slab 1 m thick, so areas are volumes
const gravity = 9.81;
const slamSpeed = 6;       // m/s

// ---------------------------------------------------------------------------------------------
// 2x2 matrices as arrays [a, b, c, d] = [[a, b], [c, d]]
// ---------------------------------------------------------------------------------------------
function mul(A, B) {
  return [A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3],
          A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3]];
}
function transpose(A) { return [A[0], A[2], A[1], A[3]]; }
function det(A) { return A[0] * A[3] - A[1] * A[2]; }
function inverse(A) { const d = det(A); return [A[3] / d, -A[1] / d, -A[2] / d, A[0] / d]; }
// Cofactor matrix J F^-T: no division, so it stays finite when the triangle is flat or inverted.
function cofactor(F) { return [F[3], -F[2], -F[1], F[0]]; }
// Rotation part of the polar decomposition F = R S (chapter 21).
function polarRotation(F) {
  const theta = Math.atan2(F[2] - F[1], F[0] + F[3]);
  const c = Math.cos(theta), s = Math.sin(theta);
  return [c, -s, s, c];
}

// Lame parameters from Young's modulus and Poisson's ratio.
function lame(E, nu) {
  return { mu: E / (2 * (1 + nu)), lambda: E * nu / ((1 + nu) * (1 - 2 * nu)) };
}

// ---------------------------------------------------------------------------------------------
// Material models. Each returns { P, psi }: first Piola-Kirchhoff stress and energy density.
// ---------------------------------------------------------------------------------------------

// Linear elasticity: small strain eps = (F + F^T)/2 - I. Not rotation invariant.
function linearP(F, mu, lam) {
  const e00 = F[0] - 1, e01 = 0.5 * (F[1] + F[2]), e11 = F[3] - 1, tr = e00 + e11;
  const P = [2 * mu * e00 + lam * tr, 2 * mu * e01, 2 * mu * e01, 2 * mu * e11 + lam * tr];
  return { P, psi: mu * (e00 * e00 + 2 * e01 * e01 + e11 * e11) + 0.5 * lam * tr * tr };
}

// St Venant-Kirchhoff: Green strain E = (F^T F - I)/2, P = F (2 mu E + lambda tr(E) I).
function stvkP(F, mu, lam) {
  const C = mul(transpose(F), F);
  const E00 = 0.5 * (C[0] - 1), E01 = 0.5 * C[1], E11 = 0.5 * (C[3] - 1), tr = E00 + E11;
  const S = [2 * mu * E00 + lam * tr, 2 * mu * E01, 2 * mu * E01, 2 * mu * E11 + lam * tr];
  return { P: mul(F, S), psi: mu * (E00 * E00 + 2 * E01 * E01 + E11 * E11) + 0.5 * lam * tr * tr };
}

// Fixed corotated (MPM course eqs. 49 and 52):
// psi = mu sum (sigma_i - 1)^2 + lambda/2 (J - 1)^2,  P = 2 mu (F - R) + lambda (J - 1) J F^-T.
function fixedCorotatedP(F, mu, lam) {
  // TODO (exercise 1): compute R = polarRotation(F), J = det(F) and cof = cofactor(F) (that is J F^-T),
  // then P = 2 mu (F - R) + lambda (J - 1) cof, entry by entry.
  // The energy density is psi = mu |F - R|^2 + lambda/2 (J - 1)^2, with |.| the sum of squared entries.
  const P = [0, 0, 0, 0];
  const psi = 0;
  return { P, psi };
}

// Neo-Hookean (MPM course eqs. 46 and 48):
// psi = mu/2 (tr(F^T F) - 2) - mu log J + lambda/2 (log J)^2,  P = mu (F - F^-T) + lambda log(J) F^-T.
// Since dJ/dF = J F^-T, this is P = mu F + g'(J) J F^-T with g(J) = -mu log J + lambda/2 (log J)^2.
// Undefined for J <= 0. Below J = 0.1 the demo continues g as a straight line, so an inverted
// triangle still feels a strong, finite push back toward positive J.
const neoHookeanMinJ = 0.1;
function neoHookeanP(F, mu, lam) {
  const J = det(F), Jc = Math.max(J, neoHookeanMinJ), logJ = Math.log(Jc), cof = cofactor(F);
  const gPrime = (lam * logJ - mu) / Jc;
  const g = -mu * logJ + 0.5 * lam * logJ * logJ + gPrime * (J - Jc);
  const P = [];
  for (let k = 0; k < 4; k++) P.push(mu * F[k] + gPrime * cof[k]);
  const trC = F[0] * F[0] + F[1] * F[1] + F[2] * F[2] + F[3] * F[3];
  return { P, psi: 0.5 * mu * (trC - 2) + g };
}

// How much stiffer than at rest the material has become under compression. The fastest wave speed
// grows with its square root, so the explicit time step must shrink with it. Neo-Hookean's
// volume term stiffens like (1 - log J) / J^2; the other models do not stiffen when compressed.
function stiffening(minJ) {
  if (model.value !== 'neo-Hookean' || minJ >= 1) return 1;
  const J = Math.max(minJ, neoHookeanMinJ);
  return (1 - Math.log(J)) / (J * J);
}

const models = {
  'linear': linearP, 'St Venant-Kirchhoff': stvkP, 'fixed corotated': fixedCorotatedP, 'neo-Hookean': neoHookeanP,
};

// ---------------------------------------------------------------------------------------------
// The mesh: an nx x ny grid of squares, two triangles each
// ---------------------------------------------------------------------------------------------
const nx = 16, ny = 8, width = 0.8, height = 0.4, cell = width / nx;
const count = (nx + 1) * (ny + 1);
const restX = new Float64Array(count), restY = new Float64Array(count);
const x = new Float64Array(count), y = new Float64Array(count);
const vx = new Float64Array(count), vy = new Float64Array(count);
const fx = new Float64Array(count), fy = new Float64Array(count);
const mass = new Float64Array(count);
const tris = [];            // vertex indices
let DmInv, restArea;        // per triangle: Dm^-1 (4 numbers) and rest area V0
for (let j = 0; j <= ny; j++) {
  for (let i = 0; i <= nx; i++) { restX[i + j * (nx + 1)] = i * cell; restY[i + j * (nx + 1)] = j * cell; }
}
for (let j = 0; j < ny; j++) {
  for (let i = 0; i < nx; i++) {
    const a = i + j * (nx + 1), b = a + 1, c = a + nx + 1, d = c + 1;
    if ((i + j) % 2 === 0) tris.push([a, b, d], [a, d, c]); else tris.push([a, b, c], [b, d, c]);
  }
}
DmInv = new Float64Array(4 * tris.length);
restArea = new Float64Array(tris.length);
for (let t = 0; t < tris.length; t++) {
  const [i0, i1, i2] = tris[t];
  const Dm = [restX[i1] - restX[i0], restX[i2] - restX[i0], restY[i1] - restY[i0], restY[i2] - restY[i0]];
  const inv = inverse(Dm);
  for (let k = 0; k < 4; k++) DmInv[4 * t + k] = inv[k];
  restArea[t] = 0.5 * Math.abs(det(Dm));
  for (const v of tris[t]) mass[v] += density * restArea[t] / 3;   // lumped mass
}
let totalRestArea = 0;
for (let t = 0; t < tris.length; t++) totalRestArea += restArea[t];

function reset() {
  for (let k = 0; k < count; k++) {
    x[k] = restX[k] - 0.5 * width; y[k] = restY[k] + 0.5;
    vx[k] = 0; vy[k] = 0;
  }
  grabbed = [];
}

function spin(omega) {
  let cx = 0, cy = 0, m = 0;
  for (let k = 0; k < count; k++) { cx += mass[k] * x[k]; cy += mass[k] * y[k]; m += mass[k]; }
  cx /= m; cy /= m;
  for (let k = 0; k < count; k++) { vx[k] = -omega * (y[k] - cy); vy[k] = omega * (x[k] - cx); }
}

function triangleF(t) {
  const [i0, i1, i2] = tris[t];
  const Ds = [x[i1] - x[i0], x[i2] - x[i0], y[i1] - y[i0], y[i2] - y[i0]];
  const o = 4 * t;
  return mul(Ds, [DmInv[o], DmInv[o + 1], DmInv[o + 2], DmInv[o + 3]]);
}

// Forces of triangle t with stress P: H = -V0 P Dm^-T. Columns of H are the forces on
// vertices 1 and 2; vertex 0 gets minus their sum, so the forces add up to zero.
function addTriangleForces(t, P) {
  const [i0, i1, i2] = tris[t];
  const o = 4 * t, V0 = restArea[t];
  const m00 = DmInv[o], m01 = DmInv[o + 1], m10 = DmInv[o + 2], m11 = DmInv[o + 3];
  // TODO (exercise 2): form H = -V0 P Dm^-T (Dm^-1 is [m00, m01, m10, m11]; mind the transpose),
  // add its first column to (fx, fy) of vertex i1, its second column to vertex i2,
  // and minus the sum of both columns to vertex i0.
}

// ---------------------------------------------------------------------------------------------
// Time stepping
// ---------------------------------------------------------------------------------------------
let grabbed = [], substeps = 1, elasticEnergy = 0, currentMinJ = 1;

function substep(h, lp) {
  const constitutive = models[model.value];
  fx.fill(0); fy.fill(0);
  elasticEnergy = 0;
  currentMinJ = Infinity;
  for (let t = 0; t < tris.length; t++) {
    const F = triangleF(t);
    currentMinJ = Math.min(currentMinJ, det(F));
    const r = constitutive(F, lp.mu, lp.lambda);
    elasticEnergy += restArea[t] * r.psi;
    addTriangleForces(t, r.P);
  }
  const g = gravityOn.value ? gravity : 0, decay = Math.exp(-damping.value * h);
  for (let k = 0; k < count; k++) {
    vx[k] = (vx[k] + h * fx[k] / mass[k]) * decay;
    vy[k] = (vy[k] + h * (fy[k] / mass[k] - g)) * decay;
  }
  // Mouse: a stiff, critically damped spring pulling the grabbed vertices to the cursor.
  const w = 40;
  for (const gb of grabbed) {
    const k = gb.k;
    vx[k] += h * (w * w * (Sim.mouse.x + gb.dx - x[k]) - 2 * w * vx[k]);
    vy[k] += h * (w * w * (Sim.mouse.y + gb.dy - y[k]) - 2 * w * vy[k]);
  }
  for (let k = 0; k < count; k++) {
    x[k] += h * vx[k]; y[k] += h * vy[k];
    if (y[k] < 0) { y[k] = 0; if (vy[k] < 0) vy[k] = 0; vx[k] *= Math.exp(-30 * h); }   // floor with friction
    if (x[k] < -1.3) { x[k] = -1.3; if (vx[k] < 0) vx[k] = 0; }
    if (x[k] > 1.3) { x[k] = 1.3; if (vx[k] > 0) vx[k] = 0; }
  }
}

function step(dt) {
  if (Sim.mouse.pressed) {
    grabbed = [];
    for (let k = 0; k < count; k++) {
      const dx = x[k] - Sim.mouse.x, dy = y[k] - Sim.mouse.y;
      if (dx * dx + dy * dy < 0.08 * 0.08) grabbed.push({ k, dx, dy });
    }
  }
  if (!Sim.mouse.down) grabbed = [];
  const lp = lame(youngs.value, poisson.value);
  // Explicit stability: the fastest wave (speed c) must not cross a triangle in one substep.
  // The substep shrinks when the material stiffens; at most 1000 substeps per frame (slow motion beyond).
  const c = Math.sqrt((lp.lambda + 2 * lp.mu) / density);
  let t = 0;
  substeps = 0;
  while (t < dt && substeps < 1000) {
    const h = Math.min(dt - t, 0.2 * cell / (c * Math.sqrt(stiffening(currentMinJ))));
    substep(h, lp);
    t += h;
    substeps++;
  }
}

// ---------------------------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------------------------
const ramp = [];
for (let i = 0; i <= 64; i++) {
  const t = i / 64 * 2 - 1, u = Math.abs(t);
  const pale = [238, 236, 228], end = t < 0 ? [47, 111, 223] : [217, 71, 58];
  const c = pale.map((p, m) => Math.round(p + (end[m] - p) * u));
  ramp.push(`rgb(${c[0]},${c[1]},${c[2]})`);
}
function colourJ(J) {
  if (!(J > 0)) return '#5b2a86';
  const t = Math.max(-1, Math.min(1, Math.log2(J) / 0.5));   // full colour at J = 0.71 or 1.41
  return ramp[Math.round((t + 1) * 32)];
}

function draw(g) {
  g.ground(0);
  let area = 0, minJ = Infinity, inverted = 0;
  for (let t = 0; t < tris.length; t++) {
    const F = triangleF(t), J = det(F);
    area += J * restArea[t];
    minJ = Math.min(minJ, J);
    if (J <= 0) inverted++;
    const [i0, i1, i2] = tris[t];
    g.polyline([x[i0], y[i0], x[i1], y[i1], x[i2], y[i2]],
      { closed: true, fill: colourJ(J), color: 'rgba(60,60,60,0.35)', width: 0.6 });
  }
  for (const gb of grabbed) g.dot(x[gb.k], y[gb.k], 3, Sim.colors.orange);
  let kinetic = 0, potential = 0;
  for (let k = 0; k < count; k++) {
    kinetic += 0.5 * mass[k] * (vx[k] * vx[k] + vy[k] * vy[k]);
    potential += mass[k] * gravity * y[k];
  }
  const lp = lame(youngs.value, poisson.value);
  Sim.show('mu, lambda (Pa)', `${Sim.format(lp.mu)}, ${Sim.format(lp.lambda)}`);
  Sim.show('substeps per frame', substeps);
  Sim.show('area / rest area', area / totalRestArea, 4);
  Sim.show('smallest J', minJ, 3);
  Sim.show('inverted triangles', inverted);
  Sim.show('kinetic energy (J)', kinetic, 2);
  Sim.show('elastic energy (J)', elasticEnergy, 2);
  Sim.plot('area ratio', Math.max(-0.5, Math.min(area / totalRestArea, 3)), { color: Sim.colors.blue, min: 0, max: 1.2 });
  Sim.plot('smallest J', Math.max(-0.5, Math.min(minJ, 3)), { color: Sim.colors.purple });
}

Sim.run({ title: 'A jelly block, four ways (starter)', view: [-1.3, -0.6, 1.3, 1.5], reset, step, draw, grid: 0.1 });
