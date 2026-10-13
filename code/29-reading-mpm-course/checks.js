// Chapter 29 - Reading the MPM course: numerical checks of the notes' equations on a random case.
//   (48)  neo-Hookean P against a finite-difference derivative of the energy (46)
//   (52)  fixed corotated P against a finite-difference derivative of (49)
//   (59)  the 2D closed form for the change of the polar rotation, d(theta), against finite differences
//   (189) the grid force -sum_p V0 P F^T grad w against a finite-difference gradient of the energy (164)
//   (196) the force differential against a finite difference of the force
//   and that the grid forces have zero total force and zero total torque.
// Matrices are 2x2, stored row-major as [a, b, c, d] = [[a, b], [c, d]].

const model = Sim.select('energy for the force test', ['fixed corotated', 'neo-Hookean'], 'fixed corotated', { reset: true });
let seed = 1;
Sim.button('new random case', () => { seed++; Sim.reset(); });

const E = 1e4, nu = 0.3;                                   // a soft jelly (Pa)
const mu = E / (2 * (1 + nu)), lam = E * nu / ((1 + nu) * (1 - 2 * nu));
const dx = 0.1, n = 11;                                    // 11 x 11 grid nodes on [0, 1]^2
const V0 = dx * dx / 4;                                    // rest volume per particle (m^2 in 2D)

// --- 2x2 matrix helpers -------------------------------------------------------------------------
const det = (F) => F[0] * F[3] - F[1] * F[2];
const invT = (F) => { const J = det(F); return [F[3] / J, -F[2] / J, -F[1] / J, F[0] / J]; };
const mul = (A, B) => [A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3], A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3]];
const T = (A) => [A[0], A[2], A[1], A[3]];
const add = (A, B, s = 1) => A.map((v, k) => v + s * B[k]);
const maxAbs = (A) => Math.max(...A.map(Math.abs));
function polarAngle(F) { return Math.atan2(F[2] - F[1], F[0] + F[3]); }
function rot(th) { const c = Math.cos(th), s = Math.sin(th); return [c, -s, s, c]; }

// --- the two energies and their stresses ---------------------------------------------------------
function psiNH(F) { const J = det(F); return mu / 2 * (F[0] ** 2 + F[1] ** 2 + F[2] ** 2 + F[3] ** 2 - 2) - mu * Math.log(J) + lam / 2 * Math.log(J) ** 2; }
function PNH(F) { const J = det(F), Fi = invT(F); return add(add(F, Fi, -1).map((v) => mu * v), Fi.map((v) => lam * Math.log(J) * v)); }
function psiFC(F) { const R = rot(polarAngle(F)), J = det(F), M = add(F, R, -1); return mu * (M[0] ** 2 + M[1] ** 2 + M[2] ** 2 + M[3] ** 2) + lam / 2 * (J - 1) ** 2; }
function PFC(F) { const R = rot(polarAngle(F)), J = det(F), Fi = invT(F); return add(add(F, R, -1).map((v) => 2 * mu * v), Fi.map((v) => lam * (J - 1) * J * v)); }
const psi = () => (model.value === 'neo-Hookean' ? psiNH : psiFC);
const P = () => (model.value === 'neo-Hookean' ? PNH : PFC);

function fdGradient(f, F, eps = 1e-6) {
  return [0, 1, 2, 3].map((k) => { const a = F.slice(), b = F.slice(); a[k] += eps; b[k] -= eps; return (f(a) - f(b)) / (2 * eps); });
}

// --- quadratic B-spline weights on the grid ------------------------------------------------------
function stencil(x, y) {
  const out = [];
  const wts = (u) => { const b = Math.floor(u - 0.5), f = u - b; return { b, w: [0.5 * (1.5 - f) ** 2, 0.75 - (f - 1) ** 2, 0.5 * (f - 0.5) ** 2], d: [(f - 1.5) / dx, -2 * (f - 1) / dx, (f - 0.5) / dx] }; };
  const sx = wts(x / dx), sy = wts(y / dx);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    out.push({ node: (sy.b + j) * n + (sx.b + i), w: sx.w[i] * sy.w[j], gx: sx.d[i] * sy.w[j], gy: sx.w[i] * sy.d[j] });
  }
  return out;
}

let parts, results, forces;

function gauss(rand) { return Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand()); }

// F_p(x_hat) = (I + sum_j (x_hat_j - x_j) grad w_jp^T) F_p^n, with u_j = x_hat_j - x_j
function deformed(p, u) {
  let G = [0, 0, 0, 0];
  for (const s of p.st) { const ux = u[2 * s.node], uy = u[2 * s.node + 1]; G = add(G, [ux * s.gx, ux * s.gy, uy * s.gx, uy * s.gy]); }
  return mul(add([1, 0, 0, 1], G), p.F);
}
function energy(u) { let e = 0; for (const p of parts) e += V0 * psi()(deformed(p, u)); return e; }
// Equation (189) generalised to moved nodes: f_i = -sum_p V0 P(F_p(x_hat)) (F_p^n)^T grad w_ip
function gridForces(u) {
  const f = new Float64Array(2 * n * n);
  for (const p of parts) {
    const tau = mul(P()(deformed(p, u)), T(p.F));
    for (const s of p.st) { f[2 * s.node] -= V0 * (tau[0] * s.gx + tau[1] * s.gy); f[2 * s.node + 1] -= V0 * (tau[2] * s.gx + tau[3] * s.gy); }
  }
  return f;
}

function reset() {
  const rand = Sim.random(seed);
  parts = [];
  for (let k = 0; k < 8; k++) {
    let F;
    do { F = [1 + 0.25 * gauss(rand), 0.25 * gauss(rand), 0.25 * gauss(rand), 1 + 0.25 * gauss(rand)]; } while (det(F) < 0.3);
    const x = 0.35 + 0.3 * rand(), y = 0.35 + 0.3 * rand();
    parts.push({ x, y, F, st: stencil(x, y) });
  }
  results = [];
  const F = parts[0].F;
  const check = (label, a, b) => {
    const err = maxAbs(add(a, b, -1)) / Math.max(maxAbs(b), 1e-30);
    results.push([label, err]);
  };
  check('(48) neo-Hookean P vs FD', PNH(F), fdGradient(psiNH, F));
  check('(52) fixed corotated P vs FD', PFC(F), fdGradient(psiFC, F));

  // (57)-(59) in 2D: d(theta) = ((R^T dF)_21 - (R^T dF)_12) / tr(S), with S = R^T F
  const dF = [gauss(rand), gauss(rand), gauss(rand), gauss(rand)], eps = 1e-6;
  const R = rot(polarAngle(F)), S = mul(T(R), F), RdF = mul(T(R), dF);
  const dthClosed = (RdF[2] - RdF[1]) / (S[0] + S[3]);
  const dthFD = (polarAngle(add(F, dF, eps)) - polarAngle(add(F, dF, -eps))) / (2 * eps);
  check('(59) 2D d(theta) of polar R vs FD', [dthClosed], [dthFD]);

  // (189): force formula vs -d(energy)/d(x_hat) at x_hat = x
  const zero = new Float64Array(2 * n * n);
  forces = gridForces(zero);
  const fd = new Float64Array(2 * n * n);
  for (let k = 0; k < 2 * n * n; k++) {
    if (forces[k] === 0) continue;
    const a = zero.slice(), b = zero.slice(); a[k] += eps; b[k] -= eps;
    fd[k] = -(energy(a) - energy(b)) / (2 * eps);
  }
  check('(189) grid force vs -grad e (FD)', Array.from(forces), Array.from(fd));

  // (196): force differential with dF_p = sum_j u_j grad w_jp^T F_p^n, vs FD of the force along u
  const u = new Float64Array(2 * n * n).map(() => gauss(rand));
  const df = new Float64Array(2 * n * n);
  for (const p of parts) {
    const dFp = add(deformed(p, u), p.F, -1);            // deformed() is linear in u, so this is exactly dF_p
    const dP = add(P()(add(p.F, dFp, eps)), P()(add(p.F, dFp, -eps)), -1).map((v) => v / (2 * eps));
    const A = mul(dP, T(p.F));
    for (const s of p.st) { df[2 * s.node] -= V0 * (A[0] * s.gx + A[1] * s.gy); df[2 * s.node + 1] -= V0 * (A[2] * s.gx + A[3] * s.gy); }
  }
  const fp = gridForces(u.map((v) => v * eps)), fm = gridForces(u.map((v) => -v * eps));
  check('(196) force differential vs FD', Array.from(df), Array.from(fp).map((v, k) => (v - fm[k]) / (2 * eps)));

  // Momentum: sum of forces and sum of torques should vanish (relative to the force scale).
  let fx = 0, fy = 0, tq = 0, scale = 0;
  for (let i = 0; i < n * n; i++) {
    const x = (i % n) * dx, y = Math.floor(i / n) * dx, a = forces[2 * i], b = forces[2 * i + 1];
    fx += a; fy += b; tq += x * b - y * a; scale += Math.hypot(a, b);
  }
  results.push(['sum of forces / sum |f|', Math.hypot(fx, fy) / scale]);
  results.push(['sum of torques / sum |f| (m)', Math.abs(tq) / scale]);
}

function draw(g) {
  for (let i = 0; i < n * n; i++) g.dot((i % n) * dx, Math.floor(i / n) * dx, 2, Sim.colors.gray);
  // Each particle drawn as its deformed little square F * [-a, a]^2.
  const a = 0.03;
  for (const p of parts) {
    const F = p.F, c = [[-a, -a], [a, -a], [a, a], [-a, a]];
    g.polyline(c.map(([X, Y]) => ({ x: p.x + F[0] * X + F[1] * Y, y: p.y + F[2] * X + F[3] * Y })), { closed: true, fill: Sim.colors.blue, alpha: 0.35, color: Sim.colors.blue, width: 1 });
  }
  let fmax = 0;
  for (let k = 0; k < forces.length; k++) fmax = Math.max(fmax, Math.abs(forces[k]));
  const s = 0.12 / Math.max(fmax, 1e-12);
  for (let i = 0; i < n * n; i++) {
    const fx = forces[2 * i], fy = forces[2 * i + 1];
    if (fx === 0 && fy === 0) continue;
    g.arrow((i % n) * dx, Math.floor(i / n) * dx, fx * s, fy * s, { color: Sim.colors.red, width: 1.5 });
  }
  for (const [label, err] of results) Sim.show(label, err.toExponential(1) + (err < 1e-5 ? '  pass' : '  FAIL'));
  Sim.show('largest grid force', fmax.toFixed(2) + ' N/m (2D)');
}

Sim.run({ title: 'Checking the MPM course equations', view: [0.1, 0.1, 0.9, 0.9], reset, step() {}, draw, grid: false, showTime: false });
