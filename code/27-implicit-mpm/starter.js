// Chapter 27 - Implicit MPM (starter)
// Backward Euler on the grid, solved as an optimisation (Gast et al. 2015; MPM course section 11.5):
//   minimise E(v) = sum_i 1/2 m_i |v_i - v*_i|^2 + Phi(x_hat(v)),   x_hat_i = x_i + h v_i,  v* = v^n + h g.
// The energy, its gradient and Hessian-vector products are done for you. Your job, in the TODO blocks:
//   1. the preconditioned conjugate gradient solver cg(),
//   2. Newton's method with a backtracking line search in newtonSolve().
// Until cg() works, nothing feels any elastic force and the beam falls off the wall.

const scene = Sim.select('scene', ['cantilever', 'jelly drop'], 'cantilever', { reset: true });
const method = Sim.select('method', ['implicit (Newton + line search)', 'semi-implicit (1 Newton step)', 'explicit'],
  'implicit (Newton + line search)');
const youngs = Sim.logSlider("Young's modulus E (Pa)", 4, 7, 1e6, { reset: true });
const stepMs = Sim.logSlider('time step h (ms)', -1, 1.53, 1000 / 60);
const maxNewton = Sim.slider('max Newton iterations', 1, 20, 8, { step: 1 });
const maxCG = Sim.slider('max CG iterations', 5, 500, 100, { step: 5 });
Sim.button('kick', kick);

// ---------------------------------------------------------------------------------------------
// Grid. Velocities are stored as flat vectors of length 2 (NN * NN): [vx0, vy0, vx1, vy1, ...].
// ---------------------------------------------------------------------------------------------
const N = 64, dx = 1 / N, invDx = N, NN = N + 1, NDOF = 2 * NN * NN, BOUND = 3;
const gravity = -9.81, rho = 1000, nu = 0.3;
const gm = new Float64Array(NN * NN);       // node mass
const vStar = new Float64Array(NDOF);       // v* = v^n + h g, the velocity without elastic forces
const vel = new Float64Array(NDOF);         // the unknown v^{n+1}
const free = new Float64Array(NDOF);        // 1 for a degree of freedom we solve for, 0 if fixed (inactive or wall)
const freeList = new Int32Array(NDOF);      // the indices of the free degrees of freedom
let nFree = 0;
const diagH = new Float64Array(NDOF);       // Jacobi preconditioner: approximate diagonal of the Hessian
const grad = new Float64Array(NDOF), dir = new Float64Array(NDOF), vTry = new Float64Array(NDOF);
const cgR = new Float64Array(NDOF), cgZ = new Float64Array(NDOF), cgP = new Float64Array(NDOF), cgAp = new Float64Array(NDOF);

// ---------------------------------------------------------------------------------------------
// Particles
// ---------------------------------------------------------------------------------------------
const MAXP = 3000;
const px = new Float64Array(MAXP), py = new Float64Array(MAXP);
const pvx = new Float64Array(MAXP), pvy = new Float64Array(MAXP);
const pC = new Float64Array(4 * MAXP), pF = new Float64Array(4 * MAXP);
const pm = new Float64Array(MAXP), pV0 = new Float64Array(MAXP);
const pcol = new Array(MAXP);
let count = 0, tipParticle = -1, tipY0 = 0;

// Per-particle stencil, rebuilt every step: the 9 nodes, their weights and q_ip = F_p^T grad w_ip.
// With q, the deformation gradient after the step is F_hat_p(v) = F_p + h sum_i v_i q_ip^T.
const sNode = new Int32Array(9 * MAXP), sW = new Float64Array(9 * MAXP);
const sQx = new Float64Array(9 * MAXP), sQy = new Float64Array(9 * MAXP);
// Per-particle state at the current Newton iterate, for Hessian-vector products.
const hF = new Float64Array(4 * MAXP);          // F_hat
const hRc = new Float64Array(MAXP), hRs = new Float64Array(MAXP), hR2 = new Float64Array(MAXP);

function addBlock(cx, cy, w, hgt, angle, vx, vy, color) {
  const s = dx / 2, c = Math.cos(angle), sn = Math.sin(angle);
  for (let a = -w / 2 + s / 2; a < w / 2; a += s) {
    for (let b = -hgt / 2 + s / 2; b < hgt / 2; b += s) {
      if (count >= MAXP) return;
      const p = count++;
      px[p] = cx + c * a - sn * b; py[p] = cy + sn * a + c * b;
      pvx[p] = vx; pvy[p] = vy;
      pC.fill(0, 4 * p, 4 * p + 4);
      pF[4 * p] = 1; pF[4 * p + 1] = 0; pF[4 * p + 2] = 0; pF[4 * p + 3] = 1;
      pV0[p] = s * s; pm[p] = rho * s * s; pcol[p] = color;
    }
  }
}

// Cantilever: clamped in the left wall, free length L = 0.5 m from the wall face at x = 3 dx, thickness 0.1 m.
const BEAM_L = 0.5, BEAM_T = 0.1, BEAM_Y = 0.6;
function reset() {
  count = 0; tipParticle = -1;
  if (scene.value === 'cantilever') {
    const x0 = 0.01, x1 = BOUND * dx + BEAM_L;
    addBlock((x0 + x1) / 2, BEAM_Y, x1 - x0, BEAM_T, 0, 0, 0, Sim.colors.blue);
    // the particle nearest the free end on the centre line, to measure the deflection
    let best = 1e9;
    for (let p = 0; p < count; p++) {
      const d = Math.abs(px[p] - (x1 - dx / 4)) + Math.abs(py[p] - BEAM_Y);
      if (d < best) { best = d; tipParticle = p; }
    }
    tipY0 = py[tipParticle];
  } else {
    addBlock(0.3, 0.2, 0.16, 0.16, 0, 0, 0, Sim.colors.red);
    addBlock(0.6, 0.6, 0.16, 0.16, Math.PI / 6, -1, 0, Sim.colors.orange);
  }
  stats.newton = 0; stats.cg = 0; stats.halvings = 0; stats.residual = 0;
}

function kick() {
  for (let p = 0; p < count; p++) if (px[p] > 0.35) pvy[p] += 1.5;
}

function lameMu() { return youngs.value / (2 * (1 + nu)); }
function lameLambda() { return youngs.value * nu / ((1 + nu) * (1 - 2 * nu)); }

// ---------------------------------------------------------------------------------------------
// Fixed corotated energy, stress and stress differential, without trigonometry.
// R is the rotation [[c, -s], [s, c]] with (c, s) proportional to (a, b) = (F00 + F11, F10 - F01).
// ---------------------------------------------------------------------------------------------
const P = new Float64Array(4);

function energyDensity(F00, F01, F10, F11, mu, la) {
  const a = F00 + F11, b = F10 - F01, r = Math.hypot(a, b) || 1;
  const c = a / r, s = b / r, J = F00 * F11 - F01 * F10;
  const e00 = F00 - c, e01 = F01 + s, e10 = F10 - s, e11 = F11 - c;          // F - R
  return mu * (e00 * e00 + e01 * e01 + e10 * e10 + e11 * e11) + 0.5 * la * (J - 1) * (J - 1);
}

// P = 2 mu (F - R) + lambda (J - 1) cof(F), cof(F) = J F^-T = [[F11, -F10], [-F01, F00]].
function stress(F00, F01, F10, F11, c, s, mu, la, P) {
  const J = F00 * F11 - F01 * F10;
  P[0] = 2 * mu * (F00 - c) + la * (J - 1) * F11;
  P[1] = 2 * mu * (F01 + s) - la * (J - 1) * F10;
  P[2] = 2 * mu * (F10 - s) - la * (J - 1) * F01;
  P[3] = 2 * mu * (F11 - c) + la * (J - 1) * F00;
}

// dP = 2 mu (dF - dR) + lambda (cof(F) : dF) cof(F) + lambda (J - 1) cof(dF)    (course eqs. 53-56)
// with dR = dtheta [[-s, -c], [c, -s]] and dtheta = (a db - b da) / (a^2 + b^2).
function stressDifferential(F00, F01, F10, F11, c, s, r2, d00, d01, d10, d11, mu, la, P) {
  const J = F00 * F11 - F01 * F10;
  const a = F00 + F11, b = F10 - F01;
  const dth = (a * (d10 - d01) - b * (d00 + d11)) / r2;
  const dJ = F11 * d00 - F10 * d01 - F01 * d10 + F00 * d11;               // cof(F) : dF
  P[0] = 2 * mu * (d00 + dth * s) + la * (dJ * F11 + (J - 1) * d11);
  P[1] = 2 * mu * (d01 + dth * c) + la * (-dJ * F10 - (J - 1) * d10);
  P[2] = 2 * mu * (d10 - dth * c) + la * (-dJ * F01 - (J - 1) * d01);
  P[3] = 2 * mu * (d11 + dth * s) + la * (dJ * F00 + (J - 1) * d00);
}

// ---------------------------------------------------------------------------------------------
// The incremental potential E(v), its gradient and Hessian-vector products. All matrix-free.
// ---------------------------------------------------------------------------------------------
let hCur = 1 / 60;

// F_hat_p(v) = F_p + h sum_i v_i q_ip^T, written into out at offset 4p.
function deformedF(v, p, out) {
  let f00 = pF[4 * p], f01 = pF[4 * p + 1], f10 = pF[4 * p + 2], f11 = pF[4 * p + 3];
  for (let k = 9 * p; k < 9 * p + 9; k++) {
    const i = 2 * sNode[k], qx = hCur * sQx[k], qy = hCur * sQy[k];
    f00 += v[i] * qx; f01 += v[i] * qy; f10 += v[i + 1] * qx; f11 += v[i + 1] * qy;
  }
  out[4 * p] = f00; out[4 * p + 1] = f01; out[4 * p + 2] = f10; out[4 * p + 3] = f11;
}

const tmpF = new Float64Array(4 * MAXP);

// E(v). Fixed and inactive degrees of freedom are part of the sum too; they never change.
function energy(v) {
  const mu = lameMu(), la = lameLambda();
  let e = 0;
  for (let n = 0; n < NN * NN; n++) {
    if (gm[n] <= 0) continue;
    const ux = v[2 * n] - vStar[2 * n], uy = v[2 * n + 1] - vStar[2 * n + 1];
    e += 0.5 * gm[n] * (ux * ux + uy * uy);
  }
  for (let p = 0; p < count; p++) {
    deformedF(v, p, tmpF);
    e += pV0[p] * energyDensity(tmpF[4 * p], tmpF[4 * p + 1], tmpF[4 * p + 2], tmpF[4 * p + 3], mu, la);
  }
  return e;
}

// E(v) and its gradient g_i = m_i (v_i - v*_i) + h sum_p V_p^0 P(F_hat_p) q_ip  (= m_i (v_i - v*_i) - h f_i).
// Also caches F_hat and its rotation for the Hessian-vector products at this v.
function energyAndGradient(v, g) {
  const mu = lameMu(), la = lameLambda();
  let e = 0;
  g.fill(0);
  for (let n = 0; n < NN * NN; n++) {
    if (gm[n] <= 0) continue;
    const ux = v[2 * n] - vStar[2 * n], uy = v[2 * n + 1] - vStar[2 * n + 1];
    e += 0.5 * gm[n] * (ux * ux + uy * uy);
    g[2 * n] = gm[n] * ux; g[2 * n + 1] = gm[n] * uy;
  }
  for (let p = 0; p < count; p++) {
    deformedF(v, p, hF);
    const F00 = hF[4 * p], F01 = hF[4 * p + 1], F10 = hF[4 * p + 2], F11 = hF[4 * p + 3];
    const a = F00 + F11, b = F10 - F01, r2 = Math.max(a * a + b * b, 1e-12), r = Math.sqrt(r2);
    hRc[p] = a / r; hRs[p] = b / r; hR2[p] = r2;
    e += pV0[p] * energyDensity(F00, F01, F10, F11, mu, la);
    stress(F00, F01, F10, F11, hRc[p], hRs[p], mu, la, P);
    const V = pV0[p] * hCur;
    for (let k = 9 * p; k < 9 * p + 9; k++) {
      const i = 2 * sNode[k], qx = sQx[k], qy = sQy[k];
      g[i] += V * (P[0] * qx + P[1] * qy);
      g[i + 1] += V * (P[2] * qx + P[3] * qy);
    }
  }
  for (let i = 0; i < NDOF; i++) g[i] *= free[i];
  return e;
}

// out = H u = M u + h^2 sum_p V_p^0 dP(F_hat_p)[dF_p] q_ip,  with dF_p = sum_j u_j q_jp^T  (course eqs. 196-197).
// u must be zero on fixed degrees of freedom; out is only meaningful on the free ones (the solver never reads
// the others, so we do not spend time clearing them).
function hessVec(u, out) {
  const mu = lameMu(), la = lameLambda(), h2 = hCur * hCur;
  out.fill(0);
  for (let f = 0; f < nFree; f++) { const i = freeList[f]; out[i] = gm[i >> 1] * u[i]; }
  for (let p = 0; p < count; p++) {
    let d00 = 0, d01 = 0, d10 = 0, d11 = 0;
    for (let k = 9 * p; k < 9 * p + 9; k++) {
      const i = 2 * sNode[k];
      d00 += u[i] * sQx[k]; d01 += u[i] * sQy[k]; d10 += u[i + 1] * sQx[k]; d11 += u[i + 1] * sQy[k];
    }
    stressDifferential(hF[4 * p], hF[4 * p + 1], hF[4 * p + 2], hF[4 * p + 3], hRc[p], hRs[p], hR2[p],
      d00, d01, d10, d11, mu, la, P);
    const V = pV0[p] * h2;
    for (let k = 9 * p; k < 9 * p + 9; k++) {
      const i = 2 * sNode[k], qx = sQx[k], qy = sQy[k];
      out[i] += V * (P[0] * qx + P[1] * qy);
      out[i + 1] += V * (P[2] * qx + P[3] * qy);
    }
  }
}

// Dot product over the free degrees of freedom only (everything else is fixed and stays zero).
function dot(a, b) {
  let s = 0;
  for (let f = 0; f < nFree; f++) { const i = freeList[f]; s += a[i] * b[i]; }
  return s;
}

// ---------------------------------------------------------------------------------------------
// Preconditioned conjugate gradient for H x = b, matrix-free: only hessVec is needed.
// Starts from x = 0. Stops when |r| < relTol |b|, after maxIter iterations, or if it meets a direction of
// negative curvature (p^T H p <= 0), which can happen because the elastic Hessian is not always positive
// definite. Returns the number of iterations.
// ---------------------------------------------------------------------------------------------
function cg(b, x, maxIter, relTol) {
  x.fill(0); cgP.fill(0);
  // TODO (exercise 1): preconditioned CG. Work only on the free entries: loop
  //   for (let f = 0; f < nFree; f++) { const i = freeList[f]; ... }
  // and use dot(a, b), hessVec(u, out) (out = H u) and the Jacobi preconditioner z_i = r_i / diagH[i].
  // Scratch vectors cgR, cgZ, cgP, cgAp are allocated for you; keep cgP zero outside the free entries.
  //   r = b, z = M^-1 r, p = z
  //   repeat: Ap = H p; if p.Ap <= 0, stop (on the first iteration, return x = z instead);
  //           alpha = (r.z) / (p.Ap); x += alpha p; r -= alpha Ap; stop if |r|^2 < relTol^2 |b|^2;
  //           z = M^-1 r; beta = (r.z)_new / (r.z)_old; p = z + beta p
  // Return the number of iterations. Until then x stays zero: no elastic forces at all.
  return 0;
}

// ---------------------------------------------------------------------------------------------
// Newton's method on E(v) with a backtracking (Armijo) line search.
// ---------------------------------------------------------------------------------------------
const stats = { newton: 0, cg: 0, halvings: 0, residual: 0 };
const NEWTON_TOL = 1e-3;            // mass-weighted RMS of g_i / m_i (m/s): the velocity error of the step

function residualNorm(g) {
  let s = 0, m = 0;
  for (let n = 0; n < NN * NN; n++) {
    if (gm[n] <= 0) continue;
    s += (g[2 * n] * g[2 * n] + g[2 * n + 1] * g[2 * n + 1]) / gm[n];
    m += gm[n];
  }
  return Math.sqrt(s / m);
}

function newtonSolve() {
  for (let i = 0; i < NDOF; i++) vel[i] = vStar[i] * free[i];      // initial guess: no elastic forces
  // TODO (exercise 2): replace this single, unguarded Newton step by Newton's method with a line search.
  //   for up to maxNewton.value iterations:
  //     E0 = energyAndGradient(vel, grad); stop if residualNorm(grad) < NEWTON_TOL; count stats.newton
  //     solve H d = -g with cg (put -grad into vTry, the result in dir; add the iterations to stats.cg)
  //     if g.d is not negative, use d = -g / diagH instead
  //     alpha = 1; while energy(vel + alpha d) > E0 + 1e-4 alpha g.d: alpha /= 2, count stats.halvings
  //     vel = vel + alpha d
  energyAndGradient(vel, grad);
  for (let i = 0; i < NDOF; i++) vTry[i] = -grad[i];
  stats.cg += cg(vTry, dir, maxCG.value, 1e-3);
  stats.newton = 1;
  for (let i = 0; i < NDOF; i++) vel[i] += dir[i];
  energyAndGradient(vel, grad);
  stats.residual = residualNorm(grad);
}

// The linearised update of course section 11.4: one Newton step from v = 0, no line search.
function semiImplicitSolve() {
  vel.fill(0);
  energyAndGradient(vel, grad);
  for (let i = 0; i < NDOF; i++) vTry[i] = -grad[i];
  stats.cg += cg(vTry, dir, maxCG.value, 1e-3);
  stats.newton = 1;
  for (let i = 0; i < NDOF; i++) vel[i] = dir[i];
  energyAndGradient(vel, grad);
  stats.residual = residualNorm(grad);
}

// Symplectic Euler (chapter 25): v = v* + h f(x^n) / m, which is v = -g(0) / m.
function explicitSolve() {
  vel.fill(0);
  energyAndGradient(vel, grad);
  for (let n = 0; n < NN * NN; n++) {
    if (gm[n] <= 0) continue;
    vel[2 * n] = -grad[2 * n] / gm[n];
    vel[2 * n + 1] = -grad[2 * n + 1] / gm[n];
  }
  stats.residual = 0;
}

// ---------------------------------------------------------------------------------------------
// One MPM step: P2G (APIC), set up the grid problem, solve for v, G2P with the F update, advection.
// ---------------------------------------------------------------------------------------------
const wx = new Float64Array(3), wy = new Float64Array(3), dwx = new Float64Array(3), dwy = new Float64Array(3);
function quadratic(f, w, dw) {
  w[0] = 0.5 * (1.5 - f) * (1.5 - f); w[1] = 0.75 - (f - 1) * (f - 1); w[2] = 0.5 * (f - 0.5) * (f - 0.5);
  dw[0] = (f - 1.5) * invDx; dw[1] = -2 * (f - 1) * invDx; dw[2] = (f - 0.5) * invDx;
}

function step(h) {
  hCur = h;
  stats.newton = 0; stats.cg = 0; stats.halvings = 0;
  gm.fill(0); vStar.fill(0);

  // P2G: mass and APIC momentum, and the per-particle stencil.
  for (let p = 0; p < count; p++) {
    const bx = Math.floor(px[p] * invDx - 0.5), by = Math.floor(py[p] * invDx - 0.5);
    const fx = px[p] * invDx - bx, fy = py[p] * invDx - by;
    quadratic(fx, wx, dwx); quadratic(fy, wy, dwy);
    const F00 = pF[4 * p], F01 = pF[4 * p + 1], F10 = pF[4 * p + 2], F11 = pF[4 * p + 3];
    const C00 = pC[4 * p], C01 = pC[4 * p + 1], C10 = pC[4 * p + 2], C11 = pC[4 * p + 3];
    let k = 9 * p;
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++, k++) {
        const node = (bx + a) * NN + (by + b);
        const w = wx[a] * wy[b], gx = dwx[a] * wy[b], gy = wx[a] * dwy[b];
        const rx = (a - fx) * dx, ry = (b - fy) * dx;
        sNode[k] = node; sW[k] = w;
        sQx[k] = F00 * gx + F10 * gy;              // q = F^T grad w
        sQy[k] = F01 * gx + F11 * gy;
        const mw = w * pm[p];
        gm[node] += mw;
        vStar[2 * node] += mw * (pvx[p] + C00 * rx + C01 * ry);
        vStar[2 * node + 1] += mw * (pvy[p] + C10 * rx + C11 * ry);
      }
    }
  }

  // Grid velocities v^n, then v* = v^n + h g. Sticky walls and inactive nodes are fixed at zero.
  nFree = 0;
  for (let i = 0; i < NN; i++) {
    for (let j = 0; j < NN; j++) {
      const n = i * NN + j;
      const wall = i < BOUND || i > N - BOUND || j < BOUND || j > N - BOUND;
      const f = gm[n] > 0 && !wall ? 1 : 0;
      free[2 * n] = f; free[2 * n + 1] = f;
      if (f) { freeList[nFree++] = 2 * n; freeList[nFree++] = 2 * n + 1; }
      if (gm[n] > 0) { vStar[2 * n] /= gm[n]; vStar[2 * n + 1] = vStar[2 * n + 1] / gm[n] + h * gravity; }
      if (!f) { vStar[2 * n] = 0; vStar[2 * n + 1] = 0; }
    }
  }

  // Jacobi preconditioner: the diagonal of M + h^2 K for linear elasticity at F_p,
  // K_ii(alpha alpha) ~ sum_p V_p^0 (mu |q|^2 + (lambda + mu) q_alpha^2).
  const mu = lameMu(), la = lameLambda(), h2 = h * h;
  for (let n = 0; n < NN * NN; n++) { diagH[2 * n] = gm[n]; diagH[2 * n + 1] = gm[n]; }
  for (let k = 0; k < 9 * count; k++) {
    const V = pV0[Math.floor(k / 9)] * h2, qx = sQx[k], qy = sQy[k], q2 = qx * qx + qy * qy, i = 2 * sNode[k];
    diagH[i] += V * (mu * q2 + (la + mu) * qx * qx);
    diagH[i + 1] += V * (mu * q2 + (la + mu) * qy * qy);
  }
  for (let i = 0; i < NDOF; i++) if (!(diagH[i] > 0)) diagH[i] = 1;

  if (method.value === 'explicit') explicitSolve();
  else if (method.value === 'semi-implicit (1 Newton step)') semiImplicitSolve();
  else newtonSolve();

  // G2P: APIC velocity and C, F <- F_hat(v) = (I + h grad v) F, advection.
  const s = 4 * invDx * invDx;
  for (let p = 0; p < count; p++) {
    const bx = Math.floor(px[p] * invDx - 0.5), by = Math.floor(py[p] * invDx - 0.5);
    const fx = px[p] * invDx - bx, fy = py[p] * invDx - by;
    let vx = 0, vy = 0, B00 = 0, B01 = 0, B10 = 0, B11 = 0, k = 9 * p;
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++, k++) {
        const i = 2 * sNode[k], w = sW[k];
        const rx = (a - fx) * dx, ry = (b - fy) * dx;
        const ux = vel[i], uy = vel[i + 1];
        vx += w * ux; vy += w * uy;
        B00 += w * ux * rx; B01 += w * ux * ry; B10 += w * uy * rx; B11 += w * uy * ry;
      }
    }
    pvx[p] = vx; pvy[p] = vy;
    pC[4 * p] = s * B00; pC[4 * p + 1] = s * B01; pC[4 * p + 2] = s * B10; pC[4 * p + 3] = s * B11;
    deformedF(vel, p, pF);
    px[p] += h * vx; py[p] += h * vy;
    const lo = dx, hi = 1 - dx;
    if (!(px[p] > lo)) px[p] = lo; else if (px[p] > hi) px[p] = hi;
    if (!(py[p] > lo)) py[p] = lo; else if (py[p] > hi) py[p] = hi;
  }
}

// ---------------------------------------------------------------------------------------------
// Drawing and readouts
// ---------------------------------------------------------------------------------------------
const FRAME = 1 / 60;
function substeps() { return Math.max(1, Math.min(60, Math.round(FRAME / (stepMs.value / 1000)))); }

function draw(g) {
  const wallColor = '#e3e0d6', b = BOUND * dx;
  g.rect(0, 0, 1, b, { fill: wallColor }); g.rect(0, 1 - b, 1, b, { fill: wallColor });
  g.rect(0, 0, b, 1, { fill: wallColor }); g.rect(1 - b, 0, b, 1, { fill: wallColor });
  let vmax = 0;
  for (let p = 0; p < count; p++) {
    const v = Math.hypot(pvx[p], pvy[p]);
    if (!(v < 1e6)) { vmax = Infinity; break; }
    vmax = Math.max(vmax, v);
  }
  g.points(px, py, count, { radius: 1.8, colors: pcol });

  const h = stepMs.value / 1000, c = Math.sqrt((lameLambda() + 2 * lameMu()) / rho);
  Sim.show('h c / dx', (h * c / dx).toFixed(1));
  Sim.show('Newton iterations', stats.newton);
  Sim.show('CG iterations', stats.cg);
  Sim.show('line search halvings', stats.halvings);
  Sim.show('residual (m/s)', stats.residual.toExponential(1));
  if (scene.value === 'cantilever' && tipParticle >= 0) {
    const Ep = youngs.value / (1 - nu * nu);   // plane strain
    const predicted = 3 * rho * 9.81 * Math.pow(BEAM_L, 4) / (2 * Ep * BEAM_T * BEAM_T);
    const defl = tipY0 - py[tipParticle];
    Sim.show('tip deflection', (defl * 100).toFixed(2) + ' cm');
    Sim.show('beam theory', (predicted * 100).toFixed(2) + ' cm');
    Sim.plot('tip deflection (cm)', isFinite(defl) ? Math.max(-20, Math.min(20, defl * 100)) : 20, { color: Sim.colors.blue });
  }
  if (!(vmax < 50)) g.text('Exploded! Press R.', 0.5, 0.5, { color: Sim.colors.red, size: 18, align: 'center' });
}

Sim.run({
  title: 'Implicit MPM: Newton, CG and a line search (starter)', view: [0, 0, 1, 1], reset, step, draw, grid: false,
  dt: () => stepMs.value / 1000, stepsPerFrame: substeps,
});
