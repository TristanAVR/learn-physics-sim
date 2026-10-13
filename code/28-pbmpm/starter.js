// Chapter 28 - Position-Based MPM (starter)
// PB-MPM (Lewin 2024): particles carry a displacement d_p (= h v over the step) and a deformation displacement
// D_p (the paper's D; in APIC terms h C_p). The transfers (P2G, grid update, G2P) and the explicit MPM
// comparison are written. Your job, in the TODO blocks:
//   1. the iteration loop of Algorithm 1 in stepPBMPM, and the integrate step,
//   2. the liquid constraint,
//   3. the elastic constraint.
// Until 1 is done nothing moves; until 2 and 3 are done everything falls like dust.

const solver = Sim.select('solver', ['PB-MPM', 'explicit MPM'], 'PB-MPM', { reset: true });
const scene = Sim.select('scene', ['dam break + jelly + sand', 'dam break', 'jelly and sand'], 'dam break + jelly + sand', { reset: true });
const stepMs = Sim.logSlider('time step h (ms)', 0.22, 1.53, 1000 / 240);
const iterations = Sim.slider('iterations K', 1, 20, 4, { step: 1 });
const liquidRelax = Sim.slider('liquid relaxation', 0, 2, 1.5, { step: 0.05 });
const elasticRelax = Sim.slider('elastic relaxation', 0, 2, 1.5, { step: 0.05 });
const viscosity = Sim.slider('liquid viscosity', 0, 1, 0.02, { step: 0.01 });
const elasticityRatio = Sim.slider('elasticity ratio (shape vs volume)', 0, 1, 1, { step: 0.05 });
const gridVolume = Sim.checkbox('liquid: correct volume from the grid', true);

// ---------------------------------------------------------------------------------------------
// Grid
// ---------------------------------------------------------------------------------------------
const N = 64, dx = 1 / N, invDx = N, NN = N + 1, BOUND = 3;
const gravity = -9.81, rho = 1000;
const gm = new Float64Array(NN * NN);           // node mass
const gdx = new Float64Array(NN * NN), gdy = new Float64Array(NN * NN);   // displacement momentum, then displacement
const gvol = new Float64Array(NN * NN);         // particle volume per node (for the liquid's grid volume)

// ---------------------------------------------------------------------------------------------
// Particles
// ---------------------------------------------------------------------------------------------
const LIQUID = 0, ELASTIC = 1, SAND = 2;
const MAXP = 5000;
const px = new Float64Array(MAXP), py = new Float64Array(MAXP);
const pdx = new Float64Array(MAXP), pdy = new Float64Array(MAXP);   // displacement d_p over one step (m)
const pD = new Float64Array(4 * MAXP);         // deformation displacement D_p (row-major 2x2, dimensionless)
const pF = new Float64Array(4 * MAXP);         // deformation gradient F_p (elastic and sand)
const pJ = new Float64Array(MAXP);             // liquid: tracked volume ratio J (the paper's L is 1/J)
const pLogJp = new Float64Array(MAXP);         // sand: stored dilation (volume correction)
const pm = new Float64Array(MAXP), pV0 = new Float64Array(MAXP);
const pmat = new Uint8Array(MAXP);
const palive = new Uint8Array(MAXP);
const pcol = new Array(MAXP);
const pvx = new Float64Array(MAXP), pvy = new Float64Array(MAXP), pC = new Float64Array(4 * MAXP);  // explicit solver
let count = 0, removed = 0, lastH = 0;

function addBlock(x0, y0, x1, y1, mat, color) {
  const s = dx / 2;
  for (let x = x0 + s / 2; x < x1; x += s) {
    for (let y = y0 + s / 2; y < y1; y += s) {
      if (count >= MAXP) return;
      const p = count++;
      px[p] = x; py[p] = y; pdx[p] = 0; pdy[p] = 0; pvx[p] = 0; pvy[p] = 0;
      pD.fill(0, 4 * p, 4 * p + 4); pC.fill(0, 4 * p, 4 * p + 4);
      pF[4 * p] = 1; pF[4 * p + 1] = 0; pF[4 * p + 2] = 0; pF[4 * p + 3] = 1;
      pJ[p] = 1; pLogJp[p] = 0;
      pV0[p] = s * s; pm[p] = rho * s * s;
      pmat[p] = mat; palive[p] = 1; pcol[p] = color;
    }
  }
}

function reset() {
  count = 0; removed = 0; lastH = 0;
  const c = Sim.colors, b = BOUND * dx;
  if (scene.value !== 'jelly and sand') addBlock(b, b, b + 0.3, b + 0.55, LIQUID, c.blue);
  if (scene.value !== 'dam break') {
    addBlock(0.55, 0.30, 0.71, 0.46, ELASTIC, c.red);
    addBlock(0.75, b, 0.91, b + 0.22, SAND, '#c99a2e');
  }
}

// ---------------------------------------------------------------------------------------------
// Small 2x2 helpers
// ---------------------------------------------------------------------------------------------
// Polar SVD F = U diag(s0, s1) V^T, as in chapter 26. U = rotation (cu, su), V = rotation (cv, sv).
const svd = { cu: 1, su: 0, cv: 1, sv: 0, s0: 1, s1: 1 };
function svd2(F00, F01, F10, F11) {
  let c = F00 + F11, s = F10 - F01;
  const r = Math.hypot(c, s);
  if (r > 1e-12) { c /= r; s /= r; } else { c = 1; s = 0; }
  const a = c * F00 + s * F10, b = c * F01 + s * F11, d = -s * F01 + c * F11;
  const e = Math.hypot(a - d, 2 * b);
  let cv = 1, sv = 0;
  if (e > 1e-12) { const c2 = (a - d) / e; cv = Math.sqrt(0.5 * (1 + c2)); sv = (b >= 0 ? 1 : -1) * Math.sqrt(0.5 * (1 - c2)); }
  svd.s0 = a * cv * cv + 2 * b * cv * sv + d * sv * sv;
  svd.s1 = a * sv * sv - 2 * b * cv * sv + d * cv * cv;
  svd.cv = cv; svd.sv = sv; svd.cu = c * cv - s * sv; svd.su = s * cv + c * sv;
}
// out = U diag(s0, s1) V^T
function compose(s0, s1, out) {
  const cu = svd.cu, su = svd.su, cv = svd.cv, sv = svd.sv;
  out[0] = cu * s0 * cv + su * s1 * sv; out[1] = cu * s0 * sv - su * s1 * cv;
  out[2] = su * s0 * cv - cu * s1 * sv; out[3] = su * s0 * sv + cu * s1 * cv;
}
const T = new Float64Array(4);      // scratch: target shape A
const Fs = new Float64Array(4);     // scratch: candidate F* = (I + D) F

// ---------------------------------------------------------------------------------------------
// 1. Constraint solve (paper Algorithm 2). Each particle edits its own D, independently of all others.
// ---------------------------------------------------------------------------------------------
function solveConstraints() {
  const visc = viscosity.value, relL = liquidRelax.value, relE = elasticRelax.value, beta = elasticityRatio.value;
  for (let p = 0; p < count; p++) {
    if (!palive[p]) continue;
    const k = 4 * p;
    let D00 = pD[k], D01 = pD[k + 1], D10 = pD[k + 2], D11 = pD[k + 3];
    if (pmat[p] === LIQUID) {
      // TODO (exercise 2): the liquid constraint.
      //   Viscosity: subtract visc times the symmetric part of D, sym(D) = (D + D^T) / 2.
      //   Volume: the step will change the volume to (1 + tr D) J, with J = pJ[p]. Find the alpha for which
      //   adding alpha I to D makes that exactly 1, and add relL * alpha to the diagonal of D.
    } else {
      const F00 = pF[k], F01 = pF[k + 1], F10 = pF[k + 2], F11 = pF[k + 3];
      // The deformation gradient this D would produce: F* = (I + D) F
      Fs[0] = (1 + D00) * F00 + D01 * F10; Fs[1] = (1 + D00) * F01 + D01 * F11;
      Fs[2] = D10 * F00 + (1 + D11) * F10; Fs[3] = D10 * F01 + (1 + D11) * F11;
      svd2(Fs[0], Fs[1], Fs[2], Fs[3]);
      let s0 = 1, s1 = 1, detLo = 0.1, detHi = 1000;
      if (pmat[p] === SAND) {
        // Sand resists compression but not stretching: the shape target keeps any stretch (s >= 1),
        // unless the particle is still dilated from earlier, in which case it does not resist at all.
        s0 = svd.s0; s1 = svd.s1;
        if (pLogJp[p] === 0) { s0 = Math.max(s0, 1); s1 = Math.max(s1, 1); }
        detHi = 1;
      }
      // TODO (exercise 3): the elastic constraint.
      //   Shape target: compose(s0, s1, T) writes U diag(s0, s1) V^T (the rotation R for elastic particles).
      //   Volume target: F* / sqrt(|det F*|), with |det F*| clamped to [detLo, detHi] and the sign of det kept.
      //   Target A = beta * shape + (1 - beta) * volume, into T.
      //   The D that would produce A exactly is A F^-1 - I (F = pF, the start-of-step F). Move D toward it:
      //   D += relE * (A F^-1 - I - D).
      if (pmat[p] === SAND) {
        // as in the PB-MPM code, sand also gets the viscous damping, which helps piles settle
        const sOff = 0.5 * (D01 + D10);
        D00 -= visc * D00; D11 -= visc * D11; D01 -= visc * sOff; D10 -= visc * sOff;
      }
    }
    pD[k] = D00; pD[k + 1] = D01; pD[k + 2] = D10; pD[k + 3] = D11;
  }
}

// ---------------------------------------------------------------------------------------------
// 2-4. P2G of mass and displacement momentum, grid update, G2P.
// ---------------------------------------------------------------------------------------------
const wx = new Float64Array(3), wy = new Float64Array(3);
function weights(f, w) {
  w[0] = 0.5 * (1.5 - f) * (1.5 - f); w[1] = 0.75 - (f - 1) * (f - 1); w[2] = 0.5 * (f - 0.5) * (f - 0.5);
}

function particleToGrid() {
  gm.fill(0); gdx.fill(0); gdy.fill(0); gvol.fill(0);
  for (let p = 0; p < count; p++) {
    if (!palive[p]) continue;
    const bx = Math.floor(px[p] * invDx - 0.5), by = Math.floor(py[p] * invDx - 0.5);
    const fx = px[p] * invDx - bx, fy = py[p] * invDx - by;
    weights(fx, wx); weights(fy, wy);
    const k = 4 * p, m = pm[p];
    const D00 = pD[k], D01 = pD[k + 1], D10 = pD[k + 2], D11 = pD[k + 3];
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const node = (bx + a) * NN + (by + b);
        const w = wx[a] * wy[b];
        const rx = (a - fx) * dx, ry = (b - fy) * dx;      // x_i - x_p
        gm[node] += w * m;
        gdx[node] += w * m * (pdx[p] + D00 * rx + D01 * ry);   // m w (d + D (x_i - x_p))
        gdy[node] += w * m * (pdy[p] + D10 * rx + D11 * ry);
        gvol[node] += w * pV0[p];
      }
    }
  }
}

// Normalise by mass, then collide the DISPLACED grid positions x_i + d_i with the walls: a node whose
// displacement would carry it out of the box is stopped (both components, as in the PB-MPM code).
function gridUpdate() {
  const lo = BOUND * dx, hi = 1 - BOUND * dx;
  for (let i = 0; i < NN; i++) {
    for (let j = 0; j < NN; j++) {
      const node = i * NN + j, m = gm[node];
      if (m <= 1e-12) { gdx[node] = 0; gdy[node] = 0; continue; }
      let ddx = gdx[node] / m, ddy = gdy[node] / m;
      const X = i * dx + ddx, Y = j * dx + ddy;
      if (X < lo || X > hi || Y < lo || Y > hi) { ddx = 0; ddy = 0; }
      gdx[node] = ddx; gdy[node] = ddy;
    }
  }
}

// d_p = sum_i w d_i,  D_p = (4 / dx^2) sum_i w d_i (x_i - x_p)^T, and the liquid's grid volume measure.
function gridToParticle() {
  const s = 4 * invDx * invDx, cellArea = dx * dx, useGrid = gridVolume.value;
  for (let p = 0; p < count; p++) {
    if (!palive[p]) continue;
    const bx = Math.floor(px[p] * invDx - 0.5), by = Math.floor(py[p] * invDx - 0.5);
    const fx = px[p] * invDx - bx, fy = py[p] * invDx - by;
    weights(fx, wx); weights(fy, wy);
    let ddx = 0, ddy = 0, B00 = 0, B01 = 0, B10 = 0, B11 = 0, fill = 0;
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const node = (bx + a) * NN + (by + b);
        const w = wx[a] * wy[b];
        const rx = (a - fx) * dx, ry = (b - fy) * dx;
        const ux = gdx[node], uy = gdy[node];
        ddx += w * ux; ddy += w * uy;
        B00 += w * ux * rx; B01 += w * ux * ry; B10 += w * uy * rx; B11 += w * uy * ry;
        fill += w * gvol[node] / cellArea;
      }
    }
    const k = 4 * p;
    pdx[p] = ddx; pdy[p] = ddy;
    pD[k] = s * B00; pD[k + 1] = s * B01; pD[k + 2] = s * B10; pD[k + 3] = s * B11;
    // How densely is the liquid packed around this particle? fill = 1 at rest density. If it is packed more
    // tightly than that, pull the tracked J toward the grid's answer 1 / fill (an objective volume measure).
    if (useGrid && pmat[p] === LIQUID && fill > 1) pJ[p] += 0.1 * (1 / fill - pJ[p]);
  }
}

// ---------------------------------------------------------------------------------------------
// 5. Integrate (once per step, after the iterations).
// ---------------------------------------------------------------------------------------------
// Drucker-Prager return mapping with volume correction (chapter 26), for the singular values in svd;
// friction angle 30 degrees, lambda / mu = 1.5 (nu = 0.3). Writes the projected values into svd.
const SIN_PHI = Math.sin(30 * Math.PI / 180), DP_ALPHA = Math.sqrt(2 / 3) * 2 * SIN_PHI / (3 - SIN_PHI);
function sandReturnMap(p) {
  const e0 = Math.log(Math.max(Math.abs(svd.s0), 1e-6)) + 0.5 * pLogJp[p];
  const e1 = Math.log(Math.max(Math.abs(svd.s1), 1e-6)) + 0.5 * pLogJp[p], tr = e0 + e1;
  if (tr >= 0) { pLogJp[p] = tr; svd.s0 = 1; svd.s1 = 1; return; }
  pLogJp[p] = 0;
  const h0 = e0 - 0.5 * tr, h1 = e1 - 0.5 * tr, hn = Math.hypot(h0, h1);
  const dGamma = hn + 2.5 * tr * DP_ALPHA;
  if (dGamma > 0 && hn > 1e-12) { svd.s0 = Math.exp(e0 - dGamma * h0 / hn); svd.s1 = Math.exp(e1 - dGamma * h1 / hn); }
  else { svd.s0 = Math.exp(e0); svd.s1 = Math.exp(e1); }
}

// After the iterations, a solid particle's new F (in Fn) goes through the safety net and plasticity.
// Returns false if the particle had to be removed.
const Fn = new Float64Array(4);
function finishSolid(p) {
  svd2(Fn[0], Fn[1], Fn[2], Fn[3]);
  // The paper deletes particles whose F is hopelessly ill-conditioned (cond F > 10^6)...
  if (!(Math.abs(svd.s0) < 1e6 * Math.abs(svd.s1) && Math.abs(svd.s1) < 1e6 * Math.abs(svd.s0))) {
    palive[p] = 0; removed++; return false;
  }
  // ...and the code clamps the singular values as a safety net.
  svd.s0 = Math.min(Math.max(svd.s0, 0.2), 1e4); svd.s1 = Math.min(Math.max(svd.s1, 0.2), 1e4);
  if (pmat[p] === SAND) sandReturnMap(p);
  compose(svd.s0, svd.s1, Fn);
  pF.set(Fn, 4 * p);
  return true;
}

function integrate(h) {
  const lo = dx, hi = 1 - dx;
  for (let p = 0; p < count; p++) {
    if (!palive[p]) continue;
    const k = 4 * p;
    const D00 = pD[k], D01 = pD[k + 1], D10 = pD[k + 2], D11 = pD[k + 3];
    if (pmat[p] === LIQUID) {
      // TODO (exercise 1): the liquid tracks only its volume: J <- J (1 + tr D), kept above 0.1.
    } else {
      // TODO (exercise 1): F <- (I + D) F into Fn[0..3], then: if (!finishSolid(p)) continue;
      //   (finishSolid applies the safety net and the sand plasticity and stores F.)
    }
    // TODO (exercise 1): move the particle by its displacement d, then add gravity for the next step,
    //   d <- d + g h^2 (gravity is the constant `gravity`, negative, in y).
    if (!(px[p] > lo)) px[p] = lo; else if (px[p] > hi) px[p] = hi;
    if (!(py[p] > lo)) py[p] = lo; else if (py[p] > hi) py[p] = hi;
  }
}

// One PB-MPM step (paper Algorithm 1).
function stepPBMPM(h) {
  // d and D are displacements over one step; if the step size changed, rescale them.
  if (lastH > 0 && h !== lastH) {
    const r = h / lastH;
    for (let p = 0; p < count; p++) { pdx[p] *= r; pdy[p] *= r; for (let i = 0; i < 4; i++) pD[4 * p + i] *= r; }
    for (let p = 0; p < count; p++) pdy[p] += gravity * (h * h - r * lastH * lastH);
  }
  if (lastH === 0) for (let p = 0; p < count; p++) pdy[p] += gravity * h * h;
  lastH = h;
  // TODO (exercise 1): Algorithm 1. Repeat iterations.value times: solveConstraints, particleToGrid,
  //   gridUpdate, gridToParticle. Then integrate(h) once.
}

// ---------------------------------------------------------------------------------------------
// Explicit MLS-MPM at the same time step, for comparison (chapters 25 and 26). Liquid: K (J-1) J,
// elastic: fixed corotated, sand: Hencky with the same Drucker-Prager return mapping.
// ---------------------------------------------------------------------------------------------
const K_LIQUID = 1e5, E_SOLID = 1e5, NU = 0.3;
const MU = E_SOLID / (2 * (1 + NU)), LA = E_SOLID * NU / ((1 + NU) * (1 - 2 * NU));
function stepExplicit(h) {
  gm.fill(0); gdx.fill(0); gdy.fill(0);
  const s = 4 * invDx * invDx;
  for (let p = 0; p < count; p++) {
    if (!palive[p]) continue;
    const bx = Math.floor(px[p] * invDx - 0.5), by = Math.floor(py[p] * invDx - 0.5);
    const fx = px[p] * invDx - bx, fy = py[p] * invDx - by;
    weights(fx, wx); weights(fy, wy);
    const k = 4 * p, m = pm[p];
    let A00, A01, A10, A11;                 // V_p^0 P F^T (the Kirchhoff stress times volume)
    if (pmat[p] === LIQUID) { A00 = A11 = pV0[p] * K_LIQUID * (pJ[p] - 1) * pJ[p]; A01 = A10 = 0; }
    else {
      svd2(pF[k], pF[k + 1], pF[k + 2], pF[k + 3]);
      let t0, t1;
      if (pmat[p] === SAND) {
        const e0 = Math.log(Math.max(Math.abs(svd.s0), 1e-6)), e1 = Math.log(Math.max(Math.abs(svd.s1), 1e-6));
        t0 = 2 * MU * e0 + LA * (e0 + e1); t1 = 2 * MU * e1 + LA * (e0 + e1);
      } else {
        const J = svd.s0 * svd.s1;
        t0 = 2 * MU * (svd.s0 - 1) * svd.s0 + LA * (J - 1) * J; t1 = 2 * MU * (svd.s1 - 1) * svd.s1 + LA * (J - 1) * J;
      }
      const c = svd.cu, sn = svd.su, V = pV0[p];
      A00 = V * (c * c * t0 + sn * sn * t1); A01 = A10 = V * c * sn * (t0 - t1); A11 = V * (sn * sn * t0 + c * c * t1);
    }
    const Q00 = m * pC[k] - h * s * A00, Q01 = m * pC[k + 1] - h * s * A01;
    const Q10 = m * pC[k + 2] - h * s * A10, Q11 = m * pC[k + 3] - h * s * A11;
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const node = (bx + a) * NN + (by + b);
        const w = wx[a] * wy[b], rx = (a - fx) * dx, ry = (b - fy) * dx;
        gm[node] += w * m;
        gdx[node] += w * (m * pvx[p] + Q00 * rx + Q01 * ry);
        gdy[node] += w * (m * pvy[p] + Q10 * rx + Q11 * ry);
      }
    }
  }
  const lo = BOUND * dx, hi = 1 - BOUND * dx;
  for (let i = 0; i < NN; i++) {
    for (let j = 0; j < NN; j++) {
      const node = i * NN + j, m = gm[node];
      if (m <= 1e-12) { gdx[node] = 0; gdy[node] = 0; continue; }
      let vx = gdx[node] / m, vy = gdy[node] / m + h * gravity;
      const X = i * dx + h * vx, Y = j * dx + h * vy;          // same wall rule, on h v
      if (X < lo || X > hi || Y < lo || Y > hi) { vx = 0; vy = 0; }
      gdx[node] = vx; gdy[node] = vy;
    }
  }
  for (let p = 0; p < count; p++) {
    if (!palive[p]) continue;
    const bx = Math.floor(px[p] * invDx - 0.5), by = Math.floor(py[p] * invDx - 0.5);
    const fx = px[p] * invDx - bx, fy = py[p] * invDx - by;
    weights(fx, wx); weights(fy, wy);
    let vx = 0, vy = 0, B00 = 0, B01 = 0, B10 = 0, B11 = 0;
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const node = (bx + a) * NN + (by + b);
        const w = wx[a] * wy[b], rx = (a - fx) * dx, ry = (b - fy) * dx;
        const ux = gdx[node], uy = gdy[node];
        vx += w * ux; vy += w * uy;
        B00 += w * ux * rx; B01 += w * ux * ry; B10 += w * uy * rx; B11 += w * uy * ry;
      }
    }
    const k = 4 * p;
    const C00 = s * B00, C01 = s * B01, C10 = s * B10, C11 = s * B11;
    pC[k] = C00; pC[k + 1] = C01; pC[k + 2] = C10; pC[k + 3] = C11;
    pvx[p] = vx; pvy[p] = vy;
    if (pmat[p] === LIQUID) pJ[p] *= (1 + h * C00) * (1 + h * C11) - h * h * C01 * C10;
    else {
      const F00 = pF[k], F01 = pF[k + 1], F10 = pF[k + 2], F11 = pF[k + 3];
      pF[k] = F00 + h * (C00 * F00 + C01 * F10); pF[k + 1] = F01 + h * (C00 * F01 + C01 * F11);
      pF[k + 2] = F10 + h * (C10 * F00 + C11 * F10); pF[k + 3] = F11 + h * (C10 * F01 + C11 * F11);
      if (pmat[p] === SAND) {
        svd2(pF[k], pF[k + 1], pF[k + 2], pF[k + 3]);
        sandReturnMap(p);
        compose(svd.s0, svd.s1, Fn);
        pF.set(Fn, k);
      }
    }
    pdx[p] = h * vx; pdy[p] = h * vy;           // kept in step with v, for the readouts
    px[p] += h * vx; py[p] += h * vy;
    const lo2 = dx, hi2 = 1 - dx;
    if (!(px[p] > lo2)) px[p] = lo2; else if (px[p] > hi2) px[p] = hi2;
    if (!(py[p] > lo2)) py[p] = lo2; else if (py[p] > hi2) py[p] = hi2;
  }
}

function step(h) {
  if (solver.value === 'PB-MPM') stepPBMPM(h); else stepExplicit(h);
}

// ---------------------------------------------------------------------------------------------
// Time stepping: one step every h seconds of simulated time, at most a few per frame.
// ---------------------------------------------------------------------------------------------
const FRAME = 1 / 60;
let acc = 0;
function stepsThisFrame() {
  const h = stepMs.value / 1000;
  acc += FRAME;
  let n = Math.floor(acc / h + 1e-9);
  acc -= n * h;
  if (n > 12) { n = 12; acc = 0; }
  return n;
}

// ---------------------------------------------------------------------------------------------
// Drawing and readouts
// ---------------------------------------------------------------------------------------------
function draw(g) {
  const wallColor = '#e3e0d6', b = BOUND * dx;
  g.rect(0, 0, 1, b, { fill: wallColor }); g.rect(0, 1 - b, 1, b, { fill: wallColor });
  g.rect(0, 0, b, 1, { fill: wallColor }); g.rect(1 - b, 0, b, 1, { fill: wallColor });
  const h = stepMs.value / 1000;
  let vmax = 0, liquidArea = 0, liquidN = 0;
  for (let p = 0; p < count; p++) {
    if (!palive[p]) continue;
    const v = Math.hypot(pdx[p], pdy[p] - (solver.value === 'PB-MPM' ? gravity * h * h : 0)) / h;
    if (!(v < 1e6)) { vmax = Infinity; break; }
    if (v > vmax) vmax = v;
    if (pmat[p] === LIQUID) { liquidArea += pJ[p]; liquidN++; }
  }
  g.points(px, py, count, { radius: 1.6, colors: pcol });
  Sim.show('particles', count - removed + (removed ? ' (' + removed + ' removed)' : ''));
  Sim.show('step', (1 / h).toFixed(0) + ' Hz, h = ' + (h * 1000).toFixed(1) + ' ms');
  Sim.show('max speed', isFinite(vmax) ? vmax.toFixed(2) + ' m/s' : 'inf');
  if (liquidN) {
    Sim.show('liquid volume (mean J)', (liquidArea / liquidN).toFixed(3));
    Sim.plot('liquid mean J', liquidArea / liquidN, { color: Sim.colors.blue, min: 0.8, max: 1.1 });
  }
  if (!(vmax < 100)) g.text('Exploded! Press R.', 0.5, 0.6, { color: Sim.colors.red, size: 18, align: 'center' });
}

Sim.run({
  plotCorner: 'top-right',
  title: 'Position-based MPM (starter)', view: [0, 0, 1, 1], reset, step, draw, grid: false,
  dt: () => stepMs.value / 1000, stepsPerFrame: stepsThisFrame,
});
