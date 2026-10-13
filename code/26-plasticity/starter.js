// Chapter 26 - Plasticity: Snow and Sand (starter)
// Explicit MLS-MPM (chapter 25's compact form) with elastoplastic materials. Each particle stores only the
// elastic part F_E of its deformation gradient. After every F update, a return mapping projects F_E back
// onto the material's elastic region, and the rest of the deformation is forgotten (it goes into F_P).
// The dough return mapping is done for you as an example. Your job, in the TODO blocks:
//   1. the snow return mapping and its hardening (Stomakhin et al. 2013),
//   2. the Drucker-Prager return mapping for sand (Klar et al. 2016).
// Until then, snow and sand are purely elastic and bounce like jelly.

const scene = Sim.select('scene', ['snowball', 'column', 'pour', 'blob drop'], 'snowball', { reset: true });
const material = Sim.select('material', ['snow', 'sand', 'dough', 'jelly (elastic)'], 'snow', { reset: true });
const thetaC = Sim.slider('snow critical compression theta_c', 0.005, 0.1, 0.025, { step: 0.001 });
const thetaS = Sim.slider('snow critical stretch theta_s', 0.001, 0.05, 0.0075, { step: 0.0005 });
const xi = Sim.slider('snow hardening xi', 0, 20, 10, { step: 0.5 });
const frictionDeg = Sim.slider('sand friction angle phi (deg)', 5, 50, 30, { step: 1 });
const doughYield = Sim.slider('dough yield stretch', 0.01, 0.3, 0.12, { step: 0.01 });
const walls = Sim.select('walls', ['sticky', 'slip'], 'sticky');

// ---------------------------------------------------------------------------------------------
// Grid
// ---------------------------------------------------------------------------------------------
const N = 64, dx = 1 / N, invDx = N, NN = N + 1, BOUND = 3;
const gravity = -9.81;
const gm = new Float64Array(NN * NN), gvx = new Float64Array(NN * NN), gvy = new Float64Array(NN * NN);

// ---------------------------------------------------------------------------------------------
// Materials (SI units). Snow from Stomakhin et al. 2013, sand from Klar et al. 2016.
// ---------------------------------------------------------------------------------------------
const SNOW = 0, SAND = 1, DOUGH = 2, JELLY = 3;
const MATERIALS = [
  { name: 'snow', E: 1.4e5, nu: 0.2, rho: 400 },
  { name: 'sand', E: 3.537e5, nu: 0.3, rho: 2200 },
  { name: 'dough', E: 6e4, nu: 0.3, rho: 1100 },
  { name: 'jelly (elastic)', E: 4e4, nu: 0.3, rho: 1000 },
];
for (const m of MATERIALS) {
  m.mu = m.E / (2 * (1 + m.nu));
  m.la = m.E * m.nu / ((1 + m.nu) * (1 - 2 * m.nu));
}

// ---------------------------------------------------------------------------------------------
// Particles
// ---------------------------------------------------------------------------------------------
const MAXP = 4000;
const px = new Float64Array(MAXP), py = new Float64Array(MAXP);
const pvx = new Float64Array(MAXP), pvy = new Float64Array(MAXP);
const pC = new Float64Array(4 * MAXP);
const pF = new Float64Array(4 * MAXP);    // ELASTIC deformation gradient F_E (row-major 2x2)
const pU = new Float64Array(2 * MAXP);    // F_E = U diag(s0, s1) V^T: cos and sin of the rotation U
const pS = new Float64Array(2 * MAXP);    // singular values s0, s1 of F_E
const pJp = new Float64Array(MAXP);       // snow: plastic volume ratio J_P = det F_P
const pLogJp = new Float64Array(MAXP);    // sand: volume correction (log of the volume gained by dilating)
const pm = new Float64Array(MAXP), pV0 = new Float64Array(MAXP);
const pmat = new Uint8Array(MAXP);
const pcol = new Array(MAXP);
let count = 0;

function addParticle(x, y, vx, vy, mat) {
  if (count >= MAXP) return;
  const p = count++;
  const s = dx / 2;
  px[p] = x; py[p] = y; pvx[p] = vx; pvy[p] = vy;
  pC[4 * p] = pC[4 * p + 1] = pC[4 * p + 2] = pC[4 * p + 3] = 0;
  pF[4 * p] = 1; pF[4 * p + 1] = 0; pF[4 * p + 2] = 0; pF[4 * p + 3] = 1;
  pU[2 * p] = 1; pU[2 * p + 1] = 0; pS[2 * p] = 1; pS[2 * p + 1] = 1;
  pJp[p] = 1; pLogJp[p] = 0;
  pV0[p] = s * s;
  pm[p] = MATERIALS[mat].rho * pV0[p];
  pmat[p] = mat;
}

function addShape(inside, x0, y0, x1, y1, vx, vy, mat) {
  const s = dx / 2;
  for (let x = x0 + s / 2; x < x1; x += s) {
    for (let y = y0 + s / 2; y < y1; y += s) {
      if (inside(x, y)) addParticle(x, y, vx, vy, mat);
    }
  }
}

let emitTime = 0, pourMaterial = SAND;
function reset() {
  count = 0;
  emitTime = 0;
  const mat = MATERIALS.findIndex((m) => m.name === material.value);
  pourMaterial = mat;
  const floor = BOUND * dx;
  if (scene.value === 'snowball') {
    const cx = 0.3, cy = 0.45, r = 0.1;
    addShape((x, y) => (x - cx) * (x - cx) + (y - cy) * (y - cy) < r * r, cx - r, cy - r, cx + r, cy + r, 3, -4, mat);
  } else if (scene.value === 'column') {
    addShape(() => true, 0.44, floor, 0.56, floor + 0.36, 0, 0, mat);
  } else if (scene.value === 'blob drop') {
    const cx = 0.5, cy = 0.45, r = 0.11;
    addShape((x, y) => (x - cx) * (x - cx) + (y - cy) * (y - cy) < r * r, cx - r, cy - r, cx + r, cy + r, 0, -1, mat);
  }
  // 'pour' starts empty; emit() adds material from the nozzle.
}

// Pour: a nozzle 3 cm wide at 0.45 m lets out a stream at 0.5 m/s until POUR_LIMIT particles exist.
const POUR_LIMIT = 1300, NOZZLE_Y = 0.45;
function emit(h) {
  if (scene.value !== 'pour' || count >= POUR_LIMIT) return;
  const s = dx / 2, speed = 0.5;
  emitTime += h;
  if (emitTime * speed < s) return;     // one new row each time the stream has moved one spacing
  emitTime -= s / speed;
  for (let x = 0.5 - 0.015 + s / 2; x < 0.5 + 0.015; x += s) addParticle(x, NOZZLE_Y, 0, -speed, pourMaterial);
}

// ---------------------------------------------------------------------------------------------
// 2x2 polar SVD: F = U diag(s0, s1) V^T with U, V rotations [[c, -s], [s, c]] (s1 < 0 if det F < 0).
// ---------------------------------------------------------------------------------------------
const svd = { cu: 1, su: 0, cv: 1, sv: 0, s0: 1, s1: 1 };
function svd2(F00, F01, F10, F11) {
  // Polar decomposition F = R S: R is the rotation with (cos, sin) along (F00 + F11, F10 - F01).
  let c = F00 + F11, s = F10 - F01;
  const r = Math.hypot(c, s);
  if (r > 1e-12) { c /= r; s /= r; } else { c = 1; s = 0; }
  const a = c * F00 + s * F10, b = c * F01 + s * F11, d = -s * F01 + c * F11;   // S = R^T F = [[a, b], [b, d]]
  // S = V diag(s0, s1) V^T, with V the rotation by phi, tan(2 phi) = 2b / (a - d). Half-angle formulas:
  const e = Math.hypot(a - d, 2 * b);
  let cv = 1, sv = 0;
  if (e > 1e-12) {
    const c2 = (a - d) / e;
    cv = Math.sqrt(0.5 * (1 + c2));
    sv = (b >= 0 ? 1 : -1) * Math.sqrt(0.5 * (1 - c2));
  }
  svd.s0 = a * cv * cv + 2 * b * cv * sv + d * sv * sv;
  svd.s1 = a * sv * sv - 2 * b * cv * sv + d * cv * cv;
  svd.cv = cv; svd.sv = sv;
  svd.cu = c * cv - s * sv;            // U = R V
  svd.su = s * cv + c * sv;
}

// Rebuild F_E = U diag(svd.s0, svd.s1) V^T into pF at offset k.
function composeF(k) {
  const cu = svd.cu, su = svd.su, cv = svd.cv, sv = svd.sv, s0 = svd.s0, s1 = svd.s1;
  pF[k] = cu * s0 * cv + su * s1 * sv;
  pF[k + 1] = cu * s0 * sv - su * s1 * cv;
  pF[k + 2] = su * s0 * cv - cu * s1 * sv;
  pF[k + 3] = su * s0 * sv + cu * s1 * cv;
}

// ---------------------------------------------------------------------------------------------
// Return mappings. On entry, svd holds the SVD of the trial F_E; each one moves the singular values
// back into the material's elastic region, updates the plastic state, and rebuilds F_E.
// ---------------------------------------------------------------------------------------------

// Snow (course section 6.5, eqs. 79-86): clamp the singular values to [1 - theta_c, 1 + theta_s].
// F = F_E F_P must not change, so J_P absorbs whatever volume J_E loses: J_P <- J_P * J_trial / J_new.
function snowReturnMap(p) {
  // TODO (exercise 1): clamp svd.s0 and svd.s1 to [1 - thetaC.value, 1 + thetaS.value], multiply pJp[p]
  // by J_trial / J_new (use the absolute value of the trial product), store the clamped values back in
  // svd.s0 and svd.s1, and rebuild F_E with composeF(4 * p). Model it on doughReturnMap below.
}

// Snow hardening (course eq. 87): mu and lambda scale with exp(xi (1 - J_P)).
function snowHardening(Jp) {
  // TODO (exercise 1): return exp(xi (1 - J_P)), with the exponent capped at 3 so the time step stays sane.
  return 1;
}

// Sand (Klar et al. 2016): Drucker-Prager cone in Hencky strain eps = log(sigma), with volume correction:
// pLogJp remembers how much the sand has dilated, and it must be compressed back before it resists.
function sandReturnMap(p) {
  const M = MATERIALS[SAND];
  const sinPhi = Math.sin(frictionDeg.value * Math.PI / 180);
  const alpha = Math.sqrt(2 / 3) * 2 * sinPhi / (3 - sinPhi);
  // TODO (exercise 2):
  //   eps_i = log|s_i| + 0.5 * pLogJp[p]  (Hencky strain plus the stored dilation), tr = eps_0 + eps_1.
  //   Case 1, tr >= 0 (expansion): s_i = 1, and pLogJp[p] = tr.
  //   Otherwise pLogJp[p] = 0, eps_hat = eps - (tr / 2) I, and
  //     dGamma = |eps_hat| + (2 lambda + 2 mu) / (2 mu) * tr * alpha   (M.la, M.mu).
  //   Case 2, dGamma <= 0 (inside the cone): s_i = exp(eps_i).
  //   Case 3, outside: s_i = exp(eps_i - dGamma * eps_hat_i / |eps_hat|).
  //   Store into svd.s0, svd.s1 and rebuild F_E with composeF(4 * p).
}

// Dough: clamp the singular values to [1/(1+y), 1+y], then rescale so the volume is unchanged.
// Shape change beyond the yield stretch is forgotten, volume change is not: the dough flows but stays dough.
function doughReturnMap(p) {
  const hi = 1 + doughYield.value, lo = 1 / hi;
  const J = svd.s0 * svd.s1;
  let s0 = Math.min(Math.max(svd.s0, lo), hi);
  let s1 = Math.min(Math.max(svd.s1, lo), hi);
  const scale = Math.sqrt(Math.abs(J / (s0 * s1)));        // in 2D, scaling both by sqrt(.) restores det
  svd.s0 = s0 * scale; svd.s1 = s1 * scale;
  composeF(4 * p);
}

// ---------------------------------------------------------------------------------------------
// Stress. For an isotropic material, the Kirchhoff stress tau = P F_E^T = U diag(tau0, tau1) U^T,
// so we only need its two principal values. Writes A = V_p^0 tau into A[0..2] (A is symmetric).
// ---------------------------------------------------------------------------------------------
const A = new Float64Array(3);
function particleStress(p) {
  const M = MATERIALS[pmat[p]];
  const s0 = pS[2 * p], s1 = pS[2 * p + 1];
  let t0, t1;
  if (pmat[p] === SAND) {
    // Hencky: tau_i = 2 mu eps_i + lambda tr(eps), eps_i = log s_i. The model Drucker-Prager is built on.
    const e0 = Math.log(Math.max(Math.abs(s0), 1e-6)), e1 = Math.log(Math.max(Math.abs(s1), 1e-6));
    t0 = 2 * M.mu * e0 + M.la * (e0 + e1);
    t1 = 2 * M.mu * e1 + M.la * (e0 + e1);
  } else {
    // Fixed corotated (chapter 25) in principal form: tau_i = 2 mu (s_i - 1) s_i + lambda (J - 1) J.
    const hard = pmat[p] === SNOW ? snowHardening(pJp[p]) : 1;
    const mu = M.mu * hard, la = M.la * hard, J = s0 * s1;
    t0 = 2 * mu * (s0 - 1) * s0 + la * (J - 1) * J;
    t1 = 2 * mu * (s1 - 1) * s1 + la * (J - 1) * J;
  }
  const c = pU[2 * p], s = pU[2 * p + 1], V = pV0[p];
  A[0] = V * (c * c * t0 + s * s * t1);       // A00
  A[1] = V * c * s * (t0 - t1);               // A01 = A10
  A[2] = V * (s * s * t0 + c * c * t1);       // A11
}

// ---------------------------------------------------------------------------------------------
// One explicit MLS-MPM step
// ---------------------------------------------------------------------------------------------
const wx = new Float64Array(3), wy = new Float64Array(3);
function weights(f, w) {
  w[0] = 0.5 * (1.5 - f) * (1.5 - f);
  w[1] = 0.75 - (f - 1) * (f - 1);
  w[2] = 0.5 * (f - 0.5) * (f - 0.5);
}

function step(h) {
  emit(h);
  gm.fill(0); gvx.fill(0); gvy.fill(0);
  const s = 4 * invDx * invDx;           // D_p^-1 for quadratic B-splines

  // P2G, MLS form: the stress enters the affine momentum (chapter 25).
  for (let p = 0; p < count; p++) {
    const bx = Math.floor(px[p] * invDx - 0.5), by = Math.floor(py[p] * invDx - 0.5);
    const fx = px[p] * invDx - bx, fy = py[p] * invDx - by;
    weights(fx, wx); weights(fy, wy);
    particleStress(p);
    const m = pm[p], k = 4 * p;
    const Q00 = m * pC[k] - h * s * A[0], Q01 = m * pC[k + 1] - h * s * A[1];
    const Q10 = m * pC[k + 2] - h * s * A[1], Q11 = m * pC[k + 3] - h * s * A[2];
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const node = (bx + a) * NN + (by + b);
        const w = wx[a] * wy[b];
        const rx = (a - fx) * dx, ry = (b - fy) * dx;
        gm[node] += w * m;
        gvx[node] += w * (m * pvx[p] + Q00 * rx + Q01 * ry);
        gvy[node] += w * (m * pvy[p] + Q10 * rx + Q11 * ry);
      }
    }
  }

  // Grid update: momentum to velocity, gravity, walls.
  const sticky = walls.value === 'sticky';
  for (let i = 0; i < NN; i++) {
    for (let j = 0; j < NN; j++) {
      const node = i * NN + j, m = gm[node];
      if (m <= 0) continue;
      let vx = gvx[node] / m, vy = gvy[node] / m + h * gravity;
      const left = i < BOUND, right = i > N - BOUND, bottom = j < BOUND, top = j > N - BOUND;
      if (left || right || bottom || top) {
        if (sticky) { vx = 0; vy = 0; }
        else {
          if ((left && vx < 0) || (right && vx > 0)) vx = 0;
          if ((bottom && vy < 0) || (top && vy > 0)) vy = 0;
        }
      }
      gvx[node] = vx; gvy[node] = vy;
    }
  }

  // G2P, the trial update F_E <- (I + h C) F_E, then the return mapping.
  for (let p = 0; p < count; p++) {
    const bx = Math.floor(px[p] * invDx - 0.5), by = Math.floor(py[p] * invDx - 0.5);
    const fx = px[p] * invDx - bx, fy = py[p] * invDx - by;
    weights(fx, wx); weights(fy, wy);
    let vx = 0, vy = 0, B00 = 0, B01 = 0, B10 = 0, B11 = 0;
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) {
        const node = (bx + a) * NN + (by + b);
        const w = wx[a] * wy[b];
        const rx = (a - fx) * dx, ry = (b - fy) * dx;
        const ux = gvx[node], uy = gvy[node];
        vx += w * ux; vy += w * uy;
        B00 += w * ux * rx; B01 += w * ux * ry; B10 += w * uy * rx; B11 += w * uy * ry;
      }
    }
    const k = 4 * p;
    const C00 = s * B00, C01 = s * B01, C10 = s * B10, C11 = s * B11;
    pC[k] = C00; pC[k + 1] = C01; pC[k + 2] = C10; pC[k + 3] = C11;
    pvx[p] = vx; pvy[p] = vy;

    // Trial elastic deformation: assume the whole step was elastic.
    const F00 = pF[k], F01 = pF[k + 1], F10 = pF[k + 2], F11 = pF[k + 3];
    pF[k] = F00 + h * (C00 * F00 + C01 * F10);
    pF[k + 1] = F01 + h * (C00 * F01 + C01 * F11);
    pF[k + 2] = F10 + h * (C10 * F00 + C11 * F10);
    pF[k + 3] = F11 + h * (C10 * F01 + C11 * F11);
    svd2(pF[k], pF[k + 1], pF[k + 2], pF[k + 3]);
    // Return mapping: what does not fit in the elastic region is forgotten.
    const mat = pmat[p];
    if (mat === SNOW) snowReturnMap(p);
    else if (mat === SAND) sandReturnMap(p);
    else if (mat === DOUGH) doughReturnMap(p);
    pU[2 * p] = svd.cu; pU[2 * p + 1] = svd.su;
    pS[2 * p] = svd.s0; pS[2 * p + 1] = svd.s1;

    px[p] += h * vx; py[p] += h * vy;
    const lo = dx, hi = 1 - dx;
    if (!(px[p] > lo)) px[p] = lo; else if (px[p] > hi) px[p] = hi;
    if (!(py[p] > lo)) py[p] = lo; else if (py[p] > hi) py[p] = hi;
  }
}

// ---------------------------------------------------------------------------------------------
// Time step chosen every frame from the wave speed and the CFL limit (chapter 25), with some margin.
// ---------------------------------------------------------------------------------------------
const FRAME = 1 / 60, MAX_SUBSTEPS = 60;
let hStep = 1e-3, nSub = 17, vmax = 0;
function chooseStep() {
  const M0 = MATERIALS[pourMaterial];
  let cmax = Math.sqrt((M0.la + 2 * M0.mu) / M0.rho);
  vmax = 0;
  for (let p = 0; p < count; p++) {
    const M = MATERIALS[pmat[p]];
    const hard = pmat[p] === SNOW ? snowHardening(pJp[p]) : 1;
    const c = Math.sqrt((M.la + 2 * M.mu) * hard / M.rho);
    if (c > cmax) cmax = c;
    const v = Math.hypot(pvx[p], pvy[p]);
    if (v > vmax) vmax = v;
  }
  const h = Math.min(0.6 * dx / cmax, 0.5 * dx / Math.max(vmax, 1e-6), 2e-3);
  nSub = Math.min(MAX_SUBSTEPS, Math.ceil(FRAME / h));
  hStep = Math.min(h, FRAME / nSub);
  return nSub;
}

// ---------------------------------------------------------------------------------------------
// Drawing and measurements
// ---------------------------------------------------------------------------------------------
const SNOW_SHADES = ['#4d6b8f', '#6b87a8', '#8ea6c2', '#b3c5d9', '#d3dde8'];
const BASE_COLORS = ['#8ea6c2', '#c99a2e', '#d08a4a', Sim.colors.red];

// Slope of a pile: height profile in columns 2 dx wide, a line fitted to each flank between 20 % and
// 80 % of the peak height, and the two angles averaged.
function pileSlopeDeg() {
  const nb = Math.floor(1 / (2 * dx));
  const top = new Float64Array(nb).fill(-1);
  const floor = BOUND * dx;
  for (let p = 0; p < count; p++) {
    const b = Math.floor(px[p] / (2 * dx));
    if (py[p] - floor > top[b]) top[b] = py[p] - floor;
  }
  let peak = 0, ib = 0;
  for (let b = 0; b < nb; b++) if (top[b] > peak) { peak = top[b]; ib = b; }
  if (peak < 0.03) return NaN;
  let sum = 0, sides = 0;
  for (const dir of [-1, 1]) {
    let n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (let b = ib; b >= 0 && b < nb; b += dir) {
      const hgt = top[b];
      if (hgt < 0) break;
      if (hgt > 0.8 * peak || hgt < 0.2 * peak) continue;
      const x = (b + 0.5) * 2 * dx;
      n++; sx += x; sy += hgt; sxx += x * x; sxy += x * hgt;
    }
    if (n >= 3) {
      const slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
      sum += Math.atan(Math.abs(slope)) * 180 / Math.PI; sides++;
    }
  }
  return sides ? sum / sides : NaN;
}

function draw(g) {
  const wallColor = '#e3e0d6', b = BOUND * dx;
  g.rect(0, 0, 1, b, { fill: wallColor });
  g.rect(0, 1 - b, 1, b, { fill: wallColor });
  g.rect(0, 0, b, 1, { fill: wallColor });
  g.rect(1 - b, 0, b, 1, { fill: wallColor });
  if (scene.value === 'pour') g.rect(0.48, NOZZLE_Y + 0.005, 0.04, 0.03, { fill: Sim.colors.gray });

  let JpMin = Infinity, JpMax = 0;
  for (let p = 0; p < count; p++) {
    const mat = pmat[p];
    if (mat === SNOW) {
      // darker = compacted (J_P < 1), lighter = loosened (J_P > 1)
      const Jp = pJp[p];
      pcol[p] = SNOW_SHADES[Math.max(0, Math.min(4, Math.round(2 + 40 * (Jp - 1))))];
      if (Jp < JpMin) JpMin = Jp;
      if (Jp > JpMax) JpMax = Jp;
    } else pcol[p] = BASE_COLORS[mat];
  }
  g.points(px, py, count, { radius: 1.8, colors: pcol });

  Sim.show('particles', count);
  Sim.show('time step h', (hStep * 1000).toFixed(3) + ' ms x ' + nSub + (nSub * hStep < FRAME * 0.99 ? ' (slow motion)' : ''));
  Sim.show('max speed', vmax.toFixed(2) + ' m/s');
  if (material.value === 'snow' && count > 0) Sim.show('J_P range', JpMin.toFixed(3) + ' .. ' + JpMax.toFixed(3));
  if (scene.value === 'column' || scene.value === 'pour') {
    const a = pileSlopeDeg();
    Sim.show('pile slope', isFinite(a) ? a.toFixed(1) + ' deg' : '-');
  }
}

Sim.run({
  title: 'Elastoplastic MPM: snow, sand, dough (starter)', view: [0, 0, 1, 1], reset, step, draw, grid: false,
  stepsPerFrame: chooseStep, dt: () => hStep,
});
