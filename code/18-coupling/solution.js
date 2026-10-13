// Chapter 18 - Coupling Two Solvers (reference solution)
// A stiff rod (an XPBD chain solved with many Gauss-Seidel iterations) is pressed into a soft membrane
// (a mattress of particles solved by its own few Gauss-Seidel sweeps). The two solvers only exchange
// contact forces and positions. Pick a coupling scheme and watch the contact force: some chatter, some
// settle onto the monolithic answer (the faint dashed membrane and the grey plot line).
// Both bodies move only vertically, so every equation stays one-dimensional and easy to read.

const mode = Sim.select('coupling', ['monolithic', 'explicit lag', 'fixed-point', 'relaxed fixed-point',
  'swapped roles', 'compliance: diagonal', 'compliance: Gershgorin'], 'compliance: Gershgorin', { reset: true });
const footSel = Sim.select('contacts under the foot', ['1', '2', '3', '5'], '3', { reset: true });
const passes = Sim.slider('coupling passes per step', 1, 20, 5, { step: 1 });
const relax = Sim.slider('relaxation theta', 0.01, 1, 0.05, { step: 0.01 });
const rodK = Sim.logSlider('rod link stiffness (N/m)', 2, 7, 1e5, { reset: true });
const pushCm = Sim.slider('push depth (cm)', 0, 8, 4, { step: 0.1 });
const memIters = Sim.slider('membrane GS sweeps per pass', 1, 60, 30, { step: 1 });

const h = 1 / 120;          // time step (s); two steps per animation frame

// Membrane: N nodes from x = -1 to 1, ends pinned at y = 0. Each node is tied to y = 0 by a soft
// spring kf (the mattress) and to its neighbours by a stiffer shear spring kc (the cover).
const N = 41, memX = new Float64Array(N);
const memMass = 0.01, kf = 20, kc = 1000;
for (let i = 0; i < N; i++) memX[i] = -1 + 2 * i / (N - 1);
const memA = memMass / (h * h) + kf + 2 * kc;   // diagonal of the membrane's implicit-Euler matrix

// Rod: nR particles in a vertical column, the top one driven by the "hand", the bottom one is the foot.
const nR = 9, rodSpacing = 0.06, rodMass = 0.05, rodLen = (nR - 1) * rodSpacing;
const monoIters = 300, rodIters = 300;

let world, ref, contacts, W, Wstatic, D, footS, footSstatic;

function makeWorld() {
  const w = {
    rodY: new Float64Array(nR), rodV: new Float64Array(nR), rodPred: new Float64Array(nR),
    rodLam: new Float64Array(nR - 1), rodStart: new Float64Array(nR),
    memY: new Float64Array(N), memV: new Float64Array(N), memPred: new Float64Array(N), memStart: new Float64Array(N),
    f: new Float64Array(contacts.length),       // contact forces (N), positive = foot pushes membrane down
    fOld: new Float64Array(contacts.length), residual: 0,
  };
  for (let i = 0; i < nR; i++) w.rodY[i] = 0.02 + rodLen - i * rodSpacing;   // foot starts 2 cm above
  return w;
}

function reset() {
  const nC = +footSel.value, mid = (N - 1) / 2;
  contacts = [];
  for (let j = 0; j < nC; j++) contacts.push(mid - Math.floor((nC - 1) / 2) + j);
  // Compliance (Green's function) of the membrane at the contact nodes: W[j][k] is how far node j moves
  // for a unit force on node k, during one implicit time step (Wstatic: for a slow, static push).
  W = membraneCompliance(memA);
  Wstatic = membraneCompliance(kf + 2 * kc);
  D = contactCompliance(W, mode.value);
  footS = rodFootCompliance(rodMass / (h * h));
  footSstatic = (nR - 1) / rodK.value;
  world = makeWorld();
  ref = makeWorld();
}

// ---------------------------------------------------------------------------------------------
// Compliance matrices
// ---------------------------------------------------------------------------------------------

// Solve the tridiagonal system with constant diagonal `diag`, off-diagonals -off, on n unknowns (Thomas).
function solveTridiag(diag, off, rhs, lastDiag) {
  const n = rhs.length, c = new Float64Array(n), d = new Float64Array(n), x = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const a = (i === n - 1 && lastDiag !== undefined) ? lastDiag : diag;
    const den = a - (i > 0 ? -off * c[i - 1] : 0);
    c[i] = -off / den;
    d[i] = (rhs[i] - (i > 0 ? -off * d[i - 1] : 0)) / den;
  }
  for (let i = n - 1; i >= 0; i--) x[i] = d[i] - (i < n - 1 ? c[i] * x[i + 1] : 0);
  return x;
}

// Push node k with a unit force and record how far every contact node moves: one column of W per contact.
function membraneCompliance(diag) {
  const nC = contacts.length, M = [];
  for (let j = 0; j < nC; j++) M.push(new Float64Array(nC));
  for (let k = 0; k < nC; k++) {
    const rhs = new Float64Array(N - 2);
    rhs[contacts[k] - 1] = 1;                       // interior unknowns are nodes 1..N-2
    const u = solveTridiag(diag, kc, rhs);
    for (let j = 0; j < nC; j++) M[j][k] = u[contacts[j] - 1];
  }
  return M;
}

// How far the foot moves for a unit force on it, with the rod top held (one step, inertia included).
function rodFootCompliance(massTerm) {
  const k = rodK.value, rhs = new Float64Array(nR - 1);
  rhs[nR - 2] = 1;
  return solveTridiag(massTerm + 2 * k, k, rhs, massTerm + k)[nR - 2];
}

// The compliance the stiff solver puts on each contact: how much the wall gives under that contact.
function contactCompliance(W, modeName) {
  const nC = W.length, D = new Float64Array(nC);
  if (modeName === 'compliance: diagonal') {
    for (let j = 0; j < nC; j++) D[j] = W[j][j];
  } else if (modeName === 'compliance: Gershgorin') {
    // Row sum of |W|: the contact also feels the membrane sinking under its neighbours' forces.
    for (let j = 0; j < nC; j++) {
      let s = 0;
      for (let k = 0; k < nC; k++) s += Math.abs(W[j][k]);
      D[j] = s;
    }
  }
  return D;
}

// ---------------------------------------------------------------------------------------------
// The two solvers
// ---------------------------------------------------------------------------------------------

// Membrane solver: Gauss-Seidel sweeps on its implicit-Euler system A y = b, warm-started from w.memY.
// forceDown[j] is the external force on contact node j. With footLevel set, the membrane instead treats
// the foot as a rigid ceiling (swapped roles) and returns the force the foot must have applied.
function solveMembrane(w, forceDown, sweeps, footLevel) {
  const y = w.memY, m = memMass / (h * h);
  const extra = new Float64Array(N), ceiling = new Float64Array(N).fill(Infinity);
  if (forceDown) for (let j = 0; j < contacts.length; j++) extra[contacts[j]] = -forceDown[j];
  if (footLevel !== undefined) for (let j = 0; j < contacts.length; j++) ceiling[contacts[j]] = footLevel;
  for (let s = 0; s < sweeps; s++) {
    for (let i = 1; i < N - 1; i++) {
      y[i] = (m * w.memPred[i] + extra[i] + kc * (y[i - 1] + y[i + 1])) / memA;
      if (y[i] > ceiling[i]) y[i] = ceiling[i];      // projected Gauss-Seidel: clamp before moving on
    }
  }
  if (footLevel !== undefined) {
    // Residual of the clamped rows = the force the foot exerted (positive = pushing down).
    const out = new Float64Array(contacts.length);
    for (let j = 0; j < contacts.length; j++) {
      const i = contacts[j];
      out[j] = Math.max(0, m * w.memPred[i] + kc * (y[i - 1] + y[i + 1]) - memA * y[i]);
    }
    return out;
  }
}

// Stiff solver: XPBD on the rod column, many Gauss-Seidel iterations, starting from the prediction.
// Each contact j treats membrane node contacts[j] as an obstacle with compliance D[j] (zero = rigid).
// fNow[j] is the force the membrane currently feels; footForceUp is used for swapped roles instead.
function solveRod(w, D, fNow, footForceUp) {
  const y = w.rodY, alphaRod = 1 / rodK.value / (h * h), wr = 1 / rodMass, nC = contacts.length;
  const swapped = footForceUp !== undefined;
  y.set(w.rodPred);
  if (swapped) y[nR - 1] += h * h * wr * footForceUp;       // a known external force goes into the prediction
  w.rodLam.fill(0);
  const lam = new Float64Array(nC);
  for (let it = 0; it < rodIters; it++) {
    for (let i = 0; i < nR - 1; i++) {
      // Vertical distance constraint C = (y_i - y_{i+1}) - spacing. Particle 0 (the top) is driven: w = 0.
      const wi = i === 0 ? 0 : wr, C = (y[i] - y[i + 1]) - rodSpacing;
      const dl = (-C - alphaRod * w.rodLam[i]) / (wi + wr + alphaRod);
      w.rodLam[i] += dl;
      y[i] += wi * dl;
      y[i + 1] -= wr * dl;
    }
    if (swapped) continue;
    for (let j = 0; j < nC; j++) {
      // The membrane gives way by D[j] metres per newton. Its unloaded level is memY + D f (where it
      // would sit with no contact force), and the contact gets XPBD compliance alpha = D[j].
      const level = w.memY[contacts[j]] + D[j] * fNow[j];
      const at = D[j] / (h * h);
      const C = y[nR - 1] - level;                 // foot must stay above the (compliant) membrane
      let dl = (-C - at * lam[j]) / (wr + at);
      if (lam[j] + dl < 0) dl = -lam[j];           // contacts only push
      lam[j] += dl;
      y[nR - 1] += wr * dl;
    }
  }
  const f = new Float64Array(nC);
  for (let j = 0; j < nC; j++) f[j] = lam[j] / (h * h);   // XPBD multiplier -> force
  return f;
}

// One loop that sees every equation: rod links, contacts (moving both sides) and membrane rows.
function solveMonolithic(w) {
  const y = w.rodY, alphaRod = 1 / rodK.value / (h * h), wr = 1 / rodMass, nC = contacts.length;
  const wm = 1 / (h * h * memA);                    // how far a membrane node moves per unit multiplier
  y.set(w.rodPred);
  w.rodLam.fill(0);
  const lam = new Float64Array(nC);
  for (let it = 0; it < monoIters; it++) {
    for (let i = 0; i < nR - 1; i++) {
      const wi = i === 0 ? 0 : wr, C = (y[i] - y[i + 1]) - rodSpacing;
      const dl = (-C - alphaRod * w.rodLam[i]) / (wi + wr + alphaRod);
      w.rodLam[i] += dl; y[i] += wi * dl; y[i + 1] -= wr * dl;
    }
    for (let j = 0; j < nC; j++) {
      const c = contacts[j], C = y[nR - 1] - w.memY[c];
      let dl = -C / (wr + wm);
      if (lam[j] + dl < 0) dl = -lam[j];
      lam[j] += dl; y[nR - 1] += wr * dl; w.memY[c] -= wm * dl;
    }
    for (let j = 0; j < nC; j++) w.f[j] = lam[j] / (h * h);
    solveMembrane(w, w.f, 1);
  }
}

// ---------------------------------------------------------------------------------------------
// Time step
// ---------------------------------------------------------------------------------------------

function stepWorld(w, scheme, topY) {
  // Predict both bodies (no gravity: the only load is the push), drive the rod top.
  w.rodStart.set(w.rodY); w.memStart.set(w.memY);
  for (let i = 0; i < nR; i++) w.rodPred[i] = w.rodY[i] + h * w.rodV[i];
  w.rodPred[0] = topY;
  for (let i = 0; i < N; i++) w.memPred[i] = (i === 0 || i === N - 1) ? 0 : w.memY[i] + h * w.memV[i];
  w.residual = 0;

  if (scheme === 'monolithic') {
    solveMonolithic(w);
  } else if (scheme === 'swapped roles') {
    // The soft solver handles the contact (rod foot = rigid ceiling); the stiff rod receives the force.
    for (let k = 0; k < passes.value; k++) {
      w.fOld.set(w.f);
      const fm = solveMembrane(w, null, memIters.value, w.rodY[nR - 1]);
      let total = 0;
      for (let j = 0; j < fm.length; j++) { w.f[j] = fm[j]; total += fm[j]; }
      solveRod(w, null, null, total);
    }
  } else {
    // Rod first (it receives positions), then the membrane (it receives forces).
    const K = scheme === 'explicit lag' ? 1 : passes.value;
    const Duse = scheme.startsWith('compliance') ? D : new Float64Array(contacts.length);
    for (let k = 0; k < K; k++) {
      w.fOld.set(w.f);
      const fRod = solveRod(w, Duse, w.f);
      const th = scheme === 'relaxed fixed-point' ? relax.value : 1;
      for (let j = 0; j < fRod.length; j++) w.f[j] = (1 - th) * w.f[j] + th * fRod[j];
      solveMembrane(w, w.f, memIters.value);
    }
  }
  for (let j = 0; j < contacts.length; j++) w.residual += Math.abs(w.f[j] - w.fOld[j]);

  for (let i = 0; i < nR; i++) w.rodV[i] = (w.rodY[i] - w.rodStart[i]) / h;
  for (let i = 0; i < N; i++) w.memV[i] = (w.memY[i] - w.memStart[i]) / h;
}

let topY = 0;
function step() {
  // The hand lowers the rod top towards its target at 5 cm/s.
  const target = rodLen - pushCm.value / 100;
  const cur = world.rodY[0];
  topY = cur + Math.max(-0.05 * h, Math.min(0.05 * h, target - cur));
  stepWorld(world, mode.value, topY);
  stepWorld(ref, 'monolithic', topY);
}

// ---------------------------------------------------------------------------------------------
// Readouts: predicted convergence factor per coupling pass
// ---------------------------------------------------------------------------------------------

// Spectral radius of (S + D)^-1 (D - W) with S = s * ones: Cholesky of P = S + D, then Jacobi
// eigenvalues of the symmetric L^-1 (D - W) L^-T.
function couplingRadius(Wm, Dv, s) {
  const n = Wm.length, P = [], Nm = [];
  for (let i = 0; i < n; i++) {
    P.push(new Float64Array(n)); Nm.push(new Float64Array(n));
    for (let j = 0; j < n; j++) { P[i][j] = s + (i === j ? Dv[i] : 0); Nm[i][j] = (i === j ? Dv[i] : 0) - Wm[i][j]; }
  }
  const L = [];
  for (let i = 0; i < n; i++) {
    L.push(new Float64Array(n));
    for (let j = 0; j <= i; j++) {
      let sum = P[i][j];
      for (let k = 0; k < j; k++) sum -= L[i][k] * L[j][k];
      L[i][j] = i === j ? Math.sqrt(sum) : sum / L[j][j];
    }
  }
  const solveL = (b) => { const x = new Float64Array(n); for (let i = 0; i < n; i++) { let s2 = b[i]; for (let k = 0; k < i; k++) s2 -= L[i][k] * x[k]; x[i] = s2 / L[i][i]; } return x; };
  const T = [];
  for (let j = 0; j < n; j++) T.push(solveL(Nm.map((row) => row[j])));        // columns of L^-1 N
  const S = [];
  for (let i = 0; i < n; i++) S.push(solveL(T.map((col) => col[i])));          // rows of L^-1 N L^-T
  return Math.max(...symmetricEigenvalues(S).map(Math.abs));
}

function symmetricEigenvalues(A) {
  const n = A.length, a = A.map((r) => Float64Array.from(r));
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] * a[p][q];
    if (off < 1e-24) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-30) continue;
      const th = 0.5 * Math.atan2(2 * a[p][q], a[q][q] - a[p][p]), c = Math.cos(th), s = Math.sin(th);
      for (let k = 0; k < n; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
      for (let k = 0; k < n; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
    }
  }
  return a.map((r, i) => r[i]);
}

// Compliance of the membrane as felt by a rigid flat foot: 1 / (1^T W^-1 1).
function footCompliance(Wm) {
  const n = Wm.length, ones = new Float64Array(n).fill(1);
  const A = Wm.map((r) => Float64Array.from(r)), b = Float64Array.from(ones);
  for (let i = 0; i < n; i++) for (let r = i + 1; r < n; r++) {
    const f = A[r][i] / A[i][i];
    for (let c = i; c < n; c++) A[r][c] -= f * A[i][c];
    b[r] -= f * b[i];
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let c = i + 1; c < n; c++) s -= A[i][c] * x[c]; x[i] = s / A[i][i]; }
  return 1 / x.reduce((p, q) => p + q, 0);
}

// ---------------------------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------------------------

function draw(g) {
  const C = Sim.colors, footHalf = 0.5 * (contacts.length - 1) * (2 / (N - 1)) + 0.02;
  // Membrane as a filled mattress, plus the monolithic reference as a dashed line.
  const pts = [];
  for (let i = 0; i < N; i++) pts.push({ x: memX[i], y: Math.max(-0.4, Math.min(0.5, world.memY[i])) });
  g.polyline(pts.concat([{ x: 1, y: -5 }, { x: -1, y: -5 }]), { fill: '#e8dcc4', closed: true });
  g.polyline(pts, { color: C.orange, width: 2.5 });
  const refPts = [];
  for (let i = 0; i < N; i++) refPts.push({ x: memX[i], y: ref.memY[i] });
  g.polyline(refPts, { color: C.ink, width: 1, dash: [4, 4], alpha: 0.6 });
  g.points(memX, world.memY.map((v) => Math.max(-0.4, Math.min(0.5, v))), N, { radius: 2, color: C.orange });

  // Rod column and its flat foot.
  const clampY = (v) => Math.max(-0.4, Math.min(0.6, v));
  const top = world.rodY[0], foot = clampY(world.rodY[nR - 1]);
  g.line(0, top, 0, foot, { color: C.blue, width: 6 });
  for (let i = 0; i < nR; i++) g.dot(0, clampY(world.rodY[i]), 3, C.ink);
  g.rect(-footHalf, foot, 2 * footHalf, 0.025, { fill: C.blue });
  g.rect(-0.06, top, 0.12, 0.04, { fill: C.gray });
  g.text('hand', 0.08, top + 0.02, { color: C.muted, size: 12 });

  // Contact forces as arrows (1 cm per newton, capped).
  let total = 0, totalRef = 0;
  for (let j = 0; j < contacts.length; j++) {
    const x = memX[contacts[j]], len = Math.min(0.3, 0.01 * world.f[j]);
    if (len > 0.002) g.arrow(x, foot, 0, -len, { color: C.red, width: 2 });
    total += world.f[j]; totalRef += ref.f[j];
  }

  const Fstatic = (pushCm.value / 100) / (footSstatic + footCompliance(Wstatic));
  const cap = Math.max(4 * Fstatic, 1);
  Sim.plot('foot force (N)', Math.min(cap, total), { color: C.red, min: 0 });
  Sim.plot('monolithic', Math.min(cap, totalRef), { color: C.gray });
  Sim.show('foot force', total.toFixed(2) + ' N');
  Sim.show('monolithic', totalRef.toFixed(2) + ' N');
  Sim.show('static prediction', Fstatic.toFixed(2) + ' N');
  Sim.show('W_jj (one step)', (W[0][0] * 1000).toFixed(3) + ' mm/N');
  Sim.show('interface change, last pass', world.residual.toExponential(1) + ' N');

  // Predicted factor by which the force error shrinks (or grows) each coupling pass.
  const r = footCompliance(W) / footS, m = mode.value;
  let rho = null;
  if (m === 'explicit lag' || m === 'fixed-point') rho = r;
  else if (m === 'relaxed fixed-point') rho = Math.abs(1 - relax.value * (1 + r));
  else if (m === 'swapped roles') rho = 1 / r;
  else if (m.startsWith('compliance')) rho = D.every((v) => v > 0) ? couplingRadius(W, D, footS) : r;
  if (rho !== null) Sim.show('predicted factor per pass', rho.toFixed(3) + (rho < 1 ? ' (converges)' : ' (diverges)'));
}

Sim.run({ title: 'A stiff rod pressed into a soft membrane', view: [-0.6, -0.25, 0.6, 0.58], reset, step, draw,
  dt: h, stepsPerFrame: 2 });
