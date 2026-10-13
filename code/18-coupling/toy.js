// Chapter 18 - Coupling Two Solvers: the two-spring toy
// A hand pushes a "rod" spring k1 down onto a "membrane" spring k2. Two separate solvers, one per spring,
// agree on the interface force f by exchanging values. The exact answer is f* = k1 k2 d / (k1 + k2).
// The contact is glued (it can pull as well as push), so the algebra stays linear and errors can grow
// without bound. Watch the plot: the red force either settles on the grey exact curve or runs away.

const scheme = Sim.select('coupling', ['explicit lag', 'fixed-point', 'relaxed fixed-point', 'compliance-augmented',
  'swapped roles'], 'explicit lag', { reset: true });
const k1 = Sim.logSlider('rod stiffness k1 (N/m)', 1, 4, 50, { reset: true });
const k2 = Sim.logSlider('membrane stiffness k2 (N/m)', 1, 4, 100, { reset: true });
const passes = Sim.slider('passes per step', 1, 20, 1, { step: 1 });
const theta = Sim.slider('relaxation theta', 0.01, 1, 0.5, { step: 0.01 });
const wScale = Sim.slider('compliance estimate / true', 0, 3, 1, { step: 0.05 });

const dMax = 0.1;                 // the hand pushes between 0 and 10 cm
let s;                            // interface state shared by the two solvers

function reset() {
  s = { t: 0, f: 0, x: 0, y: 0, d: 0, ratio: NaN };
}

// Solver 1, the rod: the hand holds its top at depth d. Given where the membrane surface is (y),
// pinning the tip there makes the rod compress by d - y, so it pushes with k1 (d - y).
function rodGivenPosition(d, y) { return k1.value * (d - y); }
// The rod given a force instead: its tip sits at d - f / k1.
function rodGivenForce(d, f) { return d - f / k1.value; }
// Solver 2, the membrane: given a force it sinks by f / k2; given a position it pushes back with k2 y.
function membraneGivenForce(f) { return f / k2.value; }
function membraneGivenPosition(x) { return k2.value * x; }

function step(dt) {
  s.t += dt;
  s.d = dMax * 0.5 * (1 - Math.cos(s.t));                       // slow push in and out, period 2 pi s
  const d = s.d, fExact = k1.value * k2.value * d / (k1.value + k2.value);
  const K = scheme.value === 'explicit lag' ? 1 : passes.value;
  const errs = [s.f - fExact];                    // error of the incoming force, against today's target
  for (let k = 0; k < K; k++) {
    const m = scheme.value;
    if (m === 'swapped roles') {
      s.x = rodGivenForce(d, s.f);                // the stiff rod receives the force ...
      s.f = membraneGivenPosition(s.x);           // ... the soft membrane receives the position
      s.y = s.x;
    } else {
      let fNew;
      if (m === 'compliance-augmented') {
        // The rod solver models the membrane as y(f) = y + What (f - f_old) and solves for f in one go:
        // f = k1 (d - y - What (f - f_old))  =>  f = k1 (d - y + What f_old) / (1 + k1 What).
        const What = wScale.value / k2.value;
        fNew = k1.value * (d - s.y + What * s.f) / (1 + k1.value * What);
      } else {
        fNew = rodGivenPosition(d, s.y);
      }
      const th = m === 'relaxed fixed-point' ? theta.value : 1;
      s.f = (1 - th) * s.f + th * fNew;
      s.y = membraneGivenForce(s.f);
      s.x = rodGivenForce(d, s.f);
    }
    errs.push(s.f - fExact);
  }
  // Measured factor: error after the last pass over the error one pass earlier.
  const prev = errs[errs.length - 2], e = errs[errs.length - 1];
  s.ratio = Math.abs(prev) > 1e-9 ? e / prev : NaN;
}

function predictedFactor() {
  const r = k1.value / k2.value, m = scheme.value;
  if (m === 'explicit lag' || m === 'fixed-point') return -r;
  if (m === 'relaxed fixed-point') return 1 - theta.value * (1 + r);
  if (m === 'compliance-augmented') {
    const What = wScale.value / k2.value, Wt = 1 / k2.value;
    return k1.value * (What - Wt) / (1 + k1.value * What);
  }
  return -1 / r;
}

function drawSpring(g, x, yTop, yBot, color, coils, amp, width) {
  const pts = [{ x, y: yTop }];
  for (let i = 0; i <= coils; i++) {
    const t = i / coils, yy = yTop - 0.03 - (yTop - yBot - 0.06) * t;
    pts.push({ x: x + (i === 0 || i === coils ? 0 : (i % 2 ? amp : -amp)), y: yy });
  }
  pts.push({ x, y: yBot });
  g.polyline(pts, { color, width });
}

function draw(g) {
  const C = Sim.colors, show = (v) => Math.max(-0.55, Math.min(0.5, v));
  const handY = 0.5 - s.d, tipY = show(-s.x), memY = show(-s.y);
  // Membrane spring from the ground up to its surface, rod spring from the hand down to its tip.
  g.ground(-0.6);
  drawSpring(g, 0, memY, -0.6, C.orange, 10, 0.08, 2);
  g.rect(-0.25, memY - 0.02, 0.5, 0.02, { fill: C.orange });
  drawSpring(g, 0, handY, tipY, C.blue, 14, 0.04, 2);
  g.rect(-0.08, tipY, 0.16, 0.02, { fill: C.blue });
  g.rect(-0.12, handY, 0.24, 0.04, { fill: C.gray });
  g.text('hand', 0.16, handY + 0.02, { color: C.muted, size: 12 });
  g.text('rod k1', 0.12, 0.5 * (handY + tipY), { color: C.blue });
  g.text('membrane k2', 0.32, -0.3, { color: C.orange });
  if (Math.abs(tipY - memY) > 0.004) g.text('mismatch', -0.55, 0.5 * (tipY + memY), { color: C.red, size: 12 });

  const fExact = k1.value * k2.value * s.d / (k1.value + k2.value);
  const fCap = 3 * k1.value * k2.value * dMax / (k1.value + k2.value);
  Sim.plot('interface force f (N)', Math.max(-fCap, Math.min(fCap, s.f)), { color: C.red });
  Sim.plot('exact', fExact, { color: C.gray });
  Sim.show('k1 / k2', (k1.value / k2.value).toFixed(3));
  Sim.show('f, exact f*', s.f.toFixed(3) + ', ' + fExact.toFixed(3) + ' N');
  Sim.show('predicted factor per pass', predictedFactor().toFixed(3));
  Sim.show('measured factor', isFinite(s.ratio) ? s.ratio.toFixed(3) : '-');
  if (!(Math.abs(s.f) < 1e6)) Sim.show('status', 'diverged: press R');
}

Sim.run({ title: 'Two springs, two solvers', view: [-0.9, -0.65, 0.9, 0.62], reset, step, draw, dt: 1 / 60, grid: false });
