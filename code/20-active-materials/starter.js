// Chapter 20 - Muscles and Travelling Waves (starter)
// Three scenes built from one ingredient, a muscle constraint whose rest length follows an activation a(t):
//  - peristaltic tube: ring muscles squeeze a soft tube in a travelling wave and push the contents along;
//  - crawling worm: segment muscles plus direction-dependent ground friction make a worm crawl;
//  - muscle under load: one muscle lifting a weight, driven by active strain or by active stress.
// Everything is XPBD with many small substeps (one constraint pass per substep).
// Your job: the four TODOs below (the muscle constraint, the travelling wave, direction-dependent friction
// and the active-stress force). Until then the muscles do nothing.

const scene = Sim.select('scene', ['peristaltic tube', 'crawling worm', 'muscle under load'], 'peristaltic tube',
  { onChange: onSceneChange });
// Wave controls (tube and worm).
const waveSpeed = Sim.slider('wave speed c (m/s)', 0.05, 1, 0.2, { step: 0.01 });
const wavelength = Sim.slider('wavelength (m)', 0.2, 1.2, 0.6, { step: 0.05 });
const strength = Sim.slider('strength (max strain)', 0, 0.9, 0.8, { step: 0.01 });
// Tube.
const contents = Sim.select('tube contents', ['ball', 'grains'], 'ball', { reset: true });
// Worm.
const muBack = Sim.slider('friction backward', 0, 2, 1.0, { step: 0.05 });
const muFwd = Sim.slider('friction forward', 0, 2, 0.1, { step: 0.05 });
// Muscle under load.
const law = Sim.select('activation law', ['active strain', 'active stress'], 'active strain', { reset: true });
const loadMass = Sim.logSlider('load mass (kg)', -2, 1, 0.5);
const forceCap = Sim.logSlider('muscle force cap (N)', 0, 3, 20);
const tauMs = Sim.slider('activation time constant (ms)', 0, 200, 50, { step: 5 });

const groups = {
  'peristaltic tube': ['wave speed c (m/s)', 'wavelength (m)', 'strength (max strain)', 'tube contents'],
  'crawling worm': ['wave speed c (m/s)', 'wavelength (m)', 'strength (max strain)', 'friction backward', 'friction forward'],
  'muscle under load': ['activation law', 'load mass (kg)', 'muscle force cap (N)', 'activation time constant (ms)'],
};
// Show only the controls that matter for the current scene (plain DOM, nothing physical here).
function showControlsFor(name) {
  const all = new Set([].concat(...Object.values(groups)));
  for (const row of document.querySelectorAll('.sim-control')) {
    const label = row.querySelector('.sim-control-name');
    if (!label || label.textContent === 'scene') continue;
    if (all.has(label.textContent)) row.style.display = groups[name].includes(label.textContent) ? '' : 'none';
  }
}
showControlsFor(scene.value);

function onSceneChange() {
  // Each scene starts with its own wavelength: the worm's body (0.5 m) holds exactly one.
  if (scene.value === 'crawling worm') wavelength.set(0.5);
  if (scene.value === 'peristaltic tube') wavelength.set(0.6);
  showControlsFor(scene.value);
  Sim.reset();
}

const substeps = 20, frameDt = 1 / 60, h = frameDt / substeps;
const g = 9.81;

// ---------------------------------------------------------------------------------------------
// Particles and constraints (flat arrays for particles, small objects for constraints)
// ---------------------------------------------------------------------------------------------
let P, dists, tethers, muscles, time, stats;

function makeParticles(n) {
  return { n: 0, x: new Float64Array(n), y: new Float64Array(n), px: new Float64Array(n), py: new Float64Array(n),
    vx: new Float64Array(n), vy: new Float64Array(n), w: new Float64Array(n), r: new Float64Array(n) };
}
function addParticle(x, y, mass, radius) {
  const i = P.n++;
  P.x[i] = x; P.y[i] = y; P.w[i] = mass > 0 ? 1 / mass : 0; P.r[i] = radius;
  return i;
}
function addDistance(a, b, alpha) {
  dists.push({ a, b, L: Math.hypot(P.x[b] - P.x[a], P.y[b] - P.y[a]), alpha });
}
function addMuscle(a, b, strainMax, alpha, cap, tensionOnly, s) {
  const L0 = Math.hypot(P.x[b] - P.x[a], P.y[b] - P.y[a]);
  muscles.push({ a, b, L0, strainMax, alpha, cap, tensionOnly, s, act: 0, lam: 0, force: 0 });
}

// ---------------------------------------------------------------------------------------------
// The muscle constraint. Active strain: the activation a in [0, 1] shortens the rest length,
// L = L0 (1 - strainMax a). Tension-only muscles pull but never push. The force is capped at m.cap.
// ---------------------------------------------------------------------------------------------
function solveMuscle(m) {
  const a = m.a, b = m.b;
  const dx = P.x[b] - P.x[a], dy = P.y[b] - P.y[a], len = Math.hypot(dx, dy);
  if (len < 1e-9) return;
  const nx = dx / len, ny = dy / len;
  // TODO (exercise 1): the muscle constraint.
  //  1. Active strain: rest = m.L0 * (1 - m.strainMax * m.act).
  //  2. C = len - rest. A tension-only muscle (m.tensionOnly) does nothing while C <= 0 (set m.force = 0).
  //  3. XPBD: dl = (-C - at * m.lam) / (wa + wb + at), with at = m.alpha / h^2.
  //  4. Cap the force: lam = clamp(m.lam + dl, -m.cap h^2, m.tensionOnly ? 0 : m.cap h^2).
  //  5. Apply d = lam - m.lam: a moves by -wa n d, b by +wb n d. Store m.lam = lam, m.force = -lam / h^2.
  m.force = 0;
}

// Active stress instead: no rest-length change, the muscle adds a contractile force along its fibre.
// Linear force-length curve: full force a*Fmax at the rest length L0, none at the shortest length.
function activeStressForce(m) {
  // TODO (exercise 4): return act * cap * fl, where fl falls linearly from 1 at len = L0 to 0 at the
  // shortest length L0 (1 - strainMax), clamped to [0, 1].
  return 0;
}

// ---------------------------------------------------------------------------------------------
// Travelling wave of activation a(s, t) = f(s - c t): one smooth contraction per wavelength.
// ---------------------------------------------------------------------------------------------
function waveActivation(s, t) {
  // TODO (exercise 2): a travelling wave a = f(s - c t). Compute the phase of s - c t within one
  // wavelength (a number in [0, 1); careful, % can return negative numbers), then return a bump that is 0 at
  // phase 0 and 1 at phase 0.5: bump = 0.5 (1 - cos(2 pi phase)), cubed to make it narrower.
  return 0;
}

// ---------------------------------------------------------------------------------------------
// Ground contact with direction-dependent friction (worm). dn is how far the contact pushed the
// particle out of the ground this substep; Coulomb friction may cancel at most mu * dn of sliding.
// ---------------------------------------------------------------------------------------------
function groundFriction(i, dn) {
  const slide = P.x[i] - P.px[i];                       // +x is forward (towards the head)
  // TODO (exercise 3): make the friction depend on the direction of sliding (forward = slide > 0).
  const mu = muBack.value;
  const maxCancel = mu * dn;
  if (Math.abs(slide) <= maxCancel) P.x[i] = P.px[i];
  else P.x[i] -= Math.sign(slide) * maxCancel;
}

// ---------------------------------------------------------------------------------------------
// Scenes
// ---------------------------------------------------------------------------------------------
const tube = { L: 1.6, R: 0.1, n: 41, wallR: 0.015, ballR: 0.06 };
const worm = { nSeg: 10, L0: 0.05, H0: 0.04, nodeR: 0.01 };
const hang = { L0: 0.3, top: 0.55 };

function reset() {
  P = makeParticles(400); dists = []; tethers = []; muscles = []; time = 0;
  stats = { content: [], wormStart: 0, track: [], maxSpeed: 0, energy: 0 };
  if (scene.value === 'peristaltic tube') buildTube();
  else if (scene.value === 'crawling worm') buildWorm();
  else buildHang();
  for (let i = 0; i < P.n; i++) { P.px[i] = P.x[i]; P.py[i] = P.y[i]; }
}

function buildTube() {
  const { L, R, n, wallR } = tube, ds = L / (n - 1);
  tube.top = []; tube.bot = [];
  for (let i = 0; i < n; i++) {
    tube.top.push(addParticle(i * ds, R, 0.02, wallR));
    tube.bot.push(addParticle(i * ds, -R, 0.02, wallR));
  }
  for (const wall of [tube.top, tube.bot]) {
    for (let i = 0; i < n; i++) {
      // The elastic tube and its surroundings hold each wall node near its rest place (50 N/m).
      tethers.push({ i: wall[i], x0: P.x[wall[i]], y0: P.y[wall[i]], alpha: 1 / 50 });
      if (i > 0) addDistance(wall[i - 1], wall[i], 1e-4);
    }
  }
  // Ring muscles: one across the tube at every node pair, tension-only, rest length = the diameter.
  for (let i = 0; i < n; i++) addMuscle(tube.top[i], tube.bot[i], 0.9, 1e-5, 10, true, i * ds);
  // Start the contents just ahead of the first contraction (at s = wavelength / 2 when t = 0).
  const rnd = Sim.random(7), x0 = 0.5 * wavelength.value + 0.12;
  if (contents.value === 'ball') {
    stats.content.push(addParticle(x0, 0, 0.05, tube.ballR));
  } else {
    // A 6 x 4 block of grains, placed without overlaps (overlapping starts explode, see chapter 32).
    for (let k = 0; k < 24; k++) {
      const gx = x0 + 0.04 * (k % 6) + 0.004 * rnd(), gy = -0.06 + 0.04 * Math.floor(k / 6);
      stats.content.push(addParticle(gx, gy, 0.002, 0.018));
    }
  }
  stats.track = [{ t: 0, x: contentX() }];
  stats.offset = 0;
}

function buildWorm() {
  const { nSeg, L0, nodeR } = worm;
  worm.nodes = [];
  // Node 0 is the head (front, +x); segment k joins node k and node k+1.
  for (let k = 0; k <= nSeg; k++) worm.nodes.push(addParticle(0.4 - k * L0, nodeR, 0.01, nodeR));
  for (let k = 0; k < nSeg; k++) {
    // A segment is a fluid-filled compartment: the longitudinal muscle shortens it, the circular muscle
    // squeezes it long again. Together they set its length, so model the pair as one two-sided muscle.
    addMuscle(worm.nodes[k], worm.nodes[k + 1], 0.5, 1e-6, 5, false, (k + 0.5) * L0);
  }
  stats.wormStart = P.x[worm.nodes[0]];
}

function buildHang() {
  const top = addParticle(0, hang.top, 0, 0.01);
  const load = addParticle(0, hang.top - hang.L0, loadMass.value, 0.03);
  hang.topI = top; hang.loadI = load;
  addMuscle(top, load, 0.3, 1e-4, forceCap.value, true, 0);
}

// What drives each muscle in each scene (the command u; activation follows it).
function command(m, t) {
  if (scene.value === 'muscle under load') return (t % 4) < 2 ? 1 : 0;    // on for 2 s, off for 2 s
  return Math.min(1, t / 0.5) * waveActivation(m.s, t);                 // fade the wave in over 0.5 s
}

function contentX() {
  let s = 0;
  for (const i of stats.content) s += P.x[i];
  return s / stats.content.length;
}

// ---------------------------------------------------------------------------------------------
// Time step
// ---------------------------------------------------------------------------------------------
function step() {
  const isTube = scene.value === 'peristaltic tube', isWorm = scene.value === 'crawling worm';
  const isHang = scene.value === 'muscle under load';
  if (isHang) {
    // Sliders act immediately on the hanging muscle.
    P.w[hang.loadI] = 1 / loadMass.value;
    muscles[0].cap = forceCap.value;
  }
  for (let sub = 0; sub < substeps; sub++) {
    time += h;
    // Activation dynamics: a relaxes towards the command with time constant tau (0 = instant).
    for (const m of muscles) {
      const u = command(m, time);
      if (isHang) {
        const tau = tauMs.value / 1000;
        m.act = tau > 0 ? m.act + (u - m.act) * (1 - Math.exp(-h / tau)) : u;
      } else {
        m.act = u;
        if (isTube || isWorm) m.strainMax = isTube ? Math.min(0.9, strength.value) : 0.5 * strength.value;
      }
      m.lam = 0;
    }
    // Predict.
    for (let i = 0; i < P.n; i++) {
      if (P.w[i] === 0) continue;
      if (!isTube) P.vy[i] -= g * h;
      P.px[i] = P.x[i]; P.py[i] = P.y[i];
    }
    if (isHang && law.value === 'active stress') {
      for (const m of muscles) {
        const F = activeStressForce(m);
        const dx = P.x[m.b] - P.x[m.a], dy = P.y[m.b] - P.y[m.a], len = Math.hypot(dx, dy) || 1;
        P.vx[m.a] += h * P.w[m.a] * F * dx / len; P.vy[m.a] += h * P.w[m.a] * F * dy / len;
        P.vx[m.b] -= h * P.w[m.b] * F * dx / len; P.vy[m.b] -= h * P.w[m.b] * F * dy / len;
      }
    }
    for (let i = 0; i < P.n; i++) {
      if (P.w[i] === 0) continue;
      P.x[i] += h * P.vx[i]; P.y[i] += h * P.vy[i];
    }
    // Constraints, one pass.
    for (const t of tethers) solveTether(t);
    for (const d of dists) solveDistance(d);
    for (const m of muscles) {
      if (isHang && law.value === 'active stress') { const a = m.act; m.act = 0; solveMuscle(m); m.act = a; }
      else solveMuscle(m);
    }
    if (isTube) solveTubeContacts();
    if (isHang) {
      // The floor under the weight (the weight hangs from its top face; the box is 2 * side tall).
      const i = hang.loadI, bottom = 2 * loadSide();
      if (P.y[i] < bottom) P.y[i] = bottom;
    }
    if (isWorm) {
      for (const i of worm.nodes) {
        const dn = P.r[i] - P.y[i];
        if (dn > 0) { P.y[i] += dn; groundFriction(i, dn); }
      }
    }
    // Velocities.
    for (let i = 0; i < P.n; i++) {
      if (P.w[i] === 0) { P.vx[i] = 0; P.vy[i] = 0; continue; }
      P.vx[i] = (P.x[i] - P.px[i]) / h; P.vy[i] = (P.y[i] - P.py[i]) / h;
    }
    if (isTube) {
      for (const i of stats.content) { const k = Math.exp(-2 * h); P.vx[i] *= k; P.vy[i] *= k; }  // viscous fluid
    }
    if (isHang) { const k = Math.exp(-10 * h); P.vx[hang.loadI] *= k; P.vy[hang.loadI] *= k; }     // viscous drag
  }
  if (isTube) wrapContents();
  if (isHang) {
    const v = Math.hypot(P.vx[hang.loadI], P.vy[hang.loadI]);
    stats.maxSpeed = Math.max(stats.maxSpeed, v);
  }
}

function solveTether(t) {
  const i = t.i, at = t.alpha / (h * h), wi = P.w[i];
  for (const axis of ['x', 'y']) {
    const C = P[axis][i] - (axis === 'x' ? t.x0 : t.y0);
    P[axis][i] += wi * (-C) / (wi + at);
  }
}

function solveDistance(d) {
  const a = d.a, b = d.b, dx = P.x[b] - P.x[a], dy = P.y[b] - P.y[a], len = Math.hypot(dx, dy);
  if (len < 1e-9) return;
  const wa = P.w[a], wb = P.w[b], at = d.alpha / (h * h), C = len - d.L;
  const dl = -C / (wa + wb + at), nx = dx / len, ny = dy / len;
  P.x[a] -= wa * nx * dl; P.y[a] -= wa * ny * dl;
  P.x[b] += wb * nx * dl; P.y[b] += wb * ny * dl;
}

function pushApart(a, b, minDist) {
  const dx = P.x[b] - P.x[a], dy = P.y[b] - P.y[a], d2 = dx * dx + dy * dy;
  if (d2 >= minDist * minDist || d2 < 1e-18) return;
  const d = Math.sqrt(d2), wa = P.w[a], wb = P.w[b];
  if (wa + wb === 0) return;
  const C = d - minDist, nx = dx / d, ny = dy / d, s = -C / (wa + wb);
  P.x[a] -= wa * nx * s; P.y[a] -= wa * ny * s;
  P.x[b] += wb * nx * s; P.y[b] += wb * ny * s;
}

// Keep particle i at least minDist from the wall segment (a, b), sharing the correction between the
// particle and the segment's two nodes by their barycentric weights.
function pushFromSegment(i, a, b, minDist) {
  const ex = P.x[b] - P.x[a], ey = P.y[b] - P.y[a], e2 = ex * ex + ey * ey;
  let t = e2 > 0 ? ((P.x[i] - P.x[a]) * ex + (P.y[i] - P.y[a]) * ey) / e2 : 0;
  t = Math.max(0, Math.min(1, t));
  const qx = P.x[a] + t * ex, qy = P.y[a] + t * ey;
  const dx = P.x[i] - qx, dy = P.y[i] - qy, d = Math.hypot(dx, dy);
  if (d >= minDist || d < 1e-12) return;
  const nx = dx / d, ny = dy / d, C = d - minDist;
  const wa = P.w[a] * (1 - t) * (1 - t), wb = P.w[b] * t * t, wi = P.w[i];
  const s = -C / (wi + wa + wb);
  P.x[i] += wi * nx * s; P.y[i] += wi * ny * s;
  P.x[a] -= P.w[a] * (1 - t) * nx * s; P.y[a] -= P.w[a] * (1 - t) * ny * s;
  P.x[b] -= P.w[b] * t * nx * s; P.y[b] -= P.w[b] * t * ny * s;
}

function solveTubeContacts() {
  const c = stats.content;
  for (const i of c) {
    // Only the few wall segments near the particle can touch it (wall nodes stay close to their rest x).
    const ds = tube.L / (tube.n - 1), k0 = Math.floor(P.x[i] / ds);
    for (const wall of [tube.top, tube.bot]) {
      for (let k = Math.max(0, k0 - 3); k <= Math.min(tube.n - 2, k0 + 3); k++) {
        pushFromSegment(i, wall[k], wall[k + 1], P.r[i] + tube.wallR);
      }
    }
  }
  for (let p = 0; p < c.length; p++) for (let q = p + 1; q < c.length; q++) pushApart(c[p], c[q], P.r[c[p]] + P.r[c[q]]);
}

// Contents near the right end jump back by a whole number of wavelengths, to the same phase of the wave,
// so the transport can be watched (and measured) for ever.
function wrapContents() {
  const lam = wavelength.value, shift = lam * Math.max(1, Math.floor((tube.L - 0.4) / lam));
  for (const i of stats.content) {
    if (P.x[i] > tube.L - 0.15) {
      P.x[i] -= shift; P.px[i] -= shift;
      stats.offset += shift / stats.content.length;
    } else if (P.x[i] < 0.05) {
      P.x[i] += shift; P.px[i] += shift;
      stats.offset -= shift / stats.content.length;
    }
  }
  stats.track.push({ t: time, x: contentX() + stats.offset });
  while (stats.track.length > 2 && stats.track[0].t < time - 3) stats.track.shift();
}

// ---------------------------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------------------------
function draw(gr) {
  const C = Sim.colors;
  if (scene.value === 'peristaltic tube') drawTube(gr, C);
  else if (scene.value === 'crawling worm') drawWorm(gr, C);
  else drawHang(gr, C);
}

function drawTube(gr, C) {
  Sim.setView([-0.1, -0.45, 1.7, 0.4]);
  // Muscles, shaded by activation.
  for (const m of muscles) {
    if (m.act > 0.02) gr.line(P.x[m.a], P.y[m.a], P.x[m.b], P.y[m.b], { color: C.red, width: 3, alpha: Math.min(1, m.act) });
  }
  for (const wall of [tube.top, tube.bot]) {
    const pts = wall.map((i) => ({ x: P.x[i], y: P.y[i] }));
    gr.polyline(pts, { color: '#b9826b', width: 7 });
  }
  for (const i of stats.content) gr.circle(P.x[i], P.y[i], P.r[i], { fill: contents.value === 'ball' ? C.blue : C.teal });
  // The activation profile a(s) along the tube, drawn below it.
  const prof = [];
  for (const m of muscles) prof.push({ x: m.s, y: -0.4 + 0.2 * m.act });
  gr.line(0, -0.4, tube.L, -0.4, { color: C.muted, width: 1 });
  gr.polyline(prof, { color: C.red, width: 2 });
  gr.text('activation a(s, t)', 0, -0.17, { color: C.red, size: 12 });
  gr.arrow(1.2, -0.17, 0.2, 0, { color: C.red, width: 1.5 });
  gr.text('c', 1.42, -0.17, { color: C.red, size: 12 });

  const tr = stats.track;
  if (tr.length > 1 && tr[tr.length - 1].t - tr[0].t > 1) {
    const v = (tr[tr.length - 1].x - tr[0].x) / (tr[tr.length - 1].t - tr[0].t);
    Sim.show('contents speed', v.toFixed(3) + ' m/s');
    Sim.show('contents speed / wave speed', (v / waveSpeed.value).toFixed(2));
    Sim.plot('speed / c', v / waveSpeed.value, { color: C.blue, min: 0, max: 1.2 });
  }
  let maxF = 0;
  for (const m of muscles) maxF = Math.max(maxF, m.force);
  Sim.show('wave period', (wavelength.value / waveSpeed.value).toFixed(2) + ' s');
  Sim.show('largest ring force', maxF.toFixed(2) + ' N');
}

// Quasi-static crawling speed for this wave shape (see the chapter): parts of the body slide forward where
// the activation is below a level q, and the friction balance mu_f * (fraction forward) = mu_b * (fraction
// backward) fixes q. The mean activation of this bump shape is 5/16.
function wormPrediction() {
  const eps = 0.5 * strength.value, F = muBack.value / Math.max(1e-9, muBack.value + muFwd.value);
  const q = Math.pow(0.5 * (1 - Math.cos(Math.PI * F)), 3);
  return eps * waveSpeed.value * (q - 5 / 16);
}

function drawWorm(gr, C) {
  const head = P.x[worm.nodes[0]];
  Sim.setView([head - 0.72, -0.1, head + 0.18, 0.25]);
  gr.ground(0);
  for (let x = Math.floor(head - 1.2); x < head + 0.6; x += 0.1) gr.line(x, 0, x, -0.02, { color: C.muted, width: 1 });
  // Each segment keeps its area: shorter means fatter (the hydrostatic skeleton).
  for (const m of muscles) {
    const x1 = P.x[m.a], x2 = P.x[m.b], len = Math.abs(x1 - x2);
    const H = worm.H0 * worm.L0 / Math.max(0.2 * worm.L0, len);
    const xm = 0.5 * (x1 + x2), ym = worm.nodeR + 0.5 * H - 0.01;
    gr.rect(Math.min(x1, x2), ym - H / 2, len, H, { fill: m.act > 0.3 ? '#c46a6a' : '#d99a8f' });
    gr.line(xm, ym - H / 2, xm, ym + H / 2, { color: '#a35d55', width: 1 });
  }
  gr.circle(head + 0.006, worm.nodeR + 0.012, 0.012, { fill: '#c46a6a' });
  gr.text('head', head + 0.03, 0.06, { color: C.muted, size: 12 });
  const dist = head - stats.wormStart;
  Sim.show('distance crawled', (dist * 100).toFixed(1) + ' cm');
  if (time > 0.5) Sim.show('mean speed', (dist / time * 100).toFixed(2) + ' cm/s');
  if (muBack.value + muFwd.value > 0) Sim.show('quasi-static prediction', (wormPrediction() * 100).toFixed(2) + ' cm/s');
  Sim.plot('head x (m)', head, { color: C.red });
}

function loadSide() { return 0.03 + 0.03 * Math.cbrt(loadMass.value); }

function drawHang(gr, C) {
  Sim.setView([-0.5, -0.05, 0.5, 0.65]);
  gr.ground(0);
  const m = muscles[0], top = hang.topI, load = hang.loadI;
  gr.rect(-0.15, hang.top, 0.3, 0.03, { fill: '#d9d6cc' });
  const len = Math.hypot(P.x[load] - P.x[top], P.y[load] - P.y[top]);
  gr.line(P.x[top], P.y[top], P.x[load], P.y[load], { color: m.act > 0.05 ? C.red : '#d99a8f', width: 6 + 10 * m.act });
  const side = loadSide();
  gr.rect(P.x[load] - side, P.y[load] - 2 * side, 2 * side, 2 * side, { fill: C.gray });
  gr.line(-0.3, hang.top - hang.L0, 0.3, hang.top - hang.L0, { color: C.muted, dash: [4, 4], width: 1 });
  gr.text('rest length', 0.12, hang.top - hang.L0 + 0.015, { color: C.muted, size: 11 });
  const W = loadMass.value * g;
  Sim.show('activation a', m.act.toFixed(2));
  Sim.show('shortening', ((hang.L0 - len) * 100).toFixed(2) + ' cm');
  Sim.show('muscle tension', (law.value === 'active stress' ? activeStressForce(m) + m.force : m.force).toFixed(2) + ' N');
  Sim.show('load weight m g', W.toFixed(2) + ' N');
  Sim.show('fastest load speed so far', stats.maxSpeed.toFixed(2) + ' m/s');
  Sim.plot('shortening (cm)', (hang.L0 - len) * 100, { color: C.red });
  Sim.plot('activation x 10', m.act * 10, { color: C.orange });
}

Sim.run({ title: 'Muscles and travelling waves (starter)', view: [-0.1, -0.45, 1.7, 0.4], reset, step, draw, dt: frameDt, grid: false });
