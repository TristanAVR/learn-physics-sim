// Chapter 31 - Hunting Explosions (reference solution)
// A deliberately fragile scene: a chain with a weight (XPBD), a box of explicit springs, and three balls.
// Plant a failure, or a mystery one, and let the watchdog find the first bad frame and the particle where
// it started. Every frame's input is recorded, so the failure can be replayed exactly. The fuzzer runs the
// scene headless with random settings and random mouse drags and counts failures by kind.

const causeSel = Sim.select('cause to plant', ['stiff box springs', 'tiny box mass', 'heavy weight',
  'overlapping balls', 'ball deep in the ground', 'coincident box corners', 'Jacobi with redundant links'],
  'stiff box springs');
Sim.button('plant cause', () => queueEvent(causeSel.value));
Sim.button('mystery 1', () => queueEvent('mystery 1'));
Sim.button('mystery 2', () => queueEvent('mystery 2'));
Sim.button('mystery 3', () => queueEvent('mystery 3'));
const watchOn = Sim.checkbox('watchdog', true);
const showVel = Sim.checkbox('draw velocities', false);
const showErr = Sim.checkbox('draw constraint errors', true);
Sim.button('continue after a trip', () => { frozen = false; report = null; });
Sim.button('replay recording', () => startReplay(Infinity));
Sim.button('replay to the frame before the trip', () => { if (firstReport) startReplay(firstReport.frame - 1); });
Sim.button('run fuzzer (200 runs)', () => { fuzz = runFuzzer(200, 300, 12345); });

const G = 9.81, frameDt = 1 / 60;
const defaults = { substeps: 4, iters: 10, kBox: 2000, weightMass: 0.5 };

// ---------------------------------------------------------------------------------------------
// The scene. Everything lives in one object S so the fuzzer can make as many as it likes.
// ---------------------------------------------------------------------------------------------
function makeScene(settings) {
  const n = 32;
  const S = {
    settings: Object.assign({}, settings), frame: 0, n: 0,
    x: new Float64Array(n), y: new Float64Array(n), vx: new Float64Array(n), vy: new Float64Array(n),
    px: new Float64Array(n), py: new Float64Array(n), w: new Float64Array(n), r: new Float64Array(n),
    body: [], dists: [], springs: [], box: [], balls: [], chain: [], jacobi: false,
    grab: -1, wasDown: false, firstNaN: null,
  };
  const add = (x, y, m, r, body) => {
    const i = S.n++;
    S.x[i] = x; S.y[i] = y; S.w[i] = m > 0 ? 1 / m : 0; S.r[i] = r; S.body.push(body);
    return i;
  };
  // Chain: pinned at the top, 10 links of 8 cm, 50 g beads, a weight at the end.
  for (let k = 0; k <= 10; k++) {
    const m = k === 0 ? 0 : (k === 10 ? settings.weightMass : 0.05);
    S.chain.push(add(0.5 + 0.08 * k, 1.7, m, k === 10 ? 0.06 : 0.015, k === 10 ? 'weight' : 'chain'));
  }
  for (let k = 0; k < 10; k++) S.dists.push({ a: S.chain[k], b: S.chain[k + 1], L: 0.08 });
  // Box: 4 corners of 0.2 kg joined by 6 explicit springs (sides and diagonals).
  const bx = -0.75, by = 0.02, s = 0.3;
  for (const [cx, cy] of [[0, 0], [s, 0], [s, s], [0, s]]) S.box.push(add(bx + cx, by + cy, 0.2, 0.02, 'box'));
  for (const [a, b] of [[0, 1], [1, 2], [2, 3], [3, 0], [0, 2], [1, 3]]) {
    const ia = S.box[a], ib = S.box[b];
    S.springs.push({ a: ia, b: ib, L: Math.hypot(S.x[ib] - S.x[ia], S.y[ib] - S.y[ia]) });
  }
  // Balls.
  for (const [x, y] of [[0.1, 0.4], [0.35, 0.8], [0.0, 1.2]]) S.balls.push(add(x, y, 0.3, 0.08, 'ball'));
  return S;
}

// Plant a failure cause into a running scene (a recorded event, so replays plant it at the same frame).
function applyEvent(S, name) {
  const B = S.box, b = S.balls;
  if (name === 'stiff box springs') S.settings.kBox *= 10;
  else if (name === 'tiny box mass') S.w[B[2]] *= 1000;
  else if (name === 'heavy weight') S.w[S.chain[10]] /= 1000;
  else if (name === 'overlapping balls') { S.x[b[1]] = S.x[b[0]] + 0.01; S.y[b[1]] = S.y[b[0]]; }
  else if (name === 'ball deep in the ground') { S.y[b[2]] = -0.3; }
  else if (name === 'coincident box corners') { S.x[B[1]] = S.x[B[0]]; S.y[B[1]] = S.y[B[0]]; }
  else if (name === 'Jacobi with redundant links') {
    S.jacobi = true;
    for (let k = 0; k + 2 < S.chain.length; k++) S.dists.push({ a: S.chain[k], b: S.chain[k + 2], L: 0.16 });
  }
  else if (name === 'mystery 1') S.settings.substeps = 1;
  else if (name === 'mystery 2') { S.dists.push({ a: S.chain[9], b: S.chain[10], L: 0.02 }); }
  else if (name === 'mystery 3') { S.x[b[2]] = S.x[B[3]]; S.y[b[2]] = S.y[B[3]] - 0.05; }
}

// ---------------------------------------------------------------------------------------------
// One frame: substeps of (explicit spring forces, integrate, project constraints, update velocities).
// input = { x, y, down } is the mouse in world coordinates. Returns nothing; S.firstNaN records the
// substep and phase where a non-finite number first appeared.
// ---------------------------------------------------------------------------------------------
function stepScene(S, input) {
  const sub = S.settings.substeps, h = frameDt / sub;
  // Mouse: grab the nearest particle on press, let go on release.
  if (input.down && !S.wasDown) {
    let best = -1, bd = 0.15;
    for (let i = 0; i < S.n; i++) {
      const d = Math.hypot(S.x[i] - input.x, S.y[i] - input.y);
      if (S.w[i] > 0 && d < bd) { bd = d; best = i; }
    }
    S.grab = best;
  }
  if (!input.down) S.grab = -1;
  S.wasDown = input.down;

  for (let s = 0; s < sub; s++) {
    // Explicit forces: gravity and the box springs (with a little damping).
    const fx = new Float64Array(S.n), fy = new Float64Array(S.n);
    for (const sp of S.springs) {
      const dx = S.x[sp.b] - S.x[sp.a], dy = S.y[sp.b] - S.y[sp.a], len = Math.hypot(dx, dy);
      const nx = dx / len, ny = dy / len;              // no guard on purpose: len = 0 gives NaN
      const relv = (S.vx[sp.b] - S.vx[sp.a]) * nx + (S.vy[sp.b] - S.vy[sp.a]) * ny;
      const f = S.settings.kBox * (len - sp.L) + 0.5 * relv;
      fx[sp.a] += f * nx; fy[sp.a] += f * ny; fx[sp.b] -= f * nx; fy[sp.b] -= f * ny;
    }
    noteNaN(S, s, 'spring forces', fx, fy);
    for (let i = 0; i < S.n; i++) {
      if (S.w[i] === 0) continue;
      S.vx[i] += h * S.w[i] * fx[i];
      S.vy[i] += h * (S.w[i] * fy[i] - G);
      S.px[i] = S.x[i]; S.py[i] = S.y[i];
      S.x[i] += h * S.vx[i]; S.y[i] += h * S.vy[i];
    }
    noteNaN(S, s, 'integration');
    // Constraints: the chain (Gauss-Seidel, or Jacobi without averaging once planted), the mouse,
    // collisions with the ground and between balls and box corners.
    for (let it = 0; it < S.settings.iters; it++) {
      if (S.jacobi) solveDistancesJacobi(S); else for (const d of S.dists) solveDistance(S, d);
    }
    if (S.grab >= 0) {
      const i = S.grab, at = 1e-3 / (h * h);
      for (const [p, t] of [['x', input.x], ['y', input.y]]) S[p][i] += S.w[i] * (t - S[p][i]) / (S.w[i] + at);
    }
    for (const a of S.balls) {
      for (const b of S.balls) if (a < b) pushApart(S, a, b);
      for (const b of S.box) pushApart(S, a, b);
    }
    for (let i = 0; i < S.n; i++) {
      const pen = S.r[i] - S.y[i];
      if (S.w[i] > 0 && pen > 0) {
        S.y[i] += pen;
        const slide = S.x[i] - S.px[i], cap = 0.5 * pen;   // Coulomb friction, mu = 0.5
        S.x[i] -= Math.abs(slide) <= cap ? slide : Math.sign(slide) * cap;
      }
    }
    noteNaN(S, s, 'constraints');
    for (let i = 0; i < S.n; i++) {
      if (S.w[i] === 0) { S.vx[i] = 0; S.vy[i] = 0; continue; }
      S.vx[i] = (S.x[i] - S.px[i]) / h; S.vy[i] = (S.y[i] - S.py[i]) / h;
    }
  }
  S.frame++;
}

function solveDistance(S, d) {
  const dx = S.x[d.b] - S.x[d.a], dy = S.y[d.b] - S.y[d.a], len = Math.hypot(dx, dy);
  const wa = S.w[d.a], wb = S.w[d.b];
  if (len < 1e-12 || wa + wb === 0) return;
  const dl = -(len - d.L) / (wa + wb), nx = dx / len, ny = dy / len;
  S.x[d.a] -= wa * nx * dl; S.y[d.a] -= wa * ny * dl;
  S.x[d.b] += wb * nx * dl; S.y[d.b] += wb * ny * dl;
}

// Jacobi: every constraint computes its correction from the same positions, and the corrections are
// summed. Without dividing by the number of constraints per particle, it overshoots.
function solveDistancesJacobi(S) {
  const cx = new Float64Array(S.n), cy = new Float64Array(S.n);
  for (const d of S.dists) {
    const dx = S.x[d.b] - S.x[d.a], dy = S.y[d.b] - S.y[d.a], len = Math.hypot(dx, dy);
    const wa = S.w[d.a], wb = S.w[d.b];
    if (len < 1e-12 || wa + wb === 0) continue;
    const dl = -(len - d.L) / (wa + wb), nx = dx / len, ny = dy / len;
    cx[d.a] -= wa * nx * dl; cy[d.a] -= wa * ny * dl;
    cx[d.b] += wb * nx * dl; cy[d.b] += wb * ny * dl;
  }
  for (let i = 0; i < S.n; i++) { S.x[i] += cx[i]; S.y[i] += cy[i]; }
}

function pushApart(S, a, b) {
  const dx = S.x[b] - S.x[a], dy = S.y[b] - S.y[a], d = Math.hypot(dx, dy), min = S.r[a] + S.r[b];
  if (d >= min || d < 1e-12) return;
  const wa = S.w[a], wb = S.w[b], C = d - min, nx = dx / d, ny = dy / d, s = -C / (wa + wb);
  S.x[a] -= wa * nx * s; S.y[a] -= wa * ny * s;
  S.x[b] += wb * nx * s; S.y[b] += wb * ny * s;
}

// ---------------------------------------------------------------------------------------------
// Monitors
// ---------------------------------------------------------------------------------------------
// Total energy (kinetic + gravity + box springs), and each particle's speed for finding the culprit.
function energy(S) {
  let E = 0;
  const speed = new Float64Array(S.n);
  for (let i = 0; i < S.n; i++) {
    if (S.w[i] === 0) continue;
    const m = 1 / S.w[i];
    speed[i] = Math.hypot(S.vx[i], S.vy[i]);
    E += 0.5 * m * speed[i] * speed[i] + m * G * S.y[i];
  }
  for (const sp of S.springs) {
    const len = Math.hypot(S.x[sp.b] - S.x[sp.a], S.y[sp.b] - S.y[sp.a]);
    E += 0.5 * S.settings.kBox * (len - sp.L) * (len - sp.L);
  }
  return { E, speed };
}

function maxChainError(S) {
  let e = 0, at = -1;
  for (let k = 0; k < 10; k++) {
    const d = S.dists[k], err = Math.abs(Math.hypot(S.x[d.b] - S.x[d.a], S.y[d.b] - S.y[d.a]) - d.L);
    if (!(err <= e)) { e = err; at = d.b; }
  }
  return { e, at };
}

// Watchdog part 1: the first particle whose position or velocity (or force, if given) is not a finite
// number, or -1.
function findNonFinite(S, fx, fy) {
  for (let i = 0; i < S.n; i++) {
    if (!Number.isFinite(S.x[i]) || !Number.isFinite(S.y[i]) || !Number.isFinite(S.vx[i]) || !Number.isFinite(S.vy[i])) return i;
    if (fx && (!Number.isFinite(fx[i]) || !Number.isFinite(fy[i]))) return i;
  }
  return -1;
}

function noteNaN(S, substep, phase, fx, fy) {
  if (S.firstNaN) return;
  const i = findNonFinite(S, fx, fy);
  if (i >= 0) S.firstNaN = { frame: S.frame, substep, phase, particle: i };
}

// Watchdog part 2: compare this frame with the last one. Returns a report or null.
// before = { E, speed } measured after this frame's input was applied, after = the same at the end.
// The culprit is the particle whose speed rose most: a light particle can be fast with little energy.
function checkFrame(S, before, after, external) {
  if (S.firstNaN) {
    const r = S.firstNaN;
    return { kind: 'NaN', frame: r.frame, particle: r.particle, detail: 'substep ' + r.substep + ', after ' + r.phase };
  }
  const tol = 1 + 0.01 * Math.abs(before.E);              // 1 J plus 1%: bigger than honest wobble
  if (!external && after.E - before.E > tol) {
    let worst = 0, gain = -Infinity;
    for (let i = 0; i < S.n; i++) if (after.speed[i] - before.speed[i] > gain) { gain = after.speed[i] - before.speed[i]; worst = i; }
    return { kind: 'energy jump', frame: S.frame - 1, particle: worst,
      detail: 'E rose by ' + (after.E - before.E).toFixed(2) + ' J in one frame' };
  }
  const ce = maxChainError(S);
  if (ce.e > 0.02) return { kind: 'stretch', frame: S.frame - 1, particle: ce.at, detail: 'chain link off by ' + (ce.e * 100).toFixed(1) + ' cm' };
  return null;
}

function stateHash(S) {
  let hsh = 0;
  for (let i = 0; i < S.n; i++) hsh = (hsh * 31 + Math.round(S.x[i] * 1e6) + 7 * Math.round(S.y[i] * 1e6)) % 1000000007;
  return hsh;
}

// ---------------------------------------------------------------------------------------------
// The live run: record every frame's input, replay on demand
// ---------------------------------------------------------------------------------------------
let world, recording, replay, pending, frozen, report, firstReport, fuzz = null;

function reset() {
  world = makeScene(defaults);
  recording = []; replay = null; pending = []; frozen = false; report = null; firstReport = null;
}

function queueEvent(name) { pending.push(name); }

function startReplay(stopFrame) {
  if (!recording.length) return;
  const rec = recording, first = firstReport;
  world = makeScene(defaults);
  replay = { frames: rec, index: 0, stopFrame };
  recording = rec; pending = []; report = null; frozen = false; firstReport = first;
  // Jump straight to the stop frame without drawing (the frames before it are known to be fine).
  if (isFinite(stopFrame)) {
    while (replay.index < Math.min(stopFrame, rec.length)) runFrame();
    frozen = true;
  }
}

function runFrame() {
  let input;
  if (replay) {
    if (replay.index >= replay.frames.length) { replay = null; frozen = true; return; }
    input = replay.frames[replay.index++];
  } else {
    input = { x: Sim.mouse.x, y: Sim.mouse.y, down: Sim.mouse.down, events: pending };
    pending = [];
    recording.push(input);
  }
  for (const e of input.events) applyEvent(world, e);
  const before = energy(world);
  stepScene(world, input);
  const after = energy(world);
  if (watchOn.value) {
    const r = checkFrame(world, before, after, world.grab >= 0);
    if (r) {
      r.hash = stateHash(world);
      if (!firstReport) firstReport = r;
      else if (r.frame === firstReport.frame) r.reproduced = r.hash === firstReport.hash;
      report = r; frozen = true;
    }
  }
}

function step() {
  if (frozen) return;
  runFrame();
}

// ---------------------------------------------------------------------------------------------
// The fuzzer: many short headless runs with random settings and random mouse drags.
// ---------------------------------------------------------------------------------------------
function randomSettings(rnd) {
  return {
    substeps: [1, 2, 4, 8][Math.floor(rnd() * 4)],
    iters: [2, 5, 10, 20][Math.floor(rnd() * 4)],
    kBox: Math.pow(10, 2.5 + 2.5 * rnd()),              // 316 to 100,000 N/m
    weightMass: 0.05 * Math.pow(10, 3 * rnd()),         // 1 to 1000 times a chain bead
  };
}

// A random mouse program: grab near a random particle at a random time and drag it somewhere.
function randomProgram(rnd, S) {
  const target = Math.floor(rnd() * S.n), start = Math.floor(rnd() * 200), len = 20 + Math.floor(rnd() * 60);
  const gx = -1 + 3 * rnd(), gy = 0.1 + 1.5 * rnd();
  const ox = S.x[target], oy = S.y[target];
  return (f) => {
    if (f < start || f > start + len) return { x: 0, y: -5, down: false, events: [] };
    const t = (f - start) / len;
    return { x: ox + (gx - ox) * t, y: oy + (gy - oy) * t, down: true, events: [] };
  };
}

function runFuzzer(runs, frames, seed) {
  const rnd = Sim.random(seed);
  const result = { runs, frames, kinds: {}, bySub: {}, grid: {} };
  for (let r = 0; r < runs; r++) {
    const settings = randomSettings(rnd);
    const S = makeScene(settings);
    const program = randomProgram(rnd, S);
    let kind = 'ok';
    for (let f = 0; f < frames && kind === 'ok'; f++) {
      const input = program(f);
      const before = energy(S);
      stepScene(S, input);
      const rep = checkFrame(S, before, energy(S), S.grab >= 0);
      if (rep) kind = rep.kind;
    }
    // Tally by kind, by substeps, and on a substeps x stiffness grid.
    result.kinds[kind] = (result.kinds[kind] || 0) + 1;
    const sk = settings.substeps;
    result.bySub[sk] = result.bySub[sk] || { n: 0, fail: 0 };
    result.bySub[sk].n++;
    if (kind !== 'ok') result.bySub[sk].fail++;
    const kb = Math.min(4, Math.floor((Math.log10(settings.kBox) - 2.5) / 0.5));
    const key = sk + ',' + kb;
    result.grid[key] = result.grid[key] || { n: 0, energy: 0 };
    result.grid[key].n++;
    if (kind === 'energy jump' || kind === 'NaN') result.grid[key].energy++;
  }
  return result;
}

// ---------------------------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------------------------
function draw(g) {
  const S = world, C = Sim.colors;
  g.ground(0);
  const ok = (v) => Number.isFinite(v) && Math.abs(v) < 50;
  // Chain, coloured by constraint error.
  for (let k = 0; k < 10; k++) {
    const d = S.dists[k];
    if (!ok(S.x[d.a]) || !ok(S.x[d.b])) continue;
    const err = Math.abs(Math.hypot(S.x[d.b] - S.x[d.a], S.y[d.b] - S.y[d.a]) - d.L);
    const col = showErr.value && err > 0.005 ? (err > 0.02 ? C.red : C.orange) : C.ink;
    g.line(S.x[d.a], S.y[d.a], S.x[d.b], S.y[d.b], { color: col, width: 2 });
  }
  g.rect(0.4, 1.7, 0.2, 0.03, { fill: '#d9d6cc' });
  for (const sp of S.springs) if (ok(S.x[sp.a]) && ok(S.x[sp.b])) g.line(S.x[sp.a], S.y[sp.a], S.x[sp.b], S.y[sp.b], { color: C.purple, width: 1.5 });
  for (let i = 0; i < S.n; i++) {
    if (!ok(S.x[i]) || !ok(S.y[i])) continue;
    const b = S.body[i];
    const col = b === 'ball' ? C.blue : b === 'box' ? C.purple : b === 'weight' ? C.gray : C.ink;
    g.circle(S.x[i], S.y[i], S.r[i], { fill: col });
    if (showVel.value) g.arrow(S.x[i], S.y[i], 0.05 * S.vx[i], 0.05 * S.vy[i], { color: C.green, width: 1.5 });
  }
  if (S.grab >= 0) g.line(S.x[S.grab], S.y[S.grab], Sim.mouse.x, Sim.mouse.y, { color: C.green, dash: [4, 4] });

  const e = energy(S).E;
  Sim.plot('total energy (J)', Number.isFinite(e) ? Math.max(-50, Math.min(200, e)) : NaN, { color: C.red });
  Sim.show('frame', S.frame);
  Sim.show('energy', Number.isFinite(e) ? e.toFixed(2) + ' J' : String(e));
  Sim.show('recording', replay ? 'replaying ' + replay.index + ' / ' + replay.frames.length : recording.length + ' frames');
  // Diagnostics: the quantities each failure cause changes.
  let wMax = 0, mMin = Infinity, mMax = 0;
  for (const i of S.box) wMax = Math.max(wMax, S.w[i]);
  for (const i of S.chain) if (S.w[i] > 0) { mMin = Math.min(mMin, 1 / S.w[i]); mMax = Math.max(mMax, 1 / S.w[i]); }
  const hs = frameDt / S.settings.substeps;
  Sim.show('box h * omega (estimate)', (hs * Math.sqrt(4 * S.settings.kBox * wMax)).toFixed(2));
  Sim.show('chain mass ratio', (mMax / mMin).toFixed(0));
  Sim.show('chain constraints', S.dists.length + (S.jacobi ? ' (Jacobi)' : ' (Gauss-Seidel)'));
  if (report) {
    const i = report.particle;
    if (ok(S.x[i]) && ok(S.y[i])) {
      g.circle(S.x[i], S.y[i], S.r[i] + 0.05, { color: C.red, width: 3 });
    }
    Sim.show('WATCHDOG', report.kind + ' at frame ' + report.frame);
    Sim.show('first particle', i + ' (' + S.body[i] + ')');
    Sim.show('detail', report.detail);
    if (report.reproduced !== undefined) Sim.show('replay reproduces it', report.reproduced ? 'yes, same state hash' : 'NO: not deterministic');
  } else if (frozen) {
    Sim.show('paused', 'at frame ' + S.frame + ': press Step (.) to advance one frame');
  }
  if (fuzz) drawFuzz(g, fuzz);
}

// The fuzzer's table goes into the info box: failures (energy jump or NaN) / runs, per cell.
function drawFuzz(g, F) {
  Sim.show('fuzzer', F.runs + ' runs x ' + F.frames + ' frames: ' +
    Object.keys(F.kinds).map((k) => k + ' ' + F.kinds[k]).join(', '));
  Sim.show('energy/NaN fails by box k', '316  1k  3.2k  10k  32k+ N/m');
  for (const sub of [1, 2, 4, 8]) {
    const cells = [];
    for (let j = 0; j < 5; j++) { const c = F.grid[sub + ',' + j]; cells.push(c ? c.energy + '/' + c.n : '-'); }
    Sim.show(sub + ' substeps', cells.join('  '));
  }
}

Sim.run({ title: 'A fragile scene and a watchdog', view: [-1.2, -0.85, 1.6, 2.2], reset, step, draw, dt: frameDt });
