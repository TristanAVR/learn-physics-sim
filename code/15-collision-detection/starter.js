// Chapter 15 - Collision Detection (starter)
// Scene 1: particles poured into a static scene described by a signed distance function (SDF), with
// particle-particle contacts found either by brute force or by a uniform spatial hash.
// Scene 2: a row of small fast balls fired at a thin wall, to see tunnelling and the speculative-margin fix.
// Your jobs, marked TODO: the capsule SDF, the particle-vs-scene response, the spatial hash, and the
// speculative contact. Until then the particles fall straight through the scene.

const POUR = 'pouring particles', BALLS = 'fast balls vs thin wall';
const mode = Sim.select('scene', [POUR, BALLS], POUR, { reset: true });
const search = Sim.select('neighbour search', ['brute force', 'spatial hash'], 'brute force');
const count = Sim.slider('particles (pouring)', 100, 2000, 600, { step: 50, reset: true });
const showCells = Sim.checkbox('draw hash cells', false);
const showField = Sim.checkbox('draw distance field', false);
const speed = Sim.slider('ball speed (m/s)', 0, 40, 16, { step: 0.5, reset: true });
const substeps = Sim.slider('substeps', 1, 16, 8, { step: 1, reset: true });
const speculative = Sim.checkbox('speculative margin |v| h', false, { reset: true });

const G = 9.81;
const R_POUR = 0.012;          // particle radius in the pouring scene (m)
const R_BALL = 0.02;           // ball radius (m)
const WALL_X = 0.8, WALL_T = 0.01, WALL_HALF_H = 0.6;   // thin wall: centre x, thickness, half height (m)
const LANES = 20;              // balls fired side by side
const MAX_N = 2000;

// Particle state in flat arrays: position, position at the start of the substep, velocity.
const x = new Float64Array(MAX_N), y = new Float64Array(MAX_N);
const px = new Float64Array(MAX_N), py = new Float64Array(MAX_N);
const vx = new Float64Array(MAX_N), vy = new Float64Array(MAX_N);
let n = 0, radius = R_POUR;

// Speculative contacts: one plane per particle, found at the start of the substep.
const specOn = new Uint8Array(MAX_N);
const specNx = new Float64Array(MAX_N), specNy = new Float64Array(MAX_N), specD = new Float64Array(MAX_N);

// Spatial hash (counting sort, after Mueller's "Ten Minute Physics", tutorial 11).
const TABLE_SIZE = 2 * MAX_N;
const cellStart = new Int32Array(TABLE_SIZE + 1);
const cellEntries = new Int32Array(MAX_N);
let cellSize = 2 * R_POUR;

let pairTests = 0, collisionMs = 0, collisionMsAvg = 0, emitted = 0, fireTime = 0;
let rng;

// ------------------------------------------------------------------------------------------------
// Signed distance functions: negative inside, positive outside, zero on the surface.
// ------------------------------------------------------------------------------------------------
function sdCircle(px, py, cx, cy, r) {
  return Math.hypot(px - cx, py - cy) - r;
}

// Half-plane through point (ax, ay) with unit normal (nx, ny) pointing into the free space.
function sdHalfPlane(px, py, ax, ay, nx, ny) {
  return (px - ax) * nx + (py - ay) * ny;
}

// Box with centre (cx, cy), half extents (hx, hy), rotated by angle (rad).
function sdBox(px, py, cx, cy, hx, hy, angle) {
  const c = Math.cos(angle), s = Math.sin(angle);
  const lx = c * (px - cx) + s * (py - cy);       // the point in the box's own frame
  const ly = -s * (px - cx) + c * (py - cy);
  const dx = Math.abs(lx) - hx, dy = Math.abs(ly) - hy;
  const outside = Math.hypot(Math.max(dx, 0), Math.max(dy, 0));
  const inside = Math.min(Math.max(dx, dy), 0);
  return outside + inside;
}

// Capsule: all points within distance r of the segment a-b.
function sdCapsule(px, py, ax, ay, bx, by, r) {
  // TODO: find the closest point q on the segment a-b (project p onto the line, clamp the parameter t
  // to [0, 1]), then return the distance from p to q minus r.
  return 1e9;
}

// The whole scene: the union of shapes is the minimum of their distances.
const SLOPE_A = [-0.45, 0.95, 0.42, 0.02, -0.3];
const SLOPE_B = [0.45, 0.6, 0.42, 0.02, 0.3];
const CAPSULE = [-0.25, 0.3, 0.15, 0.22, 0.05];
const BUMPER = [0.62, 0.0, 0.14];
function sceneSdf(px, py) {
  if (mode.value === BALLS) return sdBox(px, py, WALL_X, 0.5, WALL_T / 2, WALL_HALF_H, 0);
  let d = sdHalfPlane(px, py, 0, 0, 0, 1);                  // floor
  d = Math.min(d, sdHalfPlane(px, py, -1, 0, 1, 0));        // left wall
  d = Math.min(d, sdHalfPlane(px, py, 1, 0, -1, 0));        // right wall
  d = Math.min(d, sdBox(px, py, SLOPE_A[0], SLOPE_A[1], SLOPE_A[2], SLOPE_A[3], SLOPE_A[4]));
  d = Math.min(d, sdBox(px, py, SLOPE_B[0], SLOPE_B[1], SLOPE_B[2], SLOPE_B[3], SLOPE_B[4]));
  d = Math.min(d, sdCapsule(px, py, CAPSULE[0], CAPSULE[1], CAPSULE[2], CAPSULE[3], CAPSULE[4]));
  d = Math.min(d, sdCircle(px, py, BUMPER[0], BUMPER[1], BUMPER[2]));
  return d;
}

// Surface normal = normalised gradient of the SDF, by central differences.
const nrm = { x: 0, y: 1 };
function sdfNormal(px, py) {
  const e = 1e-4;
  const gx = sceneSdf(px + e, py) - sceneSdf(px - e, py);
  const gy = sceneSdf(px, py + e) - sceneSdf(px, py - e);
  const l = Math.hypot(gx, gy) || 1;
  nrm.x = gx / l; nrm.y = gy / l;
  return nrm;
}

// Particle against the static scene: if the particle overlaps the surface, push it out along the normal
// by the penetration depth.
function collideSdf(i) {
  // TODO: d = sceneSdf(x[i], y[i]) - radius is the gap between the particle's surface and the scene.
  // If d < 0, move the particle by -d along the normal from sdfNormal(x[i], y[i]).
}

// ------------------------------------------------------------------------------------------------
// Particle-particle contacts
// ------------------------------------------------------------------------------------------------
function solvePair(i, j) {
  pairTests++;
  const dx = x[j] - x[i], dy = y[j] - y[i];
  const d2 = dx * dx + dy * dy, D = 2 * radius;
  if (d2 >= D * D || d2 === 0) return;
  const d = Math.sqrt(d2), c = 0.5 * (D - d) / d;    // equal masses: each moves half the overlap
  x[i] -= dx * c; y[i] -= dy * c;
  x[j] += dx * c; y[j] += dy * c;
}

function solvePairsBruteForce() {
  for (let i = 0; i < n; i++)
    for (let j = i + 1; j < n; j++) solvePair(i, j);
}

function hashCell(ix, iy) {
  const h = Math.imul(ix, 92837111) ^ Math.imul(iy, 689287499);   // 32-bit integer multiply
  return Math.abs(h) % TABLE_SIZE;
}

function buildHash() {
  // TODO: counting sort of the n particles into TABLE_SIZE buckets, bucket = hashCell(cell x, cell y)
  // with cell = Math.floor(position / cellSize).
  //   1. cellStart.fill(0), then count the particles in each bucket.
  //   2. Running sum, so cellStart[h] holds the END of bucket h; set cellStart[TABLE_SIZE] = n.
  //   3. For each particle: cellStart[h]--, then cellEntries[cellStart[h]] = i.
  // Afterwards the particles of bucket h are cellEntries[cellStart[h]] ... cellEntries[cellStart[h + 1] - 1].
}

// Visit every particle stored in the 3x3 block of cells around particle i.
function forEachHashNeighbour(i, visit) {
  // TODO: find particle i's cell (ix, iy); for each of the 9 cells (ix + dx, iy + dy), dx and dy in
  // -1..1, look up its bucket and call visit(j) for every particle j stored there.
}

function solvePairsHash() {
  buildHash();
  for (let i = 0; i < n; i++) {
    forEachHashNeighbour(i, (j) => { if (j > i) solvePair(i, j); });
  }
}

// ------------------------------------------------------------------------------------------------
// Tunnelling defence: a speculative contact. Before moving, look for a surface within the distance the
// particle can travel this substep (|v| h). If there is one, remember its plane and keep the particle on
// the near side of that plane, wherever the prediction lands.
// ------------------------------------------------------------------------------------------------
function findSpeculativeContact(i, h) {
  specOn[i] = 0;
  // TODO: margin = |v| h. If the gap sceneSdf(x[i], y[i]) - radius is smaller than the margin, store a
  // contact: specOn[i] = 1, the normal (specNx, specNy) from sdfNormal, and the distance specD[i].
}

function solveSpeculativeContact(i) {
  if (!specOn[i]) return;
  // TODO: the signed gap to the stored plane is C = n . (x - p) + specD[i] - radius, where p = (px[i], py[i])
  // is where the particle started the substep. If C < 0, move the particle by -C along n.
}

// ------------------------------------------------------------------------------------------------
// Scene setup and the step
// ------------------------------------------------------------------------------------------------
function reset() {
  rng = Sim.random(7);
  n = 0; emitted = 0; collisionMsAvg = 0;
  fieldCache = null;
  if (mode.value === POUR) {
    radius = R_POUR;
    cellSize = 2 * radius;
    Sim.setView([-1.05, -0.05, 1.05, 1.45]);
  } else {
    radius = R_BALL;
    cellSize = 2 * radius;
    Sim.setView([-0.45, -0.15, 1.55, 1.15]);
    fireBalls();
  }
}

function fireBalls() {
  n = 0; fireTime = 0;
  const s = speed.value / 60 / substeps.value;      // distance moved per substep
  for (let k = 0; k < LANES; k++) {
    // evenly spread starting offsets, so the substep positions fall at every phase relative to the wall
    x[n] = -((k + 0.5) / LANES) * s; y[n] = 0.02 + k * 0.05;
    vx[n] = speed.value; vy[n] = 0;
    n++;
  }
}

function emit() {
  if (emitted >= count.value) return;
  for (let k = 0; k < 6 && emitted < count.value; k++) {
    x[n] = -0.75 + 0.003 * rng(); y[n] = 1.25 + k * 2.1 * R_POUR;
    vx[n] = 0.8; vy[n] = -0.2;
    n++; emitted++;
  }
}

function step(dt) {
  const ns = substeps.value, h = dt / ns, pouring = mode.value === POUR;
  if (pouring) emit();
  else if ((fireTime += dt) > 1.5) fireBalls();      // fire a new volley every 1.5 s
  pairTests = 0; collisionMs = 0;
  for (let s = 0; s < ns; s++) {
    for (let i = 0; i < n; i++) {
      px[i] = x[i]; py[i] = y[i];
      if (pouring) vy[i] -= G * h;
      if (speculative.value) findSpeculativeContact(i, h);
      x[i] += vx[i] * h; y[i] += vy[i] * h;
    }
    if (pouring) {
      const t0 = performance.now();
      if (search.value === 'spatial hash') solvePairsHash(); else solvePairsBruteForce();
      collisionMs += performance.now() - t0;
    }
    for (let i = 0; i < n; i++) {
      if (speculative.value) solveSpeculativeContact(i);
      collideSdf(i);
    }
    for (let i = 0; i < n; i++) {
      vx[i] = (x[i] - px[i]) / h; vy[i] = (y[i] - py[i]) / h;
      if (!pouring && x[i] > 1.45) { x[i] = 1.45; vx[i] = 0; }   // park the balls that got through
    }
  }
  collisionMsAvg += 0.05 * (collisionMs - collisionMsAvg);
}

// ------------------------------------------------------------------------------------------------
// Drawing
// ------------------------------------------------------------------------------------------------
let fieldCache = null;
function drawField(g) {
  const ctx = g.ctx, W = ctx.canvas.width, H = ctx.canvas.height, cell = 4;
  if (!fieldCache || fieldCache.W !== W || fieldCache.H !== H) {
    const w = Math.ceil(W / cell), hgt = Math.ceil(H / cell);
    const off = document.createElement('canvas');
    off.width = w; off.height = hgt;
    const octx = off.getContext('2d'), img = octx.createImageData(w, hgt);
    const s0 = g.toScreen(0, 0), s1 = g.toScreen(1, 0), scale = s1.x - s0.x;
    for (let j = 0; j < hgt; j++) {
      for (let i = 0; i < w; i++) {
        const wx = ((i + 0.5) * cell - s0.x) / scale, wy = (s0.y - (j + 0.5) * cell) / scale;
        const d = sceneSdf(wx, wy);
        const band = 0.85 + 0.15 * Math.cos(d / 0.1 * 2 * Math.PI);  // contour bands every 10 cm
        const k = 4 * (j * w + i);
        if (d > 0) { img.data[k] = 120 * band; img.data[k + 1] = 160 * band; img.data[k + 2] = 230 * band; }
        else { img.data[k] = 235 * band; img.data[k + 1] = 150 * band; img.data[k + 2] = 70 * band; }
        img.data[k + 3] = 60;
      }
    }
    octx.putImageData(img, 0, 0);
    fieldCache = { W, H, off };
  }
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(fieldCache.off, 0, 0, W, H);
}

function boxCorners(b) {
  const c = Math.cos(b[4]), s = Math.sin(b[4]), pts = [];
  for (const [sx, sy] of [[-1, -1], [1, -1], [1, 1], [-1, 1]]) {
    const lx = sx * b[2], ly = sy * b[3];
    pts.push({ x: b[0] + c * lx - s * ly, y: b[1] + s * lx + c * ly });
  }
  return pts;
}

function drawCapsule(g, c, style) {
  const ux = c[2] - c[0], uy = c[3] - c[1], l = Math.hypot(ux, uy);
  const nx = -uy / l * c[4], ny = ux / l * c[4];
  g.polyline([{ x: c[0] + nx, y: c[1] + ny }, { x: c[2] + nx, y: c[3] + ny },
    { x: c[2] - nx, y: c[3] - ny }, { x: c[0] - nx, y: c[1] - ny }], { fill: style.fill, closed: true });
  g.circle(c[0], c[1], c[4], style);
  g.circle(c[2], c[3], c[4], style);
}

function drawParticles(g, color, highlight) {
  const ctx = g.ctx, rp = radius * g.pxPerUnit * (window.devicePixelRatio || 1);
  ctx.globalAlpha = 1;
  ctx.fillStyle = color;
  ctx.beginPath();
  for (let i = 0; i < n; i++) {
    const s = g.toScreen(x[i], y[i]);
    ctx.moveTo(s.x + rp, s.y);
    ctx.arc(s.x, s.y, rp, 0, 2 * Math.PI);
  }
  ctx.fill();
  if (highlight) for (const [i, c] of highlight) g.circle(x[i], y[i], radius, { fill: c });
}

function drawPouring(g) {
  const solid = { fill: '#cfcabb' };
  g.ground(0);
  g.rect(-1.3, -0.1, 0.3, 1.7, solid);
  g.rect(1.0, -0.1, 0.3, 1.7, solid);
  g.polyline(boxCorners(SLOPE_A), { fill: solid.fill, closed: true });
  g.polyline(boxCorners(SLOPE_B), { fill: solid.fill, closed: true });
  drawCapsule(g, CAPSULE, solid);
  g.circle(BUMPER[0], BUMPER[1], BUMPER[2], solid);

  let highlight = null;
  const useHash = search.value === 'spatial hash';
  if (useHash && showCells.value && n > 0) {
    buildHash();
    const seen = new Set();
    for (let i = 0; i < n; i++) {
      const ix = Math.floor(x[i] / cellSize), iy = Math.floor(y[i] / cellSize), key = ix * 100003 + iy;
      if (seen.has(key)) continue;
      seen.add(key);
      g.rect(ix * cellSize, iy * cellSize, cellSize, cellSize, { stroke: Sim.colors.teal, width: 0.6, alpha: 0.7 });
    }
    // the particle nearest the mouse: its 3x3 query block, the candidates found there, and real contacts
    let best = -1, bestD = 0.08;
    for (let i = 0; i < n; i++) {
      const d = Math.hypot(x[i] - Sim.mouse.x, y[i] - Sim.mouse.y);
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best >= 0) {
      const ix = Math.floor(x[best] / cellSize), iy = Math.floor(y[best] / cellSize);
      g.rect((ix - 1) * cellSize, (iy - 1) * cellSize, 3 * cellSize, 3 * cellSize, { stroke: Sim.colors.purple, width: 2 });
      highlight = [];
      forEachHashNeighbour(best, (j) => {
        if (j === best) return;
        const touching = Math.hypot(x[j] - x[best], y[j] - y[best]) < 2 * radius * 1.05;
        highlight.push([j, touching ? Sim.colors.red : Sim.colors.orange]);
      });
      highlight.push([best, Sim.colors.purple]);
    }
  }
  drawParticles(g, Sim.colors.blue, highlight);
  Sim.show('particles', n);
  Sim.show('pair tests per substep', Math.round(pairTests / substeps.value));
  if (search.value === 'brute force') Sim.show('n(n-1)/2', n * (n - 1) / 2);
  Sim.show('contact time', collisionMsAvg.toFixed(2) + ' ms/frame');
  if (useHash && showCells.value) Sim.show('mouse', 'hover a particle: purple query, orange candidates, red contacts');
}

function drawBalls(g) {
  g.rect(WALL_X - WALL_T / 2, 0.5 - WALL_HALF_H, WALL_T, 2 * WALL_HALF_H, { fill: Sim.colors.ink });
  let through = 0;
  for (let i = 0; i < n; i++) {
    const passed = x[i] > WALL_X;
    if (passed) through++;
    // the positions the simulation actually visits: one per substep
    const sstep = speed.value / 60 / substeps.value;
    for (let k = 1; k <= 4 && vx[i] > 0; k++) g.dot(x[i] - k * sstep, y[i], 1.5, Sim.colors.gray);
    g.circle(x[i], y[i], radius, { fill: passed ? Sim.colors.red : Sim.colors.blue });
  }
  // Closed form (see the chapter): with s = |v| h and L = r + t/2, the fraction of evenly spread balls
  // that tunnel is 0 for s <= L, 2 - 2L/s for L < s < 2L, and 1 for s >= 2L = 2r + t.
  const h = 1 / 60 / substeps.value, s = speed.value * h, L = radius + WALL_T / 2;
  const predicted = s <= L ? 0 : s < 2 * L ? 2 - 2 * L / s : 1;
  Sim.show('substep h', (h * 1000).toFixed(2) + ' ms');
  Sim.show('distance per substep |v| h', (s * 100).toFixed(2) + ' cm');
  Sim.show('first tunnelling (r + t/2) / h', (L / h).toFixed(2) + ' m/s');
  Sim.show('all tunnel (2r + t) / h', (2 * L / h).toFixed(2) + ' m/s');
  Sim.show('predicted to tunnel', speculative.value ? '0 (margin on)' : (predicted * LANES).toFixed(1) + ' of ' + LANES);
  Sim.show('tunnelled', through + ' of ' + LANES);
}

function draw(g) {
  if (showField.value) drawField(g);
  if (mode.value === POUR) drawPouring(g); else drawBalls(g);
}

Sim.run({ title: 'Collision detection (starter)', view: [-1.05, -0.05, 1.05, 1.45], reset, step, draw, grid: false });
