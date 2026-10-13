// Chapter 17 - Cloth and Thin Shells (starter)
// A square of cloth as a grid of particles with XPBD stretch, shear and bending constraints, drawn with a
// small hand-made orthographic 3D projection. Drag the cloth with the mouse; turn on tearing and pull.
// Your jobs, marked TODO: the XPBD distance constraint, the shear and second-neighbour constraint lists,
// and tearing. Until solveDistance works, the cloth simply falls apart.
//   curtain:          starts horizontal, pinned at two corners (or along one edge), and swings down
//   drape over sphere: falls onto a ball standing on the floor

const CURTAIN = 'curtain', DRAPE = 'drape over sphere';
const scene = Sim.select('scene', [CURTAIN, DRAPE], CURTAIN, { reset: true });
const pinning = Sim.select('pinning (curtain)', ['two top corners', 'whole top edge'], 'two top corners', { reset: true });
const bending = Sim.select('bending', ['none', 'second neighbours', 'dihedral angle'], 'dihedral angle');
const stretchC = Sim.logSlider('stretch compliance (m/N)', -8, -2, 1e-7);
const shearC = Sim.logSlider('shear compliance (m/N)', -8, -1, 1e-3);
const bendC = Sim.logSlider('bending compliance (1/(N m))', -2, 6, 1e3);
const wind = Sim.slider('wind towards you (m/s)', 0, 6, 0, { step: 0.1 });
const tearing = Sim.checkbox('tearing', false);
const tearAt = Sim.slider('tear at stretch ratio', 1.05, 2, 1.3, { step: 0.01 });
const yawDeg = Sim.slider('view yaw (deg)', -90, 90, 25, { step: 1 });
const pitchDeg = Sim.slider('view pitch (deg)', -10, 60, 20, { step: 1 });

const G = 9.81;
const N = 25;                  // particles per side
const SPACING = 0.04;          // m between neighbours, so the cloth is 0.96 m square
const AREAL_DENSITY = 0.2;     // kg/m^2, a medium cotton
const SUBSTEPS = 15;
const DRAG = 1.0;              // linear air drag on each triangle, N per (m^2 m/s)
const SPHERE = { x: 0, y: 0.35, z: 0, r: 0.25 };
const THICK = 0.008;           // collision thickness of the cloth (m)
const MU = 0.5;                // friction against the sphere and the floor

const NP = N * N;
const x = new Float64Array(NP), y = new Float64Array(NP), z = new Float64Array(NP);
const px = new Float64Array(NP), py = new Float64Array(NP), pz = new Float64Array(NP);
const vx = new Float64Array(NP), vy = new Float64Array(NP), vz = new Float64Array(NP);
const w = new Float64Array(NP);
const id = (i, j) => j * N + i;          // column i, row j

let stretch, shear, bend2, hinges, tris, edgeMap, grabbed = -1, grabDepth = 0, grabW = 0, tornCount = 0;

// ------------------------------------------------------------------------------------------------
// Constraints
// ------------------------------------------------------------------------------------------------
function addDistance(list, a, b) {
  const c = { a, b, rest: dist(a, b), lambda: 0, alive: true };
  list.push(c);
  edgeMap.set(edgeKey(a, b), c);
  return c;
}
function edgeKey(a, b) { return a < b ? a * NP + b : b * NP + a; }
function dist(a, b) { return Math.hypot(x[a] - x[b], y[a] - y[b], z[a] - z[b]); }

// XPBD distance constraint C = |xa - xb| - rest, with compliance alpha (m/N); alphaTilde = alpha / h^2.
function solveDistance(c, alphaTilde) {
  if (!c.alive) return;
  const a = c.a, b = c.b, wsum = w[a] + w[b];
  if (wsum === 0) return;
  // TODO: the XPBD update from chapter 08, in 3D.
  //   C = |xa - xb| - c.rest, gradient wrt a is n = (xa - xb) / |xa - xb| (and -n wrt b)
  //   dlambda = (-C - alphaTilde * c.lambda) / (wsum + alphaTilde);  c.lambda += dlambda
  //   xa += w[a] * dlambda * n,  xb -= w[b] * dlambda * n
}

// Dihedral-angle bending across the edge e0-e1 shared by the triangles (e0, e1, w0) and (e1, e0, w1).
// The signed angle between the two triangle normals and its gradient follow Bridson et al. (2003).
const gx = new Float64Array(4), gy = new Float64Array(4), gz = new Float64Array(4);
function dihedralAngle(hg) {
  const i3 = hg.e0, i4 = hg.e1, i1 = hg.w0, i2 = hg.w1;
  const ex = x[i4] - x[i3], ey = y[i4] - y[i3], ez = z[i4] - z[i3];
  const el = Math.hypot(ex, ey, ez);
  // N1 = (x1 - x3) x (x1 - x4),  N2 = (x2 - x4) x (x2 - x3)
  const ax = x[i1] - x[i3], ay = y[i1] - y[i3], az = z[i1] - z[i3];
  const bx = x[i1] - x[i4], by = y[i1] - y[i4], bz = z[i1] - z[i4];
  const cx = x[i2] - x[i4], cy = y[i2] - y[i4], cz = z[i2] - z[i4];
  const dx = x[i2] - x[i3], dy = y[i2] - y[i3], dz = z[i2] - z[i3];
  const n1x = ay * bz - az * by, n1y = az * bx - ax * bz, n1z = ax * by - ay * bx;
  const n2x = cy * dz - cz * dy, n2y = cz * dx - cx * dz, n2z = cx * dy - cy * dx;
  const n1s = n1x * n1x + n1y * n1y + n1z * n1z, n2s = n2x * n2x + n2y * n2y + n2z * n2z;
  if (el < 1e-9 || n1s < 1e-14 || n2s < 1e-14) return null;
  // skip wrecked hinges (a triangle more than four times its rest area, e.g. next to a tear)
  if (hg.restN1 && (n1s > 16 * hg.restN1 || n2s > 16 * hg.restN2)) return null;
  hg.n1s = n1s; hg.n2s = n2s;
  // signed angle: cos from n1.n2, sin from (n1 x n2).e
  const l1 = Math.sqrt(n1s), l2 = Math.sqrt(n2s);
  const cosT = (n1x * n2x + n1y * n2y + n1z * n2z) / (l1 * l2);
  const crx = n1y * n2z - n1z * n2y, cry = n1z * n2x - n1x * n2z, crz = n1x * n2y - n1y * n2x;
  const sinT = (crx * ex + cry * ey + crz * ez) / (l1 * l2 * el);
  const theta = Math.atan2(sinT, cosT);
  // Gradients. Bridson et al. call these vectors the bending modes u1..u4; with our sign convention for
  // theta, d theta / d x_i = -u_i.
  const f1 = el / n1s, f2 = el / n2s;
  const a3 = (bx * ex + by * ey + bz * ez) / el, b3 = (cx * ex + cy * ey + cz * ez) / el;   // (x1-x4).E/|E|, (x2-x4).E/|E|
  const a4 = (ax * ex + ay * ey + az * ez) / el, b4 = (dx * ex + dy * ey + dz * ez) / el;   // (x1-x3).E/|E|, (x2-x3).E/|E|
  gx[0] = -f1 * n1x; gy[0] = -f1 * n1y; gz[0] = -f1 * n1z;                                   // -u1
  gx[1] = -f2 * n2x; gy[1] = -f2 * n2y; gz[1] = -f2 * n2z;                                   // -u2
  gx[2] = -(a3 * n1x / n1s + b3 * n2x / n2s);                                                // -u3
  gy[2] = -(a3 * n1y / n1s + b3 * n2y / n2s);
  gz[2] = -(a3 * n1z / n1s + b3 * n2z / n2s);
  gx[3] = a4 * n1x / n1s + b4 * n2x / n2s;                                                   // -u4
  gy[3] = a4 * n1y / n1s + b4 * n2y / n2s;
  gz[3] = a4 * n1z / n1s + b4 * n2z / n2s;
  return theta;
}

function solveDihedral(hg, alphaTilde) {
  if (!hg.alive) return;
  const theta = dihedralAngle(hg);
  if (theta === null) return;
  const idx = [hg.w0, hg.w1, hg.e0, hg.e1];
  let C = theta - hg.rest;
  if (C > Math.PI) C -= 2 * Math.PI; else if (C < -Math.PI) C += 2 * Math.PI;
  let denom = alphaTilde;
  for (let k = 0; k < 4; k++) denom += w[idx[k]] * (gx[k] * gx[k] + gy[k] * gy[k] + gz[k] * gz[k]);
  if (denom < 1e-12) return;
  const dlambda = (-C - alphaTilde * hg.lambda) / denom;
  hg.lambda += dlambda;
  for (let k = 0; k < 4; k++) {
    const p = idx[k];
    x[p] += w[p] * gx[k] * dlambda; y[p] += w[p] * gy[k] * dlambda; z[p] += w[p] * gz[k] * dlambda;
  }
}

// ------------------------------------------------------------------------------------------------
// Building the cloth
// ------------------------------------------------------------------------------------------------
function reset() {
  stretch = []; shear = []; bend2 = []; hinges = []; tris = []; edgeMap = new Map();
  grabbed = -1; tornCount = 0;
  const rng = Sim.random(3);
  const mass = AREAL_DENSITY * (N - 1) * SPACING * (N - 1) * SPACING / NP;
  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      const p = id(i, j), u = (i - (N - 1) / 2) * SPACING, v = j * SPACING;
      // both scenes start flat and horizontal; the curtain then swings down from its pinned edge
      x[p] = u; z[p] = v - (N - 1) * SPACING / 2;
      y[p] = (scene.value === CURTAIN ? 1.15 : 0.9) + 0.001 * (rng() - 0.5);
      vx[p] = vy[p] = vz[p] = 0;
      w[p] = 1 / mass;
    }
  }
  if (scene.value === CURTAIN) {
    if (pinning.value === 'two top corners') { w[id(0, 0)] = 0; w[id(N - 1, 0)] = 0; }
    else for (let i = 0; i < N; i++) w[id(i, 0)] = 0;
  }
  // stretch: the structural edges of the grid
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    if (i < N - 1) addDistance(stretch, id(i, j), id(i + 1, j));
    if (j < N - 1) addDistance(stretch, id(i, j), id(i, j + 1));
  }
  // shear: both diagonals of every square
  for (let j = 0; j < N - 1; j++) for (let i = 0; i < N - 1; i++) {
    // TODO: add the two diagonals of the square whose corner is (i, j) to the shear list, with
    // addDistance(shear, p, q), using id(i, j) for the particle in column i, row j.
  }
  // bending, cheap version: distance constraints that skip one particle
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    // TODO: when i < N - 2, call addBend2(a, b, c) for the three particles (i, j), (i + 1, j), (i + 2, j);
    // when j < N - 2, do the same along the column.
  }
  buildTrianglesAndHinges();
}

// A second-neighbour constraint a-c spans the two structural edges a-b and b-c; it dies when either tears.
function addBend2(a, b, c) {
  const k = { a, b: c, rest: dist(a, c), lambda: 0, alive: true, via1: edgeKey(a, b), via2: edgeKey(b, c) };
  bend2.push(k);
}

function buildTrianglesAndHinges() {
  // two triangles per square, split along the (i, j)-(i+1, j+1) diagonal
  for (let j = 0; j < N - 1; j++) for (let i = 0; i < N - 1; i++) {
    const p00 = id(i, j), p10 = id(i + 1, j), p11 = id(i + 1, j + 1), p01 = id(i, j + 1);
    tris.push({ a: p00, b: p11, c: p10, alive: true });
    tris.push({ a: p00, b: p01, c: p11, alive: true });
  }
  // every edge shared by two triangles becomes a hinge
  const owner = new Map();
  tris.forEach((t, ti) => {
    for (const [e0, e1, opp] of [[t.a, t.b, t.c], [t.b, t.c, t.a], [t.c, t.a, t.b]]) {
      const key = edgeKey(e0, e1);
      if (!owner.has(key)) { owner.set(key, { e0, e1, opp, ti }); continue; }
      const o = owner.get(key);
      const hg = { e0: o.e0, e1: o.e1, w0: o.opp, w1: opp, t0: o.ti, t1: ti, rest: 0, lambda: 0, alive: true };
      const th = dihedralAngle(hg);
      hg.rest = th === null ? 0 : th;
      hg.restN1 = hg.n1s; hg.restN2 = hg.n2s;
      hinges.push(hg);
    }
  });
}

function edgeAlive(a, b) {
  const c = edgeMap.get(edgeKey(a, b));
  return !c || c.alive;
}

// ------------------------------------------------------------------------------------------------
// Tearing: a constraint stretched beyond the threshold breaks, and whatever depended on it goes too.
// ------------------------------------------------------------------------------------------------
function tear() {
  for (const list of [stretch, shear]) {
    for (const c of list) {
      // TODO: if c is alive and its current length is more than tearAt.value times its rest length,
      // mark it dead (c.alive = false) and count it in tornCount.
    }
  }
  for (const k of bend2) if (k.alive) {
    const e1 = edgeMap.get(k.via1), e2 = edgeMap.get(k.via2);
    if ((e1 && !e1.alive) || (e2 && !e2.alive)) k.alive = false;
  }
  for (const t of tris) if (t.alive) t.alive = edgeAlive(t.a, t.b) && edgeAlive(t.b, t.c) && edgeAlive(t.c, t.a);
  for (const hg of hinges) if (hg.alive) hg.alive = tris[hg.t0].alive && tris[hg.t1].alive;
}

// ------------------------------------------------------------------------------------------------
// Collisions (with positional friction from chapter 16), wind and drag
// ------------------------------------------------------------------------------------------------
function frictionalPush(p, nx, ny, nz, depth) {
  x[p] += depth * nx; y[p] += depth * ny; z[p] += depth * nz;
  const dx = x[p] - px[p], dy = y[p] - py[p], dz = z[p] - pz[p];
  const dn = dx * nx + dy * ny + dz * nz;
  const tx = dx - dn * nx, ty = dy - dn * ny, tz = dz - dn * nz, lt = Math.hypot(tx, ty, tz);
  const k = lt < MU * depth ? 1 : Math.min(MU * depth / lt, 1);
  x[p] -= k * tx; y[p] -= k * ty; z[p] -= k * tz;
}

function collide() {
  for (let p = 0; p < NP; p++) {
    if (w[p] === 0) continue;
    if (scene.value === DRAPE) {
      const dx = x[p] - SPHERE.x, dy = y[p] - SPHERE.y, dz = z[p] - SPHERE.z, d = Math.hypot(dx, dy, dz);
      const r = SPHERE.r + THICK;
      if (d < r && d > 0) frictionalPush(p, dx / d, dy / d, dz / d, r - d);
    }
    if (y[p] < THICK) frictionalPush(p, 0, 1, 0, THICK - y[p]);
  }
}

function applyAir(h) {
  // linear drag on each triangle along its normal, relative to the wind; a third to each corner
  for (const t of tris) {
    if (!t.alive) continue;
    const a = t.a, b = t.b, c = t.c;
    const ux = x[b] - x[a], uy = y[b] - y[a], uz = z[b] - z[a];
    const qx = x[c] - x[a], qy = y[c] - y[a], qz = z[c] - z[a];
    let nx = uy * qz - uz * qy, ny = uz * qx - ux * qz, nz = ux * qy - uy * qx;
    const n2 = Math.hypot(nx, ny, nz);
    if (n2 < 1e-12) continue;
    const area = 0.5 * n2;
    nx /= n2; ny /= n2; nz /= n2;
    const rvx = (vx[a] + vx[b] + vx[c]) / 3, rvy = (vy[a] + vy[b] + vy[c]) / 3, rvz = (vz[a] + vz[b] + vz[c]) / 3 - wind.value;
    const vn = rvx * nx + rvy * ny + rvz * nz;
    for (const p of [a, b, c]) {
      // velocity change h * w * F / 3 with F = -DRAG * area * vn, capped so drag can at most stop the motion
      const k = Math.min(1, h * w[p] * DRAG * area / 3);
      vx[p] -= k * vn * nx; vy[p] -= k * vn * ny; vz[p] -= k * vn * nz;
    }
  }
}

// ------------------------------------------------------------------------------------------------
// Projection (orthographic): rotate by yaw about y, then by pitch about x. Screen = (X, Y), depth = Z.
// ------------------------------------------------------------------------------------------------
function view() {
  const yw = yawDeg.value * Math.PI / 180, pt = pitchDeg.value * Math.PI / 180;
  return { cy: Math.cos(yw), sy: Math.sin(yw), cp: Math.cos(pt), sp: Math.sin(pt) };
}
function project(V, X, Y, Z) {
  const x1 = V.cy * X + V.sy * Z, z1 = -V.sy * X + V.cy * Z;
  return { x: x1, y: V.cp * Y - V.sp * z1, z: V.sp * Y + V.cp * z1 };
}
function unproject(V, sx, sy, depth) {
  const Y = V.cp * sy + V.sp * depth, z1 = -V.sp * sy + V.cp * depth;
  return { x: V.cy * sx - V.sy * z1, y: Y, z: V.sy * sx + V.cy * z1 };
}

function handleMouse(V) {
  const m = Sim.mouse;
  if (m.pressed && grabbed < 0) {
    let best = -1, bestD = 0.06, bestZ = -Infinity;
    for (let p = 0; p < NP; p++) {
      const s = project(V, x[p], y[p], z[p]), d = Math.hypot(s.x - m.x, s.y - m.y);
      if (d < bestD - 0.01 || (d < bestD + 0.01 && d < 0.06 && s.z > bestZ)) { best = p; bestD = d; bestZ = s.z; }
    }
    if (best >= 0) { grabbed = best; grabDepth = project(V, x[best], y[best], z[best]).z; grabW = w[best]; w[best] = 0; }
  }
  if (!m.down && grabbed >= 0) { w[grabbed] = grabW; grabbed = -1; }
}

// ------------------------------------------------------------------------------------------------
// The step: XPBD with substeps, one pass over each constraint set per substep
// ------------------------------------------------------------------------------------------------
function step(dt) {
  const V = view();
  handleMouse(V);
  const h = dt / SUBSTEPS, h2 = h * h;
  const target = grabbed >= 0 ? unproject(V, Sim.mouse.x, Sim.mouse.y, grabDepth) : null;
  for (let s = 0; s < SUBSTEPS; s++) {
    applyAir(h);
    for (let p = 0; p < NP; p++) {
      px[p] = x[p]; py[p] = y[p]; pz[p] = z[p];
      if (w[p] === 0) continue;
      vy[p] -= G * h;
      x[p] += h * vx[p]; y[p] += h * vy[p]; z[p] += h * vz[p];
    }
    if (target) { x[grabbed] = target.x; y[grabbed] = target.y; z[grabbed] = target.z; }
    for (const list of [stretch, shear, bend2, hinges]) for (const c of list) c.lambda = 0;
    for (const c of stretch) solveDistance(c, stretchC.value / h2);
    for (const c of shear) solveDistance(c, shearC.value / h2);
    // the second-neighbour springs get the bending compliance times SPACING^2, which turns 1/(N m) into
    // m/N; a rough match that keeps the slider meaningful for both methods
    if (bending.value === 'second neighbours') for (const c of bend2) solveDistance(c, bendC.value * SPACING * SPACING / h2);
    if (bending.value === 'dihedral angle') for (const hg of hinges) solveDihedral(hg, bendC.value / h2);
    collide();
    for (let p = 0; p < NP; p++) {
      vx[p] = (x[p] - px[p]) / h; vy[p] = (y[p] - py[p]) / h; vz[p] = (z[p] - pz[p]) / h;
    }
    if (tearing.value) tear();
  }
}

// ------------------------------------------------------------------------------------------------
// Drawing: painter's algorithm over the triangles, shaded by how much they face the light
// ------------------------------------------------------------------------------------------------
function draw(g) {
  const V = view();
  const sp = new Array(NP);
  for (let p = 0; p < NP; p++) sp[p] = project(V, x[p], y[p], z[p]);
  const items = [];
  if (scene.value === DRAPE) {
    const f = [[-0.8, -0.8], [0.8, -0.8], [0.8, 0.8], [-0.8, 0.8]].map(([a, b]) => project(V, a, 0, b));
    g.polyline(f, { fill: '#e3e0d6', closed: true });
    const c = project(V, SPHERE.x, SPHERE.y, SPHERE.z);
    items.push({ z: c.z, sphere: c });
  }
  const drawn = new Uint8Array(NP);
  for (const t of tris) {
    if (!t.alive) continue;
    // skip triangles stretched to more than three times their size: they are no longer cloth
    if (dist(t.a, t.b) > 3 * SPACING * Math.SQRT2 || dist(t.b, t.c) > 3 * SPACING * Math.SQRT2 || dist(t.c, t.a) > 3 * SPACING * Math.SQRT2) continue;
    drawn[t.a] = drawn[t.b] = drawn[t.c] = 1;
    const A = sp[t.a], B = sp[t.b], C = sp[t.c];
    const ux = B.x - A.x, uy = B.y - A.y, uz = B.z - A.z, qx = C.x - A.x, qy = C.y - A.y, qz = C.z - A.z;
    let nx = uy * qz - uz * qy, ny = uz * qx - ux * qz, nz = ux * qy - uy * qx;
    const l = Math.hypot(nx, ny, nz) || 1;
    nx /= l; ny /= l; nz /= l;
    const light = Math.abs(-0.55 * nx + 0.6 * ny + 0.58 * nz);
    items.push({ z: (A.z + B.z + C.z) / 3, A, B, C, front: nz > 0, light });
  }
  items.sort((p, q) => p.z - q.z);
  for (const it of items) {
    if (it.sphere) {
      g.circle(it.sphere.x, it.sphere.y, SPHERE.r, { fill: '#c9c4b5' });
      g.circle(it.sphere.x - 0.06, it.sphere.y + 0.07, SPHERE.r * 0.45, { fill: '#d8d4c7' });
      continue;
    }
    const k = 0.45 + 0.55 * it.light;
    const col = it.front ? `rgb(${Math.round(70 * k)},${Math.round(120 * k)},${Math.round(225 * k)})`
      : `rgb(${Math.round(230 * k)},${Math.round(140 * k)},${Math.round(60 * k)})`;
    g.polyline([it.A, it.B, it.C], { fill: col, stroke: col, width: 0.6, closed: true });
  }
  for (let p = 0; p < NP; p++) if (!drawn[p]) g.dot(sp[p].x, sp[p].y, 1.5, Sim.colors.blue);   // loose particles
  for (let p = 0; p < NP; p++) if (w[p] === 0 && p !== grabbed) g.dot(sp[p].x, sp[p].y, 4, Sim.colors.ink);
  if (grabbed >= 0) g.dot(sp[grabbed].x, sp[grabbed].y, 5, Sim.colors.red);

  let maxStrain = 0;
  for (const c of stretch) if (c.alive) maxStrain = Math.max(maxStrain, dist(c.a, c.b) / c.rest - 1);
  Sim.show('particles', NP);
  Sim.show('largest stretch', (100 * maxStrain).toFixed(2) + ' %');
  if (tearing.value) Sim.show('torn constraints', tornCount);
  // a thin plate with these two compliances would be t = sqrt(12 (1 - nu^2) alpha_stretch / alpha_bend) thick
  if (bending.value === 'dihedral angle') {
    const t = Math.sqrt(12 * (1 - 0.3 * 0.3) * stretchC.value / bendC.value);
    Sim.show('equivalent shell thickness', t < 0.01 ? (1000 * t).toFixed(3) + ' mm' : (100 * t).toFixed(2) + ' cm');
  }
  Sim.show('mouse', 'drag the cloth');
}

Sim.run({ title: 'Cloth (starter)', view: [-1.0, -0.35, 1.0, 1.3], reset, step, draw, grid: false });
