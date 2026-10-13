// Chapter 16 - Contact and Friction (starter)
// Position-based contact with Coulomb friction, in three scenes:
//   tilting slope:      a block on a slope you tilt; it starts to slide when tan(theta) = mu_s
//   stick-slip:         a block on the floor pulled through a spring whose far end moves slowly
//   rope over cylinder: a particle rope over a fixed cylinder with a weight on each end; it holds while
//                       the heavier weight is at most exp(mu_s * wrap angle) times the lighter one
// Your jobs, marked TODO: positional friction (static and kinetic) in solveContact, and, as a challenge,
// restitution in applyRestitution. Until then every contact is frictionless and nothing bounces.

const SLOPE = 'tilting slope', STICK = 'stick-slip', ROPE = 'rope over cylinder';
const mode = Sim.select('scene', [SLOPE, STICK, ROPE], SLOPE, { reset: true });
const two = (v) => v.toFixed(2);
const muS = Sim.slider('static friction mu_s', 0, 1.2, 0.5, { step: 0.01, format: two });
const muK = Sim.slider('kinetic friction mu_k', 0, 1.2, 0.4, { step: 0.01, format: two });
const angleDeg = Sim.slider('slope angle (deg)', 0, 60, 15, { step: 0.1, format: (v) => v.toFixed(1) });
const restitution = Sim.slider('restitution e', 0, 1, 0, { step: 0.05, format: two });
const wrap = Sim.select('wrap angle (rope)', ['half turn (pi)', 'one and a half turns (3 pi)'], 'half turn (pi)', { reset: true });
const massL = Sim.slider('left weight (kg)', 0.05, 3, 0.5, { step: 0.01, format: two });
const massR = Sim.slider('right weight (kg)', 0.05, 3, 0.5, { step: 0.01, format: two });
Sim.button('slope: tilt slowly from 0 until it slips', () => { if (mode.value === SLOPE) { angleDeg.set(0); Sim.reset(); ramping = true; } });
Sim.button('slope: drop the block', () => { if (mode.value === SLOPE) dropBlock(); });
Sim.button('rope: raise the right weight until it slips', () => { if (mode.value === ROPE) { Sim.reset(); ramping = true; } });

const G = 9.81;
const SUBSTEPS = 20;
const ROPE_ITERATIONS = 16;         // the rope is light and its weights heavy, so it needs a few passes
const ROPE_PARTICLE_MASS = 0.001;   // kg per rope particle (0.05 kg per metre)
const MAX_N = 200;

// Particles: position, position at the start of the substep, velocity, inverse mass.
const x = new Float64Array(MAX_N), y = new Float64Array(MAX_N);
const px = new Float64Array(MAX_N), py = new Float64Array(MAX_N);
const vx = new Float64Array(MAX_N), vy = new Float64Array(MAX_N);
const w = new Float64Array(MAX_N);
// Contacts found this substep: normal and the normal velocity before the substep (for restitution).
const inContact = new Uint8Array(MAX_N);
const cnx = new Float64Array(MAX_N), cny = new Float64Array(MAX_N), vnBefore = new Float64Array(MAX_N);
let n = 0;

// Rope: distance constraints between neighbours, and the tension each one carried in the last substep.
const restLen = 0.02;
let tension = new Float64Array(MAX_N), leftTop = 0, rightTop = 0;

// Scene state
const BLOCK_HW = 0.1, BLOCK_HH = 0.06;          // block half width and half height (m)
const CYL_R = 0.2, ROPE_R = 0.008;               // cylinder radius, rope thickness radius (m)
const SPRING_K = 100, SPRING_REST = 0.4, DRIVER_SPEED = 0.02;   // stick-slip spring (N/m, m, m/s)
let ramping = false, slipValue = null, driverX = 0, springForce = 0, peakF = 0, troughF = Infinity;
let sliding = false, aAlong = 0, vAlongPrev = 0, maxHeight = 0, dropped = false, landed = false;

// ------------------------------------------------------------------------------------------------
// The contact: particle i against a static surface with unit normal (nx, ny), penetrating by depth > 0.
// ------------------------------------------------------------------------------------------------
function solveContact(i, nx, ny, depth) {
  // normal part: push the particle out of the surface
  x[i] += depth * nx;
  y[i] += depth * ny;
  // tangential motion of the particle over this substep (the surface does not move)
  const dx = x[i] - px[i], dy = y[i] - py[i];
  const dn = dx * nx + dy * ny;
  const tx = dx - dn * nx, ty = dy - dn * ny;
  const lt = Math.hypot(tx, ty);
  // TODO: positional friction. The normal correction was depth; the tangential motion is (tx, ty), of
  // length lt.
  //   static:  if lt < mu_s * depth, cancel all of the tangential motion.
  //   kinetic: otherwise remove only mu_k * depth of it (never more than lt itself).
  // Use Math.min(muK.value, muS.value) for mu_k: kinetic friction never exceeds static.
  inContact[i] = 1; cnx[i] = nx; cny[i] = ny;
}

// Restitution, as a velocity pass after the PBD velocity update (Mueller et al. 2020, section 3.5).
function applyRestitution(i, h) {
  // TODO (challenge): vn = current normal velocity along (cnx[i], cny[i]); vnBefore[i] = normal velocity
  // before the substep (negative when approaching). Use e = restitution.value, or e = 0 when
  // |vnBefore| <= 2 g h (resting contact). Change the normal velocity from vn to max(-e * vnBefore, 0).
}

// Distance constraint between particles a and b (zero compliance). Returns the multiplier lambda.
function solveDistance(a, b, rest) {
  const dx = x[b] - x[a], dy = y[b] - y[a], d = Math.hypot(dx, dy);
  const ws = w[a] + w[b];
  if (d === 0 || ws === 0) return 0;
  const C = d - rest, lambda = -C / ws, nx = dx / d, ny = dy / d;
  x[a] -= w[a] * lambda * nx; y[a] -= w[a] * lambda * ny;
  x[b] += w[b] * lambda * nx; y[b] += w[b] * lambda * ny;
  return lambda;
}

// ------------------------------------------------------------------------------------------------
// Scenes
// ------------------------------------------------------------------------------------------------
function slopeFrame() {
  const t = angleDeg.value * Math.PI / 180;
  return { tx: Math.cos(t), ty: Math.sin(t), nx: -Math.sin(t), ny: Math.cos(t) };
}

function dropBlock() {
  const f = slopeFrame();
  const along = x[0] * f.tx + y[0] * f.ty;
  x[0] = along * f.tx + (BLOCK_HH + 0.3) * f.nx; y[0] = along * f.ty + (BLOCK_HH + 0.3) * f.ny;
  vx[0] = 0; vy[0] = 0; maxHeight = 0; dropped = true; landed = false;
}

function reset() {
  n = 0; slipValue = null; sliding = false; aAlong = 0; vAlongPrev = 0;
  peakF = 0; troughF = Infinity; dropped = false; maxHeight = 0;
  ramping = false;
  if (mode.value === SLOPE) {
    Sim.setView([-1.0, -0.8, 1.0, 0.9]);
    const f = slopeFrame();
    x[0] = BLOCK_HH * f.nx; y[0] = BLOCK_HH * f.ny; vx[0] = 0; vy[0] = 0; w[0] = 1;
    n = 1;
  } else if (mode.value === STICK) {
    Sim.setView([-0.3, -0.3, 1.3, 0.6]);
    x[0] = 0; y[0] = BLOCK_HH; vx[0] = 0; vy[0] = 0; w[0] = 1;
    driverX = SPRING_REST;
    n = 1;
  } else {
    Sim.setView([-0.75, -0.85, 0.75, 0.35]);
    buildRope();
  }
}

function buildRope() {
  const r = CYL_R + ROPE_R, hang = 0.4;
  const wrapAngle = wrap.value.startsWith('half') ? Math.PI : 3 * Math.PI;
  const nHang = Math.round(hang / restLen), nArc = Math.round(wrapAngle * r / restLen);
  const arcStep = wrapAngle / nArc;
  n = 0;
  // left hanging part, from the weight upwards
  for (let k = 0; k < nHang; k++) { x[n] = -r; y[n] = -hang + k * restLen; n++; }
  leftTop = n - 1;
  // the wrapped part, clockwise from the left side over the top
  for (let k = 0; k <= nArc; k++) {
    const phi = Math.PI - k * arcStep;
    x[n] = r * Math.cos(phi); y[n] = r * Math.sin(phi); n++;
  }
  rightTop = n - 1;
  // right hanging part, downwards to the weight
  for (let k = 1; k <= nHang; k++) { x[n] = r; y[n] = -k * restLen; n++; }
  // the arc spacing differs very slightly from the straight spacing; use each link's own rest length
  ropeRest = new Float64Array(n - 1);
  for (let k = 0; k < n - 1; k++) ropeRest[k] = Math.hypot(x[k + 1] - x[k], y[k + 1] - y[k]);
  for (let k = 0; k < n; k++) { vx[k] = 0; vy[k] = 0; w[k] = 1 / ROPE_PARTICLE_MASS; }
  setWeights();
}
let ropeRest = new Float64Array(1);

function setWeights() {
  w[0] = 1 / massL.value;
  w[n - 1] = 1 / massR.value;
}

function collide(i) {
  if (mode.value === ROPE) {
    const dx = x[i], dy = y[i], d = Math.hypot(dx, dy), r = CYL_R + ROPE_R;
    if (d < r && d > 0) solveContact(i, dx / d, dy / d, r - d);
  } else if (mode.value === SLOPE) {
    const f = slopeFrame();
    const along = x[i] * f.tx + y[i] * f.ty;
    const gap = x[i] * f.nx + y[i] * f.ny - BLOCK_HH;
    if (gap < 0 && Math.abs(along) < 1.2) solveContact(i, f.nx, f.ny, -gap);
  } else {
    const gap = y[i] - BLOCK_HH;
    if (gap < 0) solveContact(i, 0, 1, -gap);
  }
}

function step(dt) {
  const h = dt / SUBSTEPS;
  if (mode.value === ROPE) {
    if (ramping) massR.set(Math.min(3, massR.value * (1 + 0.1 * dt)));   // +10 % per second
    setWeights();
  }
  if (mode.value === SLOPE && ramping) angleDeg.set(Math.min(60, angleDeg.value + 3 * dt));   // 3 degrees per second
  if (mode.value === STICK) driverX += DRIVER_SPEED * dt;

  for (let s = 0; s < SUBSTEPS; s++) {
    for (let i = 0; i < n; i++) {
      inContact[i] = 0;
      vnBefore[i] = 0;
      if (w[i] === 0) continue;
      px[i] = x[i]; py[i] = y[i];
      let fx = 0;
      if (mode.value === STICK) { springForce = SPRING_K * (driverX - x[i] - SPRING_REST); fx = springForce; }
      vx[i] += h * fx * w[i];
      vy[i] -= h * G;
      x[i] += h * vx[i]; y[i] += h * vy[i];
    }
    if (mode.value === ROPE) {
      for (let k = 0; k < n - 1; k++) tension[k] = 0;
      for (let it = 0; it < ROPE_ITERATIONS; it++)
        for (let k = 0; k < n - 1; k++) tension[k] -= solveDistance(k, k + 1, ropeRest[k]) / (h * h);
    }
    for (let i = 0; i < n; i++) collide(i);
    if (dropped && inContact[0]) landed = true;
    for (let i = 0; i < n; i++) {
      // velocity before the update, for the restitution pass: from the previous substep's motion
      const ovx = vx[i], ovy = vy[i];
      vx[i] = (x[i] - px[i]) / h; vy[i] = (y[i] - py[i]) / h;
      if (inContact[i]) {
        vnBefore[i] = ovx * cnx[i] + ovy * cny[i];
        applyRestitution(i, h);
      }
    }
  }
  measure(dt);
}

function measure(dt) {
  if (mode.value === SLOPE) {
    const f = slopeFrame();
    const vAlong = vx[0] * f.tx + vy[0] * f.ty;
    aAlong += 0.1 * ((vAlong - vAlongPrev) / dt - aAlong);
    vAlongPrev = vAlong;
    sliding = Math.abs(vAlong) > 0.005;
    if (ramping && sliding) { ramping = false; slipValue = angleDeg.value; }
    if (landed) maxHeight = Math.max(maxHeight, x[0] * f.nx + y[0] * f.ny - BLOCK_HH);   // highest point after the first impact
  } else if (mode.value === STICK) {
    const nowSliding = Math.abs(vx[0]) > 1e-3;
    if (nowSliding && !sliding) peakF = springForce;      // the force at which it broke free
    if (!nowSliding && sliding) troughF = springForce;    // the force at which it stuck again
    sliding = nowSliding;
  } else {
    // The rope slips as a whole when the wrap starts moving on the LIGHT side; before that, only the
    // elastic stretch near the heavy side creeps.
    const k = massL.value <= massR.value ? leftTop + 1 : rightTop - 1;
    sliding = Math.hypot(vx[k], vy[k]) > 0.002;
    if (ramping && sliding) { ramping = false; slipValue = ropeTensionRatio(); }
  }
}

// ------------------------------------------------------------------------------------------------
// Drawing
// ------------------------------------------------------------------------------------------------
function boxPoly(cx, cy, hw, hh, tx, ty) {
  const nx = -ty, ny = tx;
  return [{ x: cx - hw * tx - hh * nx, y: cy - hw * ty - hh * ny }, { x: cx + hw * tx - hh * nx, y: cy + hw * ty - hh * ny },
    { x: cx + hw * tx + hh * nx, y: cy + hw * ty + hh * ny }, { x: cx - hw * tx + hh * nx, y: cy - hw * ty + hh * ny }];
}

function drawSlope(g) {
  const f = slopeFrame(), L = 1.2, T = 0.05;
  g.polyline([{ x: -L * f.tx, y: -L * f.ty }, { x: L * f.tx, y: L * f.ty },
    { x: L * f.tx - T * f.nx, y: L * f.ty - T * f.ny }, { x: -L * f.tx - T * f.nx, y: -L * f.ty - T * f.ny }],
  { fill: '#cfcabb', closed: true });
  g.circle(0, -0.03, 0.025, { fill: Sim.colors.gray });
  g.polyline(boxPoly(x[0], y[0], BLOCK_HW, BLOCK_HH, f.tx, f.ty), { fill: sliding ? Sim.colors.red : Sim.colors.blue, closed: true });
  // gravity split into its normal and tangential parts, drawn from the block
  const s = 0.03;
  g.arrow(x[0], y[0], 0, -G * s, { color: Sim.colors.ink, width: 1.5 });
  const gn = -G * Math.cos(angleDeg.value * Math.PI / 180), gt = -G * Math.sin(angleDeg.value * Math.PI / 180);
  g.arrow(x[0], y[0], gt * f.tx * s, gt * f.ty * s, { color: Sim.colors.orange, width: 1.5 });
  g.arrow(x[0], y[0], gn * f.nx * s, gn * f.ny * s, { color: Sim.colors.teal, width: 1.5 });

  const th = angleDeg.value * Math.PI / 180;
  Sim.show('angle', angleDeg.value.toFixed(1) + ' deg');
  Sim.show('tan(angle)', Math.tan(th).toFixed(3));
  Sim.show('mu_s', muS.value.toFixed(2) + '   (slips at ' + (Math.atan(muS.value) * 180 / Math.PI).toFixed(1) + ' deg)');
  Sim.show('block', sliding ? 'sliding' : 'sticking');
  if (slipValue !== null) Sim.show('slipped at', slipValue.toFixed(1) + ' deg');
  if (sliding && !dropped) {
    Sim.show('acceleration along slope', Math.abs(aAlong).toFixed(2) + ' m/s^2');
    Sim.show('g (sin - mu_k cos)', (G * (Math.sin(th) - muK.value * Math.cos(th))).toFixed(2) + ' m/s^2');
  }
  if (dropped) {
    Sim.show('dropped from', '0.300 m');
    Sim.show('highest bounce', maxHeight.toFixed(3) + ' m');
    Sim.show('e^2 x 0.3 m', (restitution.value * restitution.value * 0.3).toFixed(3) + ' m');
  }
}

function drawStick(g) {
  g.ground(0);
  g.polyline(boxPoly(x[0], y[0], BLOCK_HW, BLOCK_HH, 1, 0), { fill: sliding ? Sim.colors.red : Sim.colors.blue, closed: true });
  // spring from the block to the driver
  const x1 = x[0] + BLOCK_HW, x2 = driverX + BLOCK_HW, pts = [], coils = 14;
  for (let k = 0; k <= coils; k++) pts.push({ x: x1 + (x2 - x1) * k / coils, y: BLOCK_HH + (k === 0 || k === coils ? 0 : (k % 2 ? 0.02 : -0.02)) });
  g.polyline(pts, { color: Sim.colors.gray, width: 1.5 });
  g.circle(x2, BLOCK_HH, 0.015, { fill: Sim.colors.ink });
  const mg = G / w[0];
  Sim.plot('spring force / mg', springForce / mg, { color: Sim.colors.purple, min: 0, max: 0.8 });
  Sim.show('spring force / mg', (springForce / mg).toFixed(3));
  if (isFinite(troughF)) {
    Sim.show('last slip at (expect mu_s)', (peakF / mg).toFixed(3));
    Sim.show('stuck again at (expect 2 mu_k - mu_s)', (troughF / mg).toFixed(3));
  }
}

// Holding ratio actually carried by the wrap: tension where the rope leaves the cylinder, heavy side
// over light side. It differs a little from the weight ratio because the hanging rope weighs something.
function ropeTensionRatio() {
  const tl = tension[leftTop - 1], tr = tension[rightTop];
  return Math.max(tl, tr) / Math.min(tl, tr);
}

function drawRope(g) {
  g.circle(0, 0, CYL_R, { fill: '#cfcabb' });
  g.circle(0, 0, 0.01, { fill: Sim.colors.gray });
  // Draw the wrapped turns spread out a little, so that one and a half turns do not hide each other.
  // (Only the picture is spread; the simulated rope lies on the cylinder.)
  const pts = [], spread = wrap.value.startsWith('half') ? 0 : 0.025;
  for (let k = 0; k < n; k++) {
    const f = k <= leftTop ? 0 : k >= rightTop ? 1 : (k - leftTop) / (rightTop - leftTop);
    const d = Math.hypot(x[k], y[k]), s = k > leftTop && k < rightTop && d > 0 ? 1 + spread * f / d : 1;
    pts.push({ x: x[k] * s + (k >= rightTop ? spread : 0), y: y[k] * s });
  }
  g.polyline(pts, { color: Sim.colors.ink, width: 2.5 });
  const sL = 0.05 * Math.sqrt(massL.value), sR = 0.05 * Math.sqrt(massR.value);
  g.rect(x[0] - sL, y[0] - 2 * sL, 2 * sL, 2 * sL, { fill: Sim.colors.blue });
  g.rect(x[n - 1] + spread - sR, y[n - 1] - 2 * sR, 2 * sR, 2 * sR, { fill: Sim.colors.red });
  const wrapAngle = wrap.value.startsWith('half') ? Math.PI : 3 * Math.PI;
  const big = Math.max(massL.value, massR.value), small = Math.min(massL.value, massR.value);
  Sim.show('weight ratio (heavy / light)', (big / small).toFixed(3));
  Sim.show('tension where it leaves, left', tension[leftTop - 1].toFixed(2) + ' N');
  Sim.show('tension where it leaves, right', tension[rightTop].toFixed(2) + ' N');
  Sim.show('tension ratio', ropeTensionRatio().toFixed(3));
  Sim.show('exp(mu_s x wrap angle)', Math.exp(muS.value * wrapAngle).toFixed(3));
  if (muK.value < muS.value) Sim.show('note', 'mu_k < mu_s: the rope creeps; expect a slip between exp(mu_k x wrap) and exp(mu_s x wrap)');
  Sim.show('rope', sliding ? 'slipping' : 'holding');
  if (slipValue !== null) Sim.show('slipped at tension ratio', slipValue.toFixed(3));
}

function draw(g) {
  if (mode.value === SLOPE) drawSlope(g);
  else if (mode.value === STICK) drawStick(g);
  else drawRope(g);
}

Sim.run({ title: 'Contact and friction (starter)', view: [-1.0, -0.8, 1.0, 0.9], reset, step, draw });
