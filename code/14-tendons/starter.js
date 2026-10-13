// Chapter 14 - Cable-Driven Bending (starter)
// A tentacle built from chapter 13's rod (rigid segments with angles, base clamped), with two cables (tendons)
// routed through guides at offset r on the left and right of every segment. A knob of radius r_p pulls one cable
// and lets out the other. Each cable is a unilateral XPBD constraint on the length of the polyline through its
// guides: it can pull (tension) but never push. Right now the cables do nothing: fill in the TODOs in cableLength
// and solveCable (exercise 1), allowedLength (slack, exercise 2) and capstanFactor (sheath friction, exercise 3).

const knob = Sim.slider('knob angle (deg)', -180, 180, 0, { step: 1 });
const pulleyMm = Sim.slider('pulley radius r_p (mm)', 5, 30, 15, { step: 0.5 });
const slackMm = Sim.slider('slack per cable (mm)', 0, 10, 0, { step: 0.1 });
const friction = Sim.slider('sheath friction mu', 0, 0.4, 0, { step: 0.01 });
const curlDeg = Sim.slider('curl the base: sheath wrap (deg)', 0, 720, 0, { step: 5 });
const cableEA = Sim.logSlider('cable stiffness EA (N)', 2, 6, 1e5);
const sweep = Sim.checkbox('sweep the knob slowly', false);
const gravityOn = Sim.checkbox('gravity (off: floating in water)', false);
Sim.button('clear graph', () => { trail = []; });

const L = 0.4;          // tentacle length (m)
const N = 16;           // segments
const l = L / N;
const r = 0.015;        // guide offset from the centreline (m)
const EI = 0.05;        // backbone bending stiffness (N m^2)
const mu = 0.2;         // mass per length (kg/m)
const sheathLength = 1.0; // cable length inside the sheath, from knob to base (m)
const gravity = 9.81;
const damping = 4.0;
const iterations = 2;

// Segment 0 is the clamped base (kinematic), pointing up; segments 1..N are free.
let x, y, th, px, py, pth, vx, vy, om, w, wI, lambdaJoint, lambdaBend, time, trail;
const lambdaCable = { left: 0, right: 0 };
const tension = { left: 0, right: 0 };
const gradX = new Float64Array(N + 1), gradY = new Float64Array(N + 1);
const guideX = new Float64Array(N + 1), guideY = new Float64Array(N + 1);
const restCableLength = N * l;   // guide-to-guide length of a straight tentacle

function reset() {
  const n = N + 1;
  x = new Float64Array(n); y = new Float64Array(n); th = new Float64Array(n);
  px = new Float64Array(n); py = new Float64Array(n); pth = new Float64Array(n);
  vx = new Float64Array(n); vy = new Float64Array(n); om = new Float64Array(n);
  w = new Float64Array(n); wI = new Float64Array(n);
  lambdaJoint = new Float64Array(n); lambdaBend = new Float64Array(n);
  const m = mu * l, inertia = m * l * l / 12;
  for (let i = 0; i < n; i++) {
    x[i] = 0; y[i] = (i - 0.5) * l; th[i] = Math.PI / 2;
    if (i > 0) { w[i] = 1 / m; wI[i] = 1 / inertia; }
  }
  time = 0;
  trail = [];
}

// ---------------------------------------------------------------------------------------------
// The rod from chapter 13: rigid stretch-shear joints and compliant bend joints.
// ---------------------------------------------------------------------------------------------
function jointLength(i) { return i === 0 ? 0.5 * l : l; }

function solveStretchShear(i) {
  const j = i + 1;
  const rax = 0.5 * l * Math.cos(th[i]), ray = 0.5 * l * Math.sin(th[i]);
  const rbx = -0.5 * l * Math.cos(th[j]), rby = -0.5 * l * Math.sin(th[j]);
  const dx = (x[j] + rbx) - (x[i] + rax), dy = (y[j] + rby) - (y[i] + ray);
  const C = Math.hypot(dx, dy);
  if (C < 1e-12) return;
  const nx = dx / C, ny = dy / C;
  const rna = rax * ny - ray * nx, rnb = rbx * ny - rby * nx;
  const wa = w[i] + wI[i] * rna * rna, wb = w[j] + wI[j] * rnb * rnb;
  const dlambda = -C / (wa + wb);
  const pX = dlambda * nx, pY = dlambda * ny;
  x[j] += w[j] * pX; y[j] += w[j] * pY; th[j] += wI[j] * (rbx * pY - rby * pX);
  x[i] -= w[i] * pX; y[i] -= w[i] * pY; th[i] -= wI[i] * (rax * pY - ray * pX);
}

function solveBend(i, h) {
  const j = i + 1, len = jointLength(i);
  const C = th[j] - th[i];
  const alphaTilde = (len / EI) / (h * h);
  const dlambda = (-C - alphaTilde * lambdaBend[i]) / (wI[i] + wI[j] + alphaTilde);
  lambdaBend[i] += dlambda;
  th[i] -= wI[i] * dlambda;
  th[j] += wI[j] * dlambda;
}

// ---------------------------------------------------------------------------------------------
// Tendons
// ---------------------------------------------------------------------------------------------
// side = +1 for the left cable, -1 for the right. Guide k sits on segment k at offset side * r along the
// segment's normal n_k = (-sin th, cos th). Guide 0 is on the base, where the cable leaves its sheath.
function computeGuides(side) {
  for (let k = 0; k <= N; k++) {
    guideX[k] = x[k] - side * r * Math.sin(th[k]);
    guideY[k] = y[k] + side * r * Math.cos(th[k]);
  }
}

function cableLength(side) {
  computeGuides(side);
  // TODO: the cable runs straight from guide to guide. Return the length of the polyline through
  // (guideX[k], guideY[k]) for k = 0..N.
  return 0;
}

// Cable pulled in at the knob: u = r_p * phi. Positive knob angle pulls the left cable, lets out the right.
function knobPull() { return (pulleyMm.value / 1000) * knob.value * Math.PI / 180; }

// The length of cable available between the base and the tip. Pulling shortens it; slack lengthens it.
function allowedLength(side) {
  const u = side * knobPull();
  // TODO (exercise 2): each cable has slackMm.value millimetres of extra length.
  return restCableLength - u;
}

// Compliance of the cable seen from the tentacle: its own stretch inside the body, plus the stretch of the cable
// inside the sheath. With friction the tension in the sheath grows towards the knob as exp(mu * s * beta / Ls),
// so the sheath part stretches (e^(mu beta) - 1) / (mu beta) times more than a frictionless cable would.
function capstanFactor() {
  const a = friction.value * curlDeg.value * Math.PI / 180;   // mu * wrap angle
  // TODO (exercise 3): return (e^a - 1) / a, which tends to 1 as a goes to 0.
  return 1;
}
function cableCompliance() {
  return (restCableLength + sheathLength * capstanFactor()) / cableEA.value;
}

// Unilateral XPBD constraint C = length - allowed <= 0. Tension only: the accumulated multiplier never goes above 0.
function solveCable(side, h) {
  const key = side > 0 ? 'left' : 'right';
  const C = cableLength(side) - allowedLength(side);  // also fills guideX/guideY
  // TODO:
  // 1. If the cable is slack (C <= 0) and has no tension yet (lambdaCable[key] === 0), return.
  // 2. Gradient of the length with respect to guide k: (unit vector of the piece arriving at k) minus
  //    (unit vector of the piece leaving k). Accumulate it in gradX[k], gradY[k] (fill them with 0 first).
  // 3. Guide k (k >= 1) is a point on segment k at offset (guideX[k] - x[k], guideY[k] - y[k]), so its generalised
  //    inverse mass is w[k] |grad_k|^2 + wI[k] (r_k x grad_k)^2. Sum over k = 1..N (the base, k = 0, is fixed).
  // 4. XPBD with alphaTilde = cableCompliance() / h^2, then clamp: the accumulated lambdaCable[key] must stay <= 0
  //    (a cable can pull but not push). Apply the clamped dlambda to x[k], y[k] and th[k].
}

// Keep h * omega <= 0.25 for the stiffest bend joint (chapter 13 used 0.5; the cable pulls on every joint at once).
const substeps = Math.ceil(Math.sqrt(2 * 24 * EI / (mu * l * l * l * l)) / 60 / 0.25);

function step(dt) {
  if (sweep.value) knob.set(Math.round(160 * Math.sin(2 * Math.PI * time / 24)));
  const h = dt / substeps, decay = 1 / (1 + damping * h);
  for (let s = 0; s < substeps; s++) {
    time += h;
    for (let i = 1; i <= N; i++) {
      if (gravityOn.value) vy[i] -= h * gravity;
      vx[i] *= decay; vy[i] *= decay; om[i] *= decay;
      px[i] = x[i]; py[i] = y[i]; pth[i] = th[i];
      x[i] += h * vx[i]; y[i] += h * vy[i]; th[i] += h * om[i];
    }
    lambdaJoint.fill(0); lambdaBend.fill(0);
    lambdaCable.left = 0; lambdaCable.right = 0;
    for (let it = 0; it < iterations; it++) {
      for (let i = 0; i < N; i++) solveStretchShear(i);
      for (let i = 0; i < N; i++) solveBend(i, h);
      solveCable(+1, h);
      solveCable(-1, h);
    }
    for (let i = 1; i <= N; i++) {
      vx[i] = (x[i] - px[i]) / h; vy[i] = (y[i] - py[i]) / h; om[i] = (th[i] - pth[i]) / h;
    }
    tension.left = -lambdaCable.left / (h * h);
    tension.right = -lambdaCable.right / (h * h);
  }
  trail.push({ knob: knob.value, tip: tipAngleDeg() });
  if (trail.length > 1500) trail.shift();
}

function tipAngleDeg() { return (th[N] - th[0]) * 180 / Math.PI; }

// The ideal constant-curvature answer: theta = (u - slack) / r, zero inside the dead band.
function idealTipDeg(knobDeg) {
  const u = Math.abs((pulleyMm.value / 1000) * knobDeg * Math.PI / 180);
  const engaged = Math.max(0, u - slackMm.value / 1000);
  return Math.sign(knobDeg) * engaged / r * 180 / Math.PI;
}

// ---------------------------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------------------------
function drawSheathAndKnob(g) {
  // The sheath leaves the base downwards and curls by the wrap angle on its way to the knob.
  const beta = curlDeg.value * Math.PI / 180, straight = 0.06, curled = 0.25, steps = 80;
  const pts = [{ x: 0, y: -l }];
  let a = -Math.PI / 2, sx = 0, sy = -l - straight;
  pts.push({ x: sx, y: sy });
  for (let i = 0; i < steps; i++) {
    a += beta / steps;
    sx += (curled / steps) * Math.cos(a); sy += (curled / steps) * Math.sin(a);
    pts.push({ x: sx, y: sy });
  }
  sx += 0.03 * Math.cos(a); sy += 0.03 * Math.sin(a);
  pts.push({ x: sx, y: sy });
  g.polyline(pts, { color: '#6d6a62', width: 7, alpha: 0.85 });
  g.polyline(pts, { color: '#d8d4ca', width: 2 });
  const kr = 0.035, phi = knob.value * Math.PI / 180;
  const kx = sx + kr * Math.cos(a), ky = sy + kr * Math.sin(a);
  g.circle(kx, ky, kr, { fill: '#e8e4da', stroke: Sim.colors.ink, width: 1.5 });
  g.line(kx, ky, kx + kr * Math.cos(phi + Math.PI / 2), ky + kr * Math.sin(phi + Math.PI / 2), { color: Sim.colors.red, width: 2.5 });
  g.text('knob', kx + kr + 0.01, ky, { color: Sim.colors.muted, size: 11 });
}

function drawGraph(g) {
  const x0 = 0.47, y0 = 0.02, wd = 0.42, ht = 0.4;
  const toX = (k) => x0 + (k + 180) / 360 * wd, toY = (t) => y0 + (Math.max(-180, Math.min(180, t)) + 180) / 360 * ht;
  g.rect(x0, y0, wd, ht, { fill: '#ffffff', stroke: Sim.colors.axis, width: 1 });
  g.line(toX(0), y0, toX(0), y0 + ht, { color: Sim.colors.grid, width: 1 });
  g.line(x0, toY(0), x0 + wd, toY(0), { color: Sim.colors.grid, width: 1 });
  const ideal = [];
  for (let k = -180; k <= 180; k += 2) ideal.push({ x: toX(k), y: toY(idealTipDeg(k)) });
  g.polyline(ideal, { color: Sim.colors.muted, width: 1.5, dash: [5, 4] });
  for (let i = 0; i < trail.length; i += 2) g.dot(toX(trail[i].knob), toY(trail[i].tip), 1.6, Sim.colors.blue);
  g.dot(toX(knob.value), toY(tipAngleDeg()), 4, Sim.colors.red);
  g.text('tip angle vs knob angle (deg)', x0 + 0.005, y0 + ht + 0.02, { color: Sim.colors.muted, size: 11 });
  g.text('dashed: ideal (u - slack) / r', x0 + 0.005, y0 - 0.02, { color: Sim.colors.muted, size: 11 });
}

function draw(g) {
  drawSheathAndKnob(g);
  g.rect(-0.05, -l, 0.1, l, { fill: '#d9d6cc' });
  // backbone and guides
  const spine = [{ x: 0, y: 0 }];
  for (let i = 1; i <= N; i++) spine.push({ x: x[i] + 0.5 * l * Math.cos(th[i]), y: y[i] + 0.5 * l * Math.sin(th[i]) });
  g.polyline(spine, { color: '#c9b8e6', width: 2 * r * g.pxPerUnit + 4, alpha: 0.55 });
  g.polyline(spine, { color: Sim.colors.purple, width: 2 });
  for (let i = 1; i <= N; i++) {
    const nx = -Math.sin(th[i]), ny = Math.cos(th[i]);
    g.line(x[i] - 1.3 * r * nx, y[i] - 1.3 * r * ny, x[i] + 1.3 * r * nx, y[i] + 1.3 * r * ny, { color: Sim.colors.ink, width: 1.5 });
  }
  // cables: solid when taut, dashed when slack
  for (const side of [+1, -1]) {
    const len = cableLength(side), taut = len > allowedLength(side) - 1e-5;
    const pts = [];
    for (let k = 0; k <= N; k++) pts.push({ x: guideX[k], y: guideY[k] });
    g.polyline(pts, { color: side > 0 ? Sim.colors.red : Sim.colors.blue, width: taut ? 2.2 : 1.2, dash: taut ? null : [3, 3] });
  }
  drawGraph(g);

  const u = knobPull();
  Sim.show('cable pull u = r_p phi', (1000 * u).toFixed(1) + ' mm');
  Sim.show('tip angle (sim)', tipAngleDeg().toFixed(1) + ' deg');
  Sim.show('ideal (u - slack) / r', idealTipDeg(knob.value).toFixed(1) + ' deg');
  const T = Math.max(tension.left, tension.right);
  Sim.show('cable tension at base', T.toFixed(2) + ' N');
  Sim.show('force at knob T e^(mu beta)', (T * Math.exp(friction.value * curlDeg.value * Math.PI / 180)).toFixed(2) + ' N');
}

Sim.run({ title: 'A tentacle with two tendons (starter)', view: [-0.45, -0.5, 0.95, 0.48], reset, step, draw, grid: 0.05 });
