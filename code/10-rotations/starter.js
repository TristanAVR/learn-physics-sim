// Chapter 10 - Rotations and Rigid Bodies (starter)
// 2D rigid bars joined by XPBD pin joints, with substeps. Drag a bar with the mouse.
// Scenes: a hanging chain of bars, a single bar released from a small angle (to time its period),
// and a free bar with no gravity (to test an off-centre impulse).
// The bars cannot rotate yet: fill in the TODOs in applyImpulse, solvePointConstraint and substep.

const scene = Sim.select('scene', ['chain of bars', 'pendulum (small swing)', 'free bar (no gravity)'], 'chain of bars', { reset: true });
const barCount = Sim.slider('bars (chain scene)', 1, 8, 4, { step: 1, reset: true });
const substeps = Sim.slider('substeps per frame', 1, 40, 10, { step: 1 });
Sim.button('kick the last bar (off-centre impulse)', () => {
  const b = bodies[bodies.length - 1];
  // 0.5 N s straight up, applied at the bar's right end.
  const end = worldPoint(b, barLength / 2, 0);
  applyImpulse(b, end.x, end.y, 0, kickImpulse);
});

const gravity = 9.81;       // m/s^2
const barLength = 0.5;      // m
const barThick = 0.06;      // m
const barMass = 1.0;        // kg
const kickImpulse = 0.5;    // N s
const dragCompliance = 1e-3; // m/N, a soft spring between the mouse and the grabbed point

let bodies, joints, drag, clock, lastAngle, lastCrossing, period;

// Moment of inertia of a solid rectangle (length L, thickness t) about its centre: m (L^2 + t^2) / 12.
function makeBar(x, y, theta) {
  const I = barMass * (barLength * barLength + barThick * barThick) / 12;
  return {
    x, y, theta, vx: 0, vy: 0, omega: 0,
    px: x, py: y, ptheta: theta,
    invMass: 1 / barMass, invInertia: 1 / I,
  };
}

function reset() {
  bodies = []; joints = []; drag = null;
  clock = 0; lastAngle = 0; lastCrossing = -1; period = NaN;
  if (scene.value === 'free bar (no gravity)') {
    bodies.push(makeBar(0, -0.6, 0));
    return;
  }
  const n = scene.value === 'chain of bars' ? barCount.value : 1;
  // Start angle: horizontal for the chain, 10 degrees from vertical for the pendulum.
  const start = scene.value === 'chain of bars' ? 0 : -Math.PI / 2 + 10 * Math.PI / 180;
  const cx = Math.cos(start), cy = Math.sin(start);
  for (let i = 0; i < n; i++) {
    const s = (i + 0.5) * barLength;
    bodies.push(makeBar(s * cx, s * cy, start));
  }
  // Joint 0 pins the left end of bar 0 to the world origin; joint i joins bar i-1's right end to bar i's left end.
  joints.push({ a: -1, ax: 0, ay: 0, b: 0, bx: -barLength / 2, by: 0, compliance: 0 });
  for (let i = 1; i < n; i++) {
    joints.push({ a: i - 1, ax: barLength / 2, ay: 0, b: i, bx: -barLength / 2, by: 0, compliance: 0 });
  }
}

// A point given in a body's local frame, in world coordinates. Body -1 is the fixed world.
function worldPoint(body, lx, ly) {
  if (!body) return { x: lx, y: ly };
  const c = Math.cos(body.theta), s = Math.sin(body.theta);
  return { x: body.x + c * lx - s * ly, y: body.y + s * lx + c * ly };
}

// Instantaneous impulse (jx, jy) at world point (px, py): changes both linear and angular velocity.
function applyImpulse(body, px, py, jx, jy) {
  body.vx += jx * body.invMass;
  body.vy += jy * body.invMass;
  // TODO (exercise 2): an impulse at a point away from the centre also changes the angular velocity,
  // by (r x J) / I, where r is the lever arm from the centre of mass to the point.
}

// XPBD positional constraint: bring world points pa (on body A) and pb (on body B) together.
// Either body may be null (fixed world). With one iteration per substep, lambda starts at zero on the
// only visit, so the alphaTilde * lambda term of the XPBD update vanishes and is left out.
function solvePointConstraint(A, pa, B, pb, alphaTilde) {
  const dx = pb.x - pa.x, dy = pb.y - pa.y;
  const c = Math.hypot(dx, dy);
  if (c < 1e-12) return;
  const nx = dx / c, ny = dy / c;
  // TODO (exercise 1): this treats each body as a point mass at its centre, so the joint cannot turn it.
  // Use the generalised inverse mass w = 1/m + (r x n)^2 / I, with r the lever arm from the body's centre
  // to its constraint point, and also rotate each body by (r x p) / I, where p = dLambda n.
  let wA = 0, wB = 0;
  if (A) wA = A.invMass;
  if (B) wB = B.invMass;
  const dLambda = -c / (wA + wB + alphaTilde);
  // The correction p = dLambda n pushes B along p and A along -p, at their points.
  if (A) { A.x -= A.invMass * dLambda * nx; A.y -= A.invMass * dLambda * ny; }
  if (B) { B.x += B.invMass * dLambda * nx; B.y += B.invMass * dLambda * ny; }
}

function substep(h) {
  const g = scene.value === 'free bar (no gravity)' ? 0 : gravity;
  for (const b of bodies) {
    b.px = b.x; b.py = b.y; b.ptheta = b.theta;
    b.vy -= h * g;
    b.x += h * b.vx; b.y += h * b.vy;
    // TODO (exercise 1): integrate the angle with the angular velocity, like the position.
  }
  for (const j of joints) {
    const A = j.a >= 0 ? bodies[j.a] : null, B = bodies[j.b];
    solvePointConstraint(A, worldPoint(A, j.ax, j.ay), B, worldPoint(B, j.bx, j.by), j.compliance / (h * h));
  }
  if (drag) {
    const B = bodies[drag.body];
    solvePointConstraint(null, { x: Sim.mouse.x, y: Sim.mouse.y }, B, worldPoint(B, drag.lx, drag.ly), dragCompliance / (h * h));
  }
  for (const b of bodies) {
    b.vx = (b.x - b.px) / h; b.vy = (b.y - b.py) / h;
    // TODO (exercise 1): the angular velocity from the change in angle over the substep.
  }
}

function step(dt) {
  handleMouse();
  const n = substeps.value;
  for (let i = 0; i < n; i++) substep(dt / n);
  timePeriod(dt);
}

function handleMouse() {
  if (Sim.mouse.pressed) {
    for (let i = 0; i < bodies.length; i++) {
      const b = bodies[i];
      // The mouse position in the bar's local frame.
      const dx = Sim.mouse.x - b.x, dy = Sim.mouse.y - b.y;
      const c = Math.cos(b.theta), s = Math.sin(b.theta);
      const lx = c * dx + s * dy, ly = -s * dx + c * dy;
      if (Math.abs(lx) <= barLength / 2 + 0.03 && Math.abs(ly) <= barThick / 2 + 0.05) { drag = { body: i, lx, ly }; break; }
    }
  }
  if (Sim.mouse.released) drag = null;
}

// Pendulum scene: time the swing from upward zero crossings of the angle from vertical.
function timePeriod(dt) {
  clock += dt;
  if (scene.value !== 'pendulum (small swing)') return;
  const b = bodies[0];
  const angle = Math.atan2(b.x, -b.y);
  if (lastAngle < 0 && angle >= 0) {
    const t = clock - dt * angle / (angle - lastAngle);   // interpolate the crossing time within the frame
    if (lastCrossing < 0 || t - lastCrossing > 0.3) {     // ignore jitter right after a crossing
      if (lastCrossing >= 0) period = t - lastCrossing;
      lastCrossing = t;
    }
  }
  lastAngle = angle;
}

function drawBar(g, b, color) {
  const pts = [];
  for (const [lx, ly] of [[-1, -1], [1, -1], [1, 1], [-1, 1]]) {
    pts.push(worldPoint(b, lx * barLength / 2, ly * barThick / 2));
  }
  g.polyline(pts, { closed: true, fill: color, stroke: Sim.colors.ink, width: 1 });
  g.dot(b.x, b.y, 2.5, Sim.colors.ink);
}

function draw(g) {
  if (scene.value === 'free bar (no gravity)') {
    const b = bodies[0];
    Sim.setView([b.x - 1.6, b.y - 1.3, b.x + 1.6, b.y + 0.9]);
  } else {
    Sim.setView([-2.2, -2.6, 2.2, 0.6]);
    g.polyline([{ x: -0.08, y: 0.08 }, { x: 0.08, y: 0.08 }, { x: 0, y: 0 }], { closed: true, fill: Sim.colors.gray });
  }
  for (let i = 0; i < bodies.length; i++) drawBar(g, bodies[i], Sim.palette[i % Sim.palette.length]);
  for (const j of joints) {
    const p = worldPoint(bodies[j.b], j.bx, j.by);
    g.circle(p.x, p.y, 0.018, { fill: '#ffffff', stroke: Sim.colors.ink, width: 1 });
  }
  if (drag) {
    const p = worldPoint(bodies[drag.body], drag.lx, drag.ly);
    g.line(p.x, p.y, Sim.mouse.x, Sim.mouse.y, { color: Sim.colors.red, dash: [4, 3] });
  }

  const I = barMass * (barLength * barLength + barThick * barThick) / 12;
  if (scene.value === 'free bar (no gravity)') {
    const b = bodies[0];
    Sim.show('speed |v| (m/s)', Math.hypot(b.vx, b.vy), 4);
    Sim.show('omega (rad/s)', b.omega, 4);
    Sim.show('after one kick, predicted', `|v| = ${(kickImpulse / barMass).toFixed(4)}, omega = ${(kickImpulse * barLength / 2 / I).toFixed(4)}`);
  } else if (scene.value === 'pendulum (small swing)') {
    // Small-angle period of a bar pivoted at its end: T = 2 pi sqrt(I_end / (m g L/2)).
    const Iend = I + barMass * (barLength / 2) * (barLength / 2);
    const T = 2 * Math.PI * Math.sqrt(Iend / (barMass * gravity * barLength / 2));
    Sim.show('measured period', isFinite(period) ? period.toFixed(4) + ' s' : 'timing...');
    Sim.show('small-angle prediction', T.toFixed(4) + ' s');
    const a0 = 10 * Math.PI / 180;
    Sim.show('corrected for 10 deg amplitude', (T * (1 + a0 * a0 / 16)).toFixed(4) + ' s');
  }
  let ke = 0, pe = 0;
  for (const b of bodies) {
    ke += 0.5 * barMass * (b.vx * b.vx + b.vy * b.vy) + 0.5 * b.omega * b.omega / b.invInertia;
    pe += barMass * gravity * b.y;
  }
  if (scene.value === 'free bar (no gravity)') pe = 0;
  Sim.plot('energy (J)', ke + pe, { color: Sim.colors.purple });
  Sim.show('drag', 'grab a bar with the mouse');
}

Sim.run({ title: 'Rigid bars and pin joints (starter)', view: [-2.2, -2.6, 2.2, 0.6], reset, step, draw });
