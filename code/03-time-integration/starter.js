// Chapter 03 - Time Integration (starter)
// Three identical mass-spring systems, each advanced with a different integrator.
// Explicit Euler is done for you. Your job: implement symplecticEuler and implicitEuler below,
// then compare the three energy curves.

const stiffness = Sim.logSlider('stiffness k (N/m)', 0, 4, 100, { reset: true });
const mass = Sim.slider('mass m (kg)', 0.1, 5, 1, { step: 0.1, reset: true });
const stepMs = Sim.slider('time step h (ms)', 1, 60, 16, { step: 1 });

const methods = [
  { name: 'explicit Euler', color: Sim.colors.red, integrate: explicitEuler },
  { name: 'symplectic Euler', color: Sim.colors.green, integrate: symplecticEuler },
  { name: 'implicit Euler', color: Sim.colors.blue, integrate: implicitEuler },
];

const x0 = 0.3; // initial stretch (m)
let bodies;

function reset() {
  bodies = methods.map(() => ({ x: x0, v: 0 }));
}

// Force of a spring of stiffness k stretched by x (rest length at x = 0).
function springForce(x) {
  return -stiffness.value * x;
}

// x(n+1) = x(n) + h v(n),  v(n+1) = v(n) + h a(n)
function explicitEuler(b, h) {
  const a = springForce(b.x) / mass.value;
  b.x = b.x + h * b.v;
  b.v = b.v + h * a;
}

// v(n+1) = v(n) + h a(n),  x(n+1) = x(n) + h v(n+1)
function symplecticEuler(b, h) {
  // TODO: update the velocity first, then the position with the NEW velocity.
}

// v(n+1) = v(n) + h a(n+1),  x(n+1) = x(n) + h v(n+1)
function implicitEuler(b, h) {
  // TODO: the acceleration is evaluated at the END of the step. For a linear spring,
  // a(n+1) = -(k/m) x(n+1) = -(k/m) (x(n) + h v(n+1)). Substitute that into the velocity
  // update and solve for v(n+1) on paper first; then write the two lines here.
}

function energy(b) {
  return 0.5 * mass.value * b.v * b.v + 0.5 * stiffness.value * b.x * b.x;
}

function step() {
  const h = stepMs.value / 1000;
  for (let i = 0; i < methods.length; i++) methods[i].integrate(bodies[i], h);
}

function drawSpring(g, x1, x2, y, color) {
  const coils = 12, amp = 0.05, pts = [];
  pts.push({ x: x1, y });
  for (let i = 0; i <= coils; i++) {
    const t = i / coils;
    pts.push({ x: x1 + 0.05 + (x2 - x1 - 0.1) * t, y: y + (i === 0 || i === coils ? 0 : (i % 2 ? amp : -amp)) });
  }
  pts.push({ x: x2, y });
  g.polyline(pts, { color, width: 1.5 });
}

function draw(g) {
  const e0 = 0.5 * stiffness.value * x0 * x0;
  const wall = -1.0, rest = 0.2;
  g.rect(wall - 0.08, -1.0, 0.08, 2.0, { fill: '#d9d6cc' });
  for (let i = 0; i < methods.length; i++) {
    const y = 0.6 - i * 0.6, b = bodies[i];
    const x = Math.max(-2, Math.min(3, rest + b.x));
    drawSpring(g, wall, x - 0.08, y, Sim.colors.gray);
    g.rect(x - 0.08, y - 0.08, 0.16, 0.16, { fill: methods[i].color });
    g.text(methods[i].name, 1.0, y, { color: methods[i].color });
    g.line(rest, y - 0.14, rest, y + 0.14, { color: Sim.colors.muted, dash: [3, 3], width: 1 });
    const e = energy(b) / e0;
    Sim.plot(methods[i].name, isFinite(e) ? Math.min(e, 3) : 3, { color: methods[i].color, min: 0 });
  }
  const omega = Math.sqrt(stiffness.value / mass.value);
  Sim.show('omega = sqrt(k/m)', omega.toFixed(2) + ' rad/s');
  Sim.show('h * omega', (stepMs.value / 1000 * omega).toFixed(3));
  Sim.show('energy plot', 'E / E0, capped at 3');
}

Sim.run({ title: 'Three integrators, one spring (starter)', view: [-1.2, -1.6, 2.2, 0.9], reset, step, draw, dt: () => stepMs.value / 1000 });
