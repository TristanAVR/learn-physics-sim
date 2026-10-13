// Chapter 04 - Time Integration: the phase-space picture
// Each dot is the state (position, velocity) of the same spring, advanced by a different integrator.
// The exact motion goes round the dashed circle forever. Explicit Euler spirals out, implicit Euler
// spirals in, and symplectic Euler stays on a slightly tilted ellipse (drawn faintly).

const stepSize = Sim.slider('time step h (s)', 0.02, 0.7, 0.25, { step: 0.01, reset: true });
const show = Sim.select('show', ['all three', 'explicit Euler', 'symplectic Euler', 'implicit Euler'], 'all three');

const omega = 2 * Math.PI * 0.5; // 0.5 Hz spring: omega = pi rad/s (k/m = omega^2)
const x0 = 1;

const methods = [
  { name: 'explicit Euler', color: Sim.colors.red, step: (s, h) => { const a = -omega * omega * s.x; s.x += h * s.v; s.v += h * a; } },
  { name: 'symplectic Euler', color: Sim.colors.green, step: (s, h) => { s.v += h * (-omega * omega * s.x); s.x += h * s.v; } },
  { name: 'implicit Euler', color: Sim.colors.blue, step: (s, h) => { s.v = (s.v - h * omega * omega * s.x) / (1 + h * h * omega * omega); s.x += h * s.v; } },
];

let states, trails, timer;

function reset() {
  states = methods.map(() => ({ x: x0, v: 0 }));
  trails = methods.map(() => [{ x: x0, y: 0 }]);
  timer = 0;
}

// The sandbox calls step 60 times a second; we take one integrator step every h seconds of real
// time, so a big h is visibly a big jump.
function step(dt) {
  timer += dt;
  const h = stepSize.value;
  while (timer >= h) {
    timer -= h;
    for (let i = 0; i < methods.length; i++) {
      methods[i].step(states[i], h);
      // Plot velocity divided by omega, so the exact motion is a circle of radius x0.
      const p = { x: states[i].x, y: states[i].v / omega };
      if (Math.hypot(p.x, p.y) < 50) trails[i].push(p);
      if (trails[i].length > 400) trails[i].shift();
    }
  }
}

function visible(i) {
  return show.value === 'all three' || show.value === methods[i].name;
}

function draw(g) {
  g.line(-2.6, 0, 2.6, 0, { color: Sim.colors.axis, width: 1 });
  g.line(0, -2.1, 0, 2.1, { color: Sim.colors.axis, width: 1 });
  g.text('position x (m)', 2.55, -0.12, { color: Sim.colors.muted, align: 'right', size: 12 });
  g.text('velocity / omega (m)', 0.06, 2.0, { color: Sim.colors.muted, size: 12 });

  const circle = [];
  for (let i = 0; i <= 120; i++) circle.push({ x: x0 * Math.cos(i / 120 * 2 * Math.PI), y: x0 * Math.sin(i / 120 * 2 * Math.PI) });
  g.polyline(circle, { color: Sim.colors.ink, width: 1.2, dash: [5, 5] });

  // Symplectic Euler's conserved "tilted energy": w^2 x^2 + v^2 - h w^2 x v = w^2 x0^2. Drawn in (x, v/w).
  const h = stepSize.value;
  if (visible(1) && h * omega < 2) {
    const ellipse = [];
    for (let i = 0; i <= 160; i++) {
      const t = i / 160 * 2 * Math.PI, cx = Math.cos(t), cy = Math.sin(t);
      // Solve x^2 + y^2 - h w x y = x0^2 along direction (cx, cy), with y = v / w.
      const r = x0 / Math.sqrt(cx * cx + cy * cy - h * omega * cx * cy);
      ellipse.push({ x: r * cx, y: r * cy });
    }
    g.polyline(ellipse, { color: Sim.colors.green, width: 1, alpha: 0.35 });
  }

  for (let i = 0; i < methods.length; i++) {
    if (!visible(i)) continue;
    const t = trails[i];
    if (t.length > 1) g.polyline(t, { color: methods[i].color, width: 1.5, alpha: 0.8 });
    for (const p of t) g.dot(p.x, p.y, 2, methods[i].color);
    const last = t[t.length - 1];
    g.circle(last.x, last.y, 0.06, { fill: methods[i].color });
  }

  Sim.show('h * omega', (h * omega).toFixed(2));
  Sim.show('explicit: energy x per step', (1 + h * h * omega * omega).toFixed(3));
  Sim.show('implicit: energy x per step', (1 / (1 + h * h * omega * omega)).toFixed(3));
  Sim.show('symplectic stable?', h * omega < 2 ? 'yes (h omega < 2)' : 'no');
}

Sim.run({ title: 'Phase space: where each integrator takes the spring', view: [-2.6, -2.1, 2.6, 2.1], reset, step, draw, grid: false });
