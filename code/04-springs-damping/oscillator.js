// Chapter 04 - Springs, Damping and Energy: one damped oscillator
// A 1 kg block on a 100 N/m spring with a damper, released from 0.3 m. Pick the damping ratio zeta.
// The trace shows the simulated displacement (solid) over the exact solution (dashed);
// the plot shows where the energy is: kinetic, spring, and their total, as fractions of the start.

const zetaSlider = Sim.slider('damping ratio zeta', 0, 2.5, 0.1, { step: 0.01, reset: true });
Sim.select('preset', ['undamped (0)', 'light (0.1)', 'medium (0.5)', 'critical (1)', 'heavy (2)'], 'light (0.1)', {
  onChange: (v) => { zetaSlider.set(+v.match(/\(([\d.]+)\)/)[1]); Sim.reset(); },
});

const k = 100, m = 1, u0 = 0.3;          // N/m, kg, m
const omega = Math.sqrt(k / m);           // 10 rad/s
const h = 0.001;                          // 1 ms steps, 17 per frame
let b, history, lastOutside;

function reset() {
  b = { x: u0, v: 0 };
  history = [];
  lastOutside = 0;
}

function c() { return 2 * zetaSlider.value * Math.sqrt(k * m); }

// Exact displacement of the damped oscillator released from rest at u0.
function exact(t) {
  const z = zetaSlider.value;
  if (z < 1 - 1e-9) {
    const wd = omega * Math.sqrt(1 - z * z);
    return u0 * Math.exp(-z * omega * t) * (Math.cos(wd * t) + z / Math.sqrt(1 - z * z) * Math.sin(wd * t));
  }
  if (z < 1 + 1e-9) return u0 * (1 + omega * t) * Math.exp(-omega * t);
  const s = Math.sqrt(z * z - 1), r1 = -omega * (z - s), r2 = -omega * (z + s);
  return u0 * (r1 * Math.exp(r2 * t) - r2 * Math.exp(r1 * t)) / (r1 - r2);
}

function step(dt) {
  // Symplectic Euler with spring and damper forces.
  const f = -k * b.x - c() * b.v;
  b.v += h * f / m;
  b.x += h * b.v;
  const t = Sim.time + dt;
  if (Math.abs(b.x) > 0.01 * u0) lastOutside = t;
  history.push({ t, x: b.x, e: exact(t) });
  if (history.length > 4500) history.shift();
}

function drawSpring(g, x1, x2, y) {
  const coils = 14, amp = 0.05, pts = [{ x: x1, y }];
  for (let i = 0; i <= coils; i++) {
    const t = i / coils;
    pts.push({ x: x1 + 0.05 + (x2 - x1 - 0.1) * t, y: y + (i === 0 || i === coils ? 0 : (i % 2 ? amp : -amp)) });
  }
  pts.push({ x: x2, y });
  g.polyline(pts, { color: Sim.colors.gray, width: 1.5 });
}

function draw(g) {
  const wall = -1.0, rest = 0, y = 1.35;
  g.rect(wall - 0.08, y - 0.3, 0.08, 0.6, { fill: '#d9d6cc' });
  drawSpring(g, wall, rest + b.x - 0.1, y + 0.05);
  // the dashpot, drawn under the spring
  g.line(wall, y - 0.12, rest + b.x - 0.1, y - 0.12, { color: Sim.colors.muted, width: 1.5 });
  g.rect(-0.75, y - 0.18, 0.35, 0.12, { stroke: Sim.colors.muted });
  g.rect(rest + b.x - 0.1, y - 0.1, 0.2, 0.2, { fill: Sim.colors.blue });
  g.line(rest, y - 0.25, rest, y + 0.25, { color: Sim.colors.muted, dash: [3, 3], width: 1 });

  // Displacement trace over the last 4.5 s: simulation solid, exact dashed.
  const tx0 = -1.1, tx1 = 1.6, ty = 0.45, sy = 1.0, span = 4.5;
  g.line(tx0, ty, tx1, ty, { color: Sim.colors.axis, width: 1 });
  g.text('displacement u(t), last 4.5 s', tx0, ty + 0.4, { color: Sim.colors.muted, size: 11 });
  if (history.length > 1) {
    const tEnd = history[history.length - 1].t, sim = [], ex = [];
    for (const s of history) {
      const px = tx1 - (tEnd - s.t) / span * (tx1 - tx0);
      sim.push({ x: px, y: ty + sy * s.x });
      ex.push({ x: px, y: ty + sy * s.e });
    }
    g.polyline(ex, { color: Sim.colors.ink, dash: [5, 4], width: 1.5 });
    g.polyline(sim, { color: Sim.colors.blue, width: 2 });
  }

  const kin = 0.5 * m * b.v * b.v, spr = 0.5 * k * b.x * b.x, e0 = 0.5 * k * u0 * u0;
  Sim.plot('kinetic', kin / e0, { color: Sim.colors.red, min: 0, max: 1 });
  Sim.plot('spring', spr / e0, { color: Sim.colors.green });
  Sim.plot('total E/E0', (kin + spr) / e0, { color: Sim.colors.ink });
  Sim.show('zeta', zetaSlider.value.toFixed(2));
  Sim.show('c = 2 zeta sqrt(k m)', c().toFixed(2) + ' N s/m');
  Sim.show('u simulated', b.x.toFixed(5) + ' m');
  Sim.show('u exact', exact(Sim.time).toFixed(5) + ' m');
  Sim.show('last time outside 1 %', lastOutside.toFixed(3) + ' s');
}

Sim.run({ title: 'One damped oscillator', view: [-1.25, -0.9, 1.75, 1.8], reset, step, draw, dt: h, stepsPerFrame: 17 });
