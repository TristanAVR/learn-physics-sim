// Chapter 01 - The Sandbox (reference solution)
// No physics yet: a dot moving on a circle, an arrow to the mouse, and a few controls,
// to show every part of the sandbox you will use in later chapters.

const speed = Sim.slider('angular speed (rad/s)', 0, 6, 1.5);
const radius = Sim.slider('radius (m)', 0.1, 1.5, 0.8, { step: 0.05 });
const trail = Sim.checkbox('draw trail', true, { reset: true });
Sim.button('jump half a turn', () => { angle += Math.PI; });

let angle;
let history;

function reset() {
  angle = 0;
  history = [];
}

// step(dt) advances the state by dt seconds. Here the "state" is just an angle.
function step(dt) {
  angle += speed.value * dt;
  if (trail.value) {
    history.push(position());
    if (history.length > 120) history.shift();
  }
}

function position() {
  return V.v(radius.value * Math.cos(angle), radius.value * Math.sin(angle));
}

// draw(g) draws the current state. It is called every frame, even while paused.
function draw(g) {
  g.circle(0, 0, radius.value, { stroke: Sim.colors.axis, width: 1, dash: [4, 4] });
  if (trail.value && history.length > 1) g.polyline(history, { color: Sim.colors.teal, width: 2, alpha: 0.5 });

  const p = position();
  g.line(0, 0, p.x, p.y, { color: Sim.colors.muted, width: 1 });
  g.circle(p.x, p.y, 0.06, { fill: Sim.colors.blue });

  // An arrow from the dot to the mouse, in world units.
  const m = V.v(Sim.mouse.x, Sim.mouse.y);
  if (Sim.mouse.inside) {
    const d = V.sub(m, p);
    g.arrow(p.x, p.y, d.x, d.y, { color: Sim.mouse.down ? Sim.colors.red : Sim.colors.orange, width: 2 });
    Sim.show('distance to mouse', V.len(d), 3);
  }

  Sim.show('angle', angle % (2 * Math.PI), 2);
  Sim.plot('x', p.x, { color: Sim.colors.blue });
  Sim.plot('y', p.y, { color: Sim.colors.red });
}

Sim.run({ title: 'Hello, sandbox', view: [-2, -1.6, 2, 1.6], reset, step, draw });
