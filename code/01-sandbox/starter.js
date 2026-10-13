// Chapter 01 - The Sandbox (starter)
// No physics yet: just get to know the sandbox. The dot sits still; make it move on a circle,
// then add the arrow to the mouse and the plots. Each TODO is a few lines.

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
  // TODO: advance the angle by the angular speed times dt.

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

  // TODO: when the mouse is over the canvas (Sim.mouse.inside), draw an arrow from the dot to the
  // mouse with g.arrow, and show its length with Sim.show.

  // TODO: plot the dot's x and y with Sim.plot.

  Sim.show('angle', angle % (2 * Math.PI), 2);
}

Sim.run({ title: 'Hello, sandbox (starter)', view: [-2, -1.6, 2, 1.6], reset, step, draw });
