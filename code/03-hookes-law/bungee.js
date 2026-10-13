// Chapter 03 - Springs and Hooke's Law: the bungee puzzle
// A 70 kg jumper steps off a bridge 40 m above the water on a 15 m cord. Choose the cord's
// stiffness so the lowest point is 2 m above the water: work it out on paper first (energy),
// then type it in and jump.

const stiffness = Sim.slider('cord stiffness k (N/m)', 20, 300, 60, { step: 0.5, reset: true });
const massS = Sim.slider('jumper mass (kg)', 40, 120, 70, { step: 1, reset: true });
Sim.button('jump', () => { Sim.reset(); jumping = true; });

const g = 9.81, bridge = 40, cordLength = 15, target = 2;
let y, v, lowest, jumping, done;

function reset() {
  y = bridge; v = 0; lowest = bridge; jumping = false; done = false;
}

// Accurate tiny steps (chapter 04 explains why tiny). The cord only pulls when stretched:
// it is a spring with Hooke's law for stretch > 0, and slack (no force) otherwise.
function step(dt) {
  if (!jumping || done) return;
  const n = 400, h = dt / n, m = massS.value;
  for (let i = 0; i < n; i++) {
    const stretch = (bridge - y) - cordLength;
    const cord = stretch > 0 ? stiffness.value * stretch : 0;
    const a = -g + cord / m;
    v += h * a;
    y += h * v;
    lowest = Math.min(lowest, y);
    if (y <= 0) { done = true; break; }
    if (v > 0 && bridge - y > cordLength) { done = true; break; }  // on the way back up: the lowest point is known
  }
}

function draw(gr) {
  gr.rect(-12, -6, 24, 6, { fill: '#9cc3e8' });
  gr.line(-12, 0, 12, 0, { color: Sim.colors.blue, width: 2 });
  gr.rect(-12, bridge, 12, 1.5, { fill: '#b8b4a8' });
  gr.line(-2, target, 6, target, { color: Sim.colors.green, dash: [6, 4], width: 1.5 });
  gr.text('target: 2 m above the water', 6.3, target, { color: Sim.colors.green, size: 12 });
  const anchor = { x: 0, y: bridge };
  const body = { x: 0.6, y };
  const len = Math.hypot(body.x - anchor.x, body.y - anchor.y);
  gr.line(anchor.x, anchor.y, body.x, body.y, { color: len > cordLength ? Sim.colors.red : Sim.colors.muted, width: 2 });
  gr.circle(body.x, body.y, 0.6, { fill: Sim.colors.orange });
  if (lowest < bridge) gr.line(-1.5, lowest, 1.5, lowest, { color: Sim.colors.orange, width: 2 });

  Sim.show('height', y.toFixed(2) + ' m');
  Sim.show('lowest so far', (lowest < bridge ? lowest.toFixed(2) : '-') + ' m');
  if (done) {
    const miss = lowest - target;
    Sim.show('result', y <= 0 ? 'SPLASH: too soft' : Math.abs(miss) < 0.25 ? 'perfect stop' : miss > 0 ? 'stopped ' + miss.toFixed(1) + ' m too high: too stiff' : 'too low: too soft');
  } else if (!jumping) Sim.show('result', 'press "jump"');
}

Sim.run({ title: 'Bungee: choose the cord', view: [-14, -6, 22, 42], reset, step, draw, grid: false });
