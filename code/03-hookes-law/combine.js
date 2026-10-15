// Chapter 03 - Springs and Hooke's Law: combining springs
// Two springs, end to end or side by side, carry a hanging mass. Next to them hangs ONE spring with the
// combined stiffness the rules predict. If the rules are right, the two masses always move together.
// The springs are massless; the mass settles with a little damping (tiny steps keep it accurate).

const arrangement = Sim.select('arrangement', ['end to end (series)', 'side by side (parallel)'], 'end to end (series)');
const k1 = Sim.slider('spring 1: k1 (N/m)', 100, 1000, 500, { step: 10 });
const k2 = Sim.slider('spring 2: k2 (N/m)', 100, 1000, 500, { step: 10 });
const mass = Sim.slider('mass (kg)', 0.2, 2, 1, { step: 0.1 });

const g = 9.81, top = 1.25, rest = 0.3;   // ceiling height and each spring's rest length (m)
let X, V;                                  // total stretch of the pair and its rate

function series() { return arrangement.value.startsWith('end'); }

// The combination rules from the text.
function combined() {
  return series() ? 1 / (1 / k1.value + 1 / k2.value) : k1.value + k2.value;
}

function reset() {
  X = 0; V = 0;
}

function step(dt) {
  const n = 100, h = dt / n, k = combined(), m = mass.value, c = 2 * 0.3 * Math.sqrt(k * m);
  for (let i = 0; i < n; i++) {
    V += h * (g - (k / m) * X - (c / m) * V);
    X += h * V;
  }
}

function drawSpring(gr, x, y0, y1, color) {
  const coils = 10, amp = 0.035, pts = [{ x, y: y0 }, { x, y: y0 - 0.02 }];
  for (let i = 1; i < 2 * coils; i++) {
    pts.push({ x: x + (i % 2 ? amp : -amp), y: y0 - 0.02 + (y1 - y0 + 0.04) * i / (2 * coils) });
  }
  pts.push({ x, y: y1 + 0.02 }, { x, y: y1 });
  gr.polyline(pts, { color, width: 1.8 });
}

function block(gr, x, y, label, fill) {
  gr.rect(x - 0.09, y - 0.12, 0.18, 0.12, { fill });
  gr.text(label, x, y - 0.06, { color: '#ffffff', size: 12, align: 'center' });
}

function labelForce(gr, x, y, newtons, side) {
  gr.text(newtons.toFixed(1) + ' N', x + (side < 0 ? -0.06 : 0.06), y, { color: Sim.colors.purple, size: 12, align: side < 0 ? 'right' : 'left' });
}

function draw(gr) {
  gr.rect(-0.75, top, 1.9, 0.04, { fill: '#cfccc2' });
  const W = mass.value * g, k = combined();
  const springForce = k * X;   // what the pair as a whole carries

  if (series()) {
    // Each spring carries the whole force, so each stretches by force / its own k.
    const x1 = springForce / k1.value, x2 = springForce / k2.value;
    const mid = top - rest - x1, end = mid - rest - x2;
    gr.line(-0.3, top - 2 * rest, 0.2, top - 2 * rest, { color: Sim.colors.muted, dash: [4, 4], width: 1 });
    gr.text('rest length of the pair', -0.32, top - 2 * rest, { color: Sim.colors.muted, size: 11, align: 'right' });
    drawSpring(gr, 0, top, mid, Sim.colors.ink);
    gr.circle(0, mid, 0.012, { fill: Sim.colors.ink });
    drawSpring(gr, 0, mid, end, Sim.colors.ink);
    labelForce(gr, -0.06, (top + mid) / 2, springForce, -1);
    labelForce(gr, -0.06, (mid + end) / 2, springForce, -1);
    gr.text('k1, stretch ' + (x1 * 100).toFixed(1) + ' cm', 0.06, (top + mid) / 2, { color: Sim.colors.red, size: 12 });
    gr.text('k2, stretch ' + (x2 * 100).toFixed(1) + ' cm', 0.06, (mid + end) / 2, { color: Sim.colors.red, size: 12 });
    block(gr, 0, end, mass.value.toFixed(1) + ' kg', Sim.colors.blue);
    Sim.show('force in each spring', springForce.toFixed(2) + ' N (the whole load)');
    Sim.show('total stretch', ((x1 + x2) * 100).toFixed(2) + ' cm = ' + (x1 * 100).toFixed(2) + ' + ' + (x2 * 100).toFixed(2));
  } else {
    // Both stretch the same X; the stiffer one carries more. The load hangs where the bar stays level.
    const end = top - rest - X, a = -0.15, b = 0.15;
    const f1 = k1.value * X, f2 = k2.value * X;
    const hang = a + (b - a) * f2 / (f1 + f2 || 1);   // balance of moments: f1 (hang - a) = f2 (b - hang)
    gr.line(-0.35, top - rest, 0.3, top - rest, { color: Sim.colors.muted, dash: [4, 4], width: 1 });
    gr.text('rest length', -0.37, top - rest, { color: Sim.colors.muted, size: 11, align: 'right' });
    drawSpring(gr, a, top, end, Sim.colors.ink);
    drawSpring(gr, b, top, end, Sim.colors.ink);
    gr.rect(a - 0.03, end - 0.02, b - a + 0.06, 0.02, { fill: Sim.colors.ink });
    gr.line(hang, end - 0.02, hang, end - 0.08, { color: Sim.colors.ink, width: 2 });
    labelForce(gr, a - 0.04, top - rest / 2, f1, -1);
    labelForce(gr, b + 0.04, top - rest / 2, f2, 1);
    block(gr, hang, end - 0.08, mass.value.toFixed(1) + ' kg', Sim.colors.blue);
    gr.text('both stretch ' + (X * 100).toFixed(1) + ' cm', -0.35, end - 0.05, { color: Sim.colors.red, size: 12, align: 'right' });
    Sim.show('forces', f1.toFixed(2) + ' N + ' + f2.toFixed(2) + ' N = ' + (f1 + f2).toFixed(2) + ' N');
    Sim.show('stretch of each', (X * 100).toFixed(2) + ' cm');
  }

  // The single equivalent spring, for comparison.
  const ex = 0.85, eend = top - rest - X;
  gr.line(ex - 0.15, top - rest, ex + 0.15, top - rest, { color: Sim.colors.muted, dash: [4, 4], width: 1 });
  drawSpring(gr, ex, top, eend, Sim.colors.green);
  block(gr, ex, eend, mass.value.toFixed(1) + ' kg', Sim.colors.green);
  gr.text('one spring of', ex, top + 0.12, { color: Sim.colors.green, size: 12, align: 'center' });
  gr.text(k.toFixed(0) + ' N/m', ex, top + 0.07, { color: Sim.colors.green, size: 12, align: 'center' });

  Sim.show('combined k', k.toFixed(1) + ' N/m');
  Sim.show('weight', W.toFixed(2) + ' N');
  Sim.show('settled stretch W / k', (W / k * 100).toFixed(2) + ' cm');
}

Sim.run({ title: 'Combining springs', view: [-1.3, -0.05, 1.3, 1.45], reset, step, draw, grid: false });
