// Chapter 03 - Springs and Hooke's Law (starter)
// A mass hangs on a spring. The four functions below are Hooke's law and what follows from it.
// The blue mass moves by YOUR formulas (no simulation at all: position as a function of time).
// The faint grey ghost is an accurate numerical simulation of the same spring. If your formulas
// are right, the two move together.

const stiffness = Sim.slider('stiffness k (N/m)', 20, 200, 50, { step: 1, reset: true });
const mass = Sim.slider('mass m (kg)', 0.1, 2, 1, { step: 0.1, reset: true });
const pull = Sim.slider('extra pull before release (cm)', 0, 30, 15, { step: 1, reset: true });

const g = 9.81;
const top = 1.2;          // where the spring hangs from (m)
const restLength = 0.5;   // natural length of the spring (m)

// Everything is measured along the stretch, which points DOWN here: x > 0 means the spring is
// longer than its rest length, and a positive force pulls the mass down.

// Hooke's law: the spring's force on the mass for a stretch x. Negative = pulling back up.
function springForce(x) {
  // TODO: Hooke's law. Use stiffness.value for k. Mind the sign.
  return 0;
}

// Hanging still, the forces cancel: spring -k x plus weight +m g is zero, so k x = m g.
function equilibriumStretch() {
  // TODO: solve k x = m g for x. Use mass.value, g and stiffness.value.
  return 0;
}

// How fast it swings, in radians per second.
function omega() {
  // TODO: the angular frequency of a mass on a spring.
  return 0;
}

// Stretch at time t after releasing it from rest, `pullDown` metres below equilibrium.
// It oscillates around the equilibrium: x(t) = x_eq + A cos(omega t).
function stretchAt(t, pullDown) {
  // TODO: start pullDown below equilibrium, swing around it at omega().
  return equilibriumStretch();
}

// --- Everything below is scaffolding: drawing, and the reference ghost. ---

let ghost;
function reset() {
  // The ghost: the true motion, simulated with many tiny steps (chapter 04 explains how).
  ghost = { x: mass.value * g / stiffness.value + pull.value / 100, v: 0 };
}

function step(dt) {
  const n = 200, h = dt / n, k = stiffness.value, m = mass.value;
  for (let i = 0; i < n; i++) {
    const a = (m * g - k * ghost.x) / m;   // weight down (+x), spring pulls back up
    ghost.v += h * a;
    ghost.x += h * ghost.v;
  }
}

function drawSpring(gr, yTop, yBottom, color, alpha) {
  const coils = 14, amp = 0.05, pts = [{ x: 0, y: yTop }];
  for (let i = 1; i < coils; i++) {
    const y = yTop + (yBottom - yTop) * i / coils;
    pts.push({ x: (i % 2 ? amp : -amp), y });
  }
  pts.push({ x: 0, y: yBottom });
  gr.polyline(pts, { color, width: 1.8, alpha });
}

function draw(gr) {
  gr.rect(-0.4, top, 0.8, 0.06, { fill: '#cfccc2' });
  const xEq = equilibriumStretch();
  const yEq = top - restLength - xEq;
  gr.line(-0.5, top - restLength, 0.5, top - restLength, { color: Sim.colors.muted, dash: [3, 4], width: 1 });
  gr.text('rest length', 0.52, top - restLength, { color: Sim.colors.muted, size: 11 });
  if (isFinite(yEq)) {
    gr.line(-0.5, yEq, 0.5, yEq, { color: Sim.colors.green, dash: [6, 4], width: 1 });
    gr.text('your equilibrium', 0.52, yEq, { color: Sim.colors.green, size: 11 });
  }

  // Ghost (reference)
  const yGhost = top - restLength - ghost.x;
  drawSpring(gr, top, yGhost + 0.08, Sim.colors.gray, 0.35);
  gr.rect(-0.08, yGhost - 0.08, 0.16, 0.16, { fill: Sim.colors.gray, alpha: 0.3 });

  // Your formulas
  const x = stretchAt(Sim.time, pull.value / 100);
  const y = top - restLength - x;
  if (isFinite(y)) {
    drawSpring(gr, top, y + 0.08, Sim.colors.ink, 1);
    gr.rect(-0.08, y - 0.08, 0.16, 0.16, { fill: Sim.colors.blue });
    // Forces are along the stretch (down = positive), so the arrows point at -F on screen.
    const F = springForce(x), W = mass.value * g, s = 0.012;
    gr.arrow(0.13, y, 0, -F * s, { color: Sim.colors.red, width: 2.5 });
    gr.arrow(0.2, y, 0, -W * s, { color: Sim.colors.purple, width: 2.5 });
    gr.text('spring F = ' + F.toFixed(1) + ' N', 0.24, y + 0.06, { color: Sim.colors.red, size: 11 });
    gr.text('weight +' + W.toFixed(1) + ' N', 0.24, y - 0.06, { color: Sim.colors.purple, size: 11 });
  }

  const w = omega();
  Sim.show('equilibrium stretch', (xEq * 100).toFixed(1) + ' cm');
  Sim.show('omega', isFinite(w) ? w.toFixed(2) + ' rad/s' : '-');
  Sim.show('period 2 pi / omega', isFinite(w) && w > 0 ? (2 * Math.PI / w).toFixed(3) + ' s' : '-');
  Sim.show('you vs ghost', isFinite(x) ? ((x - ghost.x) * 100).toFixed(2) + ' cm' : '-');
  Sim.plot('your stretch (cm)', isFinite(x) ? x * 100 : 0, { color: Sim.colors.blue });
  Sim.plot('ghost (cm)', ghost.x * 100, { color: Sim.colors.gray });
}

Sim.run({ title: "Hooke's law: hang it, predict it (starter)", view: [-1.1, -0.9, 1.6, 1.35], reset, step, draw, plotCorner: 'top-right' });
