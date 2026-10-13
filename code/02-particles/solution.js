// Chapter 02 - Particles and Newton (reference solution)
// One particle (or many) moved by Newton's second law. Every scene uses the same three functions:
// force(p) adds up the forces, integrate(p, h) advances the state, collideGround(p) bounces off y = 0.

const scene = Sim.select('scene', ['projectile', 'bouncing ball', 'terminal velocity', 'fountain'], 'projectile', { reset: true });
const speed = Sim.slider('launch speed (m/s)', 1, 20, 10, { step: 0.5, reset: true });
const angleDeg = Sim.slider('launch angle (deg)', 5, 85, 45, { step: 1, reset: true });
const dragModel = Sim.select('drag', ['none', 'linear', 'quadratic'], 'none', { reset: true });
const cLin = Sim.slider('linear drag c (N s/m)', 0, 2, 0.5, { step: 0.01 });
const cQuad = Sim.slider('quadratic drag c2 (kg/m)', 0, 0.2, 0.05, { step: 0.001 });
const restitution = Sim.slider('restitution e', 0, 1, 0.8, { step: 0.01 });
const stepMs = Sim.slider('time step h (ms)', 1, 50, 2, { step: 1, reset: true });

const mass = 0.5;      // kg, the same for every particle
const gravity = 9.81;  // m/s^2, pointing down (-y)

let particles, trail, landing, peaks, lastVy, emitted, stepsTaken, rng;

function reset() {
  rng = Sim.random(7);
  trail = [];
  landing = null;
  peaks = [];
  lastVy = 0;
  emitted = 0;
  stepsTaken = 0;
  const a = angleDeg.value * Math.PI / 180;
  if (scene.value === 'projectile') {
    particles = [{ x: V.v(0, 0), v: V.v(speed.value * Math.cos(a), speed.value * Math.sin(a)) }];
    const s = Math.max(speed.value * speed.value / gravity, 2);
    Sim.setView([-0.06 * s, -0.12 * s, 1.12 * s, 0.62 * s]);
  } else if (scene.value === 'bouncing ball') {
    particles = [{ x: V.v(0, 2), v: V.v(0, 0) }];
    Sim.setView([-2.2, -0.25, 2.2, 2.3]);
  } else if (scene.value === 'terminal velocity') {
    particles = [{ x: V.v(0, 0), v: V.v(0, 0) }];
  } else {
    particles = [];
    Sim.setView([-7, -0.6, 7, 8]);
  }
}

// Total force on particle p: gravity plus the selected drag. Forces add as vectors.
function force(p) {
  let f = V.v(0, -mass * gravity);                                         // weight, m g, downward
  if (dragModel.value === 'linear') {
    f = V.addScaled(f, p.v, -cLin.value);                                  // F = -c v
  } else if (dragModel.value === 'quadratic') {
    f = V.addScaled(f, p.v, -cQuad.value * V.len(p.v));                    // F = -c2 |v| v
  }
  return f;
}

// Advance one time step h: position by velocity, velocity by acceleration (a = F / m).
function integrate(p, h) {
  const a = V.scale(force(p), 1 / mass);
  p.x = V.addScaled(p.x, p.v, h);
  p.v = V.addScaled(p.v, a, h);
}

// Ground at y = 0: push the particle back to the surface and reflect its vertical velocity,
// keeping the fraction e (the coefficient of restitution).
function collideGround(p) {
  if (p.x.y < 0) {
    p.x.y = 0;
    if (p.v.y < 0) p.v.y = -restitution.value * p.v.y;
  }
}

function emit() {
  const a = Math.PI / 2 + (rng() - 0.5) * 0.6;
  const s = speed.value * (0.85 + 0.3 * rng());
  particles.push({ x: V.v(0, 0.01), v: V.v(s * Math.cos(a), s * Math.sin(a)), color: Sim.palette[emitted % Sim.palette.length] });
  emitted++;
  if (particles.length > 400) particles.shift();
}

function step(h) {
  if (scene.value === 'fountain') {
    if (stepsTaken++ % Math.max(1, Math.round(0.01 / h)) === 0) emit();   // about 100 per second
    for (const p of particles) { integrate(p, h); collideGround(p); }
    return;
  }
  const p = particles[0];
  if (landing) return;                       // the projectile stops where it first lands
  const old = V.copy(p.x);
  integrate(p, h);
  if (scene.value === 'projectile') {
    if (old.y >= 0 && p.x.y < 0) {
      const s = old.y / (old.y - p.x.y);     // interpolate the crossing of y = 0 inside this step
      landing = { x: old.x + s * (p.x.x - old.x), t: Sim.time + s * h };
      p.x = V.v(landing.x, 0);
    }
  } else if (scene.value === 'bouncing ball') {
    collideGround(p);
    if (lastVy > 0 && p.v.y <= 0) peaks.push(p.x.y);   // top of an arc
    lastVy = p.v.y;
  }
}

function drawParabola(g, v0, a) {
  const pts = [], vx = v0 * Math.cos(a), vy = v0 * Math.sin(a), T = 2 * vy / gravity;
  for (let i = 0; i <= 100; i++) {
    const t = T * i / 100;
    pts.push({ x: vx * t, y: vy * t - 0.5 * gravity * t * t });
  }
  g.polyline(pts, { color: Sim.colors.muted, dash: [6, 5], width: 1.5 });
  return vx * T;
}

// Exact trajectory with linear drag: tau = m / c, terminal speed vt = m g / c.
function drawLinearDragCurve(g, v0, a) {
  const c = cLin.value;
  if (c <= 0) return;
  const tau = mass / c, vt = mass * gravity / c, vx = v0 * Math.cos(a), vy = v0 * Math.sin(a), pts = [];
  for (let i = 0; i <= 400; i++) {
    const t = 0.01 * i, e = 1 - Math.exp(-t / tau);
    const y = (vy + vt) * tau * e - vt * t;
    pts.push({ x: vx * tau * e, y });
    if (y < 0) break;
  }
  g.polyline(pts, { color: Sim.colors.orange, dash: [3, 4], width: 1.5 });
}

function terminalSpeed() {
  if (dragModel.value === 'linear' && cLin.value > 0) return mass * gravity / cLin.value;
  if (dragModel.value === 'quadratic' && cQuad.value > 0) return Math.sqrt(mass * gravity / cQuad.value);
  return Infinity;
}

function draw(g) {
  const a = angleDeg.value * Math.PI / 180;
  if (scene.value === 'fountain') {
    g.ground(0);
    for (const p of particles) g.dot(p.x.x, p.x.y, 3, p.color);
    Sim.show('particles', particles.length);
    return;
  }
  const p = particles[0];

  if (scene.value === 'projectile') {
    g.ground(0);
    if (Sim.frame % 2 === 0 && !landing) trail.push(V.copy(p.x));
    const R = drawParabola(g, speed.value, a);
    if (dragModel.value === 'linear') drawLinearDragCurve(g, speed.value, a);
    g.line(R, -0.3, R, 0.3, { color: Sim.colors.muted, width: 2 });
    g.text('exact landing', R, 0.5, { color: Sim.colors.muted, align: 'center', size: 11 });
    for (const q of trail) g.dot(q.x, q.y, 2, Sim.colors.blue);
    g.dot(p.x.x, p.x.y, 6, Sim.colors.blue);
    if (!landing) g.arrow(p.x.x, p.x.y, p.v.x * 0.1, p.v.y * 0.1, { color: Sim.colors.red, width: 2 });
    Sim.show('exact range R', R.toFixed(4) + ' m');
    Sim.show('exact flight time', (2 * speed.value * Math.sin(a) / gravity).toFixed(4) + ' s');
    if (landing) {
      Sim.show('simulated landing x', landing.x.toFixed(4) + ' m');
      Sim.show('simulated flight time', landing.t.toFixed(4) + ' s');
      Sim.show('difference', (landing.x - R).toFixed(4) + ' m');
    } else {
      Sim.show('x, y', p.x.x.toFixed(3) + ', ' + p.x.y.toFixed(3) + ' m');
    }
  } else if (scene.value === 'bouncing ball') {
    g.ground(0);
    g.dot(p.x.x, p.x.y, 8, Sim.colors.blue);
    for (let i = 0; i < Math.min(peaks.length, 6); i++) {
      g.line(-1.6, peaks[i], 1.6, peaks[i], { color: Sim.colors.muted, dash: [4, 4], width: 1 });
    }
    const E = 0.5 * mass * V.lenSq(p.v) + mass * gravity * p.x.y;
    Sim.show('energy E', E.toFixed(4) + ' J');
    Sim.show('bounces', peaks.length);
    if (peaks.length > 0) Sim.show('last peak height', peaks[peaks.length - 1].toFixed(4) + ' m');
    if (peaks.length > 1) Sim.show('ratio of last two peaks', (peaks[peaks.length - 1] / peaks[peaks.length - 2]).toFixed(4));
    Sim.show('predicted ratio e^2', (restitution.value * restitution.value).toFixed(4));
    Sim.plot('height (m)', p.x.y, { color: Sim.colors.blue, min: 0 });
  } else {
    // Terminal velocity: the camera follows the particle; the marks every 5 m show the fall.
    Sim.setView([p.x.x - 6, p.x.y - 4, p.x.x + 6, p.x.y + 4]);
    for (let y = Math.ceil((p.x.y - 5) / 5) * 5; y < p.x.y + 5; y += 5) {
      g.line(-1.5, y, 1.5, y, { color: Sim.colors.muted, width: 1 });
      g.text(y.toFixed(0) + ' m', 1.7, y, { color: Sim.colors.muted, size: 11 });
    }
    g.dot(p.x.x, p.x.y, 8, Sim.colors.blue);
    g.arrow(p.x.x - 0.6, p.x.y, 0, p.v.y * 0.15, { color: Sim.colors.red, width: 2 });
    const vt = terminalSpeed();
    Sim.show('speed |v|', V.len(p.v).toFixed(4) + ' m/s');
    Sim.show('predicted v_t', isFinite(vt) ? vt.toFixed(4) + ' m/s' : 'none (no drag)');
    Sim.plot('speed (m/s)', V.len(p.v), { color: Sim.colors.blue, min: 0 });
    if (isFinite(vt)) Sim.plot('v_t', vt, { color: Sim.colors.muted });
  }
}

Sim.run({
  title: 'Particles and Newton',
  view: [-1, -1, 12, 6],
  reset, step, draw,
  dt: () => stepMs.value / 1000,
  stepsPerFrame: () => Math.max(1, Math.round(1000 / 60 / stepMs.value)),
});
