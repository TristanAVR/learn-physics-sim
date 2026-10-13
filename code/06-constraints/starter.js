// Chapter 06 - Constraints and Multipliers (starter)
// A pendulum two ways at the same time step: a stiff spring (red) and an exact distance constraint enforced by
// projecting positions (blue). Also: a double pendulum, and a rope that can go slack next to a rigid rod.
// projectDistance is empty, so the blue bob falls straight off the screen. Fill in the TODOs.

const scene = Sim.select('scene', ['pendulum: spring vs constraint', 'double pendulum', 'rope vs rod'], 'pendulum: spring vs constraint', { reset: true });
const stepMs = Sim.slider('time step h (ms)', 1, 40, 16, { step: 1 });
const stiffness = Sim.logSlider('spring stiffness k (N/m)', 2, 6, 10000);
const angleDeg = Sim.slider('release angle (deg)', 5, 170, 60, { step: 1, reset: true });
const launch = Sim.slider('launch speed, rope scene (m/s)', 0, 8, 5.86, { step: 0.01, reset: true });

const g = 9.81, L = 1.0, mass = 1.0;
let pivots, bodies, spring, links, stats, trails;

function particle(x, y, m) {
  return { x: V.v(x, y), prev: V.v(x, y), v: V.v(0, 0), m, w: m > 0 ? 1 / m : 0 };
}

function newStats() {
  return { lastCross: -1, periods: [], maxErr: 0, bottomT: null, bottomPred: null, slackAt: null, lastTension: 0 };
}

function reset() {
  const a = angleDeg.value * Math.PI / 180;
  trails = [[], []];
  if (scene.value === 'pendulum: spring vs constraint') {
    pivots = [particle(-1.3, 2.3, 0), particle(1.3, 2.3, 0)];
    spring = particle(-1.3 + L * Math.sin(a), 2.3 - L * Math.cos(a), mass);
    bodies = [particle(1.3 + L * Math.sin(a), 2.3 - L * Math.cos(a), mass)];
    links = [{ a: bodies[0], b: pivots[1], L, unilateral: false }];
    stats = [newStats(), newStats()];
    Sim.setView([-2.8, 0.9, 2.8, 3.5]);
  } else if (scene.value === 'double pendulum') {
    pivots = [particle(0, 2.3, 0)];
    const b1 = particle(L * Math.sin(a), 2.3 - L * Math.cos(a), mass);
    const b2 = particle(b1.x.x + L * Math.sin(a + 0.5), b1.x.y - L * Math.cos(a + 0.5), mass);
    bodies = [b1, b2];
    links = [{ a: b1, b: pivots[0], L, unilateral: false }, { a: b2, b: b1, L, unilateral: false }];
    stats = [newStats(), newStats()];
    Sim.setView([-2.6, 0, 2.6, 4.4]);
  } else {
    // Both bobs start at the bottom moving sideways: one on a rigid rod, one on a rope.
    pivots = [particle(0, 2.3, 0)];
    const rod = particle(0, 1.3, mass), rope = particle(0, 1.3, mass);
    rod.v = V.v(launch.value, 0); rope.v = V.v(launch.value, 0);
    bodies = [rod, rope];
    links = [{ a: rod, b: pivots[0], L, unilateral: false }, { a: rope, b: pivots[0], L, unilateral: true }];
    stats = [newStats(), newStats()];
    Sim.setView([-1.8, 0.9, 1.8, 3.5]);
  }
}

// Enforce |xa - xb| = L (or |xa - xb| <= L when unilateral) by moving the two particles along the gradient of
// C = |xa - xb| - L, each in proportion to its inverse mass. Returns the position multiplier lambda (kg m);
// lambda / h^2 is the constraint force magnitude (the tension, for a rod) during this step.
function projectDistance(c) {
  const a = c.a, b = c.b;
  // TODO 1: C = |xa - xb| - c.L and n = (xa - xb) / |xa - xb| (the gradient of C with respect to xa).
  //         lambda = -C / (wa + wb). Move a by +wa * lambda * n and b by -wb * lambda * n. Return lambda.
  //         Skip (return 0) if wa + wb is 0.
  // TODO 4: if c.unilateral and C <= 0, the rope is slack: return 0 without moving anything.
  return 0;
}

function stepSpring(p, pivot, h) {
  // Symplectic Euler with a stiff spring holding the bob at distance L from the pivot.
  const d = V.sub(p.x, pivot.x), len = V.len(d), n = V.scale(d, 1 / len);
  const f = V.addScaled(V.v(0, -p.m * g), n, -stiffness.value * (len - L));
  p.v = V.addScaled(p.v, f, h / p.m);
  p.x = V.addScaled(p.x, p.v, h);
}

function stepConstrained(h) {
  // 1. Move every free particle as if there were no constraint.
  for (const p of bodies) {
    p.prev = V.copy(p.x);
    p.v = V.addScaled(p.v, V.v(0, -g), h);
    p.x = V.addScaled(p.x, p.v, h);
  }
  // 2. Project the constraints. One link needs one pass; coupled links need several (see chapter 07).
  let lambdas = [];
  if (scene.value === 'double pendulum') {
    lambdas = [0, 0];
    for (let it = 0; it < 20; it++) {
      // TODO 3: project both links (links[0]: pivot to bob 1, links[1]: bob 1 to bob 2) and add each
      //         returned lambda to lambdas[0] and lambdas[1].
    }
  } else {
    for (const c of links) lambdas.push(projectDistance(c));
  }
  // 3. The velocity is whatever moved the particle from where it was to where it ended up.
  for (const p of bodies) p.v = V.scale(V.sub(p.x, p.prev), 1 / h);
  return lambdas;
}

// Times successive left-to-right passes through the bottom; also records the tension there.
function track(st, p, pivot, prevX, t, h, tension) {
  const rel0 = prevX - pivot.x.x, rel1 = p.x.x - pivot.x.x;
  if (rel0 < 0 && rel1 >= 0 && p.x.y < pivot.x.y) {
    const tc = t + h * (-rel0 / (rel1 - rel0));
    if (st.lastCross >= 0) st.periods.push(tc - st.lastCross);
    st.lastCross = tc;
    st.bottomT = tension;
    st.bottomPred = p.m * g + p.m * V.lenSq(p.v) / L;
  }
  st.lastTension = tension;
  st.maxErr = Math.max(st.maxErr, Math.abs(V.dist(p.x, pivot.x) - L));
}

function step(h) {
  const t = Sim.time;
  if (scene.value === 'pendulum: spring vs constraint') {
    const sx = spring.x.x, cx = bodies[0].x.x;
    stepSpring(spring, pivots[0], h);
    const lam = stepConstrained(h);
    track(stats[0], spring, pivots[0], sx, t, h, stiffness.value * (V.dist(spring.x, pivots[0].x) - L));
    track(stats[1], bodies[0], pivots[1], cx, t, h, -lam[0] / (h * h));
  } else if (scene.value === 'double pendulum') {
    stepConstrained(h);
  } else {
    const lam = stepConstrained(h);
    for (let i = 0; i < 2; i++) {
      stats[i].lastTension = -lam[i] / (h * h);
      if (stats[i].slackAt === null && Sim.time > 0.05 && stats[i].lastTension === 0) {
        const r = V.sub(bodies[i].x, pivots[0].x);
        stats[i].slackAt = Math.atan2(r.x, -r.y) * 180 / Math.PI;   // angle from straight down
      }
    }
  }
  if (Sim.frame % 2 === 0) {
    for (let i = 0; i < bodies.length && i < 2; i++) {
      trails[i].push(V.copy(bodies[i].x));
      if (trails[i].length > 300) trails[i].shift();
    }
  }
}

function exactPeriod(theta0) {
  // T = 2 pi sqrt(L/g) / AGM(1, cos(theta0 / 2)), the exact large-angle period.
  let a = 1, b = Math.cos(theta0 / 2);
  for (let i = 0; i < 20; i++) { const an = (a + b) / 2; b = Math.sqrt(a * b); a = an; }
  return 2 * Math.PI * Math.sqrt(L / g) / a;
}

function lastPeriod(st) { return st.periods.length ? st.periods[st.periods.length - 1].toFixed(4) + ' s' : '...'; }

function drawBob(gr, pivot, p, color, dashed) {
  gr.line(pivot.x.x, pivot.x.y, p.x.x, p.x.y, { color: Sim.colors.ink, width: 2, dash: dashed ? [4, 4] : undefined });
  gr.dot(pivot.x.x, pivot.x.y, 4, Sim.colors.ink);
  gr.dot(p.x.x, p.x.y, 10, color);
}

function draw(gr) {
  const h = stepMs.value / 1000;
  const colors = [Sim.colors.blue, Sim.colors.green];
  for (let i = 0; i < trails.length; i++) {
    for (const q of trails[i]) gr.dot(q.x, q.y, 1.5, colors[i]);
  }
  if (scene.value === 'pendulum: spring vs constraint') {
    for (const pv of pivots) gr.circle(pv.x.x, pv.x.y, L, { stroke: Sim.colors.axis, dash: [3, 5], width: 1 });
    drawBob(gr, pivots[0], spring, Sim.colors.red);
    drawBob(gr, pivots[1], bodies[0], Sim.colors.blue);
    gr.text('stiff spring', -1.3, 3.45, { color: Sim.colors.red, align: 'center' });
    gr.text('constraint projection', 1.3, 3.45, { color: Sim.colors.blue, align: 'center' });
    const lenS = V.dist(spring.x, pivots[0].x) - L, lenC = V.dist(bodies[0].x, pivots[1].x) - L;
    const a0 = angleDeg.value * Math.PI / 180;
    Sim.show('h * sqrt(k/m) (spring)', (h * Math.sqrt(stiffness.value / mass)).toFixed(3));
    Sim.show('length error: spring / constraint', (lenS * 1000).toFixed(3) + ' / ' + (lenC * 1000).toFixed(6) + ' mm');
    Sim.show('tension now: spring / constraint', stats[0].lastTension.toFixed(2) + ' / ' + stats[1].lastTension.toFixed(2) + ' N');
    if (stats[1].bottomT !== null) {
      Sim.show('tension at bottom (constraint)', stats[1].bottomT.toFixed(3) + ' N');
      Sim.show('m g + m v^2 / L there', stats[1].bottomPred.toFixed(3) + ' N');
    }
    Sim.show('period: spring / constraint', lastPeriod(stats[0]) + ' / ' + lastPeriod(stats[1]));
    Sim.show('small-angle 2 pi sqrt(L/g)', (2 * Math.PI * Math.sqrt(L / g)).toFixed(4) + ' s');
    Sim.show('exact for this angle', exactPeriod(a0).toFixed(4) + ' s');
    Sim.plot('spring stretch (mm)', Math.max(-50, Math.min(50, lenS * 1000)), { color: Sim.colors.red });
    Sim.plot('constraint error (mm)', lenC * 1000, { color: Sim.colors.blue });
  } else if (scene.value === 'double pendulum') {
    drawBob(gr, pivots[0], bodies[0], Sim.colors.blue);
    drawBob(gr, bodies[0], bodies[1], Sim.colors.green);
    let E = 0;
    for (const p of bodies) E += 0.5 * p.m * V.lenSq(p.v) + p.m * g * p.x.y;
    Sim.show('link errors', ((V.dist(bodies[0].x, pivots[0].x) - L) * 1000).toFixed(4) + ', ' + ((V.dist(bodies[1].x, bodies[0].x) - L) * 1000).toFixed(4) + ' mm');
    Sim.show('energy', E.toFixed(3) + ' J');
    Sim.plot('energy (J)', E, { color: Sim.colors.ink });
  } else {
    gr.circle(pivots[0].x.x, pivots[0].x.y, L, { stroke: Sim.colors.axis, dash: [3, 5], width: 1 });
    drawBob(gr, pivots[0], bodies[0], Sim.colors.blue);
    const slack = V.dist(bodies[1].x, pivots[0].x) < L - 1e-6;
    drawBob(gr, pivots[0], bodies[1], Sim.colors.green, slack);
    const v0 = launch.value, c = (2 * g * L - v0 * v0) / (3 * g * L);
    Sim.show('rod (blue) tension', stats[0].lastTension.toFixed(2) + ' N');
    Sim.show('rope (green) tension', stats[1].lastTension.toFixed(2) + ' N' + (slack ? '  slack' : ''));
    Sim.show('predicted slack angle', Math.abs(c) <= 1 ? (Math.acos(c) * 180 / Math.PI).toFixed(1) + ' deg' : 'never');
    Sim.show('rope went slack at', stats[1].slackAt === null ? '...' : Math.abs(stats[1].slackAt).toFixed(1) + ' deg');
  }
}

Sim.run({
  title: 'Constraints and multipliers (starter)',
  view: [-2.8, 0.9, 2.8, 3.5],
  reset, step, draw,
  dt: () => stepMs.value / 1000,
  stepsPerFrame: () => Math.max(1, Math.round(1000 / 60 / stepMs.value)),
  grid: 0.25,
});
