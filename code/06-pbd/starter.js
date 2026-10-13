// Chapter 06 - Position-Based Dynamics (starter)
// A rope of particles joined by distance constraints, solved the PBD way (Mueller et al. 2007):
// predict positions, project the constraints a few times, then derive velocities. Drag the rope with the mouse.
// distanceCorrection returns nothing yet, so the rope falls apart. Fill in the TODOs.

const scene = Sim.select('scene', ['rope, one end pinned', 'rope, both ends pinned'], 'rope, one end pinned', { reset: true });
const solver = Sim.select('solver', ['Gauss-Seidel', 'Jacobi'], 'Gauss-Seidel');
const iterations = Sim.slider('iterations per step', 1, 100, 10, { step: 1 });
const stiffnessK = Sim.slider('PBD stiffness k (0 to 1)', 0.01, 1, 1, { step: 0.01 });
const correctK = Sim.checkbox("use k' = 1 - (1 - k)^(1/n)", false);
const omega = Sim.slider('Jacobi over-relaxation w', 1, 2, 1, { step: 0.05 });
const substeps = Sim.slider('steps per frame', 1, 8, 1, { step: 1 });

const N = 30, segment = 0.08, mass = 0.05, gravity = 9.81;   // a 2.32 m rope of 1.5 kg
let particles, constraints, dragged, dx, dy, counts;

function reset() {
  particles = [];
  constraints = [];
  dragged = -1;
  const both = scene.value === 'rope, both ends pinned';
  const span = both ? 2.0 : (N - 1) * segment;
  for (let i = 0; i < N; i++) {
    // Start straight and horizontal; with both ends pinned, the ends are closer than the rope is long.
    const x = -span / 2 + span * i / (N - 1);
    const pinned = i === 0 || (both && i === N - 1);
    particles.push({ x: V.v(x, 2.0), prev: V.v(x, 2.0), v: V.v(0, 0), w: pinned ? 0 : 1 / mass, pinned });
  }
  for (let i = 0; i < N - 1; i++) constraints.push({ a: i, b: i + 1, rest: segment });
  dx = new Float64Array(N); dy = new Float64Array(N); counts = new Float64Array(N);
  Sim.setView([-1.9, -0.6, 1.9, 2.4]);
}

function stiffnessPerIteration() {
  const k = stiffnessK.value, n = iterations.value;
  return correctK.value ? 1 - Math.pow(1 - k, 1 / n) : k;
}

// Corrections that move particles a and b to satisfy |xa - xb| = rest, weighted by inverse mass and scaled by k.
// Returns them as { ax, ay, bx, by }, without applying them.
function distanceCorrection(c, k) {
  const pa = particles[c.a], pb = particles[c.b];
  // TODO 1: with d = xa - xb, C = |d| - rest and n = d / |d|, the projection moves
  //         a by -k * wa * C / (wa + wb) * n and b by +k * wb * C / (wa + wb) * n.
  //         Return those two vectors as { ax, ay, bx, by }; return null if wa + wb is 0.
  return null;
}

// Gauss-Seidel: project one constraint at a time, in place; later constraints see earlier corrections.
function solveGaussSeidel(k) {
  for (const c of constraints) {
    const corr = distanceCorrection(c, k);
    if (!corr) continue;
    const pa = particles[c.a], pb = particles[c.b];
    pa.x = V.v(pa.x.x + corr.ax, pa.x.y + corr.ay);
    pb.x = V.v(pb.x.x + corr.bx, pb.x.y + corr.by);
  }
}

// Jacobi: every constraint computes its correction from the same positions; corrections are summed per
// particle, averaged over the number of constraints touching it, scaled by w and applied together.
function solveJacobi(k) {
  dx.fill(0); dy.fill(0); counts.fill(0);
  // TODO 3: for every constraint, compute distanceCorrection from the CURRENT positions (do not apply it),
  //         add it to dx, dy of both particles and count the constraint in counts[] for both.
  //         Then move every particle with counts > 0 by omega.value * (dx, dy) / counts.
}

function step(dt) {
  const h = dt / substeps.value;
  const k = stiffnessPerIteration();
  for (let s = 0; s < substeps.value; s++) {
    // 1. Predict: apply external forces to the velocity, move freely.
    for (const p of particles) {
      p.prev = V.copy(p.x);
      if (p.w === 0) continue;
      p.v = V.v(p.v.x, p.v.y - gravity * h);
      p.x = V.addScaled(p.x, p.v, h);
    }
    if (dragged >= 0) particles[dragged].x = V.v(Sim.mouse.x, Sim.mouse.y);
    // 2. Project the constraints, a fixed number of times.
    for (let it = 0; it < iterations.value; it++) {
      if (solver.value === 'Gauss-Seidel') solveGaussSeidel(k); else solveJacobi(k);
    }
    // 3. Velocities from the change in position.
    for (const p of particles) p.v = V.scale(V.sub(p.x, p.prev), 1 / h);
  }
}

function handleMouse() {
  const ms = Sim.mouse;
  if (ms.pressed) {
    let best = -1, bestD = 0.25;
    for (let i = 0; i < N; i++) {
      const d = V.dist(particles[i].x, ms);
      if (!particles[i].pinned && d < bestD) { best = i; bestD = d; }
    }
    if (best >= 0) { dragged = best; particles[best].w = 0; }
  }
  if (dragged >= 0 && !ms.down) { particles[dragged].w = 1 / mass; dragged = -1; }
}

function ropeLength() {
  let L = 0;
  for (const c of constraints) L += V.dist(particles[c.a].x, particles[c.b].x);
  return L;
}

function draw(g) {
  handleMouse();
  const pts = particles.map((p) => p.x);
  g.polyline(pts, { color: Sim.colors.ink, width: 2 });
  for (let i = 0; i < N; i++) {
    const p = particles[i];
    g.dot(p.x.x, p.x.y, p.pinned ? 5 : 3, p.pinned ? Sim.colors.ink : i === dragged ? Sim.colors.orange : Sim.colors.blue);
  }
  const rest = (N - 1) * segment, ratio = ropeLength() / rest;
  const top = particles[0].x.y;
  Sim.show('h', (1000 / 60 / substeps.value).toFixed(2) + ' ms');
  Sim.show('k per iteration', stiffnessPerIteration().toFixed(4));
  Sim.show('length / rest length', ratio.toFixed(4));
  if (scene.value === 'rope, one end pinned') {
    const drop = top - particles[N - 1].x.y;
    Sim.show('end below pin', drop.toFixed(4) + ' m (rest length ' + rest.toFixed(2) + ')');
    Sim.plot('stretch (%)', (ratio - 1) * 100, { color: Sim.colors.red, min: 0 });
  } else {
    const sag = top - 0.5 * (particles[N / 2 - 1].x.y + particles[N / 2].x.y);   // N is even
    Sim.show('sag at middle', sag.toFixed(4) + ' m');
    Sim.plot('stretch (%)', (ratio - 1) * 100, { color: Sim.colors.red, min: 0 });
    Sim.plot('sag (cm)', sag * 100, { color: Sim.colors.blue });
  }
}

Sim.run({ title: 'Position-based dynamics: a rope (starter)', view: [-1.9, -0.6, 1.9, 2.4], reset, step, draw, dt: 1 / 60 });
