// Chapter 04 - Springs, Damping and Energy (starter)
// Particles joined by spring-dampers, advanced with symplectic Euler. Drag a particle with the mouse.
// Scenes: a soft square (spring grid with diagonals), a hanging chain, and a single spring for measurements.
// springForce returns zero, so the particles fall as loose dust. Fill in the TODOs.

const scene = Sim.select('scene', ['soft square', 'hanging chain', 'single spring'], 'soft square', { reset: true });
const stiffness = Sim.logSlider('spring stiffness k (N/m)', 1, 4.5, 300);
const damping = Sim.logSlider('damper c (N s/m)', -3, 2, 0.05);
const substeps = Sim.slider('substeps per frame', 1, 40, 10, { step: 1 });
const diagonals = Sim.checkbox('diagonal springs (square)', true, { reset: true });
const maxwell = Sim.checkbox('Maxwell springs (they flow)', false, { reset: true });
const relaxTime = Sim.slider('Maxwell relaxation time tau (s)', 0.1, 5, 1, { step: 0.1 });
const holdStretch = Sim.checkbox('hold 20 cm stretch (single spring)', false, { reset: true });

const gravity = 9.81;
let particles, springs, dragged, lastMouse, throwVel, startU, minU, omegaMax, holdTime;

function addParticle(x, y, mass, pinned) {
  particles.push({ x: V.v(x, y), v: V.v(0, 0), f: V.v(0, 0), m: mass, w: pinned ? 0 : 1 / mass });
  return particles.length - 1;
}

function addSpring(a, b) {
  const rest = V.dist(particles[a].x, particles[b].x);
  springs.push({ a, b, rest0: rest, rest });
}

function reset() {
  particles = []; springs = [];
  dragged = -1; lastMouse = null; throwVel = V.v(0, 0); holdTime = 0;
  if (scene.value === 'soft square') {
    const n = 8, size = 1.0, mass = 1.0 / (n * n), s = size / (n - 1);   // 1 kg in total
    const angle = 0.1, cx = -0.3, cy = 1.4;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const p = V.rotate(V.v(j * s - size / 2, i * s - size / 2), angle);
        addParticle(cx + p.x, cy + p.y, mass, false);
      }
    }
    const id = (i, j) => i * n + j;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (j < n - 1) addSpring(id(i, j), id(i, j + 1));
        if (i < n - 1) addSpring(id(i, j), id(i + 1, j));
        if (diagonals.value && i < n - 1 && j < n - 1) {
          addSpring(id(i, j), id(i + 1, j + 1));
          addSpring(id(i, j + 1), id(i + 1, j));
        }
      }
    }
    Sim.setView([-1.6, -0.3, 1.6, 2.5]);
  } else if (scene.value === 'hanging chain') {
    const n = 20, link = 0.1, mass = 0.025;            // 20 links of 10 cm, 0.5 kg in total
    addParticle(-1.0, 2.0, mass, true);
    for (let i = 1; i <= n; i++) {
      addParticle(-1.0 + i * link, 2.0, mass, false);
      addSpring(i - 1, i);
    }
    Sim.setView([-2.3, -0.3, 2.3, 2.5]);
  } else {
    // One 1 kg mass below a fixed anchor. Start 20 cm below the equilibrium point, at rest.
    const rest = 0.8, mass = 1.0;
    addParticle(0, 2.2, 0, true);
    const stretch = holdStretch.value ? 0.2 : mass * gravity / stiffness.value + 0.2;
    addParticle(0, 2.2 - rest - stretch, mass, holdStretch.value);
    springs.push({ a: 0, b: 1, rest0: rest, rest });
    Sim.setView([-1.5, -0.3, 1.5, 2.5]);
  }
  startU = displacement();
  minU = startU;
  modeX = null;
  omegaMax = fastestOmega(300);
}

// Force on particle a from the spring-damper (a, b). Particle b receives the opposite force.
// Hooke along the spring axis, plus a damper on the rate of stretching (relative velocity along the axis).
function springForce(s) {
  const pa = particles[s.a], pb = particles[s.b];
  // TODO 1: unit vector n from a to b, the stretch (current length - s.rest), and the rate of stretching
  //         (relative velocity pb.v - pa.v projected on n). Tension = k * stretch + c * rate,
  //         with k = stiffness.value and c = damping.value. Return the force on a: tension * n.
  return V.v(0, 0);
}

// Maxwell element: a dashpot in series with each spring. The dashpot lets the rest length creep
// toward the current length with time constant tau = c_M / k (solved exactly for one step).
function relaxRestLength(s, h) {
  // TODO 4: move s.rest toward the current length of the spring, as rest += (len - rest) * h / tau,
  //         or exactly: rest = len + (rest - len) * exp(-h / tau), with tau = relaxTime.value.
}

function step(dt) {
  const h = dt / substeps.value;
  for (let k = 0; k < substeps.value; k++) {
    for (const p of particles) p.f = V.v(0, -p.m * gravity);
    for (const s of springs) {
      const f = springForce(s);
      particles[s.a].f = V.add(particles[s.a].f, f);
      particles[s.b].f = V.sub(particles[s.b].f, f);
    }
    // Symplectic Euler: velocity first, then position with the new velocity.
    for (let i = 0; i < particles.length; i++) {
      const p = particles[i];
      if (p.w === 0 || i === dragged) continue;
      p.v = V.addScaled(p.v, p.f, h * p.w);
      p.x = V.addScaled(p.x, p.v, h);
      if (p.x.y < 0) { p.x.y = 0; p.v = V.v(0, Math.max(0, p.v.y)); }   // sticky, inelastic ground
    }
    if (maxwell.value) for (const s of springs) relaxRestLength(s, h);
  }
  if (holdStretch.value) holdTime += dt;
  const u = displacement();
  if (u < minU) minU = u;
}

// Single-spring scene: displacement of the mass below its equilibrium point (m).
function displacement() {
  if (scene.value !== 'single spring') return 0;
  const p = particles[1], eq = particles[0].x.y - springs[0].rest - p.m * gravity / stiffness.value;
  return eq - p.x.y;
}

function handleMouse(dt) {
  const ms = Sim.mouse;
  if (ms.pressed) {
    let best = -1, bestD = 0.3;
    for (let i = 0; i < particles.length; i++) {
      const d = V.dist(particles[i].x, ms);
      if (particles[i].w > 0 && d < bestD) { best = i; bestD = d; }
    }
    dragged = best;
  }
  if (dragged >= 0) {
    const p = particles[dragged];
    if (lastMouse) throwVel = V.scale(V.sub(V.v(ms.x, ms.y), lastMouse), 1 / dt);
    p.x = V.v(ms.x, Math.max(0, ms.y));
    p.v = V.v(0, 0);
    if (ms.released || !ms.down) { p.v = throwVel; dragged = -1; }
  }
  lastMouse = V.v(ms.x, ms.y);
}

// Largest natural frequency of the network around its current shape: power iteration on M^-1 K, where K is
// the stiffness matrix linearised at the current positions (axial k, plus tension / length sideways).
let modeX = null, modeY = null;
function fastestOmega(iterations) {
  const N = particles.length;
  if (!modeX || modeX.length !== N) {
    modeX = new Float64Array(N); modeY = new Float64Array(N);
    const rng = Sim.random(3);
    for (let i = 0; i < N; i++) { modeX[i] = rng() - 0.5; modeY[i] = rng() - 0.5; }
  }
  const fx = new Float64Array(N), fy = new Float64Array(N);
  let lambda = 0;
  for (let it = 0; it < iterations; it++) {
    fx.fill(0); fy.fill(0);
    for (const s of springs) {
      const d = V.sub(particles[s.b].x, particles[s.a].x), len = V.len(d), n = V.scale(d, 1 / len);
      const dux = modeX[s.b] - modeX[s.a], duy = modeY[s.b] - modeY[s.a];
      const axial = dux * n.x + duy * n.y;
      const side = stiffness.value * (len - s.rest) / len;          // geometric stiffness T / L
      const tx = stiffness.value * axial * n.x + side * (dux - axial * n.x);
      const ty = stiffness.value * axial * n.y + side * (duy - axial * n.y);
      fx[s.a] += tx; fy[s.a] += ty; fx[s.b] -= tx; fy[s.b] -= ty;  // this is -K u
    }
    let dot = 0, un = 0, norm = 0;
    for (let i = 0; i < N; i++) {
      const w = i === dragged ? 0 : particles[i].w;
      fx[i] *= -w; fy[i] *= -w;                                     // M^-1 K u
      dot += modeX[i] * fx[i] + modeY[i] * fy[i];
      un += modeX[i] * modeX[i] + modeY[i] * modeY[i];
      norm += fx[i] * fx[i] + fy[i] * fy[i];
    }
    norm = Math.sqrt(norm);
    if (norm === 0) return 0;
    lambda = dot / un;
    for (let i = 0; i < N; i++) { modeX[i] = fx[i] / norm; modeY[i] = fy[i] / norm; }
  }
  return Math.sqrt(Math.max(lambda, 0));
}

function energies() {
  let kin = 0, spring = 0, grav = 0;
  for (const p of particles) {
    if (p.w === 0) continue;
    kin += 0.5 * p.m * V.lenSq(p.v);
    grav += p.m * gravity * p.x.y;
  }
  for (const s of springs) {
    const e = V.dist(particles[s.a].x, particles[s.b].x) - s.rest;
    spring += 0.5 * stiffness.value * e * e;
  }
  return { kin, spring, grav };
}

function draw(g) {
  handleMouse(1 / 60);
  g.ground(0);
  for (const s of springs) {
    const a = particles[s.a].x, b = particles[s.b].x;
    const strain = (V.dist(a, b) - s.rest0) / s.rest0;
    const color = strain > 0.02 ? Sim.colors.red : strain < -0.02 ? Sim.colors.blue : Sim.colors.gray;
    g.line(a.x, a.y, b.x, b.y, { color, width: 1.5 });
  }
  for (let i = 0; i < particles.length; i++) {
    const p = particles[i];
    const r = scene.value === 'single spring' ? 9 : 3.5;
    g.dot(p.x.x, p.x.y, r, p.w === 0 ? Sim.colors.ink : i === dragged ? Sim.colors.orange : Sim.colors.blue);
  }
  if (dragged >= 0) g.line(particles[dragged].x.x, particles[dragged].x.y, Sim.mouse.x, Sim.mouse.y, { color: Sim.colors.orange, dash: [3, 3] });

  if (!Sim.paused) omegaMax = fastestOmega(20);
  const h = 1 / 60 / substeps.value;
  Sim.show('h (substep)', (h * 1000).toFixed(3) + ' ms');
  Sim.show('omega_max', omegaMax.toFixed(1) + ' rad/s');
  Sim.show('h * omega_max', (h * omegaMax).toFixed(3) + (h * omegaMax > 2 ? '  (> 2: unstable)' : ''));

  const E = energies();
  if (scene.value === 'single spring') {
    const p = particles[1], k = stiffness.value, m = p.m;
    const s = springs[0], len = V.dist(particles[0].x, p.x);
    const force = k * (len - s.rest);
    Sim.show('zeta = c / (2 sqrt(k m))', (damping.value / (2 * Math.sqrt(k * m))).toFixed(3));
    Sim.show('critical c = 2 sqrt(k m)', (2 * Math.sqrt(k * m)).toFixed(3) + ' N s/m');
    if (holdStretch.value) {
      Sim.show('spring force', force.toFixed(3) + ' N');
      Sim.show('force / initial', (force / (k * 0.2)).toFixed(4));
      if (maxwell.value) Sim.show('exp(-t / tau)', Math.exp(-holdTime / relaxTime.value).toFixed(4));
      Sim.plot('spring force (N)', force, { color: Sim.colors.red, min: 0 });
    } else {
      Sim.show('displacement u', displacement().toFixed(4) + ' m');
      Sim.show('overshoot', (Math.max(0, -minU) / startU * 100).toFixed(2) + ' %');
      Sim.plot('u (m)', displacement(), { color: Sim.colors.blue });
    }
  } else {
    Sim.show('kinetic', E.kin.toFixed(3) + ' J');
    Sim.show('spring', E.spring.toFixed(3) + ' J');
    Sim.show('gravity', E.grav.toFixed(3) + ' J');
    Sim.plot('total E (J)', E.kin + E.spring + E.grav, { color: Sim.colors.ink, min: 0 });
    Sim.plot('kinetic', E.kin, { color: Sim.colors.red });
    Sim.plot('spring', E.spring, { color: Sim.colors.green });
  }
}

Sim.run({ title: 'Springs, damping and energy (starter)', view: [-1.6, -0.3, 1.6, 2.5], reset, step, draw, dt: 1 / 60 });
