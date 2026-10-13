// Chapter 10 - Rotations and Rigid Bodies: the tumbling box (starter)
// A box with three different principal moments of inertia, spinning in empty space with no torque.
// Spin it about its smallest or largest axis and it keeps spinning; spin it about the intermediate axis and
// any tiny perturbation grows until the box flips over, again and again (the "tennis racket" effect).
// Orientation is a unit quaternion. The 3D view is a fixed orthographic camera, drawn with the 2D sandbox.
// The box does not move yet: fill in the two TODOs in integrate().

const axisChoice = Sim.select('spin about', ['smallest I (x)', 'intermediate I (y)', 'largest I (z)'], 'intermediate I (y)', { reset: true });
const spinRate = Sim.slider('spin rate (rad/s)', 1, 10, 4, { step: 0.1, reset: true });
const perturbation = Sim.logSlider('perturbation', -4, -1, 1e-3, { reset: true });
const substeps = Sim.slider('substeps per frame', 1, 50, 20, { step: 1 });
const gyroscopic = Sim.checkbox('gyroscopic term', true);

// Box of 1.0 x 0.6 x 0.2 m (half extents below), 1 kg. Principal moments about the box axes:
// I_x = m (b^2 + c^2) / 3 for half extents (a, b, c), and so on.
const half = [0.5, 0.3, 0.1];
const mass = 1;
const inertia = [
  mass * (half[1] * half[1] + half[2] * half[2]) / 3,
  mass * (half[0] * half[0] + half[2] * half[2]) / 3,
  mass * (half[0] * half[0] + half[1] * half[1]) / 3,
];

let q, omegaW, E0, L0, start, t, growth;

// ---------------------------------------------------------------------------------------------
// Quaternions as { w, x, y, z }, vectors as [x, y, z].
// ---------------------------------------------------------------------------------------------
function qmul(a, b) {
  return {
    w: a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
    x: a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
    y: a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
    z: a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
  };
}
function qnormalize(a) {
  const l = Math.hypot(a.w, a.x, a.y, a.z);
  return { w: a.w / l, x: a.x / l, y: a.y / l, z: a.z / l };
}
function qconj(a) { return { w: a.w, x: -a.x, y: -a.y, z: -a.z }; }
// Rotate vector v by unit quaternion a: v' = a v a*, written out without building the products.
function qrotate(a, v) {
  const tx = 2 * (a.y * v[2] - a.z * v[1]), ty = 2 * (a.z * v[0] - a.x * v[2]), tz = 2 * (a.x * v[1] - a.y * v[0]);
  return [
    v[0] + a.w * tx + (a.y * tz - a.z * ty),
    v[1] + a.w * ty + (a.z * tx - a.x * tz),
    v[2] + a.w * tz + (a.x * ty - a.y * tx),
  ];
}
function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }

function reset() {
  q = { w: 1, x: 0, y: 0, z: 0 };
  const k = axisChoice.value.startsWith('smallest') ? 0 : axisChoice.value.startsWith('intermediate') ? 1 : 2;
  const w = [0, 0, 0];
  w[k] = spinRate.value;
  // A small nudge about both other axes.
  w[(k + 1) % 3] = perturbation.value * spinRate.value;
  w[(k + 2) % 3] = perturbation.value * spinRate.value;
  omegaW = w;
  const wb = bodyOmega();
  E0 = kineticEnergy(wb);
  L0 = angularMomentum(wb);
  start = perpendicular(wb);
  t = 0;
  growth = { t10: NaN, t100: NaN };
}

function bodyOmega() { return qrotate(qconj(q), omegaW); }
function kineticEnergy(wb) { return 0.5 * (inertia[0] * wb[0] * wb[0] + inertia[1] * wb[1] * wb[1] + inertia[2] * wb[2] * wb[2]); }
function angularMomentum(wb) { return Math.hypot(inertia[0] * wb[0], inertia[1] * wb[1], inertia[2] * wb[2]); }
function spinIndex() { return axisChoice.value.startsWith('smallest') ? 0 : axisChoice.value.startsWith('intermediate') ? 1 : 2; }
function perpendicular(wb) { const k = spinIndex(); return Math.hypot(wb[(k + 1) % 3], wb[(k + 2) % 3]); }

// One substep of torque-free rigid rotation.
function integrate(h) {
  // Euler's equations in the body frame: I dw/dt = tau - w x (I w), with tau = 0 here.
  let wb = bodyOmega();
  if (gyroscopic.value) {
    // TODO (exercise 3, part 2): the gyroscopic term. Compute I w (componentwise, the body frame is
    // principal), then g = w x (I w), and update each component: w_k <- w_k - h g_k / I_k.
  }
  omegaW = qrotate(q, wb);
  // TODO (exercise 3, part 1): integrate the orientation: q <- normalize(q + h/2 (0, omega) q), with omega
  // in world coordinates (omegaW). Use qmul and qnormalize.
}

function step(dt) {
  const n = substeps.value;
  for (let i = 0; i < n; i++) integrate(dt / n);
  t += dt;
  // Measure the growth rate of the perturbation while it is still small.
  const p = perpendicular(bodyOmega());
  if (isNaN(growth.t10) && p > 10 * start) growth.t10 = t;
  if (isNaN(growth.t100) && p > 100 * start) growth.t100 = t;
}

// ---------------------------------------------------------------------------------------------
// Drawing: a fixed orthographic camera. camera() maps a world vector to [screen x, screen y, depth].
// ---------------------------------------------------------------------------------------------
const yaw = -0.6, pitch = 0.45;
function camera(v) {
  const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const x1 = cy * v[0] + sy * v[2], z1 = -sy * v[0] + cy * v[2];
  const y2 = cp * v[1] - sp * z1, z2 = sp * v[1] + cp * z1;
  return [x1, y2, z2];
}

const faceColors = [Sim.colors.red, Sim.colors.green, Sim.colors.blue];

function shade(hex, f) {
  const r = parseInt(hex.slice(1, 3), 16), g = parseInt(hex.slice(3, 5), 16), b = parseInt(hex.slice(5, 7), 16);
  const m = (c) => Math.round(255 - (255 - c) * f);
  return `rgb(${m(r)}, ${m(g)}, ${m(b)})`;
}

function draw(g) {
  // Faces: for each axis k and sign s, the four corners in cyclic order. Draw the ones facing the camera.
  for (let k = 0; k < 3; k++) {
    for (const s of [-1, 1]) {
      const normal = [0, 0, 0]; normal[k] = s;
      const nc = camera(qrotate(q, normal));
      if (nc[2] <= 0) continue;
      const u = (k + 1) % 3, v = (k + 2) % 3, pts = [];
      for (const [su, sv] of [[-1, -1], [1, -1], [1, 1], [-1, 1]]) {
        const corner = [0, 0, 0];
        corner[k] = s * half[k]; corner[u] = su * half[u]; corner[v] = sv * half[v];
        const c = camera(qrotate(q, corner));
        pts.push({ x: c[0], y: c[1] });
      }
      g.polyline(pts, { closed: true, fill: shade(faceColors[k], 0.35 + 0.65 * nc[2]), stroke: Sim.colors.ink, width: 1 });
    }
  }
  // The body's spin axis, and the (fixed) direction of the angular momentum.
  const k = spinIndex(), axis = [0, 0, 0];
  axis[k] = 0.75;
  const a1 = camera(qrotate(q, axis)), a0 = camera(qrotate(q, axis.map((x) => -x)));
  g.line(a0[0], a0[1], a1[0], a1[1], { color: Sim.colors.ink, width: 2.5 });
  g.dot(a1[0], a1[1], 4, faceColors[k]);
  const wb = bodyOmega();
  const Lb = [inertia[0] * wb[0], inertia[1] * wb[1], inertia[2] * wb[2]];
  const Lw = qrotate(q, Lb), Ll = Math.hypot(Lw[0], Lw[1], Lw[2]);
  const Lc = camera([0.95 * Lw[0] / Ll, 0.95 * Lw[1] / Ll, 0.95 * Lw[2] / Ll]);
  g.line(-Lc[0], -Lc[1], Lc[0], Lc[1], { color: Sim.colors.muted, dash: [5, 4], width: 1 });
  g.text('L', Lc[0] + 0.03, Lc[1], { color: Sim.colors.muted });

  const W = spinRate.value;
  Sim.plot('wx/W', wb[0] / W, { color: Sim.colors.red, min: -1.05, max: 1.05 });
  Sim.plot('wy/W', wb[1] / W, { color: Sim.colors.green });
  Sim.plot('wz/W', wb[2] / W, { color: Sim.colors.blue });
  Sim.show('I_x, I_y, I_z (kg m^2)', inertia.map((x) => x.toFixed(4)).join(', '));
  Sim.show('energy / start', (kineticEnergy(wb) / E0).toFixed(5));
  Sim.show('|L| / start', (angularMomentum(wb) / L0).toFixed(5));
  const [I1, I2, I3] = inertia;
  if (k === 1) {
    const sigma = W * Math.sqrt((I2 - I1) * (I3 - I2) / (I1 * I3));
    Sim.show('predicted growth rate', sigma.toFixed(3) + ' 1/s');
    const meas = Math.log(10) / (growth.t100 - growth.t10);
    Sim.show('measured growth rate', isFinite(meas) ? meas.toFixed(3) + ' 1/s' : 'measuring...');
  } else {
    // Stable axes: the perturbation wobbles with angular frequency W sqrt((Ik - Ij)(Ik - Il) / (Ij Il)).
    const j = (k + 1) % 3, l = (k + 2) % 3;
    const f = W * Math.sqrt((inertia[k] - inertia[j]) * (inertia[k] - inertia[l]) / (inertia[j] * inertia[l]));
    Sim.show('stable: wobble frequency', f.toFixed(3) + ' rad/s');
  }
}

Sim.run({ title: 'The tumbling box (starter)', view: [-1.4, -1.7, 2.4, 1.3], reset, step, draw, grid: false });
