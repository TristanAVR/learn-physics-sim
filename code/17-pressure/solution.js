// Chapter 17 - Pressure and Inflated Bodies (reference solution)
// A 2D balloon: a closed loop of particles held by XPBD distance constraints (the wall), filled either with
// an ideal gas (pressure P = nRT / A, applied as a force on every edge) or kept at a fixed area by a
// global XPBD area constraint. Squeeze it with the plunger (slider or mouse), or drop it.
// 2D units: pressure is force per length of wall (N/m), tension is a force (N), area is in m^2.

const GAS = 'ideal gas (P A = const)', AREA = 'area constraint';
const model = Sim.select('gas model', [GAS, AREA], GAS, { reset: true });
const P0 = Sim.slider('initial pressure P0 (N/m)', 50, 2000, 500, { step: 10, reset: true });
const R0 = Sim.slider('balloon radius R0 (m)', 0.15, 0.45, 0.3, { step: 0.01, reset: true });
const wallY = Sim.logSlider('wall stretch modulus Y (N)', 2.5, 3.5, 2000, { reset: true });
const gravity = Sim.checkbox('gravity', true);
const areaC = Sim.logSlider('area compliance (m^3/N)', -9, -3, 1e-9);
const plunger = Sim.slider('plunger height (m)', 0.15, 1.5, 1.4, { step: 0.01 });
Sim.button('drop from 1 m', () => {
  for (let i = 0; i < N; i++) { y[i] += 1.0; vx[i] = 0; vy[i] = 0; }
  plunger.set(1.5); plungerY = 1.5;
});

const G = 9.81;
const N = 60;                 // wall particles
const WALL_MASS = 1.0;        // kg in total (a rubber wall, per metre of depth)
const SUBSTEPS = 20;
const ITERATIONS = 20;        // passes over the wall (and area) constraints per substep: a stiff hoop converges slowly

const x = new Float64Array(N), y = new Float64Array(N);
const px = new Float64Array(N), py = new Float64Array(N);
const vx = new Float64Array(N), vy = new Float64Array(N);
const w = new Float64Array(N);
const edgeLambda = new Float64Array(N);   // multiplier of the wall constraint i -> i+1, this substep
let restLen = 0, nRT = 0, area0 = 0, areaLambda = 0, plungerForce = 0, plungerContacts = 0;
let plungerY = 1.4;           // where the plunger is; it follows the slider at a limited speed

// ------------------------------------------------------------------------------------------------
// Geometry of the closed outline (particles in counter-clockwise order)
// ------------------------------------------------------------------------------------------------
function polygonArea() {
  // shoelace formula: A = 1/2 sum (x_i y_{i+1} - x_{i+1} y_i)
  let a = 0;
  for (let i = 0; i < N; i++) {
    const j = (i + 1) % N;
    a += x[i] * y[j] - x[j] * y[i];
  }
  return 0.5 * a;
}

function currentPressure() {
  if (model.value === GAS) return nRT / polygonArea();
  return areaLambda;          // the area constraint's multiplier, already divided by h^2: see solveArea
}

// ------------------------------------------------------------------------------------------------
// Ideal gas: pressure force on every edge, P times the edge length along the outward normal,
// shared equally by its two particles.
// ------------------------------------------------------------------------------------------------
function applyPressure(h) {
  const P = nRT / polygonArea();
  for (let i = 0; i < N; i++) {
    const j = (i + 1) % N;
    const ex = x[j] - x[i], ey = y[j] - y[i];
    // outward normal times length, for a counter-clockwise outline: (ey, -ex)
    const fx = 0.5 * P * ey, fy = -0.5 * P * ex;
    vx[i] += h * w[i] * fx; vy[i] += h * w[i] * fy;
    vx[j] += h * w[j] * fx; vy[j] += h * w[j] * fy;
  }
}

// ------------------------------------------------------------------------------------------------
// Area constraint C = A - A0 with compliance alpha (m^3/N). The gradient of the shoelace area with
// respect to particle i is 1/2 (y_{i+1} - y_{i-1}, x_{i-1} - x_{i+1}): the outward normals of the two
// neighbouring edges times half their lengths.
// ------------------------------------------------------------------------------------------------
const gradX = new Float64Array(N), gradY = new Float64Array(N);
function solveArea(alphaTilde, h) {
  const C = polygonArea() - area0;
  let denom = alphaTilde;
  for (let i = 0; i < N; i++) {
    const next = (i + 1) % N, prev = (i + N - 1) % N;
    gradX[i] = 0.5 * (y[next] - y[prev]);
    gradY[i] = 0.5 * (x[prev] - x[next]);
    denom += w[i] * (gradX[i] * gradX[i] + gradY[i] * gradY[i]);
  }
  const dlambda = (-C - alphaTilde * areaLambdaRaw) / denom;
  areaLambdaRaw += dlambda;
  for (let i = 0; i < N; i++) {
    x[i] += w[i] * gradX[i] * dlambda;
    y[i] += w[i] * gradY[i] * dlambda;
  }
  areaLambda = areaLambdaRaw / (h * h);   // force per unit of area gradient = the pressure
}
let areaLambdaRaw = 0;

// ------------------------------------------------------------------------------------------------
// The wall: XPBD distance constraints between neighbours
// ------------------------------------------------------------------------------------------------
function solveEdge(i, alphaTilde) {
  const j = (i + 1) % N;
  const dx = x[j] - x[i], dy = y[j] - y[i], d = Math.hypot(dx, dy);
  if (d < 1e-12) return;
  const C = d - restLen;
  const dl = (-C - alphaTilde * edgeLambda[i]) / (w[i] + w[j] + alphaTilde);
  edgeLambda[i] += dl;
  const nx = dx / d, ny = dy / d;
  x[i] -= w[i] * dl * nx; y[i] -= w[i] * dl * ny;
  x[j] += w[j] * dl * nx; y[j] += w[j] * dl * ny;
}

// Even edges first, then odd ones. Sweeping round the ring in one direction would drag the particles
// along the sweep a little every pass, and a pre-stressed ring would slowly start to spin.
function solveWall(alphaTilde) {
  for (let i = 0; i < N; i += 2) solveEdge(i, alphaTilde);
  for (let i = 1; i < N; i += 2) solveEdge(i, alphaTilde);
}

// ------------------------------------------------------------------------------------------------
// Scene
// ------------------------------------------------------------------------------------------------
function reset() {
  const r = R0.value;
  for (let i = 0; i < N; i++) {
    const a = 2 * Math.PI * i / N - Math.PI / 2;
    x[i] = r * Math.cos(a); y[i] = r + r * Math.sin(a);       // resting on the ground
    vx[i] = 0; vy[i] = 0; w[i] = N / WALL_MASS;
  }
  // The wall is pre-stretched so that, at radius R0 and pressure P0, its tension is P0 R0 (Laplace).
  const chord = 2 * r * Math.sin(Math.PI / N), T0 = P0.value * r;
  restLen = chord / (1 + T0 / wallY.value);
  area0 = polygonArea();
  nRT = P0.value * area0;              // ideal gas: P A = nRT, fixed by the initial state
  areaLambdaRaw = 0; areaLambda = P0.value;
  plungerY = plunger.value;
  Sim.setView([-1.0, -0.1, 1.0, 1.55]);
}

function step(dt) {
  const h = dt / SUBSTEPS, h2 = h * h;
  if (Sim.mouse.down && Sim.mouse.inside) plunger.set(Math.max(0.15, Math.min(1.5, Sim.mouse.y)));
  const PLUNGER_SPEED = 0.5;   // m/s; a plunger that jumped would push the top of the wall through the bottom
  const wallAlpha = restLen / wallY.value;           // compliance of one wall segment (m/N)
  plungerForce = 0;
  for (let s = 0; s < SUBSTEPS; s++) {
    if (model.value === GAS) applyPressure(h);
    for (let i = 0; i < N; i++) {
      px[i] = x[i]; py[i] = y[i];
      if (gravity.value) vy[i] -= G * h;
      x[i] += h * vx[i]; y[i] += h * vy[i];
    }
    edgeLambda.fill(0);
    areaLambdaRaw = 0;
    for (let it = 0; it < ITERATIONS; it++) {
      solveWall(wallAlpha / h2);
      if (model.value === AREA) solveArea(areaC.value / h2, h);
    }
    // contacts: the ground below, the plunger above (both frictionless)
    plungerY += Math.max(-PLUNGER_SPEED * h, Math.min(PLUNGER_SPEED * h, plunger.value - plungerY));
    plungerContacts = 0;
    for (let i = 0; i < N; i++) {
      if (y[i] < 0) y[i] = 0;
      if (y[i] > plungerY) {
        plungerForce += (y[i] - plungerY) / (w[i] * h2) / SUBSTEPS;   // m d / h^2, averaged over substeps
        y[i] = plungerY;
        plungerContacts++;
      }
    }
    for (let i = 0; i < N; i++) { vx[i] = (x[i] - px[i]) / h; vy[i] = (y[i] - py[i]) / h; }
  }
  lastH = h;
}
let lastH = 1 / 60 / SUBSTEPS;

// ------------------------------------------------------------------------------------------------
// Drawing and measurements
// ------------------------------------------------------------------------------------------------
function draw(g) {
  g.ground(0);
  const pts = [];
  for (let i = 0; i < N; i++) pts.push({ x: x[i], y: y[i] });
  g.polyline(pts, { fill: '#f3c9b8', closed: true });
  g.polyline(pts, { color: Sim.colors.red, width: 2, closed: true });
  g.rect(-1.2, plungerY, 2.4, 0.04, { fill: Sim.colors.gray });
  g.line(0, plungerY + 0.04, 0, 1.7, { color: Sim.colors.gray, width: 4 });

  const A = polygonArea(), P = currentPressure();
  // centroid and mean radius
  let cx = 0, cy = 0;
  for (let i = 0; i < N; i++) { cx += x[i]; cy += y[i]; }
  cx /= N; cy /= N;
  let r = 0, T = 0, xmin = Infinity, xmax = -Infinity;
  for (let i = 0; i < N; i++) {
    r += Math.hypot(x[i] - cx, y[i] - cy);
    T += -edgeLambda[i] / (lastH * lastH);          // tension = multiplier / h^2 (positive when stretched)
    xmin = Math.min(xmin, x[i]); xmax = Math.max(xmax, x[i]);
  }
  r /= N; T /= N;
  Sim.show('area A / A0', (A / area0).toFixed(4));
  Sim.show('pressure P', P.toFixed(1) + ' N/m');
  Sim.show('mean wall tension T', T.toFixed(2) + ' N');
  if (plungerContacts === 0) {
    Sim.show('P x mean radius (Laplace: = T)', (P * r).toFixed(2) + ' N');
  } else {
    // squeezed between two plates, the free sides are half-circles spanning the gap
    Sim.show('T / P (Laplace: radius of the sides)', (T / P).toFixed(4) + ' m');
    Sim.show('half the gap', (plungerY / 2).toFixed(4) + ' m');
    Sim.show('plunger force', plungerForce.toFixed(2) + ' N');
    Sim.show('P x width - 2T', (P * (xmax - xmin) - 2 * T).toFixed(2) + ' N');
  }
  Sim.plot('pressure / P0', P / P0.value, { color: Sim.colors.purple });
  Sim.plot('area / A0', A / area0, { color: Sim.colors.teal });
}

Sim.run({ title: 'Pressure', view: [-1.0, -0.1, 1.0, 1.55], reset, step, draw });
