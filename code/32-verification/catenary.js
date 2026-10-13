// Chapter 32 - Verifying a Simulator: a catenary test for a PBD chain
// A chain of N equal links (total length 2 m) hangs between two pins 1.6 m apart. At rest it should lie on
// the catenary y = a cosh(x / a) + c, with a fixed by the length: 2 a sinh(s / (2 a)) = L.
// The test lets the chain settle from a V shape, then measures the largest vertical distance between the
// nodes and the catenary. The table runs the same test headless for four chain resolutions.

const nSel = Sim.select('links N', ['10', '20', '40', '80'], '20', { reset: true });
const iters = Sim.logSlider('Gauss-Seidel iterations', 0, 2.5, 50, { step: 0.05 });
const tolPct = Sim.slider('tolerance (% of sag)', 0.1, 5, 1, { step: 0.1 });
Sim.button('rerun the table', () => { table = runTable(); });

const span = 1.6, length = 2.0, substeps = 4, h = 1 / 60 / substeps, G = 9.81;
const a = catenaryParameter(span, length);
const sag = a * (Math.cosh(span / (2 * a)) - 1);

// Solve 2 a sinh(s / 2a) = L for a by bisection (the left side falls monotonically from infinity to s).
function catenaryParameter(s, L) {
  let lo = 0.01, hi = 100;
  for (let k = 0; k < 200; k++) {
    const mid = 0.5 * (lo + hi);
    if (2 * mid * Math.sinh(s / (2 * mid)) > L) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}
function catenaryY(x) { return a * Math.cosh(x / a) - a * Math.cosh(span / (2 * a)); }

// A chain of N links between pins at (-s/2, 0) and (s/2, 0), starting as a V of the right length.
function makeChain(N) {
  const n = N + 1, c = { N, x: new Float64Array(n), y: new Float64Array(n), px: new Float64Array(n),
    py: new Float64Array(n), vx: new Float64Array(n), vy: new Float64Array(n), L: length / N };
  const depth = Math.sqrt((length / 2) * (length / 2) - (span / 2) * (span / 2));
  for (let i = 0; i < n; i++) {
    const t = i / N;
    c.x[i] = -span / 2 + span * t;
    c.y[i] = -depth * (1 - Math.abs(2 * t - 1));
  }
  return c;
}

// One substep of PBD: gravity, predict, Gauss-Seidel distance projections, velocities, light damping.
// All interior nodes have the same mass, so the inverse masses cancel out of the corrections.
function substep(c, nIter) {
  const n = c.N + 1;
  for (let i = 1; i < n - 1; i++) {
    c.vy[i] -= G * h;
    c.px[i] = c.x[i]; c.py[i] = c.y[i];
    c.x[i] += h * c.vx[i]; c.y[i] += h * c.vy[i];
  }
  for (let it = 0; it < nIter; it++) {
    for (let k = 0; k < c.N; k++) {
      const a0 = k, b0 = k + 1;
      const wa = a0 === 0 ? 0 : 1, wb = b0 === n - 1 ? 0 : 1;
      const dx = c.x[b0] - c.x[a0], dy = c.y[b0] - c.y[a0], len = Math.hypot(dx, dy);
      const corr = (len - c.L) / ((wa + wb) * len);
      c.x[a0] += wa * dx * corr; c.y[a0] += wa * dy * corr;
      c.x[b0] -= wb * dx * corr; c.y[b0] -= wb * dy * corr;
    }
  }
  for (let i = 1; i < n - 1; i++) {
    c.vx[i] = 0.995 * (c.x[i] - c.px[i]) / h;
    c.vy[i] = 0.995 * (c.y[i] - c.py[i]) / h;
  }
}

// The test: settle for 15 s of simulated time, then compare with the catenary.
function measure(c) {
  let dev = 0, stretch = 0;
  for (let i = 0; i <= c.N; i++) dev = Math.max(dev, Math.abs(c.y[i] - catenaryY(c.x[i])));
  for (let k = 0; k < c.N; k++) stretch = Math.max(stretch, Math.hypot(c.x[k + 1] - c.x[k], c.y[k + 1] - c.y[k]) / c.L - 1);
  return { dev, rel: dev / sag, stretch };
}
function catenaryTest(N, nIter) {
  const c = makeChain(N);
  for (let s = 0; s < 15 * 60 * substeps; s++) substep(c, nIter);
  return measure(c);
}

function runTable() {
  const it = Math.round(iters.value);
  const rows = [10, 20, 40, 80].map((N) => Object.assign({ N }, catenaryTest(N, it)));
  rows.iterations = it;
  return rows;
}

let chain, table = null, simTime;
function reset() {
  chain = makeChain(+nSel.value);
  simTime = 0;
  if (!table) table = runTable();          // takes a moment: four chains, 15 simulated seconds each
}

function step() {
  for (let s = 0; s < substeps; s++) substep(chain, Math.round(iters.value));
  simTime += 1 / 60;
}

function draw(g) {
  const C = Sim.colors;
  // The catenary (dashed) and the chain.
  const curve = [];
  for (let k = 0; k <= 100; k++) { const x = -span / 2 + span * k / 100; curve.push({ x, y: catenaryY(x) }); }
  g.polyline(curve, { color: C.blue, width: 1.5, dash: [6, 4] });
  const pts = [];
  for (let i = 0; i <= chain.N; i++) pts.push({ x: chain.x[i], y: chain.y[i] });
  g.polyline(pts, { color: C.ink, width: 2 });
  g.points(chain.x, chain.y, chain.N + 1, { radius: 2.5, color: C.red });
  for (const sx of [-span / 2, span / 2]) g.rect(sx - 0.04, 0, 0.08, 0.05, { fill: '#d9d6cc' });

  const m = measure(chain);
  Sim.show('catenary a', a.toFixed(4) + ' m');
  Sim.show('sag', (sag * 100).toFixed(2) + ' cm');
  Sim.show('live: max deviation', (m.dev * 1000).toFixed(2) + ' mm (' + (100 * m.rel).toFixed(2) + '% of sag)');
  Sim.show('live: max link stretch', (100 * m.stretch).toFixed(3) + ' %');
  Sim.show('settled test (15 s), iterations', table.iterations);
  for (const r of table) {
    const pass = 100 * r.rel <= tolPct.value;
    Sim.show((pass ? 'PASS' : 'FAIL') + ' N = ' + r.N, (r.dev * 1000).toFixed(2) + ' mm, ' + (100 * r.rel).toFixed(2) +
      '% of sag, stretch ' + (100 * r.stretch).toFixed(3) + '%');
  }
}

Sim.run({ title: 'Catenary test', view: [-0.95, -0.7, 0.95, 0.8], reset, step, draw, dt: 1 / 60 });
