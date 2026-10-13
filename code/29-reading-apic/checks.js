// Chapter 29 - Reading APIC: a runnable check of the paper's claims.
// Two identical spinning elastic disks, simulated with explicit MPM. The left one uses PIC transfers,
// the right one APIC (equations (8)-(10) of the paper). No gravity, no walls: total angular momentum should
// stay constant. APIC keeps it to round-off; PIC loses it every step.
// The readouts also evaluate equation (9), the moment matrix D_p, by brute force at the mouse position for
// linear, quadratic and cubic weights.

const omega0 = Sim.slider('spin (rad/s)', 0.5, 4, 2, { step: 0.1, reset: true });
const substeps = Sim.slider('substeps per frame', 5, 20, 10, { step: 1 });

const dx = 1 / 24;                  // grid spacing (m)
const nx = 48, ny = 24;             // cells; the domain is [0, 2] x [0, 1] m
const nodes = (nx + 1) * (ny + 1);
const rho = 1000;                   // density (kg/m^3)
const E = 5e4, nu = 0.3;            // Young's modulus (Pa), Poisson's ratio
const mu = E / (2 * (1 + nu)), lam = E * nu / ((1 + nu) * (1 - 2 * nu));   // Lame parameters
const radius = 0.3;

let sims;   // [PIC, APIC]

function makeDisk(cx, cy, apic) {
  const xs = [], ys = [];
  const spacing = dx / 2;           // 2 x 2 particles per cell
  for (let x = cx - radius; x <= cx + radius; x += spacing) {
    for (let y = cy - radius; y <= cy + radius; y += spacing) {
      const px = x + spacing / 2, py = y + spacing / 2;
      if ((px - cx) ** 2 + (py - cy) ** 2 < radius * radius) { xs.push(px); ys.push(py); }
    }
  }
  const n = xs.length;
  const s = {
    apic, n, cx, cy,
    x: Float64Array.from(xs), y: Float64Array.from(ys),
    vx: new Float64Array(n), vy: new Float64Array(n),
    B: new Float64Array(4 * n),     // APIC affine state, row-major [Bxx, Bxy, Byx, Byy]
    F: new Float64Array(4 * n),     // deformation gradient
    m: (spacing * spacing) * rho, V0: spacing * spacing,
    gm: new Float64Array(nodes), gvx: new Float64Array(nodes), gvy: new Float64Array(nodes),
    colors: new Array(n), L0: 0, L: 0,
  };
  const w = omega0.value;
  for (let p = 0; p < n; p++) {
    const rx = s.x[p] - cx, ry = s.y[p] - cy;
    s.vx[p] = -w * ry; s.vy[p] = w * rx;              // rigid rotation v = omega x r
    if (s.apic) {                                     // the matching affine state: C = omega [[0,-1],[1,0]], B = C D
      const d = dx * dx / 4;
      s.B[4 * p + 1] = -w * d; s.B[4 * p + 2] = w * d;
    }
    s.F[4 * p] = 1; s.F[4 * p + 3] = 1;
    s.colors[p] = (rx > 0 && Math.abs(ry) < dx * 0.75) ? Sim.colors.red : (apic ? Sim.colors.blue : Sim.colors.orange);
  }
  return s;
}

function reset() {
  sims = [makeDisk(0.5, 0.5, false), makeDisk(1.5, 0.5, true)];
  for (const s of sims) { p2g(s); s.L0 = s.L; }
}

// 1D quadratic B-spline weights and derivatives (d/dx_p) for a particle at u = x / dx (in cells).
const qw = new Float64Array(3), qd = new Float64Array(3);
function quadratic(u) {
  const base = Math.floor(u - 0.5), f = u - base;
  qw[0] = 0.5 * (1.5 - f) ** 2; qw[1] = 0.75 - (f - 1) ** 2; qw[2] = 0.5 * (f - 0.5) ** 2;
  qd[0] = (f - 1.5) / dx;       qd[1] = -2 * (f - 1) / dx;    qd[2] = (f - 0.5) / dx;
  return base;
}

// Particle to grid. PIC: equation (1). APIC: equation (8) with D_p^-1 = 4 / dx^2.
// Also adds the explicit elastic force f_i = -sum_p V0 P F^T grad w_ip (fixed corotated), so the step is MPM.
function p2g(s, h = 0) {
  s.gm.fill(0); s.gvx.fill(0); s.gvy.fill(0);
  const fx = new Float64Array(nodes), fy = new Float64Array(nodes);
  const Dinv = 4 / (dx * dx);
  for (let p = 0; p < s.n; p++) {
    const bx = quadratic(s.x[p] / dx), wx = [qw[0], qw[1], qw[2]], dwx = [qd[0], qd[1], qd[2]];
    const by = quadratic(s.y[p] / dx), wy = [qw[0], qw[1], qw[2]], dwy = [qd[0], qd[1], qd[2]];
    const k = 4 * p;
    const Cxx = s.B[k] * Dinv, Cxy = s.B[k + 1] * Dinv, Cyx = s.B[k + 2] * Dinv, Cyy = s.B[k + 3] * Dinv;
    // Stress term V0 P F^T for fixed corotated: P = 2 mu (F - R) + lam (J - 1) J F^-T, so P F^T = 2 mu (F - R) F^T + lam (J - 1) J I
    const a = s.F[k], b = s.F[k + 1], c = s.F[k + 2], d = s.F[k + 3];
    const th = Math.atan2(c - b, a + d), co = Math.cos(th), si = Math.sin(th);   // polar rotation R of F
    const J = a * d - b * c;
    const Ma = a - co, Mb = b + si, Mc = c - si, Md = d - co;                    // F - R
    const Sxx = s.V0 * (2 * mu * (Ma * a + Mb * b) + lam * (J - 1) * J);
    const Sxy = s.V0 * (2 * mu * (Ma * c + Mb * d));
    const Syx = s.V0 * (2 * mu * (Mc * a + Md * b));
    const Syy = s.V0 * (2 * mu * (Mc * c + Md * d) + lam * (J - 1) * J);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const gi = bx + i, gj = by + j, id = gj * (nx + 1) + gi;
        const w = wx[i] * wy[j];
        const rx = gi * dx - s.x[p], ry = gj * dx - s.y[p];
        let ux = s.vx[p], uy = s.vy[p];
        if (s.apic) { ux += Cxx * rx + Cxy * ry; uy += Cyx * rx + Cyy * ry; }
        s.gm[id] += w * s.m;
        s.gvx[id] += w * s.m * ux;
        s.gvy[id] += w * s.m * uy;
        const gx = dwx[i] * wy[j], gy = wx[i] * dwy[j];                         // grad w_ip
        fx[id] -= Sxx * gx + Sxy * gy;
        fy[id] -= Syx * gx + Syy * gy;
      }
    }
  }
  // Angular momentum on the grid right after P2G: equation (3), which equals equation (11) for APIC.
  let L = 0;
  for (let id = 0; id < nodes; id++) {
    if (s.gm[id] === 0) continue;
    const xi = (id % (nx + 1)) * dx, yi = Math.floor(id / (nx + 1)) * dx;
    L += xi * s.gvy[id] - yi * s.gvx[id];
  }
  s.L = L;
  // Grid velocity and symplectic Euler update.
  for (let id = 0; id < nodes; id++) {
    if (s.gm[id] === 0) continue;
    s.gvx[id] = (s.gvx[id] + h * fx[id]) / s.gm[id];
    s.gvy[id] = (s.gvy[id] + h * fy[id]) / s.gm[id];
  }
}

// Grid to particle. PIC: equation (2). APIC: also equation (10). Then the F update and advection.
function g2p(s, h) {
  for (let p = 0; p < s.n; p++) {
    const bx = quadratic(s.x[p] / dx), wx = [qw[0], qw[1], qw[2]], dwx = [qd[0], qd[1], qd[2]];
    const by = quadratic(s.y[p] / dx), wy = [qw[0], qw[1], qw[2]], dwy = [qd[0], qd[1], qd[2]];
    let vx = 0, vy = 0, Bxx = 0, Bxy = 0, Byx = 0, Byy = 0, Gxx = 0, Gxy = 0, Gyx = 0, Gyy = 0;
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const gi = bx + i, gj = by + j, id = gj * (nx + 1) + gi;
        const w = wx[i] * wy[j], ux = s.gvx[id], uy = s.gvy[id];
        const rx = gi * dx - s.x[p], ry = gj * dx - s.y[p];
        vx += w * ux; vy += w * uy;
        Bxx += w * ux * rx; Bxy += w * ux * ry; Byx += w * uy * rx; Byy += w * uy * ry;
        const gx = dwx[i] * wy[j], gy = wx[i] * dwy[j];
        Gxx += ux * gx; Gxy += ux * gy; Gyx += uy * gx; Gyy += uy * gy;         // velocity gradient
      }
    }
    const k = 4 * p;
    s.vx[p] = vx; s.vy[p] = vy;
    if (s.apic) { s.B[k] = Bxx; s.B[k + 1] = Bxy; s.B[k + 2] = Byx; s.B[k + 3] = Byy; }
    // F <- (I + h grad v) F
    const a = s.F[k], b = s.F[k + 1], c = s.F[k + 2], d = s.F[k + 3];
    s.F[k] = (1 + h * Gxx) * a + h * Gxy * c;  s.F[k + 1] = (1 + h * Gxx) * b + h * Gxy * d;
    s.F[k + 2] = h * Gyx * a + (1 + h * Gyy) * c; s.F[k + 3] = h * Gyx * b + (1 + h * Gyy) * d;
    s.x[p] += h * vx; s.y[p] += h * vy;
  }
}

function step(dt) {
  const n = substeps.value, h = dt / n;
  for (let k = 0; k < n; k++) for (const s of sims) { p2g(s, h); g2p(s, h); }
}

// Brute-force D_p / dx^2 from equation (9) for a given 1D weight function.
function momentMatrix(kind, px, py) {
  const stencil = (u) => {
    if (kind === 'linear') { const b = Math.floor(u), s = u - b; return { b, w: [1 - s, s] }; }
    if (kind === 'quadratic') { const b = quadratic(u); return { b, w: [qw[0], qw[1], qw[2]] }; }
    const b = Math.floor(u) - 1, s = u - Math.floor(u);
    return { b, w: [(1 - s) ** 3 / 6, 0.5 * s ** 3 - s * s + 2 / 3, 0.5 * (1 - s) ** 3 - (1 - s) ** 2 + 2 / 3, s ** 3 / 6] };
  };
  const ux = px / dx, uy = py / dx, sx = stencil(ux), sy = stencil(uy);
  let Dxx = 0, Dxy = 0, Dyy = 0;
  for (let i = 0; i < sx.w.length; i++) {
    for (let j = 0; j < sy.w.length; j++) {
      const w = sx.w[i] * sy.w[j], rx = sx.b + i - ux, ry = sy.b + j - uy;
      Dxx += w * rx * rx; Dxy += w * rx * ry; Dyy += w * ry * ry;
    }
  }
  return [Dxx, Dxy, Dyy];
}

function draw(g) {
  g.rect(0, 0, 2, 1, { stroke: Sim.colors.muted, width: 1 });
  for (const s of sims) {
    g.points(s.x, s.y, s.n, { radius: 1.6, colors: s.colors });
    g.text(s.apic ? 'APIC' : 'PIC', s.cx, 0.93, { align: 'center', size: 14, color: Sim.colors.ink });
    Sim.plot(s.apic ? 'APIC L/L0' : 'PIC L/L0', s.L / s.L0, { color: s.apic ? Sim.colors.blue : Sim.colors.orange, min: 0, max: 1.05 });
  }
  Sim.show('APIC: (L - L0) / L0', (sims[1].L - sims[1].L0) / sims[1].L0);
  Sim.show('PIC:  (L - L0) / L0', (sims[0].L - sims[0].L0) / sims[0].L0);

  // Probe: D_p at the mouse, in units of dx^2.
  const px = Math.min(Math.max(Sim.mouse.x, 3 * dx), 2 - 3 * dx), py = Math.min(Math.max(Sim.mouse.y, 3 * dx), 1 - 3 * dx);
  for (const kind of ['linear', 'quadratic', 'cubic']) {
    const [a, b, c] = momentMatrix(kind, px, py);
    Sim.show(`D/dx^2 ${kind}`, `[${a.toFixed(3)}, ${b.toFixed(3)}; ${b.toFixed(3)}, ${c.toFixed(3)}]`);
  }
  // Show the quadratic stencil around the probe, dot size = weight.
  const bx = quadratic(px / dx), wx = [qw[0], qw[1], qw[2]];
  const by = quadratic(py / dx), wy = [qw[0], qw[1], qw[2]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) g.dot((bx + i) * dx, (by + j) * dx, 2 + 10 * wx[i] * wy[j], Sim.colors.green);
  g.dot(px, py, 3, Sim.colors.ink);
}

Sim.run({ title: 'PIC vs APIC: angular momentum and D_p', view: [-0.05, -0.12, 2.05, 1.05], reset, step, draw, grid: dx });
