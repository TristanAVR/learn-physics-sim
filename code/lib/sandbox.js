/*
 * sandbox.js - a tiny 2D playground for writing physics simulations.
 *
 * Load it with a classic <script> tag (no build step, works from file:// and GitHub Pages).
 * Your script defines reset/step/draw functions and calls Sim.run({...}).
 *
 *   const g0 = Sim.slider('gravity', 0, 20, 9.81);
 *   let p;
 *   function reset() { p = { x: 0, y: 1, vx: 1, vy: 0 }; }
 *   function step(dt) { p.vy -= g0.value * dt; p.x += p.vx * dt; p.y += p.vy * dt; }
 *   function draw(g) { g.circle(p.x, p.y, 0.05, { fill: Sim.colors.blue }); }
 *   Sim.run({ title: 'Falling', view: [-1, -0.2, 3, 1.5], reset, step, draw });
 *
 * World coordinates: metres, y points UP. `view` is [xmin, ymin, xmax, ymax]; it is fitted into the
 * canvas keeping the aspect ratio. Keys: Space play/pause, Period step one frame, R reset.
 */
(function () {
  'use strict';

  const Sim = {};
  window.Sim = Sim;

  // ---------------------------------------------------------------------------------------------
  // Small 2D vector helpers on {x, y} objects. They return new objects; fine for small scenes.
  // For thousands of particles use flat Float64Arrays instead (see the MPM chapters).
  // ---------------------------------------------------------------------------------------------
  const V = {
    v: (x = 0, y = 0) => ({ x, y }),
    copy: (a) => ({ x: a.x, y: a.y }),
    add: (a, b) => ({ x: a.x + b.x, y: a.y + b.y }),
    sub: (a, b) => ({ x: a.x - b.x, y: a.y - b.y }),
    scale: (a, s) => ({ x: a.x * s, y: a.y * s }),
    addScaled: (a, b, s) => ({ x: a.x + b.x * s, y: a.y + b.y * s }),
    dot: (a, b) => a.x * b.x + a.y * b.y,
    cross: (a, b) => a.x * b.y - a.y * b.x,
    len: (a) => Math.hypot(a.x, a.y),
    lenSq: (a) => a.x * a.x + a.y * a.y,
    dist: (a, b) => Math.hypot(a.x - b.x, a.y - b.y),
    norm: (a) => { const l = Math.hypot(a.x, a.y); return l > 0 ? { x: a.x / l, y: a.y / l } : { x: 0, y: 0 }; },
    perp: (a) => ({ x: -a.y, y: a.x }),
    lerp: (a, b, t) => ({ x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t }),
    rotate: (a, angle) => { const c = Math.cos(angle), s = Math.sin(angle); return { x: c * a.x - s * a.y, y: s * a.x + c * a.y }; },
  };
  Sim.V = V;
  window.V = V;

  Sim.colors = {
    bg: '#fbfaf7', grid: '#ecebe6', axis: '#d6d4cc', ink: '#2b2b2b', muted: '#8a877d',
    blue: '#2f6fdf', red: '#d9473a', green: '#2e9e5b', orange: '#e58a1f', purple: '#7d4fc4',
    teal: '#1e9aa6', yellow: '#d8b22a', gray: '#9a978e',
  };
  const palette = ['#2f6fdf', '#d9473a', '#2e9e5b', '#e58a1f', '#7d4fc4', '#1e9aa6', '#d8b22a'];
  Sim.palette = palette;

  /** Seeded random numbers, so a scene looks the same every reset. */
  Sim.random = function (seed = 1) {
    let s = seed >>> 0 || 1;
    return function () {
      s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
      return s / 4294967296;
    };
  };

  // ---------------------------------------------------------------------------------------------
  // DOM: header bar, control panel, canvas, overlays
  // ---------------------------------------------------------------------------------------------
  const params = new URLSearchParams(location.search);
  const embedded = params.has('embed');
  let built = false;
  let root, header, titleEl, panel, canvas, ctx, infoEl, errorEl, playBtn;
  const controls = [];

  function build() {
    if (built) return;
    built = true;
    document.body.classList.add('sim-body');
    if (embedded) document.body.classList.add('sim-embedded');
    root = el('div', 'sim-root', document.body);
    header = el('div', 'sim-header', root);
    titleEl = el('div', 'sim-title', header);
    const buttons = el('div', 'sim-buttons', header);
    playBtn = button(buttons, 'Pause', () => setPaused(!state.paused), 'Space');
    button(buttons, 'Step', () => { setPaused(true); state.stepOnce = true; }, 'Period (.)');
    button(buttons, 'Reset', () => doReset(), 'R');
    const main = el('div', 'sim-main', root);
    const stage = el('div', 'sim-stage', main);
    canvas = el('canvas', 'sim-canvas', stage);
    infoEl = el('div', 'sim-info', stage);
    errorEl = el('pre', 'sim-error', stage);
    panel = el('div', 'sim-panel', main);
    ctx = canvas.getContext('2d');
    window.addEventListener('resize', () => resize());
    window.addEventListener('keydown', onKey);
    canvas.addEventListener('pointerdown', onPointer);
    canvas.addEventListener('pointermove', onPointer);
    window.addEventListener('pointerup', onPointer);
    canvas.addEventListener('contextmenu', (e) => e.preventDefault());
  }

  function el(tag, cls, parent) {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (parent) parent.appendChild(e);
    return e;
  }

  function button(parent, label, fn, key) {
    const b = el('button', 'sim-button', parent);
    b.textContent = label;
    if (key) b.title = key;
    b.addEventListener('click', (e) => { fn(); e.currentTarget.blur(); });
    return b;
  }

  // ---------------------------------------------------------------------------------------------
  // Controls. Call these at the top level of your script, before Sim.run.
  // Each returns an object whose .value is always current. Pass { reset: true } to reset the
  // scene when the control changes.
  // ---------------------------------------------------------------------------------------------
  Sim.slider = function (label, min, max, value, opts = {}) {
    build();
    const step = opts.step !== undefined ? opts.step : (max - min) / 1000;
    const row = el('label', 'sim-control', panel);
    const name = el('span', 'sim-control-name', row);
    name.textContent = label;
    const out = el('span', 'sim-control-value', row);
    const input = el('input', '', row);
    input.type = 'range';
    input.min = min; input.max = max; input.step = step; input.value = value;
    const c = { value: +value };
    const fmt = opts.format || ((v) => formatNumber(v));
    out.textContent = fmt(c.value);
    input.addEventListener('input', () => {
      c.value = +input.value;
      out.textContent = fmt(c.value);
      if (opts.reset) doReset();
      if (opts.onChange) opts.onChange(c.value);
    });
    c.set = (v) => { c.value = +v; input.value = v; out.textContent = fmt(c.value); };
    controls.push(c);
    return c;
  };

  /** A slider over powers of ten, for stiffness or compliance spanning many decades. */
  Sim.logSlider = function (label, minExp, maxExp, value, opts = {}) {
    build();
    const row = el('label', 'sim-control', panel);
    const name = el('span', 'sim-control-name', row);
    name.textContent = label;
    const out = el('span', 'sim-control-value', row);
    const input = el('input', '', row);
    input.type = 'range';
    input.min = minExp; input.max = maxExp; input.step = opts.step || 0.01;
    input.value = Math.log10(value);
    const c = { value: +value };
    out.textContent = formatNumber(c.value);
    input.addEventListener('input', () => {
      c.value = Math.pow(10, +input.value);
      out.textContent = formatNumber(c.value);
      if (opts.reset) doReset();
      if (opts.onChange) opts.onChange(c.value);
    });
    c.set = (v) => { c.value = +v; input.value = Math.log10(v); out.textContent = formatNumber(c.value); };
    controls.push(c);
    return c;
  };

  Sim.checkbox = function (label, value, opts = {}) {
    build();
    const row = el('label', 'sim-control sim-check', panel);
    const input = el('input', '', row);
    input.type = 'checkbox';
    input.checked = !!value;
    const name = el('span', 'sim-control-name', row);
    name.textContent = label;
    const c = { value: !!value };
    input.addEventListener('change', () => {
      c.value = input.checked;
      if (opts.reset) doReset();
      if (opts.onChange) opts.onChange(c.value);
    });
    c.set = (v) => { c.value = !!v; input.checked = !!v; };
    controls.push(c);
    return c;
  };

  Sim.select = function (label, options, value, opts = {}) {
    build();
    const row = el('label', 'sim-control', panel);
    const name = el('span', 'sim-control-name', row);
    name.textContent = label;
    const input = el('select', '', row);
    for (const o of options) {
      const opt = el('option', '', input);
      opt.value = o; opt.textContent = o;
    }
    input.value = value;
    const c = { value };
    input.addEventListener('change', () => {
      c.value = input.value;
      if (opts.reset) doReset();
      if (opts.onChange) opts.onChange(c.value);
    });
    c.set = (v) => { c.value = v; input.value = v; };
    controls.push(c);
    return c;
  };

  Sim.button = function (label, fn) {
    build();
    const b = button(panel, label, fn);
    b.classList.add('sim-panel-button');
    return b;
  };

  // ---------------------------------------------------------------------------------------------
  // Readouts. Call both from draw(). Sim.show prints a value in the top-left info box this frame;
  // Sim.plot records one sample per simulated frame (nothing while paused) and draws a scrolling
  // graph in the bottom-right corner. Plots are cleared on reset.
  // ---------------------------------------------------------------------------------------------
  let shown = [];
  Sim.show = function (label, value, digits) {
    const text = typeof value === 'number' ? formatNumber(value, digits) : String(value);
    shown.push([label, text]);
  };

  const plots = new Map();
  let recording = true;
  Sim.plot = function (name, value, opts = {}) {
    if (!recording) return;
    let p = plots.get(name);
    if (!p) {
      p = { name, values: [], color: opts.color || palette[plots.size % palette.length], opts };
      plots.set(name, p);
    }
    p.opts = Object.assign(p.opts, opts);
    p.values.push(value);
    if (p.values.length > 600) p.values.shift();
  };

  function formatNumber(v, digits) {
    if (!isFinite(v)) return String(v);
    if (digits !== undefined) return v.toFixed(digits);
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e5)) return v.toExponential(2);
    if (Number.isInteger(v)) return String(v);
    return v.toFixed(a < 1 ? 4 : a < 100 ? 3 : 1);
  }
  Sim.format = formatNumber;

  // ---------------------------------------------------------------------------------------------
  // Mouse, in world coordinates. Sim.mouse.pressed / released are true for one frame.
  // ---------------------------------------------------------------------------------------------
  const mouse = { x: 0, y: 0, down: false, pressed: false, released: false, button: 0, inside: false };
  Sim.mouse = mouse;
  let pendingPress = false, pendingRelease = false;

  function onPointer(e) {
    const r = canvas.getBoundingClientRect();
    const px = (e.clientX - r.left) * (canvas.width / r.width);
    const py = (e.clientY - r.top) * (canvas.height / r.height);
    const w = toWorld(px, py);
    mouse.x = w.x; mouse.y = w.y;
    mouse.inside = px >= 0 && py >= 0 && px <= canvas.width && py <= canvas.height;
    if (e.type === 'pointerdown') {
      mouse.down = true; mouse.button = e.button; pendingPress = true;
      canvas.setPointerCapture && canvas.setPointerCapture(e.pointerId);
    } else if (e.type === 'pointerup' && mouse.down) {
      mouse.down = false; pendingRelease = true;
    }
  }

  function onKey(e) {
    if (e.target && (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT')) return;
    if (e.code === 'Space') { setPaused(!state.paused); e.preventDefault(); }
    else if (e.key === '.') { setPaused(true); state.stepOnce = true; }
    else if (e.key === 'r' || e.key === 'R') doReset();
    if (state.opts && state.opts.onKey) state.opts.onKey(e);
  }

  // ---------------------------------------------------------------------------------------------
  // View transform and the drawing API (world coordinates, y up)
  // ---------------------------------------------------------------------------------------------
  const view = { xmin: -1, ymin: -1, xmax: 1, ymax: 1, scale: 100, ox: 0, oy: 0 };
  let dpr = 1;

  function fitView() {
    const w = canvas.width, h = canvas.height;
    const vw = view.xmax - view.xmin, vh = view.ymax - view.ymin;
    view.scale = Math.min(w / vw, h / vh);
    view.ox = (w - vw * view.scale) / 2 - view.xmin * view.scale;
    view.oy = (h + vh * view.scale) / 2 + view.ymin * view.scale;
  }
  const toScreenX = (x) => view.ox + x * view.scale;
  const toScreenY = (y) => view.oy - y * view.scale;
  function toWorld(px, py) { return { x: (px - view.ox) / view.scale, y: (view.oy - py) / view.scale }; }

  /** Change the visible region at any time: [xmin, ymin, xmax, ymax]. */
  Sim.setView = function (v) {
    view.xmin = v[0]; view.ymin = v[1]; view.xmax = v[2]; view.ymax = v[3];
    if (canvas) fitView();
  };

  function style(s, fallbackFill, fallbackStroke) {
    s = s || {};
    return {
      fill: s.fill !== undefined ? s.fill : fallbackFill,
      stroke: s.stroke !== undefined ? s.stroke : (s.color !== undefined ? s.color : fallbackStroke),
      width: (s.width || 1.5) * dpr,
      alpha: s.alpha !== undefined ? s.alpha : 1,
      dash: s.dash,
    };
  }

  function apply(st) {
    ctx.globalAlpha = st.alpha;
    if (st.fill) ctx.fillStyle = st.fill;
    if (st.stroke) { ctx.strokeStyle = st.stroke; ctx.lineWidth = st.width; }
    ctx.setLineDash(st.dash ? st.dash.map((d) => d * dpr) : []);
  }

  const g = {
    /** Pixels per world unit at the current zoom. */
    get pxPerUnit() { return view.scale / dpr; },
    /** Circle of world radius r. Style: { fill, stroke, width, alpha, dash }. */
    circle(x, y, r, s) {
      const st = style(s, s && (s.stroke || s.color) ? null : Sim.colors.blue, null);
      apply(st);
      ctx.beginPath();
      ctx.arc(toScreenX(x), toScreenY(y), Math.max(r * view.scale, 0.5), 0, Math.PI * 2);
      if (st.fill) ctx.fill();
      if (st.stroke) ctx.stroke();
      ctx.globalAlpha = 1;
    },
    /** Dot with a fixed radius in pixels, whatever the zoom. */
    dot(x, y, px = 4, color = Sim.colors.ink) {
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(toScreenX(x), toScreenY(y), px * dpr, 0, Math.PI * 2);
      ctx.fill();
    },
    line(x1, y1, x2, y2, s) {
      const st = style(s, null, Sim.colors.ink);
      apply(st);
      ctx.beginPath();
      ctx.moveTo(toScreenX(x1), toScreenY(y1));
      ctx.lineTo(toScreenX(x2), toScreenY(y2));
      ctx.stroke();
      ctx.globalAlpha = 1;
    },
    /** points: array of {x, y}, or a flat array [x0, y0, x1, y1, ...]. Style adds { closed }. */
    polyline(points, s) {
      const st = style(s, null, Sim.colors.ink);
      if (s && s.fill && !s.stroke && !s.color) st.stroke = null;
      apply(st);
      ctx.beginPath();
      const flat = typeof points[0] === 'number';
      const n = flat ? points.length / 2 : points.length;
      for (let i = 0; i < n; i++) {
        const x = flat ? points[2 * i] : points[i].x, y = flat ? points[2 * i + 1] : points[i].y;
        if (i === 0) ctx.moveTo(toScreenX(x), toScreenY(y)); else ctx.lineTo(toScreenX(x), toScreenY(y));
      }
      if (s && s.closed) ctx.closePath();
      if (st.fill) ctx.fill();
      if (st.stroke) ctx.stroke();
      ctx.globalAlpha = 1;
    },
    rect(x, y, w, h, s) {
      const st = style(s, s && (s.stroke || s.color) ? null : Sim.colors.gray, null);
      apply(st);
      const sx = toScreenX(x), sy = toScreenY(y + h);
      if (st.fill) ctx.fillRect(sx, sy, w * view.scale, h * view.scale);
      if (st.stroke) ctx.strokeRect(sx, sy, w * view.scale, h * view.scale);
      ctx.globalAlpha = 1;
    },
    /** Arrow from (x, y) along (dx, dy) in world units. */
    arrow(x, y, dx, dy, s) {
      const st = style(s, null, Sim.colors.red);
      apply(st);
      const x1 = toScreenX(x), y1 = toScreenY(y), x2 = toScreenX(x + dx), y2 = toScreenY(y + dy);
      const len = Math.hypot(x2 - x1, y2 - y1);
      ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
      if (len > 2) {
        const head = Math.min(10 * dpr, len * 0.4), ux = (x2 - x1) / len, uy = (y2 - y1) / len;
        ctx.fillStyle = st.stroke;
        ctx.beginPath();
        ctx.moveTo(x2, y2);
        ctx.lineTo(x2 - head * ux + head * 0.45 * uy, y2 - head * uy - head * 0.45 * ux);
        ctx.lineTo(x2 - head * ux - head * 0.45 * uy, y2 - head * uy + head * 0.45 * ux);
        ctx.closePath(); ctx.fill();
      }
      ctx.globalAlpha = 1;
    },
    /** Text at a world position. Style: { color, size (px), align: 'left'|'center'|'right', baseline }. */
    text(str, x, y, s = {}) {
      ctx.fillStyle = s.color || Sim.colors.ink;
      ctx.font = `${(s.size || 13) * dpr}px system-ui, sans-serif`;
      ctx.textAlign = s.align || 'left';
      ctx.textBaseline = s.baseline || 'middle';
      ctx.fillText(str, toScreenX(x), toScreenY(y));
    },
    /** Many particles at once from flat arrays: xs[i], ys[i] for i < count. Fast. */
    points(xs, ys, count, s = {}) {
      const r = (s.radius || 2) * dpr;
      const colors = s.colors;
      ctx.fillStyle = s.color || Sim.colors.blue;
      let last = null;
      for (let i = 0; i < count; i++) {
        if (colors) {
          const c = colors[i];
          if (c !== last) { ctx.fillStyle = c; last = c; }
        }
        ctx.fillRect(toScreenX(xs[i]) - r, toScreenY(ys[i]) - r, 2 * r, 2 * r);
      }
    },
    /** Background grid lines every `spacing` world units. */
    grid(spacing = 0.1, color = Sim.colors.grid) {
      const tl = toWorld(0, 0), br = toWorld(canvas.width, canvas.height);
      ctx.strokeStyle = color; ctx.lineWidth = 1; ctx.setLineDash([]);
      ctx.beginPath();
      if ((br.x - tl.x) / spacing > 400) return;
      for (let x = Math.ceil(tl.x / spacing) * spacing; x <= br.x; x += spacing) {
        const sx = Math.round(toScreenX(x)) + 0.5; ctx.moveTo(sx, 0); ctx.lineTo(sx, canvas.height);
      }
      for (let y = Math.ceil(br.y / spacing) * spacing; y <= tl.y; y += spacing) {
        const sy = Math.round(toScreenY(y)) + 0.5; ctx.moveTo(0, sy); ctx.lineTo(canvas.width, sy);
      }
      ctx.stroke();
    },
    /** Ground or wall as a hatched half-plane below y (default) for a quick floor. */
    ground(y = 0, s = {}) {
      const tl = toWorld(0, 0), br = toWorld(canvas.width, canvas.height);
      this.rect(tl.x, Math.min(br.y, y) - 10, br.x - tl.x, y - Math.min(br.y, y) + 10, { fill: s.fill || '#e3e0d6' });
      this.line(tl.x, y, br.x, y, { color: s.color || Sim.colors.muted, width: 1.5 });
    },
    /** Raw canvas access for anything else. Coordinates are in device pixels. */
    get ctx() { return ctx; },
    toScreen(x, y) { return { x: toScreenX(x), y: toScreenY(y) }; },
    /** Device pixels to world coordinates. */
    toWorld(px, py) { return toWorld(px, py); },
  };
  Sim.graphics = g;

  // ---------------------------------------------------------------------------------------------
  // The loop
  // ---------------------------------------------------------------------------------------------
  const state = { opts: null, paused: false, stepOnce: false, time: 0, frame: 0, lastError: null };
  Object.defineProperty(Sim, 'time', { get: () => state.time });
  Object.defineProperty(Sim, 'frame', { get: () => state.frame });
  Object.defineProperty(Sim, 'paused', { get: () => state.paused });

  function setPaused(p) {
    state.paused = p;
    if (playBtn) playBtn.textContent = p ? 'Play' : 'Pause';
  }

  function doReset() {
    state.time = 0;
    state.frame = 0;
    plots.clear();
    clearError();
    if (state.opts && state.opts.reset) guard(() => state.opts.reset());
    recording = true;
    renderOnce();
  }
  Sim.reset = doReset;

  function guard(fn) {
    try { fn(); return true; }
    catch (err) { showError(err); return false; }
  }

  function showError(err) {
    state.lastError = err;
    setPaused(true);
    errorEl.textContent = (err && err.stack) ? err.stack : String(err);
    errorEl.style.display = 'block';
    console.error(err);
  }
  function clearError() { state.lastError = null; if (errorEl) errorEl.style.display = 'none'; }

  window.addEventListener('error', (e) => { if (built) showError(e.error || e.message); });

  function resize(render = true) {
    recording = false;
    dpr = window.devicePixelRatio || 1;
    const r = canvas.getBoundingClientRect();
    canvas.width = Math.max(1, Math.round(r.width * dpr));
    canvas.height = Math.max(1, Math.round(r.height * dpr));
    fitView();
    if (render !== false) renderOnce();
  }

  function renderOnce() {
    if (!canvas || !state.opts) return;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.globalAlpha = 1;
    ctx.fillStyle = state.opts.background || Sim.colors.bg;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    if (state.opts.grid !== false) g.grid(state.opts.grid || gridSpacing());
    shown = [];
    if (state.opts.draw) guard(() => state.opts.draw(g));
    drawPlots();
    drawInfo();
  }

  function gridSpacing() {
    const vw = view.xmax - view.xmin;
    const raw = vw / 20, p = Math.pow(10, Math.floor(Math.log10(raw)));
    return raw / p < 2 ? p : raw / p < 5 ? 2 * p : 5 * p;
  }

  function drawInfo() {
    const lines = [];
    if (state.opts.showTime !== false) lines.push(['t', state.time.toFixed(2) + ' s']);
    for (const s of shown) lines.push(s);
    infoEl.innerHTML = lines.map(([k, v]) => `<div><span>${escapeHtml(k)}</span><b>${escapeHtml(v)}</b></div>`).join('');
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  }

  function drawPlots() {
    if (plots.size === 0) return;
    const w = Math.min(canvas.width * 0.45, 420 * dpr), h = Math.min(canvas.height * 0.28, 140 * dpr);
    const corner = (state.opts && state.opts.plotCorner) || 'bottom-right';
    const x0 = corner.endsWith('left') ? 10 * dpr : canvas.width - w - 10 * dpr;
    const y0 = corner.startsWith('top') ? 10 * dpr : canvas.height - h - 10 * dpr;
    ctx.globalAlpha = 0.92;
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(x0, y0, w, h);
    ctx.globalAlpha = 1;
    ctx.strokeStyle = Sim.colors.axis; ctx.lineWidth = 1; ctx.setLineDash([]);
    ctx.strokeRect(x0 + 0.5, y0 + 0.5, w, h);
    let lo = Infinity, hi = -Infinity;
    for (const p of plots.values()) {
      if (p.opts.min !== undefined) lo = Math.min(lo, p.opts.min);
      if (p.opts.max !== undefined) hi = Math.max(hi, p.opts.max);
      for (const v of p.values) if (isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    }
    if (!isFinite(lo) || !isFinite(hi)) return;
    if (hi - lo < 1e-12) { hi += 1e-6; lo -= 1e-6; }
    const pad = 16 * dpr;
    const sx = (w - 8 * dpr) / 600, sy = (h - pad - 6 * dpr) / (hi - lo);
    ctx.font = `${10 * dpr}px system-ui, sans-serif`;
    ctx.textBaseline = 'top'; ctx.textAlign = 'left';
    let lx = x0 + 6 * dpr;
    for (const p of plots.values()) {
      ctx.strokeStyle = p.color; ctx.lineWidth = 1.5 * dpr;
      ctx.beginPath();
      let started = false;
      for (let i = 0; i < p.values.length; i++) {
        const v = p.values[i];
        if (!isFinite(v)) { started = false; continue; }
        const px = x0 + 4 * dpr + i * sx, py = y0 + h - 4 * dpr - (v - lo) * sy;
        if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.fillStyle = p.color;
      const last = p.values[p.values.length - 1];
      const label = `${p.name} ${formatNumber(last)}`;
      ctx.fillText(label, lx, y0 + 3 * dpr);
      lx += ctx.measureText(label).width + 12 * dpr;
    }
    ctx.fillStyle = Sim.colors.muted;
    ctx.textAlign = 'right';
    ctx.fillText(formatNumber(hi), x0 + w - 4 * dpr, y0 + 3 * dpr);
    ctx.textBaseline = 'bottom';
    ctx.fillText(formatNumber(lo), x0 + w - 4 * dpr, y0 + h - 2 * dpr);
  }

  /**
   * Advance n frames synchronously, as if n animation frames had run (ignores pause), then draw once.
   * For tests and headless checks; returns false if a step threw.
   */
  Sim.advance = function (n = 1) {
    const o = state.opts;
    if (!o) return false;
    const steps = typeof o.stepsPerFrame === 'function' ? o.stepsPerFrame() : (o.stepsPerFrame || 1);
    for (let f = 0; f < n && !state.lastError; f++) {
      recording = true;
      const dt = typeof o.dt === 'function' ? o.dt() : o.dt;
      for (let i = 0; i < steps; i++) {
        if (!guard(() => o.step(dt))) break;
        state.time += dt;
      }
      state.frame++;
      if (f < n - 1) { shown = []; if (o.draw) guard(() => o.draw(g)); }
    }
    renderOnce();
    return !state.lastError;
  };
  Object.defineProperty(Sim, 'error', { get: () => state.lastError });

  function frame() {
    const o = state.opts;
    mouse.pressed = pendingPress; mouse.released = pendingRelease;
    pendingPress = false; pendingRelease = false;
    recording = false;
    if ((!state.paused || state.stepOnce) && !state.lastError) {
      recording = true;
      const steps = typeof o.stepsPerFrame === 'function' ? o.stepsPerFrame() : (o.stepsPerFrame || 1);
      const dt = typeof o.dt === 'function' ? o.dt() : o.dt;
      for (let i = 0; i < steps; i++) {
        if (!guard(() => o.step(dt))) break;
        state.time += dt;
      }
      state.frame++;
      state.stepOnce = false;
    }
    renderOnce();
    requestAnimationFrame(frame);
  }

  /**
   * Start the simulation.
   *   title          shown in the header
   *   view           [xmin, ymin, xmax, ymax] in world units (metres)
   *   reset()        build the initial state (also called on R and by { reset: true } controls)
   *   step(dt)       advance the state by dt seconds
   *   draw(g)        draw the state with the graphics API
   *   dt             time step in seconds (default 1/60), or a function returning it (e.g. from a slider)
   *   stepsPerFrame  how many step(dt) calls per animation frame (default 1), or a function
   *   paused         start paused (default false)
   *   grid           false for no grid, or a spacing in world units
   *   plotCorner     where Sim.plot draws: 'bottom-right' (default), 'bottom-left', 'top-right', 'top-left'
   */
  Sim.run = function (opts) {
    build();
    state.opts = Object.assign({ dt: 1 / 60, stepsPerFrame: 1 }, opts);
    titleEl.textContent = opts.title || document.title || 'Simulation';
    if (opts.title) document.title = opts.title;
    if (opts.view) Sim.setView(opts.view);
    if (!panel.children.length) panel.classList.add('sim-panel-empty');
    setPaused(!!opts.paused);
    resize(false);
    doReset();
    requestAnimationFrame(frame);
  };
})();
