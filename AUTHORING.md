# Authoring guide

How every chapter of Learn Physics Simulation is written. Read it fully before writing or editing a chapter, then read
the model chapter, `chapters/04-time-integration.html` with `code/04-time-integration/`, and match it.

## 1. What the site is

A LearnOpenGL-style course on physics simulation, from a falling particle to the Material Point Method. Each chapter
introduces one idea (the problem it solves, the physics, the maths, small code snippets), shows it running, and then
asks the reader to implement it themselves in a starter file, with hints and a reference solution. By the end, the
reader can read three papers unaided:

- *The Affine Particle-In-Cell Method*, Jiang, Schroeder, Selle, Teran, Stomakhin, SIGGRAPH 2015 (ACM TOG 34(4) 51).
- *The Material Point Method for Simulating Continuum Materials*, Jiang, Schroeder, Teran, Stomakhin, Selle, SIGGRAPH
  2016 course notes.
- *A Position Based Material Point Method*, Chris Lewin (EA SEED), SIGGRAPH 2024 talk; open-source WebGPU code at
  https://github.com/electronicarts/pbmpm (BSD-3).

## 2. Confidentiality and scope (hard rules)

- **Keep every example general.** Do not refer to any employer, client, product or industry application. The methods
  are general; use everyday examples: ropes, hoses, cables, chains, tentacles, cable-driven continuum robots, fishing
  rods, cloth, balloons, jelly, rubber, snow, sand, water, honey, earthworms and peristaltic pumps.
- Cite papers and public resources freely. Do not copy text from the papers: summarise in your own words, use their
  equation numbers and section names, and quote at most a short phrase. Do not reproduce figures.
- Write your own code. The PB-MPM repository may be read to get the algorithm right, but do not paste it.

## 3. Voice

- Second person, friendly, precise, unhurried. Explain *why* before *how*. Assume a capable programmer who has
  forgotten most of their maths; introduce every symbol the first time it appears.
- Concrete before abstract: start from something the reader can picture or play with, then the diagram or the
  equation, then the code.
- Short paragraphs. No marketing language, no filler, no "simply" or "just" before something hard.
- British or American spelling is fine; be consistent within a chapter (the model chapter uses British).
- Numbers carry units. Worked examples use real numbers the reader can check.

## 4. Chapter structure

Each chapter is one HTML file in `chapters/`, named by its id from `assets/chapters.js`. Copy the `<head>`, the body
attributes and the page skeleton from the model chapter exactly; change only the title and `data-chapter`.

1. `<h1>` title and a one-line `<p class="subtitle">`.
2. Opening: the problem this chapter solves, in two or three short paragraphs, ending with what the reader will be
   able to do.
3. The idea, built up in `<h2>` sections: physics, maths, small code snippets (`<pre><code class="language-js">`).
   Derive the important results; state the rest with a pointer.
4. At least one embedded demo of the finished technique, placed where it helps understanding (usually early, so the
   reader plays before the maths):
   `<div class="demo" data-src="<id>/solution.js" data-height="460" data-caption="..."></div>`.
5. `<h2>Your turn</h2>` with two to four exercises (`<div class="exercise"><h4>...</h4>...</div>`; the CSS prefixes
   "Your turn:"), each followed by `<details class="hint">` hints and a final `<details>` solution that shows the key
   code and links to `../code/<id>/solution.js`. The first exercise is always "make the starter work". Later ones ask
   the reader to predict a number before measuring it, extend the code, or break it on purpose.
6. Optional boxes: `<div class="note">` (key insight), `<div class="warning">` (common mistake), and
   `<div class="paper"><h4>...</h4>...</div>` (Paper corner: where one of the three papers says the same thing, with
   section and equation numbers and a translation of their notation into ours). MPM chapters should have several.
7. `<h2>Summary</h2>` as a bullet list, then `<h3>Questions to come back to</h3>`: three to five recall questions,
   including one that needs a calculation.

Length: 2,500-4,500 words of prose per chapter (the reading chapters may be longer). Maths uses KaTeX: `$...$` inline,
`$$...$$` display. Inside HTML, escape `<` and `&` in text (`&lt;`, `&amp;`), and avoid a literal `$` other than as a
maths delimiter.

## 4b. Presenting maths

Every formula must be followable by a reader who meets it for the first time. These rules come from the learning
research on worked examples, cognitive load and multiple representations; `chapters/04-time-integration.html` is the
reference implementation.

1. **Notation before formulas (pre-training).** Introduce the symbols on their own, before the first equation that
   combines them: a short "the notation" section with a picture where possible (a timeline, a diagram) and a table of
   symbol, how to read it, units. Never introduce a symbol and use it in a complicated expression in the same sentence.
2. **One idea at a time, then combine (isolated, then interacting elements).** Show each piece separately before the
   formula that uses them together.
3. **Numbers before symbols.** Work a concrete example with easy numbers (k = 1, m = 1, h = 0.5) that the reader can
   check on a calculator before, or right next to, the general algebra. Use a table when several methods are compared.
4. **Compare methods side by side, not one after another.** When a chapter has alternatives (explicit/symplectic/
   implicit, PIC/FLIP/APIC, Jacobi/Gauss-Seidel), put them in one table, aligned so the single difference is obvious.
5. **Label terms where they are, not in a paragraph below (no split attention).** Use `\underbrace{...}_{\text{...}}`
   for the meaning of a term, and for multi-step derivations an `aligned` block with a short justification per line
   (`&& \text{multiply out the bracket}`). Every line does one thing.
6. **Colour roles, consistently, never alone.** `\now{...}` (blue) is known at the start of a step, `\nxt{...}`
   (orange) is new, computed by the step; in prose use `<span class="now">` / `<span class="nxt">`. Introduce the
   meaning with words and a legend; colour is a second cue, never the only one.
7. **Show every algebra step a reader would have to do on paper.** "It can be shown" is not allowed for anything the
   chapter relies on. Optional derivations go in `<details class="deeper">` ("Going deeper"), which keeps the depth
   without blocking the main path.
8. **A picture for every key result.** Prefer a demo or figure that makes the result visible (phase space for
   integrators, a spectral radius as a shrinking error) and point the prose at what to look for.
9. **Prediction first, check often.** Open with a guess the reader commits to, and close the loop later. After each
   major idea add a short `<div class="check">` question with a collapsible answer.
10. **Several routes to hard results.** If the standard proof needs tools the reader may lack (eigenvalues), give an
    elementary route first (expand and cancel; a discriminant) and the standard route as "Going deeper".

## 5. Code

Every chapter with an exercise has `code/<id>/starter.js` and `code/<id>/solution.js`; a chapter may add other demo
files (e.g. `code/<id>/comparison.js`). Each file is a classic script (no modules, no imports) using the sandbox
global `Sim`. It runs at `code/run.html?src=<id>/<file>.js`.

- **The starter runs without errors** and shows something, with clearly marked `// TODO:` blocks for the reader. It
  contains all the scaffolding (UI, drawing, setup) so the reader writes only the physics the chapter teaches.
- **The solution** is the starter with the TODOs filled in, nothing else different (diffing them must show only the
  exercise). Keep both files in sync.
- Plain JavaScript that ports easily to C#, C++ or Python: functions and plain objects or typed arrays, no clever
  idioms, no external libraries. Comment the physics, not the syntax.
- SI units in the code (metres, kilograms, seconds). Gravity 9.81 m/s^2 downward (y is up in the sandbox).
- Deterministic: use `Sim.random(seed)` instead of `Math.random()`.
- Performance: the demo should run at 60 fps on a laptop. For particles in the thousands use flat `Float64Array`s
  and `g.points(...)`. Keep MPM grids around 63-128 cells across and particle counts to a few thousand.

### Sandbox API (`code/lib/sandbox.js`)

```js
// Controls: call at the top level before Sim.run. Each returns { value } (always current).
const k   = Sim.slider('stiffness k (N/m)', 0, 1000, 100, { step: 1, reset: false });
const c   = Sim.logSlider('compliance', -8, -2, 1e-5, { reset: true });   // exponents of ten
const on  = Sim.checkbox('gravity', true, { reset: false });
const m   = Sim.select('method', ['PBD', 'XPBD'], 'XPBD', { reset: true });
Sim.button('kick', () => { ... });
// { reset: true } resets the scene when the control changes; { onChange: fn } is also available.

// Readouts: call from draw().
Sim.show('label', value, digits?);     // top-left info box, this frame
Sim.plot('energy', value, { color, min, max });  // scrolling graph, one sample per simulated frame

// Mouse, in world coordinates; pressed/released are true for one frame.
Sim.mouse.x, Sim.mouse.y, Sim.mouse.down, Sim.mouse.pressed, Sim.mouse.released

// Vectors on {x, y}: V.v(x, y), V.add, V.sub, V.scale, V.addScaled(a, b, s), V.dot, V.cross (scalar),
// V.len, V.lenSq, V.dist, V.norm, V.perp, V.lerp, V.rotate(a, angle), V.copy

// Misc: Sim.random(seed) -> () => [0, 1); Sim.colors.{blue, red, green, orange, purple, teal, yellow, gray, ink,
// muted}; Sim.palette (array); Sim.time; Sim.frame; Sim.setView([xmin, ymin, xmax, ymax]); Sim.reset(); Sim.format(v)

// Drawing (world coordinates, y up). Styles: { fill, stroke | color, width (px), alpha, dash: [px, px] }.
g.circle(x, y, r, style); g.dot(x, y, radiusPx, color); g.line(x1, y1, x2, y2, style);
g.polyline(points /* [{x,y}] or flat [x0,y0,...] */, { ...style, closed, fill });
g.rect(x, y, w, h, style); g.arrow(x, y, dx, dy, style); g.text(str, x, y, { color, size, align, baseline });
g.points(xs, ys, count, { radius /* px */, color, colors /* per-point array */ }); g.grid(spacing); g.ground(y);
g.pxPerUnit; g.ctx (raw 2D context, device pixels); g.toScreen(x, y); g.toWorld(px, py)

// Run. dt and stepsPerFrame may be numbers or functions (e.g. read from a slider).
Sim.run({ title, view: [xmin, ymin, xmax, ymax], reset, step(dt), draw(g), dt: 1/60, stepsPerFrame: 1,
          paused: false, grid: true | false | spacing, background, showTime: true, onKey(e),
          plotCorner: 'bottom-right' | 'bottom-left' | 'top-right' | 'top-left' });
// Every control has .set(value) to change it from code.
```

Keys: Space play/pause, `.` step one frame, R reset. Errors in `step`/`draw` pause the run and are shown on the canvas.

## 6. Notation (use it everywhere, so chapters and papers line up)

| Symbol | Meaning |
|---|---|
| $h$ | time step (s). Substep: $h_s = h / n$. |
| $x_i, v_i, m_i, w_i = 1/m_i$ | position, velocity, mass, inverse mass of particle $i$ (bold or not; be consistent in a chapter) |
| $C(\mathbf{x})$, $\nabla C$ | constraint function and its gradient |
| $\lambda$, $\alpha$, $\tilde\alpha = \alpha/h^2$ | Lagrange multiplier, compliance (1/stiffness), time-scaled compliance |
| $k$, $c$, $\omega = \sqrt{k/m}$, $\zeta$ | stiffness, damping coefficient, angular frequency, damping ratio |
| $E$, $\nu$, $\mu$, $\lambda$ (Lamé) | Young's modulus, Poisson's ratio, Lamé parameters; say "Lamé $\lambda$" when it could clash with the multiplier |
| $EI$, $GJ$, $EA$ | bending, torsional, axial stiffness of a rod |
| $\mu_s$, $\mu_k$ | static and kinetic friction coefficients |
| MPM: particle $p$, grid node $i$ | following the papers |
| $x_p, v_p, m_p, V_p^0$ | particle position, velocity, mass, initial volume |
| $\mathbf{F}_p$, $J = \det\mathbf{F}$ | deformation gradient, volume ratio |
| $\mathbf{C}_p$, $\mathbf{B}_p$, $\mathbf{D}_p$ | APIC affine velocity matrix $\mathbf{C}_p = \mathbf{B}_p\mathbf{D}_p^{-1}$; $\mathbf{D}_p = \tfrac14\Delta x^2\mathbf{I}$ (quadratic) or $\tfrac13\Delta x^2\mathbf{I}$ (cubic) |
| $w_{ip}$, $\nabla w_{ip}$, $\Delta x$ | interpolation weight and gradient, grid spacing |
| $\mathbf{P}$, $\boldsymbol\sigma$, $\Psi$ | first Piola-Kirchhoff stress, Cauchy stress, energy density |
| PB-MPM: $\mathbf{d}_p$, $\mathbf{D}^{\text{def}}_p$ | particle displacement over the step ($= h\,v$), and the affine "deformation displacement" (the paper writes $\mathbf{D}$; call it the *deformation displacement* and warn about the clash with APIC's $\mathbf{D}_p$) |

## 7. Checking your work

- Open `http://localhost:8765/chapters/<id>.html` and `http://localhost:8765/code/run.html?src=<id>/<file>.js`
  (a local server serves the repo root on port 8765). If you have the browser tools, open your **own** tab
  (`tabs_create`), pass its `tabId` on every call, and never touch other tabs. Check: no console errors
  (`read_console_messages`), the demo looks right in a screenshot, the starter runs, the maths renders. Close your tab
  when done.
- Check every number in the prose by computing it (the Python on this machine works: `python3 -c "..."`).
- Check the physics of the demo against a closed form where one exists (a period, a deflection, a ratio), and say so
  in the chapter.
- Re-read the confidentiality rules in section 2 before finishing.
