# Learn Physics Simulation

A hands-on course in writing physics simulations, from a falling particle to the Material Point Method, in the style of
[LearnOpenGL](https://learnopengl.com): each chapter introduces one idea, shows it running in the page, and then hands
you a starter file to implement it yourself.

By the end you can read three papers unaided: *The Affine Particle-In-Cell Method* (Jiang et al., SIGGRAPH 2015),
*The Material Point Method for Simulating Continuum Materials* (Jiang et al., SIGGRAPH 2016 course) and
*A Position Based Material Point Method* (Lewin, SIGGRAPH 2024).

## Reading it

Open `index.html`. Everything is static HTML and plain JavaScript; there is no build step.

Some browsers restrict pages opened from disk, so serving the folder is the most reliable:

```
python -m http.server 8000
```

then open http://localhost:8000.

## Doing the exercises

Each chapter's code lives in `code/<chapter>/`: `starter.js` (runs, with `TODO`s for you) and `solution.js`. Run a script
with `code/run.html?src=<chapter>/starter.js`, edit it, and refresh. The small library every script uses is
`code/lib/sandbox.js`; chapter 01 explains it.

## Publishing on GitHub Pages

Push the repository to GitHub, then in the repository settings choose **Pages > Deploy from a branch**, branch `main`,
folder `/ (root)`. The `.nojekyll` file tells Pages to serve the files as they are.

## Layout

| Path | What |
|---|---|
| `index.html` | home page and table of contents |
| `chapters/` | one HTML page per chapter |
| `code/lib/` | the sandbox library |
| `code/<chapter>/` | starter and solution for each chapter |
| `code/run.html` | runs a script: `run.html?src=<chapter>/<file>.js` |
| `assets/` | site styles and scripts; `chapters.js` is the table of contents |
| `AUTHORING.md` | how chapters are written |

## The papers

The papers discussed belong to their authors and publishers; this site summarises and explains them and links to the
originals.
