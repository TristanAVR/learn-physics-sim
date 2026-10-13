// Builds the sidebar and previous/next links, embeds demos, renders maths and highlights code.
(function () {
  'use strict';

  const body = document.body;
  const current = body.dataset.chapter || '';
  const base = body.dataset.root || '.';
  const chapters = window.CHAPTERS || [];

  function chapterHref(id) { return `${base}/chapters/${id}.html`; }

  function buildSidebar() {
    const nav = document.getElementById('sidebar');
    if (!nav) return;
    let html = `<a class="brand" href="${base}/index.html">Learn Physics Simulation</a>`;
    let part = null;
    for (const c of chapters) {
      if (c.part !== part) {
        if (part !== null) html += '</ul>';
        part = c.part;
        html += `<div class="part">${part}</div><ul>`;
      }
      const n = c.id.slice(0, 2);
      html += `<li${c.id === current ? ' class="current"' : ''}><a href="${chapterHref(c.id)}"><span class="num">${n}</span>${c.title}</a></li>`;
    }
    html += '</ul>';
    nav.innerHTML = html;
    const cur = nav.querySelector('.current');
    if (cur) cur.scrollIntoView({ block: 'center' });
  }

  function buildPager() {
    const i = chapters.findIndex((c) => c.id === current);
    if (i < 0) return;
    const pager = document.createElement('nav');
    pager.className = 'pager';
    const prev = chapters[i - 1], next = chapters[i + 1];
    pager.innerHTML =
      (prev ? `<a class="prev" href="${chapterHref(prev.id)}">&larr; ${prev.title}</a>` : '<span></span>') +
      (next ? `<a class="next" href="${chapterHref(next.id)}">${next.title} &rarr;</a>` : '<span></span>');
    const article = document.querySelector('article');
    if (article) article.appendChild(pager);
  }

  function buildToc() {
    const list = document.getElementById('toc');
    if (!list) return;
    let html = '';
    let part = null;
    for (const c of chapters) {
      if (c.part !== part) {
        if (part !== null) html += '</ol>';
        part = c.part;
        html += `<h3>${part}</h3><ol>`;
      }
      html += `<li><a href="${chapterHref(c.id)}"><span class="num">${c.id.slice(0, 2)}</span> ${c.title}</a></li>`;
    }
    html += '</ol>';
    list.innerHTML = html;
  }

  // <div class="demo" data-src="04-time-integration/solution.js" data-height="420"></div>
  function embedDemos() {
    for (const d of document.querySelectorAll('.demo[data-src]')) {
      const src = d.dataset.src;
      const h = d.dataset.height || 440;
      const url = `${base}/code/run.html?embed=1&src=${encodeURIComponent(src)}`;
      const caption = d.dataset.caption ? `<div class="demo-caption">${d.dataset.caption}</div>` : '';
      d.innerHTML =
        `<div class="demo-frame" style="height:${h}px"><button class="demo-start">&#9654; Run the demo</button></div>` +
        `<div class="demo-links">${caption}<a href="${url.replace('embed=1&', '')}" target="_blank" rel="noopener">open full screen</a></div>`;
      d.querySelector('.demo-start').addEventListener('click', (e) => {
        const frame = e.currentTarget.parentElement;
        frame.innerHTML = `<iframe src="${url}" title="demo" loading="lazy"></iframe>`;
      });
    }
  }

  function renderMath() {
    if (!window.renderMathInElement) return;
    window.renderMathInElement(document.querySelector('article') || document.body, {
      delimiters: [
        { left: '$$', right: '$$', display: true },
        { left: '\\[', right: '\\]', display: true },
        { left: '\\(', right: '\\)', display: false },
        { left: '$', right: '$', display: false },
      ],
      ignoredTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code'],
      throwOnError: false,
      // Colour roles used across chapters: \now{..} is known at the start of a step, \nxt{..} is new
      // (computed by the step). Always pair colour with words or labels; never rely on colour alone.
      macros: {
        // Hex colours without '#': inside a macro, '#3' would mean "parameter 3".
        '\\now': '\\textcolor{3d7be0}{#1}',
        '\\nxt': '\\textcolor{e07b14}{#1}',
      },
    });
  }

  function highlight() {
    if (!window.hljs) return;
    for (const b of document.querySelectorAll('pre code')) window.hljs.highlightElement(b);
  }

  function sidebarToggle() {
    const t = document.getElementById('menu-toggle');
    if (t) t.addEventListener('click', () => body.classList.toggle('menu-open'));
  }

  // Puzzles. Markup (see AUTHORING.md):
  //   <div class="puzzle" id="p-scale" data-answer="2452" data-tol="0.02" data-unit="N/m">
  //     <h4>title</h4> <p>scenario ...</p>
  //     <details class="hint"><summary>Hint 1</summary>...</details> ...
  //     <div class="worked">worked solution</div>
  //   </div>
  // Multiple choice: data-choices="Yes|No" data-answer="No". The worked solution opens when the answer is right,
  // or on request after a first attempt. Solved puzzles are remembered in this browser only.
  function storageGet(k) { try { return localStorage.getItem(k); } catch (e) { return null; } }
  function storageSet(k, v) { try { localStorage.setItem(k, v); } catch (e) { /* private mode: no memory */ } }

  function parseNumber(s) {
    const cleaned = String(s).replace(/,/g, '').replace(/−/g, '-').replace(/\s+/g, ' ').trim();
    const m = cleaned.match(/[-+]?(\d+\.?\d*|\.\d+)(e[-+]?\d+)?/i);
    return m ? parseFloat(m[0]) : NaN;
  }

  function judge(value, answer, tol) {
    if (!isFinite(value)) return { ok: false, invalid: true, msg: 'Enter a number first.' };
    if (answer === 0) return { ok: Math.abs(value) <= tol, msg: '' };
    const ratio = value / answer;
    if (Math.abs(ratio - 1) <= tol) return { ok: true, msg: '' };
    if (ratio < 0) return { ok: false, msg: 'Not yet: check the sign.' };
    const decades = Math.log10(ratio);
    if (Math.abs(decades - Math.round(decades)) < 0.05 && Math.round(decades) !== 0) {
      return { ok: false, msg: `Not yet, but you are off by a factor of ${Math.pow(10, Math.abs(Math.round(decades)))}: check your units (cm or m? g or kg? ms or s?).` };
    }
    return { ok: false, msg: ratio > 1 ? 'Not yet: your answer is too high.' : 'Not yet: your answer is too low.' };
  }

  function setupPuzzles() {
    const pageKey = 'puzzle:' + (current || location.pathname) + ':';
    let n = 0;
    for (const pz of document.querySelectorAll('.puzzle')) {
      const id = pageKey + (pz.id || 'p' + n++);
      const worked = pz.querySelector('.worked');
      const firstHint = pz.querySelector('details');
      const form = document.createElement('div');
      form.className = 'puzzle-form';
      const choices = pz.dataset.choices ? pz.dataset.choices.split('|') : null;
      let input = null;
      if (choices) {
        form.innerHTML = choices.map((c, i) =>
          `<label class="puzzle-choice"><input type="radio" name="${id}" value="${c.replace(/"/g, '&quot;')}"> ${c}</label>`).join('');
      } else {
        form.innerHTML = `<input type="text" inputmode="decimal" autocomplete="off" aria-label="Your answer"> <span class="puzzle-unit">${pz.dataset.unit || ''}</span>`;
        input = form.querySelector('input');
      }
      const check = document.createElement('button');
      check.className = 'puzzle-check';
      check.textContent = 'Check';
      form.appendChild(check);
      const reveal = document.createElement('button');
      reveal.className = 'puzzle-reveal';
      reveal.textContent = 'Show the worked solution';
      reveal.hidden = true;
      form.appendChild(reveal);
      const feedback = document.createElement('div');
      feedback.className = 'puzzle-feedback';
      feedback.setAttribute('aria-live', 'polite');
      if (firstHint) pz.insertBefore(form, firstHint); else if (worked) pz.insertBefore(form, worked); else pz.appendChild(form);
      form.after(feedback);
      if (worked) worked.hidden = true;

      const solved = () => {
        pz.classList.add('solved');
        feedback.textContent = 'Correct.' + (worked ? ' Compare your reasoning with the worked solution.' : '');
        feedback.className = 'puzzle-feedback ok';
        if (worked) worked.hidden = false;
        reveal.hidden = true;
        storageSet(id, '1');
      };
      check.addEventListener('click', () => {
        let res;
        if (choices) {
          const picked = form.querySelector('input:checked');
          if (!picked) { feedback.textContent = 'Pick an answer first.'; feedback.className = 'puzzle-feedback'; return; }
          res = { ok: picked.value === pz.dataset.answer, msg: 'Not yet. Try a hint, or think about it once more.' };
        } else {
          res = judge(parseNumber(input.value), parseFloat(pz.dataset.answer), parseFloat(pz.dataset.tol || '0.02'));
        }
        if (res.ok) { solved(); return; }
        feedback.textContent = res.msg;
        feedback.className = 'puzzle-feedback no';
        if (worked && !res.invalid) reveal.hidden = false;
      });
      if (input) {
        input.addEventListener('keydown', (e) => { if (e.key === 'Enter') check.click(); });
        input.addEventListener('input', () => { feedback.textContent = ''; feedback.className = 'puzzle-feedback'; });
      }
      reveal.addEventListener('click', () => { if (worked) worked.hidden = false; reveal.hidden = true; });
      if (storageGet(id) === '1') {
        pz.classList.add('solved');
        feedback.textContent = 'You solved this one before. Try it again from memory.';
        feedback.className = 'puzzle-feedback';
      }
    }
  }

  document.addEventListener('DOMContentLoaded', () => {
    buildSidebar();
    buildToc();
    buildPager();
    embedDemos();
    setupPuzzles();
    renderMath();
    highlight();
    sidebarToggle();
  });
})();
