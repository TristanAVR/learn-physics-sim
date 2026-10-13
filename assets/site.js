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

  // <div class="demo" data-src="03-time-integration/solution.js" data-height="420"></div>
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

  document.addEventListener('DOMContentLoaded', () => {
    buildSidebar();
    buildToc();
    buildPager();
    embedDemos();
    renderMath();
    highlight();
    sidebarToggle();
  });
})();
