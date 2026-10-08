/* Page behaviour: theme, navigation, the name sequence, list filters,
   copy buttons and quick navigation. Everything here is an enhancement;
   the page reads fine without it. */
(function () {
  'use strict';

  var root = document.documentElement;
  var $ = function (sel, ctx) { return (ctx || document).querySelector(sel); };
  var $$ = function (sel, ctx) { return Array.prototype.slice.call((ctx || document).querySelectorAll(sel)); };

  /* ---------- Theme ---------- */
  var themeBtn = $('#theme-toggle');
  function setTheme(theme, persist) {
    root.dataset.theme = theme;
    if (persist) { try { localStorage.setItem('theme', theme); } catch (e) {} }
    if (themeBtn) themeBtn.setAttribute('aria-label', theme === 'dark' ? 'Switch to light theme' : 'Switch to dark theme');
    document.dispatchEvent(new CustomEvent('themechange', { detail: theme }));
  }
  function toggleTheme() { setTheme(root.dataset.theme === 'dark' ? 'light' : 'dark', true); }
  if (themeBtn) {
    themeBtn.addEventListener('click', toggleTheme);
    themeBtn.setAttribute('aria-label', root.dataset.theme === 'dark' ? 'Switch to light theme' : 'Switch to dark theme');
  }
  if (window.matchMedia) {
    var mq = matchMedia('(prefers-color-scheme: dark)');
    var follow = function (e) {
      var saved = null;
      try { saved = localStorage.getItem('theme'); } catch (err) {}
      if (!saved) setTheme(e.matches ? 'dark' : 'light', false);
    };
    if (mq.addEventListener) mq.addEventListener('change', follow);
  }

  /* ---------- Navigation ---------- */
  var nav = $('#nav');
  var navLinks = $('#nav-links');
  var menuBtn = $('#menu-toggle');
  var progress = $('#progress-bar');

  function closeMenu() {
    if (!navLinks) return;
    navLinks.classList.remove('is-open');
    if (menuBtn) menuBtn.setAttribute('aria-expanded', 'false');
  }
  if (menuBtn && navLinks) {
    menuBtn.addEventListener('click', function () {
      var open = navLinks.classList.toggle('is-open');
      menuBtn.setAttribute('aria-expanded', String(open));
    });
    navLinks.addEventListener('click', function (e) { if (e.target.closest('a')) closeMenu(); });
    document.addEventListener('keydown', function (e) { if (e.key === 'Escape') closeMenu(); });
  }

  var ticking = false;
  function onScroll() {
    ticking = false;
    var y = window.scrollY || 0;
    if (nav) nav.classList.toggle('is-stuck', y > 8);
    if (progress) {
      var max = document.documentElement.scrollHeight - window.innerHeight;
      progress.style.transform = 'scaleX(' + (max > 0 ? Math.min(1, y / max) : 0) + ')';
    }
  }
  window.addEventListener('scroll', function () {
    if (!ticking) { ticking = true; requestAnimationFrame(onScroll); }
  }, { passive: true });
  window.addEventListener('resize', onScroll);
  onScroll();

  // highlight the link for the section currently in view
  var linkFor = {};
  $$('#nav-links a').forEach(function (a) { linkFor[a.getAttribute('href').slice(1)] = a; });
  if ('IntersectionObserver' in window) {
    var spy = new IntersectionObserver(function (entries) {
      entries.forEach(function (entry) {
        if (!entry.isIntersecting) return;
        var id = entry.target.id;
        Object.keys(linkFor).forEach(function (key) {
          if (key === id) linkFor[key].setAttribute('aria-current', 'true');
          else linkFor[key].removeAttribute('aria-current');
        });
      });
    }, { rootMargin: '-35% 0px -60% 0px' });
    $$('main section[id]').forEach(function (s) { spy.observe(s); });
  }

  /* ---------- Name as a peptide ---------- */
  var AA = {
    A: ['Ala', 'Alanine', 'hydrophobic'],
    V: ['Val', 'Valine', 'hydrophobic'],
    I: ['Ile', 'Isoleucine', 'hydrophobic'],
    N: ['Asn', 'Asparagine', 'polar'],
    S: ['Ser', 'Serine', 'polar'],
    H: ['His', 'Histidine', 'basic'],
    K: ['Lys', 'Lysine', 'basic'],
    R: ['Arg', 'Arginine', 'basic'],
    D: ['Asp', 'Aspartate', 'acidic'],
    O: ['Pyl', 'Pyrrolysine', 'special']
  };
  var seq = $('#seq');
  var readout = $('#seq-readout');
  if (seq && readout) {
    var tiles = $$('.res', seq);
    var restText = readout.textContent;
    var sticky = null;

    var say = function (parts) {
      readout.textContent = '';
      parts.forEach(function (part, i) {
        if (i === 0) {
          var b = document.createElement('b');
          b.textContent = part;
          readout.appendChild(b);
        } else {
          readout.appendChild(document.createTextNode(' · ' + part));
        }
      });
    };
    var rest = function () { readout.textContent = restText; };

    var describe = function (tile) {
      var info = AA[tile.dataset.aa];
      var note = tile.dataset.aa === 'O' ? 'the 22nd genetically encoded amino acid' : info[2];
      say([tile.dataset.aa + tile.dataset.n, info[0], info[1], note]);
    };

    tiles.forEach(function (tile, i) {
      tile.style.setProperty('--i', i);
      tile.dataset.n = i + 1;
      tile.addEventListener('pointerenter', function () { describe(tile); });
      tile.addEventListener('pointerleave', function () { if (sticky) describe(sticky); else rest(); });
      tile.addEventListener('click', function () {
        if (sticky) sticky.classList.remove('is-on');
        sticky = sticky === tile ? null : tile;
        if (sticky) { sticky.classList.add('is-on'); describe(sticky); } else rest();
      });
    });

    // hovering a legend entry picks out that class of residue
    $$('.seq-legend li').forEach(function (li) {
      var cls = $('.swatch', li).dataset.class;
      li.addEventListener('pointerenter', function () {
        var letters = [], count = 0;
        tiles.forEach(function (tile) {
          var on = AA[tile.dataset.aa][2] === cls;
          tile.classList.toggle('is-class', on);
          if (on) { count++; if (letters.indexOf(tile.dataset.aa) < 0) letters.push(tile.dataset.aa); }
        });
        seq.classList.add('has-class');
        say([count + (count === 1 ? ' residue' : ' residues'), li.textContent.trim(), letters.join(', ')]);
      });
      li.addEventListener('pointerleave', function () {
        seq.classList.remove('has-class');
        tiles.forEach(function (tile) { tile.classList.remove('is-class'); });
        if (sticky) describe(sticky); else rest();
      });
    });
  }

  /* ---------- Filters (software, publications, skills) ---------- */
  var UNITS = { software: ['tool', 'tools'], pubs: ['publication', 'publications'], skills: ['skill', 'skills'] };

  $$('.filterbar[data-filter]').forEach(function (bar) {
    var name = bar.dataset.filter;
    var lists = $$('[data-filter-list="' + name + '"]');
    var nested = lists.some(function (l) { return l.hasAttribute('data-filter-nested'); });
    var items = [];
    lists.forEach(function (list) {
      $$(nested ? '.skill-list li' : ':scope > li', list).forEach(function (el) {
        items.push({ el: el, text: el.textContent.toLowerCase(), tags: (el.dataset.tags || '').split(' ') });
      });
    });
    var chips = $$('.chip', bar);
    var input = $('[data-search]', bar);
    var count = $('[data-count]', bar);
    var empty = $('[data-empty="' + name + '"]');
    var unit = UNITS[name] || ['item', 'items'];
    var value = 'all';

    function apply() {
      var q = input ? input.value.trim().toLowerCase() : '';
      var shown = 0;
      items.forEach(function (item) {
        var ok = (value === 'all' || item.tags.indexOf(value) >= 0) && (!q || item.text.indexOf(q) >= 0);
        item.el.hidden = !ok;
        if (nested) item.el.classList.toggle('is-match', ok && !!q);
        if (ok) shown++;
      });
      if (nested) {
        $$('.skill-group', lists[0]).forEach(function (group) {
          group.hidden = !$$('.skill-list li', group).some(function (li) { return !li.hidden; });
        });
      } else {
        lists.forEach(function (list) {
          var any = $$(':scope > li', list).some(function (li) { return !li.hidden; });
          var heading = list.previousElementSibling;
          if (heading && heading.classList.contains('sub')) heading.hidden = !any;
        });
      }
      if (count) {
        count.textContent = shown === items.length
          ? items.length + ' ' + unit[1]
          : shown + ' of ' + items.length + ' ' + unit[1];
      }
      if (empty) empty.hidden = shown > 0;
    }

    chips.forEach(function (chip) {
      chip.addEventListener('click', function () {
        value = chip.dataset.value;
        chips.forEach(function (c) { c.setAttribute('aria-pressed', String(c === chip)); });
        apply();
      });
    });
    if (input) input.addEventListener('input', apply);
    apply();
  });

  /* ---------- Copy buttons ---------- */
  function copy(text, btn) {
    var done = function () {
      var span = $('span', btn);
      if (!span) return;
      var old = span.textContent;
      span.textContent = 'Copied';
      setTimeout(function () { span.textContent = old; }, 1600);
    };
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(done, function () {});
      return;
    }
    var ta = document.createElement('textarea');
    ta.value = text;
    ta.setAttribute('readonly', '');
    ta.style.position = 'fixed';
    ta.style.opacity = '0';
    document.body.appendChild(ta);
    ta.select();
    try { if (document.execCommand('copy')) done(); } catch (e) {}
    document.body.removeChild(ta);
  }
  $$('[data-copy]').forEach(function (btn) {
    btn.addEventListener('click', function () {
      var src = $(btn.dataset.copy);
      if (src) copy(src.textContent.replace(/^#.*\n/gm, '').replace(/\n{2,}/g, '\n').trim() + '\n', btn);
    });
  });
  $$('[data-copy-text]').forEach(function (btn) {
    btn.addEventListener('click', function () { copy(btn.dataset.copyText, btn); });
  });

  /* ---------- Experience timeline: keep "Present" current ---------- */
  var xpList = $('.xp-list');
  if (xpList) {
    var now = new Date();
    var years = Math.max(7, now.getFullYear() - 2020 + 1);
    var months = years * 12;
    xpList.style.setProperty('--months', months);
    var axis = $('.xp-axis ol', xpList);
    if (axis) {
      for (var y = 2020 + axis.children.length; y < 2020 + years; y++) {
        var li = document.createElement('li');
        li.textContent = y;
        axis.appendChild(li);
      }
    }
    $$('.xp-track', xpList).forEach(function (t) { t.style.backgroundSize = 'calc(100% / ' + years + ') 100%'; });
    $$('[data-present]', xpList).forEach(function (bar) {
      bar.style.setProperty('--e', Math.min(months, (now.getFullYear() - 2020) * 12 + now.getMonth() + 1));
    });
  }

  /* ---------- Quick navigation ---------- */
  var dialog = $('#palette');
  var openBtn = $('#palette-open');
  if (dialog && typeof dialog.showModal === 'function') {
    var list = $('#palette-list');
    var field = $('#palette-input');
    var go = function (hash) { return function () { var t = $(hash); if (t) t.scrollIntoView(); }; };
    var open = function (url) { return function () { window.open(url, '_blank', 'noopener'); }; };
    var commands = [
      { label: 'About', kind: 'Section', run: go('#about') },
      { label: 'Selected work', kind: 'Section', run: go('#work') },
      { label: 'The PhIP-Seq platform', kind: 'Project', run: go('#work-platform') },
      { label: 'PhIP-Seq antibody profiling', kind: 'Project', run: go('#work-phipseq') },
      { label: 'EpitopeFindeR 2', kind: 'Project', run: go('#work-epitope') },
      { label: 'HAScan', kind: 'Project', run: go('#work-hascan') },
      { label: 'InFlux neutralization pipeline', kind: 'Project', run: go('#work-influx') },
      { label: 'Influenza antigenic cartography', kind: 'Project', run: go('#work-cartography') },
      { label: 'drydock', kind: 'Project', run: go('#work-drydock') },
      { label: 'Molecular dynamics and free-energy landscapes', kind: 'Project', run: go('#work-md') },
      { label: 'Software', kind: 'Section', run: go('#software') },
      { label: 'Experience', kind: 'Section', run: go('#experience') },
      { label: 'Publications', kind: 'Section', run: go('#publications') },
      { label: 'Skills', kind: 'Section', run: go('#skills') },
      { label: 'Education & recognition', kind: 'Section', run: go('#background') },
      { label: 'Get in touch', kind: 'Section', run: go('#contact') },
      { label: 'Download CV', kind: 'Action', run: function () { var a = $('a[download]'); if (a) a.click(); } },
      { label: 'Email avikarkada@gmail.com', kind: 'Action', run: function () { location.href = 'mailto:avikarkada@gmail.com'; } },
      { label: 'Switch colour theme', kind: 'Action', run: toggleTheme },
      { label: 'GitHub', kind: 'Link', run: open('https://github.com/avinashkarkada') },
      { label: 'Google Scholar', kind: 'Link', run: open('https://scholar.google.com/citations?user=KitPr8YAAAAJ&hl=en') },
      { label: 'LinkedIn', kind: 'Link', run: open('https://www.linkedin.com/in/avinash-karkada/') },
      { label: 'ORCID', kind: 'Link', run: open('https://orcid.org/0000-0001-6604-3351') }
    ];
    var matches = [];
    var active = 0;

    var render = function () {
      var q = field.value.trim().toLowerCase();
      matches = commands.filter(function (c) { return !q || (c.label + ' ' + c.kind).toLowerCase().indexOf(q) >= 0; });
      active = Math.min(active, Math.max(0, matches.length - 1));
      list.textContent = '';
      if (!matches.length) {
        var none = document.createElement('li');
        none.className = 'none';
        none.textContent = 'No matches';
        list.appendChild(none);
        field.removeAttribute('aria-activedescendant');
        return;
      }
      matches.forEach(function (c, i) {
        var li = document.createElement('li');
        li.id = 'palette-opt-' + i;
        li.setAttribute('role', 'option');
        li.setAttribute('aria-selected', String(i === active));
        li.appendChild(document.createTextNode(c.label));
        var kind = document.createElement('small');
        kind.textContent = c.kind;
        li.appendChild(kind);
        li.addEventListener('click', function () { run(i); });
        li.addEventListener('pointermove', function () { if (active !== i) { active = i; mark(); } });
        list.appendChild(li);
      });
      field.setAttribute('aria-activedescendant', 'palette-opt-' + active);
    };
    var mark = function () {
      $$('li[role="option"]', list).forEach(function (li, i) {
        li.setAttribute('aria-selected', String(i === active));
        if (i === active) li.scrollIntoView({ block: 'nearest' });
      });
      field.setAttribute('aria-activedescendant', 'palette-opt-' + active);
    };
    var run = function (i) {
      var cmd = matches[i];
      dialog.close();
      if (cmd) cmd.run();
    };
    var show = function () {
      if (dialog.open) return;
      closeMenu();
      field.value = '';
      active = 0;
      render();
      dialog.showModal();
      field.focus();
    };

    field.addEventListener('input', function () { active = 0; render(); });
    field.addEventListener('keydown', function (e) {
      if (e.key === 'ArrowDown') { e.preventDefault(); if (matches.length) { active = (active + 1) % matches.length; mark(); } }
      else if (e.key === 'ArrowUp') { e.preventDefault(); if (matches.length) { active = (active - 1 + matches.length) % matches.length; mark(); } }
      else if (e.key === 'Enter') { e.preventDefault(); run(active); }
    });
    dialog.addEventListener('click', function (e) { if (e.target === dialog) dialog.close(); });
    if (openBtn) openBtn.addEventListener('click', show);
    document.addEventListener('keydown', function (e) {
      var el = document.activeElement;
      var typing = el && (el.tagName === 'INPUT' || el.tagName === 'TEXTAREA' || el.isContentEditable);
      if ((e.key === 'k' || e.key === 'K') && (e.metaKey || e.ctrlKey)) { e.preventDefault(); show(); }
      else if (e.key === '/' && !typing && !e.metaKey && !e.ctrlKey && !e.altKey) { e.preventDefault(); show(); }
    });
  } else if (openBtn) {
    openBtn.hidden = true;
  }

  /* ---------- Reveal on scroll ---------- */
  var revealEls = $$('.sec-head, .sec-intro, .arc, .about-grid > *, .work-text, .work .fig, .sw-grid > .sw, .bg-grid > *');
  if ('IntersectionObserver' in window && !(window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches)) {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (entry) {
        if (entry.isIntersecting) { entry.target.classList.add('is-in'); io.unobserve(entry.target); }
      });
    }, { rootMargin: '0px 0px -8% 0px' });
    revealEls.forEach(function (el) { el.classList.add('reveal'); io.observe(el); });
  }

  var yearEl = $('#year');
  if (yearEl) yearEl.textContent = new Date().getFullYear();
})();
