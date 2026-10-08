/* Small shared helpers for the figures: seeded randomness, SVG creation,
   a tooltip, and an animation loop that only runs while its figure is on screen. */
(function () {
  'use strict';

  var K = window.Kit = {};
  var SVG_NS = 'http://www.w3.org/2000/svg';

  K.reduced = !!(window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches);

  K.clamp = function (v, lo, hi) { return v < lo ? lo : v > hi ? hi : v; };
  K.lerp = function (a, b, t) { return a + (b - a) * t; };

  // mulberry32: small, fast, good enough for simulated demo data
  K.rng = function (seed) {
    var t = seed >>> 0;
    return function () {
      t += 0x6D2B79F5;
      var r = Math.imul(t ^ (t >>> 15), 1 | t);
      r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
      return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
  };

  // standard normal draw (Box-Muller) from a uniform generator
  K.gauss = function (rnd) {
    var u = 0, v = 0;
    while (u === 0) u = rnd();
    while (v === 0) v = rnd();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  };

  K.seed = function () { return (Math.random() * 0xFFFFFFFF) >>> 0; };

  K.css = function (name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  };

  K.svg = function (tag, attrs, parent) {
    var el = document.createElementNS(SVG_NS, tag);
    if (attrs) for (var k in attrs) if (attrs[k] != null) el.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(el);
    return el;
  };

  K.text = function (parent, x, y, str, attrs) {
    var el = K.svg('text', attrs || null, parent);
    el.setAttribute('x', x);
    el.setAttribute('y', y);
    el.textContent = str;
    return el;
  };

  // pointer position in an SVG's viewBox coordinates
  K.svgPoint = function (svg, evt) {
    var m = svg.getScreenCTM();
    if (!m) return { x: 0, y: 0 };
    var p = new DOMPoint(evt.clientX, evt.clientY).matrixTransform(m.inverse());
    return { x: p.x, y: p.y };
  };

  K.fmt = function (n) { return Math.round(n).toLocaleString('en-US'); };

  /* Tooltip. Content always goes in through textContent. */
  var tipEl = null;
  K.tip = {
    show: function (x, y, title, sub) {
      tipEl = tipEl || document.getElementById('tip');
      if (!tipEl) return;
      tipEl.textContent = '';
      var b = document.createElement('b');
      b.textContent = title;
      tipEl.appendChild(b);
      if (sub) {
        var s = document.createElement('span');
        s.textContent = sub;
        tipEl.appendChild(s);
      }
      tipEl.hidden = false;
      var r = tipEl.getBoundingClientRect();
      var left = K.clamp(x + 14, 8, window.innerWidth - r.width - 8);
      var top = y - r.height - 12;
      if (top < 8) top = y + 18;
      tipEl.style.transform = 'translate(' + Math.round(left) + 'px,' + Math.round(top) + 'px)';
    },
    hide: function () {
      tipEl = tipEl || document.getElementById('tip');
      if (tipEl) tipEl.hidden = true;
    }
  };

  /* Reports whether an element is on screen in a visible tab. */
  K.visible = function (el, cb) {
    var onScreen = false, last = null;
    function emit() {
      var v = onScreen && !document.hidden;
      if (v !== last) { last = v; cb(v); }
    }
    if ('IntersectionObserver' in window) {
      new IntersectionObserver(function (entries) {
        onScreen = entries[entries.length - 1].isIntersecting;
        emit();
      }, { threshold: 0.05 }).observe(el);
    } else {
      onScreen = true;
    }
    document.addEventListener('visibilitychange', emit);
    emit();
  };

  /* requestAnimationFrame loop tied to an element. step(dt) gets milliseconds. */
  K.loop = function (el, step) {
    var want = false, vis = false, raf = 0, last = 0;
    function frame(t) {
      raf = 0;
      var dt = last ? Math.min(50, t - last) : 16;
      last = t;
      step(dt);
      if (want && vis && !raf) raf = requestAnimationFrame(frame);
    }
    function sync() {
      if (want && vis && !raf) { last = 0; raf = requestAnimationFrame(frame); }
      else if ((!want || !vis) && raf) { cancelAnimationFrame(raf); raf = 0; }
    }
    K.visible(el, function (v) { vis = v; sync(); });
    return {
      start: function () { want = true; sync(); },
      stop: function () { want = false; sync(); },
      isRunning: function () { return want; },
      isVisible: function () { return vis; }
    };
  };

  /* Runs fn once, the first time el scrolls into view. */
  K.onceVisible = function (el, fn) {
    if (!('IntersectionObserver' in window)) { fn(); return; }
    var io = new IntersectionObserver(function (entries) {
      if (entries.some(function (e) { return e.isIntersecting; })) { io.disconnect(); fn(); }
    }, { threshold: 0.35 });
    io.observe(el);
  };

  K.onTheme = function (fn) { document.addEventListener('themechange', fn); };

  /* Calls cb(isNarrow) now, and again whenever the element's box crosses `limit` pixels wide.
     Charts use it to switch between a wide and a compact layout rather than shrinking text. */
  K.compact = function (el, limit, cb) {
    var narrow = null;
    function check() {
      var w = el.getBoundingClientRect().width;
      if (!w) return;
      var next = w < limit;
      if (next !== narrow) { narrow = next; cb(narrow); }
    }
    if ('ResizeObserver' in window) new ResizeObserver(check).observe(el);
    else window.addEventListener('resize', check);
    check();
    if (narrow === null) { narrow = false; cb(false); }
  };

  /* Play / pause button with an icon and an optional text label. */
  K.setToggle = function (btn, playing, labelPlay, labelPause) {
    if (!btn) return;
    var use = btn.querySelector('use');
    if (use) use.setAttribute('href', playing ? '#i-pause' : '#i-play');
    var span = btn.querySelector('span');
    var label = playing ? labelPause : labelPlay;
    if (span) span.textContent = label;
    else btn.setAttribute('aria-label', label);
  };

  // canvas sized to its box at device resolution; calls draw after each resize
  K.fitCanvas = function (canvas, draw) {
    var state = { w: 0, h: 0, dpr: 1 };
    function resize() {
      var r = canvas.getBoundingClientRect();
      if (!r.width || !r.height) return;
      state.dpr = Math.min(window.devicePixelRatio || 1, 2);
      state.w = r.width;
      state.h = r.height;
      canvas.width = Math.round(r.width * state.dpr);
      canvas.height = Math.round(r.height * state.dpr);
      draw();
    }
    if ('ResizeObserver' in window) new ResizeObserver(resize).observe(canvas);
    else window.addEventListener('resize', resize);
    resize();
    return state;
  };
})();
