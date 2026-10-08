/* Fig: what HAScan reads out. Antibody binding at each position along HA, across a
   panel of strains, for one serum sample. Stem positions are conserved, so they
   react across the whole panel. Head positions are strain-specific, so they react
   only on strains close to what that donor has met. Switching donor moves the head
   pattern and leaves the stem band where it is. Simulated donors. */
(function () {
  'use strict';

  var K = window.Kit;
  var svg = document.getElementById('hascan-plot');
  var fig = document.getElementById('fig-hascan');
  if (!K || !svg || !fig) return;

  var POSITIONS = 60, HEAD = 40;           // columns 0..39 are HA1 head, 40..59 HA2 stem
  var SUBTYPES = [
    { name: 'H1N1', years: [1977, 1983, 1988, 1995, 1999, 2003, 2007, 2009, 2013, 2017, 2021, 2024] },
    { name: 'H3N2', years: [1975, 1980, 1987, 1993, 1998, 2004, 2008, 2012, 2016, 2019, 2022, 2024] }
  ];
  var STRAINS = [];
  SUBTYPES.forEach(function (st, g) {
    st.years.forEach(function (y) { STRAINS.push({ subtype: st.name, group: g, year: y }); });
  });
  var ROWS = STRAINS.length;

  // a donor is a list of past exposures: which subtype, and roughly when
  var DONORS = [
    { id: 'A', exposures: [{ group: 0, year: 1982 }, { group: 1, year: 1999 }, { group: 0, year: 2010 }] },
    { id: 'B', exposures: [{ group: 1, year: 1988 }, { group: 1, year: 2016 }, { group: 0, year: 2021 }] }
  ];
  var donor = DONORS[0];

  /* ---- which positions are conserved and which drift ---- */
  var site = (function () {
    var rnd = K.rng(31), out = [];
    for (var p = 0; p < POSITIONS; p++) {
      if (p >= HEAD) {
        // stem: mostly conserved, a few positions simply not immunogenic
        out.push({ conserved: true, strength: rnd() < 0.45 ? 0.55 + 0.4 * rnd() : 0.05 });
      } else {
        // head: antibody targets that differ between strains
        out.push({ conserved: false, strength: rnd() < 0.5 ? 0.5 + 0.5 * rnd() : 0.06, drift: 6 + 10 * rnd() });
      }
    }
    return out;
  })();

  // reactivity of one donor's serum against one strain at one position
  function reactivity(d, s, p, rnd) {
    var sp = site[p], signal;
    if (sp.conserved) {
      signal = sp.strength * (0.75 + 0.25 * rnd());
    } else {
      var best = 0;
      d.exposures.forEach(function (e) {
        if (e.group !== s.group) return;
        // the further the strain is from a strain the donor met, the weaker the match
        best = Math.max(best, Math.exp(-Math.pow((s.year - e.year) / sp.drift, 2)));
      });
      signal = sp.strength * best;
    }
    return K.clamp(signal + 0.06 * (rnd() - 0.5), 0, 1);
  }

  function matrix(d) {
    var rnd = K.rng(7 + d.id.charCodeAt(0)), m = [], counts = [];
    for (var r = 0; r < ROWS; r++) {
      var row = new Float32Array(POSITIONS), n = 0;
      for (var p = 0; p < POSITIONS; p++) {
        row[p] = reactivity(d, STRAINS[r], p, rnd);
        if (p < HEAD && row[p] > 0.35) n++;   // head sites only: the strain-specific part
      }
      m.push(row);
      counts.push(n);
    }
    return { m: m, counts: counts };
  }
  var data = matrix(donor);

  /* ---- colour ---- */
  var RAMP_LIGHT = [[248, 247, 242], [205, 226, 251], [134, 182, 239], [57, 135, 229], [24, 79, 149]];
  var RAMP_DARK = [[26, 29, 31], [28, 92, 171], [57, 135, 229], [134, 182, 239], [205, 226, 251]];
  function colour(v) {
    var ramp = document.documentElement.dataset.theme === 'dark' ? RAMP_DARK : RAMP_LIGHT;
    var t = K.clamp(v, 0, 1) * (ramp.length - 1), i = Math.min(ramp.length - 2, Math.floor(t)), f = t - i;
    var a = ramp[i], b = ramp[i + 1];
    return 'rgb(' + Math.round(a[0] + (b[0] - a[0]) * f) + ',' + Math.round(a[1] + (b[1] - a[1]) * f) + ',' + Math.round(a[2] + (b[2] - a[2]) * f) + ')';
  }

  /* ---- chart ---- */
  var W, H, L, R, T, CELL, ROW_H, narrow = false;
  var GROUP_GAP = 16;
  var rowY = function (r) { return T + r * ROW_H + (STRAINS[r].group ? GROUP_GAP : 0); };
  var cells = [], bars = [], hover = null;
  var readout = fig.querySelector('[data-readout]');
  var donorBtns = Array.prototype.slice.call(fig.querySelectorAll('[data-donor]'));

  function build(isNarrow) {
    narrow = isNarrow;
    W = narrow ? 360 : 640;
    L = narrow ? 40 : 58; R = narrow ? 46 : 74; T = narrow ? 40 : 46;
    ROW_H = narrow ? 9 : 11;
    CELL = (W - L - R) / POSITIONS;
    H = T + ROWS * ROW_H + GROUP_GAP + (narrow ? 24 : 28);
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.textContent = '';

    // region brackets across the top
    var regions = [[0, HEAD, narrow ? 'HA1 head' : 'HA1 head · variable'], [HEAD, POSITIONS, narrow ? 'HA2 stem' : 'HA2 stem · conserved']];
    var axis = K.svg('g', { 'class': 'axis' }, svg);
    regions.forEach(function (reg) {
      var x1 = L + reg[0] * CELL, x2 = L + reg[1] * CELL;
      K.svg('path', { d: 'M' + (x1 + 1) + ' ' + (T - 10) + 'v4H' + (x2 - 1) + 'v-4' }, axis);
      K.text(svg, (x1 + x2) / 2, T - 15, reg[2], { 'text-anchor': 'middle' });
    });
    K.text(svg, L, H - 6, 'position along HA');

    cells = []; bars = [];
    var grid = K.svg('g', null, svg);
    var barLayer = K.svg('g', null, svg);
    STRAINS.forEach(function (s, r) {
      var y = rowY(r);
      if (r === 0 || STRAINS[r - 1].group !== s.group) K.text(svg, 0, y - 4, s.subtype, { 'class': 'lbl-strong' });
      if (!narrow || r % 2 === 0) K.text(svg, L - 6, y + ROW_H - 2.5, String(s.year), { 'text-anchor': 'end' });
      var row = [];
      for (var p = 0; p < POSITIONS; p++) {
        row.push(K.svg('rect', { 'class': 'cell', x: L + p * CELL, y: y, width: CELL + 0.6, height: ROW_H - 1.4 }, grid));
      }
      cells.push(row);
      bars.push(K.svg('rect', { 'class': 'rowbar', x: W - R + 8, y: y + 1, height: ROW_H - 3.4, rx: 1.5 }, barLayer));
    });
    K.text(svg, W - 2, T - 15, narrow ? 'head' : 'head sites', { 'text-anchor': 'end' });

    var pad = K.svg('rect', { 'class': 'hit', x: L, y: T, width: POSITIONS * CELL, height: ROWS * ROW_H + GROUP_GAP }, svg);
    pad.addEventListener('pointermove', function (e) {
      var pt = K.svgPoint(svg, e);
      var p = K.clamp(Math.floor((pt.x - L) / CELL), 0, POSITIONS - 1);
      var r = 0;
      while (r < ROWS - 1 && pt.y > rowY(r) + ROW_H) r++;
      hover = { r: r, p: p };
      paint();
      var s = STRAINS[r], v = data.m[r][p];
      K.tip.show(e.clientX, e.clientY,
        v > 0.35 ? 'Bound' : 'Not bound',
        s.subtype + ' ' + s.year + ' · position ' + (p + 1) + ' · ' + (p < HEAD ? 'HA1 head' : 'HA2 stem'));
    });
    pad.addEventListener('pointerleave', function () { hover = null; K.tip.hide(); paint(); });
  }

  function paint() {
    var maxCount = Math.max.apply(null, data.counts) || 1;
    for (var r = 0; r < ROWS; r++) {
      for (var p = 0; p < POSITIONS; p++) {
        cells[r][p].style.fill = colour(data.m[r][p]);
      }
      bars[r].setAttribute('width', Math.max(0.8, data.counts[r] / maxCount * (R - 14)));
      bars[r].classList.toggle('is-on', !!hover && hover.r === r);
    }
    if (!readout) return;
    readout.textContent = '';
    var stem = 0, head = 0, n = 0;
    for (r = 0; r < ROWS; r++) {
      for (p = 0; p < POSITIONS; p++) {
        if (data.m[r][p] <= 0.35) continue;
        if (p < HEAD) head++; else stem++;
        n++;
      }
    }
    var reactive = data.counts.filter(function (c) { return c >= 3; }).length;
    var b = document.createElement('b');
    b.textContent = 'Donor ' + donor.id;
    readout.appendChild(b);
    readout.appendChild(document.createTextNode(
      ' · head binding on ' + reactive + ' of ' + ROWS + ' strains · stem binding on the whole panel'));
  }

  donorBtns.forEach(function (btn) {
    btn.addEventListener('click', function () {
      donor = DONORS.filter(function (d) { return d.id === btn.dataset.donor; })[0] || DONORS[0];
      donorBtns.forEach(function (b) { b.setAttribute('aria-pressed', String(b === btn)); });
      data = matrix(donor);
      paint();
    });
  });

  K.compact(svg.parentNode, 440, function (isNarrow) { build(isNarrow); paint(); });
  K.onTheme(paint);

  // switch donor once, so the moving head pattern is visible without clicking
  if (!K.reduced) {
    K.onceVisible(fig, function () {
      setTimeout(function () {
        if (donor !== DONORS[0]) return;
        var btn = fig.querySelector('[data-donor="B"]');
        if (btn) btn.click();
      }, 2400);
    });
  }
})();
