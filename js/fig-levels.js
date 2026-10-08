/* Fig. 2: differential antibody reactivity at two levels, on a simulated cohort.
   160 cases and 160 controls, 30 species with 12 peptides each. A person exposed
   to a species reacts to some of its peptides. For three planted species, cases
   are more often exposed. Each peptide is tested on its own (Fisher's exact
   test), then each species on its count of reactive peptides per person
   (Mann-Whitney test), with Benjamini-Hochberg false-discovery control at each
   level. Switching level moves every peptide's point onto its species. */
(function () {
  'use strict';

  var K = window.Kit;
  var svg = document.getElementById('levels-plot');
  var fig = document.getElementById('fig-levels');
  if (!K || !svg || !fig) return;

  /* ---- simulation and tests ---- */
  var SPECIES = 30, PEPTIDES = 12, GROUP = 160, PLANTED = 3, LIFT = 0.28, ALPHA = 0.05;

  var logFact = [0];
  for (var f = 1; f <= 2 * GROUP + 2; f++) logFact.push(logFact[f - 1] + Math.log(f));

  // two-sided Fisher exact test: a reactive of GROUP in one arm, b of GROUP in the other
  function fisher(a, b) {
    var k = a + b, lo = Math.max(0, k - GROUP), hi = Math.min(k, GROUP);
    var base = 2 * logFact[GROUP] - logFact[2 * GROUP] + logFact[k] + logFact[2 * GROUP - k];
    function logProb(x) { return base - logFact[x] - logFact[GROUP - x] - logFact[k - x] - logFact[GROUP - k + x]; }
    var observed = logProb(a), p = 0;
    for (var x = lo; x <= hi; x++) { var v = logProb(x); if (v <= observed + 1e-7) p += Math.exp(v); }
    return Math.min(1, p);
  }

  // complementary error function, relative error below 1.2e-7 (Chebyshev fit from Numerical Recipes)
  function erfc(x) {
    var z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    var r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 +
      t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }

  // two-sided Mann-Whitney test, normal approximation with tie and continuity corrections
  function mannWhitney(xs, ys) {
    var all = [], i;
    for (i = 0; i < xs.length; i++) all.push([xs[i], 0]);
    for (i = 0; i < ys.length; i++) all.push([ys[i], 1]);
    all.sort(function (a, b) { return a[0] - b[0]; });
    var n = all.length, rankSum = 0, ties = 0;
    for (i = 0; i < n;) {
      var j = i;
      while (j < n && all[j][0] === all[i][0]) j++;
      var rank = (i + j + 1) / 2, t = j - i;
      ties += t * t * t - t;
      for (var k = i; k < j; k++) if (all[k][1] === 0) rankSum += rank;
      i = j;
    }
    var n1 = xs.length, n2 = ys.length;
    var u = rankSum - n1 * (n1 + 1) / 2, mean = n1 * n2 / 2;
    var variance = n1 * n2 / 12 * ((n + 1) - ties / (n * (n - 1)));
    if (variance <= 0) return 1;
    var z = (Math.abs(u - mean) - 0.5) / Math.sqrt(variance);
    return z <= 0 ? 1 : erfc(z / Math.SQRT2);
  }

  // Benjamini-Hochberg: which tests pass, and the p-value cut-off that decided it
  function benjaminiHochberg(ps) {
    var m = ps.length, cut = -1, k;
    var order = ps.map(function (p, i) { return i; }).sort(function (a, b) { return ps[a] - ps[b]; });
    for (k = 0; k < m; k++) if (ps[order[k]] <= ALPHA * (k + 1) / m) cut = k;
    var pass = [];
    for (k = 0; k < m; k++) pass.push(false);
    for (k = 0; k <= cut; k++) pass[order[k]] = true;
    return { pass: pass, threshold: ALPHA * Math.max(1, cut + 1) / m };
  }

  function simulate(seed) {
    var rnd = K.rng(seed), planted = [], s, p, g, i;
    while (planted.length < PLANTED) {
      s = Math.floor(rnd() * SPECIES);
      if (planted.indexOf(s) < 0) planted.push(s);
    }
    var pep = { x: [], p: [], rate: [] }, sp = { x: [], p: [] };
    for (s = 0; s < SPECIES; s++) {
      var exposure = 0.25 + 0.4 * rnd(), isPlanted = planted.indexOf(s) >= 0;
      var response = [], hits = [], counts = [new Uint8Array(GROUP), new Uint8Array(GROUP)];
      for (p = 0; p < PEPTIDES; p++) { response.push(0.15 + 0.45 * rnd()); hits.push([0, 0]); }
      for (g = 0; g < 2; g++) {                 // 0 = cases, 1 = controls
        for (i = 0; i < GROUP; i++) {
          var exposed = rnd() < (g === 0 && isPlanted ? Math.min(0.95, exposure + LIFT) : exposure);
          for (p = 0; p < PEPTIDES; p++) {
            if (rnd() < (exposed ? response[p] : 0.02)) { hits[p][g]++; counts[g][i]++; }
          }
        }
      }
      var mean = 0;
      for (p = 0; p < PEPTIDES; p++) {
        var diff = (hits[p][0] - hits[p][1]) / GROUP * 100;
        pep.x.push(diff);
        pep.p.push(fisher(hits[p][0], hits[p][1]));
        pep.rate.push([hits[p][0] / GROUP * 100, hits[p][1] / GROUP * 100]);
        mean += diff / PEPTIDES;
      }
      sp.x.push(mean);                          // a species sits at the mean of its peptides
      sp.p.push(mannWhitney(Array.prototype.slice.call(counts[0]), Array.prototype.slice.call(counts[1])));
    }
    pep.fdr = benjaminiHochberg(pep.p);
    sp.fdr = benjaminiHochberg(sp.p);
    return { planted: planted, pep: pep, sp: sp };
  }
  /* ---- chart ---- */

  var N = SPECIES * PEPTIDES, X_MAX = 32, Y_MAX = 8;
  var W, H, L, R, T, B, narrow = false;
  var xOf = function (v) { return L + (K.clamp(v, -X_MAX, X_MAX) + X_MAX) / (2 * X_MAX) * (W - L - R); };
  var yOf = function (p) { return B - Math.min(Y_MAX, -Math.log10(Math.max(p, 1e-300))) / Y_MAX * (B - T); };

  var data = null, dots = [], rule, ruleLabel, hover = -1;
  var level = 0, mix = 0;                       // 0 = peptide level, 1 = species level; mix eases between them
  var readout = fig.querySelector('[data-readout]');
  var levelBtns = Array.prototype.slice.call(fig.querySelectorAll('[data-level]'));

  function build(isNarrow) {
    narrow = isNarrow;
    W = narrow ? 360 : 640; H = narrow ? 330 : 360;
    L = narrow ? 30 : 44; R = narrow ? 8 : 14; T = 28; B = H - 50;
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.textContent = '';

    var grid = K.svg('g', { 'class': 'grid' }, svg);
    [2, 4, 6, 8].forEach(function (v) {
      var y = B - v / Y_MAX * (B - T);
      K.svg('line', { x1: L, x2: W - R, y1: y, y2: y }, grid);
      K.text(svg, L - 7, y + 4, String(v), { 'text-anchor': 'end' });
    });
    K.text(svg, L - 7, B + 4, '0', { 'text-anchor': 'end' });
    var axis = K.svg('g', { 'class': 'axis' }, svg);
    K.svg('line', { x1: L, x2: W - R, y1: B, y2: B }, axis);
    K.svg('line', { x1: xOf(0), x2: xOf(0), y1: T, y2: B }, axis);
    [-30, -20, -10, 0, 10, 20, 30].forEach(function (v) {
      if (narrow && Math.abs(v) === 10) return;
      K.svg('line', { x1: xOf(v), x2: xOf(v), y1: B, y2: B + 5 }, axis);
      K.text(svg, xOf(v), B + 20, (v > 0 ? '+' : '') + v, { 'text-anchor': 'middle' });
    });
    K.text(svg, L, 12, '−log10 p', { 'class': 'lbl' });
    K.text(svg, W - R, H - 8, narrow ? 'cases minus controls, % points' : 'reactivity in cases minus controls, percentage points', { 'text-anchor': 'end' });

    rule = K.svg('line', { 'class': 'rule', x1: L, x2: W - R }, svg);
    var layer = K.svg('g', null, svg);
    ruleLabel = K.text(svg, W - R, 0, '5% FDR', { 'text-anchor': 'end', 'class': 'lbl-strong' });
    dots = [];
    for (var i = 0; i < N; i++) dots.push(K.svg('circle', { 'class': 'dot' }, layer));
    if (data) {
      // planted signals are drawn last so they sit on top
      dots.forEach(function (d, j) { if (data.planted.indexOf(Math.floor(j / PEPTIDES)) >= 0) layer.appendChild(d); });
    }

    var pad = K.svg('rect', { 'class': 'hit', x: L, y: T - 10, width: W - L - R, height: B - T + 10 }, svg);
    pad.addEventListener('pointermove', function (e) {
      if (e.pointerType === 'touch') return;
      var pt = K.svgPoint(svg, e), best = -1, bestD = 16 * 16;
      for (var j = 0; j < N; j++) {
        var c = position(j), dx = c.x - pt.x, dy = c.y - pt.y, d = dx * dx + dy * dy;
        if (d < bestD) { bestD = d; best = j; }
      }
      if (best !== hover) { hover = best; render(); }
      if (best < 0) { K.tip.hide(); return; }
      var s = Math.floor(best / PEPTIDES);
      if (level === 0) {
        K.tip.show(e.clientX, e.clientY, 'p = ' + pretty(data.pep.p[best]) + (data.pep.fdr.pass[best] ? ' · passes' : ''),
          'peptide ' + (best % PEPTIDES + 1) + ' of species ' + (s + 1) + ' · cases ' + Math.round(data.pep.rate[best][0]) + '%, controls ' + Math.round(data.pep.rate[best][1]) + '%');
      } else {
        K.tip.show(e.clientX, e.clientY, 'p = ' + pretty(data.sp.p[s]) + (data.sp.fdr.pass[s] ? ' · passes' : ''),
          'species ' + (s + 1) + ' · ' + PEPTIDES + ' peptides');
      }
    });
    pad.addEventListener('pointerleave', function () { hover = -1; K.tip.hide(); render(); });
  }

  function pretty(p) { return p >= 0.001 ? p.toFixed(3) : p.toExponential(1); }

  function ease(t) { return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2; }

  // where peptide i is drawn right now: its own spot, its species' spot, or somewhere between
  function position(i) {
    var s = Math.floor(i / PEPTIDES), e = ease(mix);
    return {
      x: K.lerp(xOf(data.pep.x[i]), xOf(data.sp.x[s]), e),
      y: K.lerp(yOf(data.pep.p[i]), yOf(data.sp.p[s]), e)
    };
  }

  function render() {
    if (!data || !dots.length) return;
    var e = ease(mix), bySpecies = mix >= 0.5;
    var small = narrow ? 2.4 : 3, large = narrow ? 5 : 6.5;
    for (var i = 0; i < N; i++) {
      var s = Math.floor(i / PEPTIDES), c = position(i), d = dots[i];
      d.setAttribute('cx', c.x.toFixed(1));
      d.setAttribute('cy', c.y.toFixed(1));
      d.setAttribute('r', K.lerp(small, large, e).toFixed(2));
      d.setAttribute('class', 'dot' +
        ((bySpecies ? data.sp.fdr.pass[s] : data.pep.fdr.pass[i]) ? ' is-on' : '') +
        (data.planted.indexOf(s) >= 0 ? ' is-planted' : '') +
        (hover >= 0 && (bySpecies ? Math.floor(hover / PEPTIDES) === s : hover === i) ? ' is-hover' : ''));
    }
    var y = K.lerp(yOf(data.pep.fdr.threshold), yOf(data.sp.fdr.threshold), e);
    rule.setAttribute('y1', y); rule.setAttribute('y2', y);
    ruleLabel.setAttribute('y', y - 6);
  }

  function describe() {
    if (!readout || !data) return;
    var hit = 0, other = 0, i;
    if (level === 0) {
      for (i = 0; i < N; i++) if (data.pep.fdr.pass[i]) { if (data.planted.indexOf(Math.floor(i / PEPTIDES)) >= 0) hit++; else other++; }
      readout.textContent = 'Peptide level: ' + hit + ' of ' + PLANTED * PEPTIDES + ' planted peptides pass';
    } else {
      for (i = 0; i < SPECIES; i++) if (data.sp.fdr.pass[i]) { if (data.planted.indexOf(i) >= 0) hit++; else other++; }
      readout.textContent = 'Species level: ' + hit + ' of ' + PLANTED + ' planted species pass';
    }
    if (other) readout.textContent += ' · ' + other + ' other';
  }

  /* ---- switching level ---- */
  var loop = K.loop(svg, function (dt) {
    var step = dt / 900;
    mix = level ? Math.min(1, mix + step) : Math.max(0, mix - step);
    render();
    if (mix === level) loop.stop();
  });

  function setLevel(next) {
    level = next;
    levelBtns.forEach(function (b) { b.setAttribute('aria-pressed', String(+b.dataset.level === level)); });
    describe();
    if (K.reduced) { mix = level; render(); } else loop.start();
  }

  levelBtns.forEach(function (btn) {
    btn.addEventListener('click', function () { setLevel(+btn.dataset.level); });
  });

  var fresh = fig.querySelector('[data-act="new"]');
  if (fresh) fresh.addEventListener('click', function () {
    data = simulate(K.seed());
    build(narrow);
    render();
    describe();
  });

  data = simulate(57);
  K.compact(svg.parentNode, 440, function (isNarrow) { build(isNarrow); render(); });
  describe();
  // show the point of the figure once: start at peptide level, then gather into species
  if (!K.reduced) K.onceVisible(fig, function () { setTimeout(function () { if (level === 0) setLevel(1); }, 1600); });
})();
