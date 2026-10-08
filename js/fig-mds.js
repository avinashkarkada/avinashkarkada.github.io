/* Fig. 5: a sequence-based antigenic map of a simulated lineage.
   A toy "antigenic site" of 34 residues drifts over 56 seasons. Pairwise
   distances between the unique sequences come from BLOSUM62, and the map is
   laid out by metric multidimensional scaling (SMACOF iterations, a few per
   frame, so the relaxation is visible). */
(function () {
  'use strict';

  var K = window.Kit;
  var svg = document.getElementById('mds-plot');
  var fig = document.getElementById('fig-mds');
  if (!K || !svg || !fig) return;

  var SITE_LENGTH = 34, SEASONS = 56, PANEL = 8;
  var EASE = 0.3, STEP_MS = 1000 / 120;

  // BLOSUM62, lower triangle, residues in the order ARNDCQEGHILKMFPSTWYV
  var BLOSUM = [4, -1, 5, -2, 0, 6, -2, -2, 1, 6, 0, -3, -3, -3, 9, -1, 1, 0, 0, -3, 5, -1, 0, 0, 2, -4, 2, 5, 0, -2, 0, -1, -3, -2, -2, 6, -2, 0, 1, -1, -3, 0, 0, -2, 8, -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, 1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, 0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4];
  function score(a, b) { return a < b ? BLOSUM[b * (b + 1) / 2 + a] : BLOSUM[a * (a + 1) / 2 + b]; }

  /* ---- simulated lineage ---- */
  var strains = [];           // { seq: residue indices, season, panel }
  (function () {
    var rnd = K.rng(14), trunk = [], seen = {}, i;
    for (i = 0; i < SITE_LENGTH; i++) trunk.push(Math.floor(rnd() * 20));
    function mutate(seq, changes) {
      seq = seq.slice();
      for (var m = 0; m < changes; m++) {
        var p = Math.floor(rnd() * SITE_LENGTH), a;
        do { a = Math.floor(rnd() * 20); } while (a === seq[p]);
        seq[p] = a;
      }
      return seq;
    }
    for (var season = 1; season <= SEASONS; season++) {
      // most seasons the circulating sequence picks up one substitution; now and then it jumps by three
      var r = rnd();
      if (r < 0.12) trunk = mutate(trunk, 3);
      else if (r < 0.62) trunk = mutate(trunk, 1);
      var pool = [trunk];
      if (rnd() < 0.5) pool.push(mutate(trunk, 1));     // a side branch that does not persist
      pool.forEach(function (seq) {
        var key = seq.join(',');
        if (!seen[key]) { seen[key] = true; strains.push({ seq: seq, season: season, panel: false }); }
      });
    }
    // a panel spread evenly through time, like the reference strains an assay would carry
    for (i = 0; i < PANEL; i++) strains[Math.round(i * (strains.length - 1) / (PANEL - 1))].panel = true;
  })();
  var N = strains.length;

  // distance between two sequences: mean over positions of how far the pair's
  // BLOSUM62 score falls below the two self-scores
  var D = [];
  (function () {
    var i, j, p;
    for (i = 0; i < N; i++) D.push(new Float64Array(N));
    for (i = 0; i < N; i++) {
      for (j = i + 1; j < N; j++) {
        var a = strains[i].seq, b = strains[j].seq, sum = 0;
        for (p = 0; p < SITE_LENGTH; p++) sum += (score(a[p], a[p]) + score(b[p], b[p])) / 2 - score(a[p], b[p]);
        D[i][j] = D[j][i] = sum / SITE_LENGTH;
      }
    }
  })();

  /* ---- MDS ---- */
  function copy(X) { return X.map(function (p) { return [p[0], p[1]]; }); }

  function centred(X) {
    var cx = 0, cy = 0;
    X.forEach(function (p) { cx += p[0]; cy += p[1]; });
    cx /= N; cy /= N;
    return X.map(function (p) { return [p[0] - cx, p[1] - cy]; });
  }

  function stress(X) {
    var num = 0, den = 0;
    for (var i = 0; i < N; i++) {
      for (var j = i + 1; j < N; j++) {
        var d = Math.hypot(X[i][0] - X[j][0], X[i][1] - X[j][1]);
        num += (d - D[i][j]) * (d - D[i][j]);
        den += D[i][j] * D[i][j];
      }
    }
    return Math.sqrt(num / den);
  }

  // one eased SMACOF (Guttman transform) step; returns the largest move
  function relax(X, pinned) {
    var next = [], moved = 0, i, j;
    for (i = 0; i < N; i++) {
      var sx = 0, sy = 0;
      for (j = 0; j < N; j++) {
        if (i === j) continue;
        var dx = X[i][0] - X[j][0], dy = X[i][1] - X[j][1];
        var b = D[i][j] / (Math.hypot(dx, dy) || 1e-9);
        sx += b * dx; sy += b * dy;
      }
      next.push([sx / N, sy / N]);
    }
    for (i = 0; i < N; i++) {
      if (i === pinned) continue;
      var mx = EASE * (next[i][0] - X[i][0]), my = EASE * (next[i][1] - X[i][1]);
      X[i][0] += mx; X[i][1] += my;
      moved = Math.max(moved, Math.abs(mx), Math.abs(my));
    }
    return moved;
  }

  // relax until nothing moves; returns how many steps that took
  function settle(X, tolerance, limit) {
    for (var k = 0; k < limit; k++) if (relax(X, -1) < tolerance) return k;
    return limit;
  }

  // rotation or reflection (as a 2x2 matrix for row vectors) that best maps X onto R
  function align(X, R) {
    X = centred(X); R = centred(R);
    var m00 = 0, m01 = 0, m10 = 0, m11 = 0;
    for (var i = 0; i < N; i++) {
      m00 += X[i][0] * R[i][0]; m01 += X[i][0] * R[i][1];
      m10 += X[i][1] * R[i][0]; m11 += X[i][1] * R[i][1];
    }
    var a = m00 + m11, b = m01 - m10, a2 = m00 - m11, b2 = m01 + m10, th;
    if (Math.hypot(a, b) >= Math.hypot(a2, b2)) {
      th = Math.atan2(b, a);
      return [Math.cos(th), Math.sin(th), -Math.sin(th), Math.cos(th)];
    }
    th = Math.atan2(b2, a2);
    return [Math.cos(th), Math.sin(th), Math.sin(th), -Math.cos(th)];
  }

  function transform(X, q) {
    return centred(X).map(function (p) { return [p[0] * q[0] + p[1] * q[2], p[0] * q[1] + p[1] * q[3]]; });
  }

  // Reference layout: start with the strains in time order along a line, then relax.
  var spread = 0;
  D.forEach(function (row) { row.forEach(function (v) { if (v > spread) spread = v; }); });
  var home = (function () {
    var rnd = K.rng(5);
    var X = strains.map(function (s) { return [(s.season / SEASONS - 0.5) * spread, 0.05 * K.gauss(rnd)]; });
    settle(X, 1e-6, 4000);
    return centred(X);
  })();
  var best = stress(home);
  var X = copy(home);
  var extent = (function () {
    var x = 0, y = 0;
    home.forEach(function (p) { x = Math.max(x, Math.abs(p[0])); y = Math.max(y, Math.abs(p[1])); });
    return { x: x, y: y };
  })();

  // A scrambled start can land in a worse local minimum, or crawl. As mapping
  // software does, try starts until one reaches the best layout in reasonable
  // time, then turn that start so the run ends in a consistent orientation.
  // The relaxation itself is not altered.
  function scrambledStart() {
    var rnd = K.rng(K.seed());
    for (var attempt = 0; attempt < 12; attempt++) {
      var start = home.map(function (p) { return [p[0] + 0.4 * spread * K.gauss(rnd), p[1] + 0.4 * spread * K.gauss(rnd)]; });
      var end = copy(start);
      var steps = settle(end, 1e-4, 520);
      if (steps < 520 && stress(end) <= best + 0.004) return transform(start, align(end, home));
    }
    return copy(home);
  }

  /* ---- chart scaffold: rebuilt when the figure switches between wide and compact ---- */
  var W, H, UNIT, OX, OY, HALF_X, HALF_Y;
  var sx = function (x) { return OX + x * UNIT; };
  var sy = function (y) { return OY - y * UNIT; };
  var layer = null, path = null;
  var readout = fig.querySelector('[data-readout]');
  var dragging = -1;
  var nodes = [];
  // sequential blue: light theme runs light to dark with time, dark theme the other way round
  var RAMP_LIGHT = [[134, 182, 239], [57, 135, 229], [28, 92, 171], [13, 54, 107]];
  var RAMP_DARK = [[24, 79, 149], [57, 135, 229], [134, 182, 239], [205, 226, 251]];

  function seasonColour(season) {
    var ramp = document.documentElement.dataset.theme === 'dark' ? RAMP_DARK : RAMP_LIGHT;
    var p = (season - 1) / (SEASONS - 1) * (ramp.length - 1), i = Math.min(ramp.length - 2, Math.floor(p)), f = p - i;
    var a = ramp[i], b = ramp[i + 1];
    return 'rgb(' + Math.round(a[0] + (b[0] - a[0]) * f) + ',' + Math.round(a[1] + (b[1] - a[1]) * f) + ',' + Math.round(a[2] + (b[2] - a[2]) * f) + ')';
  }

  function paint() {
    nodes.forEach(function (g, i) { g.lastChild.style.fill = seasonColour(strains[i].season); });
  }

  function build(narrow) {
    W = narrow ? 360 : 560; H = narrow ? 300 : 400;
    var margin = narrow ? 26 : 40;
    UNIT = Math.min((W / 2 - margin) / extent.x, (H / 2 - margin) / extent.y);
    OX = W / 2; OY = H / 2;
    HALF_X = (W / 2 - 8) / UNIT; HALF_Y = (H / 2 - 8) / UNIT;
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.textContent = '';
    dragging = -1;
    nodes = [];

    var grid = K.svg('g', { 'class': 'grid' }, svg);
    var nx = Math.floor(W / 2 / UNIT), ny = Math.floor(H / 2 / UNIT);
    for (var gx = -nx; gx <= nx; gx++) K.svg('line', { x1: sx(gx), x2: sx(gx), y1: 0, y2: H }, grid);
    for (var gy = -ny; gy <= ny; gy++) K.svg('line', { x1: 0, x2: W, y1: sy(gy), y2: sy(gy) }, grid);
    K.svg('rect', { 'class': 'frame', x: 0.5, y: 0.5, width: W - 1, height: H - 1, rx: 6 }, svg);
    K.text(svg, W - 10, H - 10, 'MDS 1', { 'text-anchor': 'end' });
    K.text(svg, 10, 18, 'MDS 2');

    path = K.svg('path', { 'class': 'drift' }, svg);
    layer = K.svg('g', null, svg);
    for (var i = 0; i < N; i++) addStrain(i);
    // panel strains sit on top so their markers are never buried
    nodes.forEach(function (g, i) { if (strains[i].panel) layer.appendChild(g); });
    paint();
  }

  function addStrain(i) {
    var s = strains[i];
    var name = 'Season ' + s.season + (s.panel ? ', panel strain' : '');
    var g = K.svg('g', { 'class': 'strain' + (s.panel ? ' is-panel' : ''), tabindex: 0, role: 'img', 'aria-label': name + '. Arrow keys move it.' }, layer);
    K.svg('circle', { 'class': 'halo', r: 12 }, g);
    if (s.panel) K.svg('path', { 'class': 'dot', d: 'M0 -7.5L6.8 4.6H-6.8Z' }, g);
    else K.svg('circle', { 'class': 'dot', r: 4.5 }, g);
    g.style.touchAction = 'none';
    nodes.push(g);

    g.addEventListener('pointerdown', function (e) {
      e.preventDefault();
      dragging = i;
      g.classList.add('is-drag');
      try { g.setPointerCapture(e.pointerId); } catch (err) {}
      K.tip.hide();
      loop.start();
    });
    g.addEventListener('pointermove', function (e) {
      if (dragging === i) {
        var p = K.svgPoint(svg, e);
        X[i][0] = K.clamp((p.x - OX) / UNIT, -HALF_X, HALF_X);
        X[i][1] = K.clamp((OY - p.y) / UNIT, -HALF_Y, HALF_Y);
      } else if (e.pointerType !== 'touch') {
        K.tip.show(e.clientX, e.clientY, 'Season ' + s.season, s.panel ? 'panel strain' : 'simulated strain');
      }
    });
    var drop = function (e) {
      if (dragging !== i) return;
      dragging = -1;
      g.classList.remove('is-drag');
      try { g.releasePointerCapture(e.pointerId); } catch (err) {}
      loop.start();
    };
    g.addEventListener('pointerup', drop);
    g.addEventListener('pointercancel', drop);
    g.addEventListener('pointerleave', function () { K.tip.hide(); });
    g.addEventListener('keydown', function (e) {
      var dx = e.key === 'ArrowLeft' ? -1 : e.key === 'ArrowRight' ? 1 : 0;
      var dy = e.key === 'ArrowDown' ? -1 : e.key === 'ArrowUp' ? 1 : 0;
      if (!dx && !dy) return;
      e.preventDefault();
      X[i][0] = K.clamp(X[i][0] + dx * 0.4, -HALF_X, HALF_X);
      X[i][1] = K.clamp(X[i][1] + dy * 0.4, -HALF_Y, HALF_Y);
      loop.start();
    });
  }

  function render(state) {
    if (!layer) return;
    var d = '';
    for (var i = 0; i < N; i++) {
      var x = sx(X[i][0]).toFixed(1), y = sy(X[i][1]).toFixed(1);
      nodes[i].setAttribute('transform', 'translate(' + x + ' ' + y + ')');
      d += (i ? 'L' : 'M') + x + ' ' + y;
    }
    path.setAttribute('d', d);          // strains are stored in time order, so this traces the drift
    if (readout) readout.textContent = 'stress ' + stress(X).toFixed(3) + ' · ' + state;
  }

  /* ---- loop: fixed 120 steps per second, stops once the map has settled ---- */
  var clock = 0;
  var loop = K.loop(svg, function (dt) {
    clock += dt;
    var moved = Infinity, steps = 0;
    while (clock >= STEP_MS && steps < 4) {
      clock -= STEP_MS; steps++;
      moved = relax(X, dragging);
      if (dragging < 0) {
        // An MDS layout is only defined up to a shift and a turn, so ease the whole
        // map back to the middle of the frame and to its reference orientation.
        var cx = 0, cy = 0, dot = 0, cross = 0, i;
        X.forEach(function (p) { cx += p[0]; cy += p[1]; });
        cx /= N; cy /= N;
        for (i = 0; i < N; i++) {
          var px = X[i][0] - cx, py = X[i][1] - cy;
          dot += px * home[i][0] + py * home[i][1];
          cross += px * home[i][1] - py * home[i][0];
        }
        var turn = Math.atan2(cross, dot) * 0.1, cos = Math.cos(turn), sin = Math.sin(turn);
        for (i = 0; i < N; i++) {
          var qx = X[i][0] - cx * 0.2, qy = X[i][1] - cy * 0.2;
          X[i][0] = qx * cos - qy * sin;
          X[i][1] = qx * sin + qy * cos;
        }
        moved = Math.max(moved, Math.abs(cx) * 0.2, Math.abs(cy) * 0.2, Math.abs(turn) * spread);
      }
    }
    if (clock > 200) clock = 0;
    var done = dragging < 0 && moved < 0.0004;
    render(dragging >= 0 ? 'holding one strain' : done ? 'settled' : 'relaxing');
    if (done) loop.stop();
  });

  function scramble() {
    X = scrambledStart();
    clock = 0;
    render('relaxing');
    loop.start();
  }

  var btn = fig.querySelector('[data-act="scramble"]');
  if (btn) btn.addEventListener('click', scramble);

  K.compact(svg.parentNode, 440, function (narrow) {
    build(narrow);
    render(loop.isRunning() ? 'relaxing' : 'settled');
  });
  K.onTheme(paint);
  if (!K.reduced) K.onceVisible(fig, scramble);
})();
