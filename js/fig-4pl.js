/* Fig. 4: four-parameter logistic fit to simulated neutralization data.
   Eight three-fold serum dilutions, two replicates each. The curve is refitted
   by Levenberg-Marquardt least squares every time a point moves. */
(function () {
  'use strict';

  var K = window.Kit;
  var svg = document.getElementById('fourpl-plot');
  var fig = document.getElementById('fig-4pl');
  if (!K || !svg || !fig) return;

  /* ---- model and fit ---- */
  // theta = [top, bottom, t50, slope]; t is log10(dilution); response falls as dilution rises
  function model(th, t) { return th[1] + (th[0] - th[1]) / (1 + Math.pow(10, th[3] * (t - th[2]))); }
  var LO = [40, -30, 0.3, 0.15], HI = [140, 60, 5.7, 8];

  function sse(th, ts, ys) {
    var s = 0;
    for (var i = 0; i < ts.length; i++) { var r = ys[i] - model(th, ts[i]); s += r * r; }
    return s;
  }

  // Gaussian elimination with partial pivoting
  function solve(A, b) {
    var n = b.length, M = A.map(function (row, i) { return row.concat(b[i]); }), i, j, k;
    for (i = 0; i < n; i++) {
      var p = i;
      for (k = i + 1; k < n; k++) if (Math.abs(M[k][i]) > Math.abs(M[p][i])) p = k;
      var tmp = M[i]; M[i] = M[p]; M[p] = tmp;
      if (Math.abs(M[i][i]) < 1e-12) return null;
      for (k = i + 1; k < n; k++) {
        var f = M[k][i] / M[i][i];
        for (j = i; j <= n; j++) M[k][j] -= f * M[i][j];
      }
    }
    var x = new Array(n);
    for (i = n - 1; i >= 0; i--) {
      var s = M[i][n];
      for (j = i + 1; j < n; j++) s -= M[i][j] * x[j];
      x[i] = s / M[i][i];
    }
    return x;
  }

  function fit4pl(ts, ys, start) {
    var th = start.slice(), cost = sse(th, ts, ys), damp = 1e-2, n = ts.length, i, j, k;
    for (var it = 0; it < 80; it++) {
      var J = [], r = [];
      for (i = 0; i < n; i++) {
        var f0 = model(th, ts[i]), row = [];
        r.push(ys[i] - f0);
        for (j = 0; j < 4; j++) {
          var h = 1e-5 * Math.max(1, Math.abs(th[j])), t2 = th.slice();
          t2[j] += h;
          row.push((model(t2, ts[i]) - f0) / h);
        }
        J.push(row);
      }
      var A = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], g = [0, 0, 0, 0];
      for (i = 0; i < n; i++) {
        for (j = 0; j < 4; j++) {
          g[j] += J[i][j] * r[i];
          for (k = 0; k < 4; k++) A[j][k] += J[i][j] * J[i][k];
        }
      }
      var improved = false;
      for (var tries = 0; tries < 8 && !improved; tries++) {
        var Ad = A.map(function (row, a) { return row.map(function (v, b) { return a === b ? v * (1 + damp) + 1e-9 : v; }); });
        var d = solve(Ad, g);
        if (!d) { damp *= 4; continue; }
        var cand = th.map(function (v, a) { return K.clamp(v + d[a], LO[a], HI[a]); });
        var c2 = sse(cand, ts, ys);
        if (c2 < cost) {
          var gain = cost - c2;
          th = cand; cost = c2; damp = Math.max(1e-7, damp * 0.3); improved = true;
          if (gain < 1e-9 * (1 + cost)) return { theta: th, sse: cost };
        } else {
          damp *= 4;
        }
      }
      if (!improved) break;
    }
    return { theta: th, sse: cost };
  }

  // log10 dilution at which the curve crosses 50%, or null if it never does
  function halfPoint(th) {
    if (!(th[0] > 50 && th[1] < 50)) return null;
    var t = th[2] + Math.log10((th[0] - th[1]) / (50 - th[1]) - 1) / th[3];
    return isFinite(t) ? t : null;
  }

  // Half-width, in log10 units, of an approximate 95% interval for that crossing:
  // the delta method on the least-squares covariance, residual variance * inverse(J'J).
  function halfWidth(ts, ys, th) {
    var n = ts.length, A = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], g = [], i, j, k;
    var step = function (j) { return 1e-5 * Math.max(1, Math.abs(th[j])); };
    for (i = 0; i < n; i++) {
      var f0 = model(th, ts[i]), row = [];
      for (j = 0; j < 4; j++) { var t2 = th.slice(); t2[j] += step(j); row.push((model(t2, ts[i]) - f0) / step(j)); }
      for (j = 0; j < 4; j++) for (k = 0; k < 4; k++) A[j][k] += row[j] * row[k];
    }
    var h0 = halfPoint(th);
    if (h0 === null) return null;
    for (j = 0; j < 4; j++) {
      var t3 = th.slice(); t3[j] += step(j);
      var h1 = halfPoint(t3);
      if (h1 === null) return null;
      g.push((h1 - h0) / step(j));
    }
    var x = solve(A, g);
    if (!x) return null;
    var variance = sse(th, ts, ys) / (n - 4) * (g[0] * x[0] + g[1] * x[1] + g[2] * x[2] + g[3] * x[3]);
    if (!(variance >= 0) || !isFinite(variance)) return null;
    return 2.179 * Math.sqrt(variance);     // t quantile for 16 points less 4 parameters
  }

  // crude starting point: the range of the data, and the dilution nearest its half-way value
  function guess(ts, ys) {
    var hi = -Infinity, lo = Infinity, i;
    for (i = 0; i < ys.length; i++) { hi = Math.max(hi, ys[i]); lo = Math.min(lo, ys[i]); }
    var half = (hi + lo) / 2, t50 = (ts[0] + ts[ts.length - 1]) / 2, best = Infinity;
    for (i = 0; i < ys.length; i++) {
      if (Math.abs(ys[i] - half) < best) { best = Math.abs(ys[i] - half); t50 = ts[i]; }
    }
    return [K.clamp(hi, LO[0], HI[0]), K.clamp(lo, LO[1], HI[1]), t50, 1];
  }

  /* ---- data ---- */
  var DILUTIONS = [20, 60, 180, 540, 1620, 4860, 14580, 43740];
  var REPS = 2, Y_MIN = -8, Y_MAX = 112;
  var ts = [], ys = [], theta = null, noise = 5, seed = 4104;

  function simulate() {
    var rnd = K.rng(seed);
    var truth = [94 + 6 * rnd(), 1 + 5 * rnd(), K.lerp(Math.log10(300), Math.log10(4000), rnd()), 0.8 + 0.9 * rnd()];
    ts = []; ys = [];
    DILUTIONS.forEach(function (d) {
      for (var rep = 0; rep < REPS; rep++) {
        var t = Math.log10(d);
        ts.push(t);
        ys.push(K.clamp(model(truth, t) + noise * K.gauss(rnd), Y_MIN, Y_MAX));
      }
    });
    theta = null;
  }

  function refit() {
    var cold = fit4pl(ts, ys, guess(ts, ys));
    var best = cold;
    if (theta) {
      var warm = fit4pl(ts, ys, theta);
      if (warm.sse < cold.sse) best = warm;
    }
    theta = best.theta;
  }

  /* ---- chart scaffold: rebuilt when the figure switches between wide and compact ---- */
  var W, H, L, R, T, B;
  var T_MIN = Math.log10(20) - 0.14, T_MAX = Math.log10(43740) + 0.14;
  var xOf = function (t) { return L + (t - T_MIN) / (T_MAX - T_MIN) * (W - L - R); };
  var yOf = function (v) { return B - (v - Y_MIN) / (Y_MAX - Y_MIN) * (B - T); };
  var vOf = function (y) { return Y_MIN + (B - y) / (B - T) * (Y_MAX - Y_MIN); };
  var marker, markH, markV, markCI, curve, markDot, markLabel, pointLayer;
  var points = [];
  var readout = fig.querySelector('[data-readout]');

  function build(narrow) {
    W = narrow ? 360 : 520; H = narrow ? 320 : 360;
    L = narrow ? 34 : 46; R = narrow ? 10 : 16; T = 28; B = H - 56;
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.textContent = '';

    var grid = K.svg('g', { 'class': 'grid' }, svg);
    [0, 25, 50, 75, 100].forEach(function (v) {
      K.svg('line', { x1: L, x2: W - R, y1: yOf(v), y2: yOf(v) }, grid);
      K.text(svg, L - 7, yOf(v) + 4, String(v), { 'text-anchor': 'end' });
    });
    var axis = K.svg('g', { 'class': 'axis' }, svg);
    K.svg('line', { x1: L, x2: W - R, y1: B, y2: B }, axis);
    DILUTIONS.forEach(function (d, i) {
      var x = xOf(Math.log10(d));
      K.svg('line', { x1: x, x2: x, y1: B, y2: B + 5 }, axis);
      if (i % 2 === 0) K.text(svg, x, B + 20, '1:' + K.fmt(d), { 'text-anchor': 'middle' });
    });
    K.text(svg, L, 12, 'neutralization (%)', { 'class': 'lbl' });
    K.text(svg, W - R, H - 8, 'serum dilution, log scale', { 'text-anchor': 'end' });

    marker = K.svg('g', null, svg);
    markH = K.svg('line', { 'class': 'rule' }, marker);
    markV = K.svg('line', { 'class': 'rule' }, marker);
    curve = K.svg('path', { 'class': 'line-s1' }, svg);
    markCI = K.svg('path', { 'class': 'ci' }, marker);
    markDot = K.svg('circle', { 'class': 'nt-dot', r: 4.5 }, marker);
    markLabel = K.text(marker, 0, 0, '', { 'class': 'lbl-strong' });
    pointLayer = K.svg('g', null, svg);
  }

  function label(i) {
    return 'Dilution 1 to ' + K.fmt(DILUTIONS[Math.floor(i / REPS)]) + ', replicate ' + (i % REPS + 1) +
      ', ' + Math.round(ys[i]) + ' percent neutralization. Use the up and down arrow keys to move it.';
  }

  function buildPoints() {
    pointLayer.textContent = '';
    points = ts.map(function (t, i) {
      var g = K.svg('g', { 'class': 'pt', tabindex: 0, role: 'slider', 'aria-orientation': 'vertical', 'aria-valuemin': Y_MIN, 'aria-valuemax': Y_MAX }, pointLayer);
      var jitter = (i % REPS === 0 ? -1 : 1) * 3.2;
      var x = xOf(t) + jitter;
      K.svg('circle', { 'class': 'halo', cx: x, r: 13 }, g);
      K.svg('circle', { 'class': 'dot', cx: x, r: 5 }, g);
      g.style.touchAction = 'none';

      g.addEventListener('pointerdown', function (e) {
        e.preventDefault();
        g.classList.add('is-drag');
        try { g.setPointerCapture(e.pointerId); } catch (err) {}
        K.tip.hide();
      });
      g.addEventListener('pointermove', function (e) {
        if (g.classList.contains('is-drag')) {
          ys[i] = K.clamp(vOf(K.svgPoint(svg, e).y), Y_MIN, Y_MAX);
          update();
        } else if (e.pointerType !== 'touch') {
          K.tip.show(e.clientX, e.clientY, Math.round(ys[i]) + '% neutralization',
            'dilution 1:' + K.fmt(DILUTIONS[Math.floor(i / REPS)]) + ' · replicate ' + (i % REPS + 1));
        }
      });
      var drop = function (e) {
        g.classList.remove('is-drag');
        try { g.releasePointerCapture(e.pointerId); } catch (err) {}
      };
      g.addEventListener('pointerup', drop);
      g.addEventListener('pointercancel', drop);
      g.addEventListener('pointerleave', function () { K.tip.hide(); });
      g.addEventListener('keydown', function (e) {
        var step = e.key === 'ArrowUp' ? 2 : e.key === 'ArrowDown' ? -2 : 0;
        if (!step) return;
        e.preventDefault();
        ys[i] = K.clamp(ys[i] + step, Y_MIN, Y_MAX);
        update();
      });
      return g;
    });
  }

  function update() {
    refit();

    points.forEach(function (g, i) {
      var y = yOf(ys[i]);
      g.firstChild.setAttribute('cy', y);
      g.lastChild.setAttribute('cy', y);
      g.setAttribute('aria-valuenow', Math.round(ys[i]));
      g.setAttribute('aria-label', label(i));
    });

    var d = '';
    for (var s = 0; s <= 90; s++) {
      var t = K.lerp(T_MIN, T_MAX, s / 90);
      d += (s ? 'L' : 'M') + xOf(t).toFixed(1) + ' ' + yOf(K.clamp(model(theta, t), Y_MIN, Y_MAX)).toFixed(1);
    }
    curve.setAttribute('d', d);

    // IC50: the dilution at which the fitted curve crosses 50%
    var t50 = halfPoint(theta);
    if (t50 !== null && !(t50 >= T_MIN && t50 <= T_MAX)) t50 = null;
    var half = t50 === null ? null : halfWidth(ts, ys, theta);
    if (t50 === null) {
      marker.setAttribute('visibility', 'hidden');
    } else {
      var x = xOf(t50), y50 = yOf(50);
      marker.setAttribute('visibility', 'visible');
      markH.setAttribute('x1', L); markH.setAttribute('x2', x); markH.setAttribute('y1', y50); markH.setAttribute('y2', y50);
      markV.setAttribute('x1', x); markV.setAttribute('x2', x); markV.setAttribute('y1', y50); markV.setAttribute('y2', B);
      markDot.setAttribute('cx', x); markDot.setAttribute('cy', y50);
      var right = x < W - R - 120;
      markLabel.setAttribute('x', right ? x + 10 : x - 10);
      markLabel.setAttribute('y', y50 - 10);
      markLabel.setAttribute('text-anchor', right ? 'start' : 'end');
      markLabel.textContent = 'IC50 1:' + K.fmt(Math.pow(10, t50));
      if (half === null) {
        markCI.setAttribute('d', '');
      } else {
        // interval drawn as a whisker along the 50% line
        var x0 = xOf(K.clamp(t50 - half, T_MIN, T_MAX)), x1 = xOf(K.clamp(t50 + half, T_MIN, T_MAX));
        markCI.setAttribute('d', 'M' + x0 + ' ' + (y50 - 5) + 'v10M' + x0 + ' ' + y50 + 'H' + x1 + 'M' + x1 + ' ' + (y50 - 5) + 'v10');
      }
    }

    if (readout) {
      readout.textContent = '';
      var b = document.createElement('b');
      b.textContent = t50 === null ? 'IC50 outside the tested range' : 'IC50 1:' + K.fmt(Math.pow(10, t50));
      readout.appendChild(b);
      var rest = '';
      if (half !== null) rest += ' · 95% CI 1:' + K.fmt(Math.pow(10, t50 - half)) + ' to 1:' + K.fmt(Math.pow(10, t50 + half));
      readout.appendChild(document.createTextNode(rest + ' · slope ' + theta[3].toFixed(2)));
    }
  }

  function regenerate() {
    simulate();
    if (!pointLayer) return;
    buildPoints();
    update();
  }

  var fresh = fig.querySelector('[data-act="new"]');
  if (fresh) fresh.addEventListener('click', function () { seed = K.seed(); regenerate(); });
  var noiseBtns = Array.prototype.slice.call(fig.querySelectorAll('[data-noise]'));
  noiseBtns.forEach(function (btn) {
    btn.addEventListener('click', function () {
      noise = parseFloat(btn.dataset.noise);
      noiseBtns.forEach(function (b) { b.setAttribute('aria-pressed', String(b === btn)); });
      regenerate();
    });
  });

  simulate();
  K.compact(svg.parentNode, 440, function (narrow) {
    build(narrow);
    buildPoints();
    update();
  });
})();
