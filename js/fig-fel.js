/* Fig. 6: a free-energy landscape recovered by sampling.
   Walkers follow overdamped Langevin dynamics on a toy 2D potential with three
   basins. Their positions are binned, and the map shows F = -kT ln(P / Pmax). */
(function () {
  'use strict';

  var K = window.Kit;
  var canvas = document.getElementById('fel-canvas');
  var fig = document.getElementById('fig-fel');
  if (!K || !canvas || !canvas.getContext || !fig) return;

  var ctx = canvas.getContext('2d');
  var X0 = -2.4, X1 = 2.4, Y0 = -1.7, Y1 = 1.7;
  var GX = 60, GY = 42;                    // histogram bins
  var WALKERS = 12, SUBSTEPS = 24, DT = 0.004, MAX_SAMPLES = 450000, F_MAX = 6;
  var WELLS = [                            // depth (kT at T = 1), centre, width
    { a: 5.0, x: -1.25, y: -0.35, s: 0.45 },
    { a: 3.3, x: 0.15, y: 0.78, s: 0.40 },
    { a: 4.2, x: 1.35, y: -0.45, s: 0.50 }
  ];
  // sequential blue ramp, light to dark
  var RAMP = [[205, 226, 251], [134, 182, 239], [57, 135, 229], [28, 92, 171], [13, 54, 107]];

  // gradient of the potential: Gaussian wells plus a steep wall near the edges
  function gradient(x, y, out) {
    var gx = 0, gy = 0;
    for (var k = 0; k < WELLS.length; k++) {
      var w = WELLS[k], dx = x - w.x, dy = y - w.y, s2 = w.s * w.s;
      var e = w.a * Math.exp(-(dx * dx + dy * dy) / (2 * s2)) / s2;
      gx += e * dx; gy += e * dy;
    }
    gx += 4 * Math.pow(x / 2.2, 7) / 2.2;
    gy += 4 * Math.pow(y / 1.5, 7) / 1.5;
    out[0] = gx; out[1] = gy;
  }

  /* ---- state ---- */
  var temp = 1, samples = 0;
  var hist = new Float32Array(GX * GY), smooth = new Float32Array(GX * GY);
  var wx = new Float32Array(WALKERS), wy = new Float32Array(WALKERS);
  var rnd = K.rng(K.seed()), g = [0, 0];
  var off = document.createElement('canvas');
  off.width = GX; off.height = GY;
  var offCtx = off.getContext('2d');
  var img = offCtx.createImageData(GX, GY);
  var readout = fig.querySelector('[data-readout]');
  var toggle = fig.querySelector('[data-act="toggle"]');
  var colors = null;

  function reset() {
    hist.fill(0);
    samples = 0;
    // every walker starts in the deepest basin, so the others have to be found
    for (var i = 0; i < WALKERS; i++) {
      wx[i] = WELLS[0].x + 0.15 * K.gauss(rnd);
      wy[i] = WELLS[0].y + 0.15 * K.gauss(rnd);
    }
  }

  function advance() {
    var amp = Math.sqrt(2 * temp * DT);
    for (var i = 0; i < WALKERS; i++) {
      var x = wx[i], y = wy[i];
      for (var s = 0; s < SUBSTEPS; s++) {
        gradient(x, y, g);
        x += -g[0] * DT + amp * K.gauss(rnd);
        y += -g[1] * DT + amp * K.gauss(rnd);
        x = K.clamp(x, X0 + 1e-3, X1 - 1e-3);
        y = K.clamp(y, Y0 + 1e-3, Y1 - 1e-3);
        var cx = Math.floor((x - X0) / (X1 - X0) * GX), cy = Math.floor((y - Y0) / (Y1 - Y0) * GY);
        hist[cy * GX + cx]++;
      }
      wx[i] = x; wy[i] = y;
    }
    samples += WALKERS * SUBSTEPS;
  }

  /* ---- drawing ---- */
  function readColors() {
    colors = {
      dark: document.documentElement.dataset.theme === 'dark',
      ink: K.css('--ink'), muted: K.css('--muted'), line: K.css('--line-2'), surface: K.css('--surface')
    };
  }

  function rampAt(t) {
    var p = K.clamp(t, 0, 1) * (RAMP.length - 1), i = Math.min(RAMP.length - 2, Math.floor(p)), f = p - i;
    var a = RAMP[i], b = RAMP[i + 1];
    return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
  }

  // free energy in kT for a (smoothed) bin count, or -1 where nothing was sampled
  function energy(count, peak) { return count > 0 ? Math.min(F_MAX, -Math.log(count / peak)) : -1; }

  function blur() {
    var peak = 0;
    for (var y = 0; y < GY; y++) {
      for (var x = 0; x < GX; x++) {
        var sum = 0, wsum = 0;
        for (var dy = -1; dy <= 1; dy++) {
          for (var dx = -1; dx <= 1; dx++) {
            var xx = x + dx, yy = y + dy;
            if (xx < 0 || yy < 0 || xx >= GX || yy >= GY) continue;
            var w = (dx === 0 ? 2 : 1) * (dy === 0 ? 2 : 1);
            sum += w * hist[yy * GX + xx]; wsum += w;
          }
        }
        var v = sum / wsum;
        smooth[y * GX + x] = v;
        if (v > peak) peak = v;
      }
    }
    return peak;
  }

  var view = K.fitCanvas(canvas, function () { draw(); });
  var PAD = { l: 30, r: 6, t: 6, b: 24 };

  function plotRect() { return { x: PAD.l, y: PAD.t, w: view.w - PAD.l - PAD.r, h: view.h - PAD.t - PAD.b }; }

  function draw() {
    if (!view || !view.w) return;
    if (!colors) readColors();
    var dpr = view.dpr, r = plotRect();
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, view.w, view.h);

    var peak = blur(), data = img.data;
    for (var y = 0; y < GY; y++) {
      for (var x = 0; x < GX; x++) {
        var c = smooth[(GY - 1 - y) * GX + x], k = (y * GX + x) * 4;   // flip: +y is up
        var f = energy(c, peak);
        if (f < 0) { data[k + 3] = 0; continue; }
        // deepest basin sits furthest from the surface colour in either theme
        var depth = 1 - f / F_MAX;
        var rgb = rampAt(colors.dark ? 1 - depth : depth);
        data[k] = rgb[0]; data[k + 1] = rgb[1]; data[k + 2] = rgb[2];
        data[k + 3] = 255 * Math.min(1, c / 1.5);
      }
    }
    offCtx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(off, r.x, r.y, r.w, r.h);

    ctx.strokeStyle = colors.line;
    ctx.lineWidth = 1;
    ctx.strokeRect(r.x + 0.5, r.y + 0.5, r.w - 1, r.h - 1);

    ctx.fillStyle = colors.muted;
    ctx.font = '11px "IBM Plex Mono", ui-monospace, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'alphabetic';
    ctx.fillText('PC1', r.x + r.w, view.h - 6);
    ctx.save();
    ctx.translate(12, r.y);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('PC2', 0, 4);
    ctx.restore();

    for (var i = 0; i < WALKERS; i++) {
      var px = r.x + (wx[i] - X0) / (X1 - X0) * r.w, py = r.y + (1 - (wy[i] - Y0) / (Y1 - Y0)) * r.h;
      ctx.beginPath();
      ctx.arc(px, py, 3.2, 0, Math.PI * 2);
      ctx.fillStyle = colors.ink;
      ctx.fill();
      ctx.lineWidth = 1.5;
      ctx.strokeStyle = colors.surface;
      ctx.stroke();
    }

    if (readout) readout.textContent = K.fmt(samples) + ' samples · ' + WALKERS + ' walkers';
  }

  /* ---- loop ---- */
  var loop = K.loop(canvas, function () {
    advance();
    draw();
    if (samples >= MAX_SAMPLES) { loop.stop(); label(); }
  });

  function label() {
    K.setToggle(toggle, loop.isRunning(), samples >= MAX_SAMPLES ? 'Run again' : samples ? 'Resume' : 'Run', 'Pause');
  }
  function start() {
    if (samples >= MAX_SAMPLES) reset();
    loop.start();
    label();
  }

  if (toggle) toggle.addEventListener('click', function () {
    if (loop.isRunning()) { loop.stop(); label(); } else start();
  });
  var resetBtn = fig.querySelector('[data-act="reset"]');
  if (resetBtn) resetBtn.addEventListener('click', function () { reset(); draw(); label(); });
  var tempBtns = Array.prototype.slice.call(fig.querySelectorAll('[data-temp]'));
  tempBtns.forEach(function (btn) {
    btn.addEventListener('click', function () {
      temp = parseFloat(btn.dataset.temp);
      tempBtns.forEach(function (b) { b.setAttribute('aria-pressed', String(b === btn)); });
      reset();
      draw();
      label();
    });
  });

  canvas.addEventListener('pointermove', function (e) {
    if (e.pointerType === 'touch') return;
    var box = canvas.getBoundingClientRect(), r = plotRect();
    var u = (e.clientX - box.left - r.x) / r.w, v = 1 - (e.clientY - box.top - r.y) / r.h;
    if (u < 0 || u >= 1 || v < 0 || v >= 1) { K.tip.hide(); return; }
    var c = smooth[Math.floor(v * GY) * GX + Math.floor(u * GX)], peak = 0;
    for (var i = 0; i < smooth.length; i++) if (smooth[i] > peak) peak = smooth[i];
    var f = energy(c, peak);
    K.tip.show(e.clientX, e.clientY,
      f < 0 ? 'Not sampled yet' : f >= F_MAX ? 'over ' + F_MAX + ' kT' : f.toFixed(1) + ' kT',
      'PC1 ' + K.lerp(X0, X1, u).toFixed(1) + ' · PC2 ' + K.lerp(Y0, Y1, v).toFixed(1));
  });
  canvas.addEventListener('pointerleave', function () { K.tip.hide(); });

  K.onTheme(function () { readColors(); draw(); });

  reset();
  if (K.reduced) {
    // no autoplay: sample quietly up front and show the finished landscape
    while (samples < MAX_SAMPLES) advance();
    draw();
    label();
  } else {
    loop.start();
    label();
  }
})();
