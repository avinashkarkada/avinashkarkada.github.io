/* Fig: continuous and discontinuous epitopes on HA1.
   The same H3 model as Fig. 1, trimmed to the HA1 head. A continuous epitope is
   one run of sequence. A discontinuous one is assembled from segments that sit far
   apart in the sequence and touch only once the chain is folded, which is why
   epitope mapping has to be done against the structure and not the sequence.
   Distances are measured from the model. */
(function () {
  'use strict';

  var K = window.Kit, data = window.HA_TRACE;
  var canvas = document.getElementById('epitope-canvas');
  var ruler = document.getElementById('epitope-ruler');
  var fig = document.getElementById('fig-epitope');
  if (!K || !data || !canvas || !canvas.getContext || !fig) return;

  var ctx = canvas.getContext('2d');
  var HA1 = data.ha1Length, OFF = data.ha1Offset, LEN = data.seq.length;
  var FIRST = 1 + OFF, LAST = HA1 + OFF;          // HA1 in H3 numbering

  var SITES = [
    { id: 'A', label: 'Continuous', segments: [[140, 146]] },
    { id: 'C', label: 'Discontinuous', segments: [[44, 54], [275, 280]] }
  ];
  var site = SITES[0];

  /* ---- coordinates: HA1 only, all three chains ---- */
  var chains = data.chains.map(function (c) {
    var pts = new Float32Array(HA1 * 3);
    for (var i = 0; i < HA1; i++) {
      pts[i * 3] = c.xyz[i * 3] * data.scale;
      pts[i * 3 + 1] = c.xyz[i * 3 + 1] * data.scale;
      pts[i * 3 + 2] = c.xyz[i * 3 + 2] * data.scale;
    }
    return pts;
  });
  var idxOf = function (h3) { return h3 - FIRST; };
  // both sites sit between residues 44 and 280, so the lower stem of HA1 is left out
  // and the drawing fills the frame instead of trailing off the bottom
  var DRAW_LO = idxOf(36), DRAW_HI = idxOf(292);

  // centre on the trimer axis, and on the midpoint of the HA1 height
  var centre = (function () {
    var sx = 0, sz = 0, lo = Infinity, hi = -Infinity, n = 0;
    chains.forEach(function (p) {
      for (var i = DRAW_LO; i <= DRAW_HI; i++) {
        sx += p[i * 3]; sz += p[i * 3 + 2];
        lo = Math.min(lo, p[i * 3 + 1]); hi = Math.max(hi, p[i * 3 + 1]);
        n++;
      }
    });
    return [sx / n, (lo + hi) / 2, sz / n];
  })();
  // extents used to fit the drawing to the canvas, whatever the rotation
  var extent = (function () {
    var xz = 0, y = 0;
    chains.forEach(function (p) {
      for (var i = DRAW_LO; i <= DRAW_HI; i++) {
        xz = Math.max(xz, Math.hypot(p[i * 3] - centre[0], p[i * 3 + 2] - centre[2]));
        y = Math.max(y, Math.abs(p[i * 3 + 1] - centre[1]));
      }
    });
    return { xz: xz, y: y };
  })();
  var radius = Math.hypot(extent.xz, extent.y);

  function pointsOf(seg) {
    var out = [];
    for (var h = seg[0]; h <= seg[1]; h++) out.push(idxOf(h));
    return out;
  }
  function centroid(seg) {
    var p = chains[0], sx = 0, sy = 0, sz = 0, n = 0;
    pointsOf(seg).forEach(function (i) { sx += p[i * 3]; sy += p[i * 3 + 1]; sz += p[i * 3 + 2]; n++; });
    return [sx / n, sy / n, sz / n];
  }
  function closestApproach(a, b) {
    var p = chains[0], best = Infinity;
    pointsOf(a).forEach(function (i) {
      pointsOf(b).forEach(function (j) {
        best = Math.min(best, Math.hypot(p[i * 3] - p[j * 3], p[i * 3 + 1] - p[j * 3 + 1], p[i * 3 + 2] - p[j * 3 + 2]));
      });
    });
    return best;
  }
  function spread(seg) {
    var p = chains[0], best = 0, pts = pointsOf(seg);
    pts.forEach(function (i) {
      pts.forEach(function (j) {
        best = Math.max(best, Math.hypot(p[i * 3] - p[j * 3], p[i * 3 + 1] - p[j * 3 + 1], p[i * 3 + 2] - p[j * 3 + 2]));
      });
    });
    return best;
  }

  /* ---- view ---- */
  var yaw = 0.4, pitch = -0.35, spin = 0.2, velocity = spin;
  var dragging = false, lastX = 0, lastY = 0, colors = null;
  var CAMERA = 460, LAYERS = 7;
  var view = K.fitCanvas(canvas, function () { draw(); });
  var screenPts = chains.map(function () { return new Float32Array(HA1 * 3); });

  function readColors() {
    colors = { quiet: K.css('--line-2'), s1: K.css('--s1'), s2: K.css('--s2'), page: K.css('--page'), ink: K.css('--ink') };
  }

  function project() {
    var cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
    var scale = Math.min(view.w / (2 * extent.xz), view.h / (2 * extent.y)) * 0.94;
    var ox = view.w / 2, oy = view.h / 2;
    for (var c = 0; c < chains.length; c++) {
      var p = chains[c], s = screenPts[c];
      for (var i = DRAW_LO; i <= DRAW_HI; i++) {
        var x = p[i * 3] - centre[0], y = p[i * 3 + 1] - centre[1], z = p[i * 3 + 2] - centre[2];
        var x1 = x * cy + z * sy, z1 = -x * sy + z * cy;
        var y2 = y * cp - z1 * sp, z2 = y * sp + z1 * cp;
        var f = CAMERA / (CAMERA - z2);
        s[i * 3] = ox + x1 * f * scale;
        s[i * 3 + 1] = oy - y2 * f * scale;
        s[i * 3 + 2] = z2;
      }
    }
  }

  function draw() {
    if (!view || !view.w) return;
    if (!colors) readColors();
    project();
    ctx.setTransform(view.dpr, 0, 0, view.dpr, 0, 0);
    ctx.clearRect(0, 0, view.w, view.h);
    ctx.lineCap = 'round';
    ctx.lineJoin = 'round';

    var depth = radius;
    var binOf = function (z) { return Math.floor(K.clamp((z + depth) / (2 * depth), 0, 0.9999) * LAYERS); };

    for (var layer = 0; layer < LAYERS; layer++) {
      var t = layer / (LAYERS - 1);
      ctx.beginPath();
      var any = false;
      for (var c = 0; c < chains.length; c++) {
        var s = screenPts[c], p = chains[c];
        for (var i = DRAW_LO + 1; i <= DRAW_HI; i++) {
          var dx = p[i * 3] - p[i * 3 - 3], dy = p[i * 3 + 1] - p[i * 3 - 2], dz = p[i * 3 + 2] - p[i * 3 - 1];
          if (dx * dx + dy * dy + dz * dz > 4.4 * 4.4) continue;
          if (binOf((s[i * 3 + 2] + s[i * 3 - 1]) / 2) !== layer) continue;
          ctx.moveTo(s[i * 3 - 3], s[i * 3 - 2]);
          ctx.lineTo(s[i * 3], s[i * 3 + 1]);
          any = true;
        }
      }
      if (!any) continue;
      ctx.globalAlpha = K.lerp(0.15, 0.5, t);
      ctx.lineWidth = K.lerp(0.9, 1.7, t);
      ctx.strokeStyle = colors.quiet;
      ctx.stroke();
    }

    // the selected segments, drawn over the faint trace
    ctx.globalAlpha = 1;
    site.segments.forEach(function (seg, k) {
      var s = screenPts[0];
      ctx.fillStyle = k === 0 ? colors.s1 : colors.s2;
      ctx.strokeStyle = colors.page;
      ctx.lineWidth = 1.6;
      pointsOf(seg).forEach(function (i) {
        ctx.beginPath();
        ctx.arc(s[i * 3], s[i * 3 + 1], 5, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      });
    });
    ctx.globalAlpha = 1;
  }

  /* ---- sequence ruler ---- */
  function drawRuler() {
    if (!ruler) return;
    var W = 560, H = 42, L = 8, R = 8, y = 18, h = 11;
    ruler.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    ruler.textContent = '';
    var xOf = function (h3) { return L + (h3 - FIRST) / (LAST - FIRST) * (W - L - R); };
    K.svg('rect', { 'class': 'ruler-track', x: L, y: y, width: W - L - R, height: h, rx: 3 }, ruler);
    site.segments.forEach(function (seg, k) {
      K.svg('rect', { 'class': 'ruler-seg s' + (k + 1), x: xOf(seg[0]), y: y - 2, width: Math.max(3, xOf(seg[1]) - xOf(seg[0])), height: h + 4, rx: 2 }, ruler);
      K.text(ruler, (xOf(seg[0]) + xOf(seg[1])) / 2, y - 6, seg[0] + '–' + seg[1], { 'text-anchor': 'middle', 'class': 'lbl-strong' });
    });
    if (site.segments.length > 1) {
      var a = site.segments[0], b = site.segments[1];
      var x1 = xOf(a[1]), x2 = xOf(b[0]), mid = (x1 + x2) / 2;
      K.svg('path', { 'class': 'ruler-gap', d: 'M' + x1 + ' ' + (y + h + 7) + 'H' + x2 }, ruler);
      K.text(ruler, mid, y + h + 20, (b[0] - a[1]) + ' residues apart', { 'text-anchor': 'middle' });
    }
    K.text(ruler, L, y + h + 20, 'HA1 ' + FIRST);
    K.text(ruler, W - R, y + h + 20, String(LAST), { 'text-anchor': 'end' });
  }

  function describe() {
    var out = fig.querySelector('[data-readout]');
    if (!out) return;
    out.textContent = '';
    var b = document.createElement('b');
    if (site.segments.length === 1) {
      b.textContent = 'One run of ' + (site.segments[0][1] - site.segments[0][0] + 1) + ' residues';
      out.appendChild(b);
      out.appendChild(document.createTextNode(' · spans ' + spread(site.segments[0]).toFixed(1) + ' Å on the surface'));
    } else {
      var a = site.segments[0], c = site.segments[1];
      b.textContent = (c[0] - a[1]) + ' residues apart in sequence';
      out.appendChild(b);
      out.appendChild(document.createTextNode(' · ' + closestApproach(a, c).toFixed(1) + ' Å apart once folded'));
    }
  }

  /* ---- animation and input ---- */
  var loop = K.loop(canvas, function (dt) {
    var s = dt / 1000;
    if (!dragging) {
      velocity += (spin - velocity) * Math.min(1, s * 3);
      yaw += velocity * s;
    }
    draw();
  });
  if (!K.reduced) loop.start(); else draw();

  canvas.addEventListener('pointerdown', function (e) {
    dragging = true; lastX = e.clientX; lastY = e.clientY; velocity = 0;
    canvas.classList.add('is-drag');
    try { canvas.setPointerCapture(e.pointerId); } catch (err) {}
    if (!loop.isRunning()) loop.start();
  });
  canvas.addEventListener('pointermove', function (e) {
    if (!dragging) return;
    yaw += (e.clientX - lastX) * 0.01;
    pitch = K.clamp(pitch + (e.clientY - lastY) * 0.006, -1.2, 1.2);
    lastX = e.clientX; lastY = e.clientY;
    if (!loop.isRunning()) draw();
  });
  function release(e) {
    if (!dragging) return;
    dragging = false;
    canvas.classList.remove('is-drag');
    try { canvas.releasePointerCapture(e.pointerId); } catch (err) {}
  }
  canvas.addEventListener('pointerup', release);
  canvas.addEventListener('pointercancel', release);
  canvas.tabIndex = 0;
  canvas.addEventListener('keydown', function (e) {
    if (e.key === 'ArrowLeft') yaw -= 0.15;
    else if (e.key === 'ArrowRight') yaw += 0.15;
    else if (e.key === 'ArrowUp') pitch = K.clamp(pitch - 0.1, -1.2, 1.2);
    else if (e.key === 'ArrowDown') pitch = K.clamp(pitch + 0.1, -1.2, 1.2);
    else return;
    e.preventDefault();
    if (!loop.isRunning()) draw();
  });

  var btns = Array.prototype.slice.call(fig.querySelectorAll('[data-epitope]'));
  btns.forEach(function (btn) {
    btn.addEventListener('click', function () {
      site = SITES.filter(function (s) { return s.id === btn.dataset.epitope; })[0] || SITES[0];
      btns.forEach(function (b) { b.setAttribute('aria-pressed', String(b === btn)); });
      drawRuler();
      describe();
      if (!loop.isRunning()) draw();
    });
  });

  K.onTheme(function () { readColors(); draw(); });
  drawRuler();
  describe();
})();
