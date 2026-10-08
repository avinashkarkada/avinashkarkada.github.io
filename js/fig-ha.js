/* Fig. 1: C-alpha trace of a full-length H3 hemagglutinin trimer (a theoretical
   model, see data-ha.js), drawn on a 2D canvas with a hand-rolled rotation and
   perspective projection. Two views: by subunit, or with the classic antigenic
   sites picked out on the head. */
(function () {
  'use strict';

  var K = window.Kit, data = window.HA_TRACE;
  var canvas = document.getElementById('ha-canvas');
  var fig = document.getElementById('fig-ha');
  if (!K || !data || !canvas || !canvas.getContext || !fig) return;

  var ctx = canvas.getContext('2d');
  var readout = fig.querySelector('[data-readout]');
  var toggle = fig.querySelector('[data-act="toggle"]');
  var THREE = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val' };

  // Wiley/Wilson antigenic sites A-E on H3, as residue ranges in mature-HA (H3) numbering
  var SITES = [
    { name: 'A', ranges: [[140, 146]] },
    { name: 'B', ranges: [[155, 160], [187, 196]] },
    { name: 'C', ranges: [[44, 54], [275, 280]] },
    { name: 'D', ranges: [[207, 212]] },
    { name: 'E', ranges: [[62, 65], [83, 83], [260, 264]] }
  ];

  /* ---- flatten the chains into typed arrays ---- */
  var LEN = data.seq.length, HA1 = data.ha1Length, n = LEN * data.chains.length;
  var P = new Float32Array(n * 3);      // model coordinates, angstroms
  var S = new Float32Array(n * 3);      // screen x, y and depth
  var unit = new Uint8Array(n);         // 0 = HA1, 1 = HA2
  var site = new Int8Array(n).fill(-1); // index into SITES, or -1
  var chainOf = [];
  var segA = [], segB = [], segT = [];  // bonded pairs and their position along the chain
  var radiusY = 0, radiusXZ = 0, tmTop = -Infinity, tmBottom = Infinity;

  var siteOfResidue = {};
  SITES.forEach(function (s, k) {
    s.ranges.forEach(function (r) { for (var h3 = r[0]; h3 <= r[1]; h3++) siteOfResidue[h3] = k; });
  });

  (function () {
    var i = 0;
    data.chains.forEach(function (c) {
      for (var j = 0; j < LEN; j++, i++) {
        var x = c.xyz[j * 3] * data.scale, y = c.xyz[j * 3 + 1] * data.scale, z = c.xyz[j * 3 + 2] * data.scale;
        P[i * 3] = x; P[i * 3 + 1] = y; P[i * 3 + 2] = z;
        unit[i] = j < HA1 ? 0 : 1;
        if (j < HA1 && siteOfResidue[j + 1 + data.ha1Offset] != null) site[i] = siteOfResidue[j + 1 + data.ha1Offset];
        chainOf.push(c.id);
        radiusY = Math.max(radiusY, Math.abs(y));
        radiusXZ = Math.max(radiusXZ, Math.sqrt(x * x + z * z));
        if (j + 1 >= data.tm[0] && j + 1 <= data.tm[1]) { tmTop = Math.max(tmTop, y); tmBottom = Math.min(tmBottom, y); }
        if (j > 0) {
          var dx = x - P[i * 3 - 3], dy = y - P[i * 3 - 2], dz = z - P[i * 3 - 1];
          // neighbouring C-alphas sit about 3.8 angstroms apart; anything longer is a gap
          if (dx * dx + dy * dy + dz * dz < 4.4 * 4.4) { segA.push(i - 1); segB.push(i); segT.push(j / LEN); }
        }
      }
    });
  })();
  var MEMBRANE_R = radiusXZ * 1.55;     // radius of the patch of membrane that is drawn

  // residue name in the numbering people use for H3: HA1 8-329, then HA2 1-221
  function label(i) {
    var j = i % LEN;
    var num = j < HA1 ? j + 1 + data.ha1Offset : j + 1 - HA1;
    return (THREE[data.seq.charAt(j)] || data.seq.charAt(j)) + ' ' + num;
  }

  /* ---- view state ---- */
  var yaw = 0.6, pitch = 0.16, spin = 0.22;   // radians, radians, radians per second
  var velocity = spin;
  var dragging = false, hovering = false, lastX = 0, lastY = 0;
  var picked = -1;
  var mode = 'units';                         // 'units' or 'sites'
  var grown = K.reduced ? 1 : 0;              // 0..1, the trace draws itself in once
  var colors = null;
  var LAYERS = 9, CAMERA = 520;
  var view, scale = 1, ox = 0, oy = 0;

  function readColors() {
    colors = {
      units: [K.css('--s1'), K.css('--s2')],
      sites: [K.css('--s1'), K.css('--s2'), K.css('--s3'), K.css('--s4'), K.css('--s5')],
      quiet: K.css('--muted'), line: K.css('--line-2'), fill: K.css('--surface-2'), page: K.css('--page'), text: K.css('--muted')
    };
  }

  function project() {
    var cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
    var pad = 12;
    scale = Math.min((view.w - 2 * pad) / (2 * MEMBRANE_R), (view.h - 2 * pad) / (2 * radiusY));
    ox = view.w / 2; oy = view.h / 2;
    for (var i = 0; i < n; i++) {
      var x = P[i * 3], y = P[i * 3 + 1], z = P[i * 3 + 2];
      var x1 = x * cy + z * sy;
      var z1 = -x * sy + z * cy;
      var y2 = y * cp - z1 * sp;
      var z2 = y * sp + z1 * cp;
      var f = CAMERA / (CAMERA - z2);
      S[i * 3] = ox + x1 * f * scale;
      S[i * 3 + 1] = oy - y2 * f * scale;
      S[i * 3 + 2] = z2;
    }
  }

  // the lipid bilayer around the transmembrane helices, as a tilted slab
  function drawMembrane() {
    var rx = MEMBRANE_R * scale, ry = Math.abs(rx * Math.sin(pitch));
    var top = oy - tmTop * Math.cos(pitch) * scale, bottom = oy - tmBottom * Math.cos(pitch) * scale;
    ctx.beginPath();
    ctx.ellipse(ox, top, rx, ry, 0, Math.PI, 0);          // upper outline
    ctx.lineTo(ox + rx, bottom);
    ctx.ellipse(ox, bottom, rx, ry, 0, 0, Math.PI);       // lower outline
    ctx.closePath();
    ctx.globalAlpha = 0.6;
    ctx.fillStyle = colors.fill;
    ctx.fill();
    ctx.globalAlpha = 1;
    ctx.lineWidth = 1;
    ctx.strokeStyle = colors.line;
    ctx.stroke();
    ctx.beginPath();                                      // rim of whichever face is turned towards the viewer
    if (pitch >= 0) ctx.ellipse(ox, top, rx, ry, 0, 0, Math.PI);
    else ctx.ellipse(ox, bottom, rx, ry, 0, Math.PI, 0);
    ctx.stroke();
    ctx.fillStyle = colors.text;
    ctx.font = '11px "IBM Plex Mono", ui-monospace, monospace';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText('membrane', ox - rx + 12, (top + bottom) / 2 + ry * 0.5);
  }

  function draw() {
    if (!view || !view.w) return;
    if (!colors) readColors();
    project();
    var dpr = view.dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, view.w, view.h);
    ctx.lineCap = 'round';
    ctx.lineJoin = 'round';
    drawMembrane();

    var depth = radiusY * Math.abs(Math.sin(pitch)) + radiusXZ;   // furthest any atom can sit from the screen plane
    var count = segA.length, bySite = mode === 'sites';
    var binOf = function (z) { return Math.floor(K.clamp((z + depth) / (2 * depth), 0, 0.9999) * LAYERS); };

    // paint back to front in depth slices so near strands sit over far ones
    for (var layer = 0; layer < LAYERS; layer++) {
      var t = layer / (LAYERS - 1);
      for (var u = 0; u < 2; u++) {
        ctx.beginPath();
        var any = false;
        for (var s = 0; s < count; s++) {
          if (segT[s] > grown) continue;
          var a = segA[s];
          if (unit[a] !== u) continue;
          var b = segB[s];
          if (binOf((S[a * 3 + 2] + S[b * 3 + 2]) / 2) !== layer) continue;
          ctx.moveTo(S[a * 3], S[a * 3 + 1]);
          ctx.lineTo(S[b * 3], S[b * 3 + 1]);
          any = true;
        }
        if (!any) continue;
        ctx.globalAlpha = bySite ? K.lerp(0.14, 0.5, t) : K.lerp(0.22, 1, t);
        ctx.lineWidth = bySite ? K.lerp(0.9, 1.8, t) : K.lerp(1, 2.5, t);
        ctx.strokeStyle = bySite ? colors.quiet : colors.units[u];
        ctx.stroke();
      }
      if (bySite && grown >= 1) {
        ctx.globalAlpha = K.lerp(0.35, 1, t);
        for (var k = 0; k < SITES.length; k++) {
          ctx.beginPath();
          for (var i = 0; i < n; i++) {
            if (site[i] !== k || binOf(S[i * 3 + 2]) !== layer) continue;
            ctx.moveTo(S[i * 3] + 3, S[i * 3 + 1]);
            ctx.arc(S[i * 3], S[i * 3 + 1], K.lerp(2, 3.6, t), 0, Math.PI * 2);
          }
          ctx.fillStyle = colors.sites[k];
          ctx.fill();
        }
      }
    }
    ctx.globalAlpha = 1;

    if (picked >= 0) {
      ctx.beginPath();
      ctx.arc(S[picked * 3], S[picked * 3 + 1], 6, 0, Math.PI * 2);
      ctx.fillStyle = bySite ? (site[picked] >= 0 ? colors.sites[site[picked]] : colors.quiet) : colors.units[unit[picked]];
      ctx.fill();
      ctx.lineWidth = 2;
      ctx.strokeStyle = colors.page;
      ctx.stroke();
    }
  }

  view = K.fitCanvas(canvas, function () { draw(); });

  function pick(mx, my) {
    var best = -1, bestZ = -Infinity;
    for (var i = 0; i < n; i++) {
      var dx = S[i * 3] - mx, dy = S[i * 3 + 1] - my;
      if (dx * dx + dy * dy < 81 && S[i * 3 + 2] > bestZ) { bestZ = S[i * 3 + 2]; best = i; }
    }
    return best;
  }

  function describe() {
    if (!readout) return;
    readout.textContent = '';
    if (picked < 0) { readout.textContent = 'Drag to rotate'; return; }
    var b = document.createElement('b');
    b.textContent = label(picked);
    readout.appendChild(b);
    var rest = ' · ' + (unit[picked] ? 'HA2' : 'HA1') + ' · chain ' + chainOf[picked];
    if (site[picked] >= 0) rest += ' · site ' + SITES[site[picked]].name;
    readout.appendChild(document.createTextNode(rest));
  }

  /* ---- animation ---- */
  var playing = !K.reduced;
  var loop = K.loop(canvas, function (dt) {
    var s = dt / 1000;
    if (grown < 1) grown = Math.min(1, grown + s / 1.6);
    if (!dragging) {
      var target = playing && !hovering ? spin : 0;
      // ease into the spin, but brake quickly when paused
      velocity += (target - velocity) * Math.min(1, s * (playing ? 3 : 10));
      yaw += velocity * s;
    }
    draw();
    // nothing left to animate: stop until the next interaction
    if (!playing && !dragging && grown >= 1 && Math.abs(velocity) < 0.002) loop.stop();
  });
  loop.start();

  function wake() { if (!loop.isRunning()) loop.start(); }

  if (toggle) {
    K.setToggle(toggle, playing, 'Start rotation', 'Pause rotation');
    toggle.addEventListener('click', function () {
      playing = !playing;
      K.setToggle(toggle, playing, 'Start rotation', 'Pause rotation');
      wake();
    });
  }

  /* ---- colour mode ---- */
  var modeBtns = Array.prototype.slice.call(fig.querySelectorAll('[data-mode]'));
  var legends = Array.prototype.slice.call(fig.querySelectorAll('[data-legend]'));
  modeBtns.forEach(function (btn) {
    btn.addEventListener('click', function () {
      mode = btn.dataset.mode;
      modeBtns.forEach(function (b) { b.setAttribute('aria-pressed', String(b === btn)); });
      legends.forEach(function (l) { l.hidden = l.dataset.legend !== mode; });
      if (loop.isRunning()) return;
      draw();
    });
  });

  /* ---- pointer ---- */
  function local(e) {
    var r = canvas.getBoundingClientRect();
    return { x: e.clientX - r.left, y: e.clientY - r.top };
  }

  canvas.addEventListener('pointerdown', function (e) {
    dragging = true;
    lastX = e.clientX; lastY = e.clientY;
    velocity = 0;
    canvas.classList.add('is-drag');
    try { canvas.setPointerCapture(e.pointerId); } catch (err) {}
    wake();
  });
  canvas.addEventListener('pointermove', function (e) {
    if (dragging) {
      var dx = e.clientX - lastX, dy = e.clientY - lastY;
      lastX = e.clientX; lastY = e.clientY;
      yaw += dx * 0.01;
      pitch = K.clamp(pitch + dy * 0.006, -0.7, 0.7);
      velocity = K.clamp(dx * 0.6, -4, 4);
      if (!loop.isRunning()) draw();
      return;
    }
    if (e.pointerType === 'touch') return;
    hovering = true;
    var p = local(e);
    var hit = pick(p.x, p.y);
    if (hit !== picked) { picked = hit; describe(); if (!loop.isRunning()) draw(); }
    wake();
  });
  function release(e) {
    if (!dragging) return;
    dragging = false;
    canvas.classList.remove('is-drag');
    try { canvas.releasePointerCapture(e.pointerId); } catch (err) {}
    wake();
  }
  canvas.addEventListener('pointerup', release);
  canvas.addEventListener('pointercancel', release);
  canvas.addEventListener('pointerleave', function () {
    hovering = false;
    if (picked >= 0) { picked = -1; describe(); }
    wake();
  });

  // arrow keys rotate when the canvas has focus
  canvas.tabIndex = 0;
  canvas.addEventListener('keydown', function (e) {
    if (e.key === 'ArrowLeft') yaw -= 0.15;
    else if (e.key === 'ArrowRight') yaw += 0.15;
    else if (e.key === 'ArrowUp') pitch = K.clamp(pitch - 0.1, -0.7, 0.7);
    else if (e.key === 'ArrowDown') pitch = K.clamp(pitch + 0.1, -0.7, 0.7);
    else return;
    e.preventDefault();
    if (!loop.isRunning()) draw();
  });

  K.onTheme(function () { readColors(); draw(); });
})();
