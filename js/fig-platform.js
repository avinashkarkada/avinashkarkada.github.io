/* Fig: the PhIP-Seq platform on the Rockfish HPC, from sequencing run to shared
   results. Stage names, order and the hit-calling thresholds follow the lab's
   pipeline documentation. Select a stage to read what it does. */
(function () {
  'use strict';

  var K = window.Kit;
  var svg = document.getElementById('platform-plot');
  var fig = document.getElementById('fig-platform');
  if (!K || !svg || !fig) return;

  var STAGES = [
    { short: 'Align', name: 'Alignment and counts',
      text: 'Demultiplexed FASTQ reads are aligned to the peptide library key, producing a peptide-by-sample count matrix and run QC.' },
    { short: 'Split', name: 'Library splitting',
      text: 'A screen usually pools several libraries, so the counts are split per library against that library’s own peptides and annotation key.' },
    { short: 'Background', name: 'Beads-only background',
      text: 'Each sample is compared against beads-only controls with edgeR, giving enrichment and fold change relative to non-specific binding.' },
    { short: 'Hits', name: 'Enrichment and hit calling',
      text: 'Enrichment p-values are FDR-adjusted, and a peptide is called a hit only when it clears the thresholds set for false-discovery rate, fold change and raw count.' },
    { short: 'Summarise', name: 'Annotation and summarisation',
      text: 'PhIPmake2 joins gene, protein and organism annotations, collapses peptides to proteins, and writes hit-masked matrices. Streamed to disk so memory stays flat whatever the screen size.' },
    { short: 'Aggregate', name: 'Taxon-level aggregation',
      text: 'ARscore and ARscape score whole organisms or peptide groups rather than single peptides, giving aggregate reactivity per taxon.' }
  ];
  var INTAKE = [
    { label: 'Raw reads', text: 'Sequencing runs land in raw-data storage, one folder per run.' },
    { label: 'Metadata', text: 'Alignment keys and sample sheets that tell the pipeline which library it is looking at.' }
  ];
  var OUTPUTS = [
    { label: 'Shared results', text: 'Processed outputs are synced out so collaborators can work from them without an HPC account.' },
    { label: 'Cold storage', text: 'Older runs move to archive storage once they are no longer in active analysis.' }
  ];

  var W, H, narrow = false, nodes = [], selected = 0;
  var info = fig.querySelector('[data-stage-info]');

  function build(isNarrow) {
    narrow = isNarrow;
    svg.textContent = '';
    nodes = [];

    var defs = K.svg('defs', null, svg);
    var marker = K.svg('marker', { id: 'plat-arrow', viewBox: '0 0 8 8', refX: 7, refY: 4, markerWidth: 6, markerHeight: 6, orient: 'auto' }, defs);
    K.svg('path', { d: 'M0 0.5 L7.5 4 L0 7.5 Z', 'class': 'dag-arrow' }, marker);

    var edges = K.svg('g', null, svg);
    var boxes = K.svg('g', null, svg);

    function chip(x, y, w, h, item, cls) {
      var g = K.svg('g', { 'class': 'plat-node ' + cls, tabindex: 0, role: 'button', 'aria-label': item.label || item.name }, boxes);
      K.svg('rect', { x: x, y: y, width: w, height: h, rx: 7 }, g);
      K.text(g, x + w / 2, y + h / 2 + 4, item.short || item.label, { 'text-anchor': 'middle' });
      var pick = function () { select(item, g); };
      g.addEventListener('click', pick);
      g.addEventListener('pointerenter', pick);
      g.addEventListener('focus', pick);
      g.addEventListener('keydown', function (e) { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); pick(); } });
      nodes.push({ g: g, item: item });
      return { x: x, y: y, w: w, h: h };
    }

    function link(a, b, vertical) {
      var d;
      if (!vertical) {
        d = 'M' + (a.x + a.w) + ' ' + (a.y + a.h / 2) + 'H' + (b.x - 3);
      } else {
        // drop out of a, run across, then down into b, so offset boxes still connect
        var ax = a.x + a.w / 2, bx = b.x + b.w / 2, mid = (a.y + a.h + b.y) / 2;
        d = Math.abs(ax - bx) < 1
          ? 'M' + ax + ' ' + (a.y + a.h) + 'V' + (b.y - 3)
          : 'M' + ax + ' ' + (a.y + a.h) + 'V' + mid + 'H' + bx + 'V' + (b.y - 3);
      }
      K.svg('path', { d: d, 'class': 'dag-edge', 'marker-end': 'url(#plat-arrow)' }, edges);
    }

    if (narrow) {
      W = 360;
      var rowH = 30, gap = 14, y = 10, half = (W - 24 - 12) / 2;
      var ins = INTAKE.map(function (it, i) { return chip(12 + i * (half + 12), y, half, rowH, it, 'is-store'); });
      y += rowH + gap;
      var st = STAGES.map(function (s, i) {
        var b = chip(12, y, W - 24, rowH, s, 'is-stage');
        if (i === 0) ins.forEach(function (a) { link(a, b, true); });
        y += rowH + gap;
        return b;
      });
      for (var i = 1; i < st.length; i++) link(st[i - 1], st[i], true);
      var outs = OUTPUTS.map(function (it, j) { return chip(12 + j * (half + 12), y, half, rowH, it, 'is-store'); });
      outs.forEach(function (b) { link(st[st.length - 1], b, true); });
      H = y + rowH + 12;
    } else {
      W = 640;
      var pad = 10, bw = (W - 2 * pad - 5 * 10) / 6, bh = 34;
      var yTop = 10, yMid = 74, yBot = 138;
      var inW = 120;
      var insW = [chip(pad, yTop, inW, 28, INTAKE[0], 'is-store'), chip(pad + inW + 14, yTop, inW, 28, INTAKE[1], 'is-store')];
      var stw = STAGES.map(function (s, i) { return chip(pad + i * (bw + 10), yMid, bw, bh, s, 'is-stage'); });
      insW.forEach(function (a) { link(a, stw[0], true); });
      for (var j = 1; j < stw.length; j++) link(stw[j - 1], stw[j], false);
      var outW = 140, ox = W - pad - outW;
      var outs2 = [chip(ox - outW - 14, yBot, outW, 28, OUTPUTS[0], 'is-store'), chip(ox, yBot, outW, 28, OUTPUTS[1], 'is-store')];
      outs2.forEach(function (b) { link(stw[stw.length - 1], b, true); });
      H = yBot + 28 + 12;
    }
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
  }

  function select(item, g) {
    selected = item;
    nodes.forEach(function (n) { n.g.classList.toggle('is-active', n.g === g); });
    if (!info) return;
    info.textContent = '';
    var name = document.createElement('p');
    name.className = 'stage-name';
    name.textContent = item.name || item.label;
    var desc = document.createElement('p');
    desc.className = 'stage-desc';
    desc.textContent = item.text;
    info.appendChild(name);
    info.appendChild(desc);
  }

  K.compact(svg.parentNode, 480, function (isNarrow) {
    build(isNarrow);
    var first = nodes.filter(function (n) { return n.item === STAGES[0]; })[0];
    if (first) select(STAGES[0], first.g);
  });
})();
