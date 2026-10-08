/* Fig. 3: rule graph of the InFlux Snakemake workflow.
   Rule names, descriptions and dependencies come from the public repository
   (workflow/Snakefile and README); edges are reduced to direct dependencies. */
(function () {
  'use strict';

  var K = window.Kit;
  var svg = document.getElementById('dag-plot');
  var fig = document.getElementById('fig-dag');
  if (!K || !svg || !fig) return;

  var SAMPLES = 8;
  var CX = 180, LEFT = 94, RIGHT = 266, WIDE = 216, NARROW = 160, NODE_H = 30, ROW = 50, TOP = 10;

  var RULES = [
    { id: 'sanitize_barcodes', row: 0, x: CX, w: WIDE,
      text: 'Reads the barcode reference file, removes BOM artifacts and formatting issues, converts sequences to uppercase, and writes barcodes_clean.fa.' },
    { id: 'anchor_barcodes', row: 1, x: CX, w: WIDE,
      text: 'Converts cleaned barcode sequences into a 5′-anchored FASTA for Cutadapt.' },
    { id: 'cutadapt_e', row: 2, x: CX, w: WIDE, perSample: true,
      text: 'Runs Cutadapt per sample for each configured error-rate label. Produces JSON reports and .done sentinel files.' },
    { id: 'qc_per_sample_eps', row: 3, x: LEFT, w: NARROW, perSample: true,
      text: 'Computes per-sample QC statistics using FASTQ read counts and Cutadapt JSON output.' },
    { id: 'aggregate_counts_e', row: 3, x: RIGHT, w: NARROW,
      text: 'Aggregates barcode counts per sample into a single counts matrix.' },
    { id: 'qc_summary_all', row: 4, x: LEFT, w: NARROW,
      text: 'Combines all per-sample QC files into a single summary TSV.' },
    { id: 'counts_to_annotated_e', row: 5, x: CX, w: WIDE,
      text: 'Merges the count matrix into the annotation sheet. Adds QC metrics and spike_in_fraction.' },
    { id: 'annotated_to_me_e', row: 6, x: CX, w: WIDE,
      text: 'Converts raw counts to molecular equivalence values using the configured spike-in columns.' },
    { id: 'titration_fit_platewise_e', row: 7, x: CX, w: WIDE,
      text: 'Performs plate-wise virus titration line fitting. Produces CSV summaries and PDF visualizations.' },
    { id: 'percent_nt_by_fraction_e', row: 8, x: CX, w: WIDE,
      text: 'Computes percent neutralization values from the line-fit outputs.' },
    { id: 'nt_fourpl_fit_e', row: 9, x: CX, w: WIDE,
      text: 'Fits 4-parameter logistic neutralization curves. Produces final CSV and PDF outputs for pmAb and serum.' }
  ];
  var EDGES = [
    ['sanitize_barcodes', 'anchor_barcodes'],
    ['anchor_barcodes', 'cutadapt_e'],
    ['cutadapt_e', 'qc_per_sample_eps'],
    ['cutadapt_e', 'aggregate_counts_e'],
    ['qc_per_sample_eps', 'qc_summary_all'],
    ['qc_summary_all', 'counts_to_annotated_e'],
    ['aggregate_counts_e', 'counts_to_annotated_e'],
    ['counts_to_annotated_e', 'annotated_to_me_e'],
    ['annotated_to_me_e', 'titration_fit_platewise_e'],
    ['titration_fit_platewise_e', 'percent_nt_by_fraction_e'],
    ['percent_nt_by_fraction_e', 'nt_fourpl_fit_e']
  ];
  // replay order; rules in the same stage run side by side
  var STAGES = [
    ['sanitize_barcodes'], ['anchor_barcodes'], ['cutadapt_e'],
    ['qc_per_sample_eps', 'aggregate_counts_e'], ['qc_summary_all'],
    ['counts_to_annotated_e'], ['annotated_to_me_e'], ['titration_fit_platewise_e'],
    ['percent_nt_by_fraction_e'], ['nt_fourpl_fit_e']
  ];

  var byId = {};
  RULES.forEach(function (r) { r.y = TOP + r.row * ROW; byId[r.id] = r; });
  var totalJobs = RULES.reduce(function (sum, r) { return sum + (r.perSample ? SAMPLES : 1); }, 0);

  /* ---- draw ---- */
  var defs = K.svg('defs', null, svg);
  var marker = K.svg('marker', { id: 'dag-arrowhead', viewBox: '0 0 8 8', refX: 7, refY: 4, markerWidth: 7, markerHeight: 7, orient: 'auto' }, defs);
  K.svg('path', { d: 'M0 0.5 L7.5 4 L0 7.5 Z', 'class': 'dag-arrow' }, marker);

  var edgeLayer = K.svg('g', null, svg);
  var edgeEls = EDGES.map(function (e) {
    var a = byId[e[0]], b = byId[e[1]];
    var x1 = a.x, y1 = a.y + NODE_H, x2 = b.x, y2 = b.y - 2, mid = (y1 + y2) / 2;
    var d = 'M' + x1 + ' ' + y1 + ' C' + x1 + ' ' + mid + ' ' + x2 + ' ' + mid + ' ' + x2 + ' ' + y2;
    return { from: e[0], to: e[1], el: K.svg('path', { d: d, 'class': 'dag-edge', 'marker-end': 'url(#dag-arrowhead)' }, edgeLayer) };
  });

  var nodeLayer = K.svg('g', null, svg);
  RULES.forEach(function (r) {
    var g = K.svg('g', { 'class': 'dag-node', tabindex: 0, role: 'button', 'aria-label': 'Rule ' + r.id }, nodeLayer);
    K.svg('rect', { x: r.x - r.w / 2, y: r.y, width: r.w, height: NODE_H, rx: 7 }, g);
    K.text(g, r.x, r.y + 19.5, r.id, { 'text-anchor': 'middle' });
    if (r.perSample && r.w === WIDE) {
      r.count = K.text(g, r.x + r.w / 2 - 8, r.y + 19.5, '×' + SAMPLES, { 'text-anchor': 'end', 'class': 'count' });
    }
    r.el = g;
    var choose = function () { select(r.id); };
    g.addEventListener('click', choose);
    g.addEventListener('pointerenter', function () { if (!running) choose(); });
    g.addEventListener('focus', choose);
    g.addEventListener('keydown', function (e) {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); choose(); }
    });
  });

  /* ---- details panel ---- */
  var info = document.getElementById('dag-info');
  var log = document.getElementById('dag-log');
  var toggle = fig.querySelector('[data-act="toggle"]');
  var running = false;
  var timers = [];

  function select(id) {
    var r = byId[id];
    RULES.forEach(function (o) { o.el.classList.toggle('is-active', o === r); });
    if (!info) return;
    info.textContent = '';
    var name = document.createElement('p');
    name.className = 'dag-rule';
    name.textContent = 'rule ' + r.id;
    var desc = document.createElement('p');
    desc.className = 'dag-desc';
    desc.textContent = r.text;
    info.appendChild(name);
    info.appendChild(desc);
  }

  function write(line) {
    if (!log) return;
    log.textContent += (log.textContent ? '\n' : '') + line;
    log.scrollTop = log.scrollHeight;
  }

  function clearState() {
    timers.forEach(clearTimeout);
    timers = [];
    RULES.forEach(function (r) { r.el.classList.remove('is-running', 'is-done'); });
    edgeEls.forEach(function (e) { e.el.classList.remove('is-done'); });
    if (log) log.textContent = '';
  }

  function later(fn, ms) { timers.push(setTimeout(fn, ms)); }

  function replay() {
    clearState();
    running = true;
    K.setToggle(toggle, true, 'Replay run', 'Stop');
    var done = 0, t = 250;
    write('Building DAG of jobs... ' + totalJobs + ' jobs, ' + SAMPLES + ' samples');

    STAGES.forEach(function (stage) {
      var longest = 0;
      stage.forEach(function (id, k) {
        var r = byId[id], jobs = r.perSample ? SAMPLES : 1, each = r.perSample ? 110 : 520;
        later(function () {
          r.el.classList.add('is-running');
          if (k === 0) select(id);
        }, t);
        later(function () {
          r.el.classList.remove('is-running');
          r.el.classList.add('is-done');
          edgeEls.forEach(function (e) { if (e.from === id) e.el.classList.add('is-done'); });
          done += jobs;
          var pct = Math.round(done / totalJobs * 100);
          write('[' + (pct < 10 ? '  ' : pct < 100 ? ' ' : '') + pct + '%] ' + id + (jobs > 1 ? '  ×' + jobs : ''));
        }, t + jobs * each);
        longest = Math.max(longest, jobs * each);
      });
      t += longest + 120;
    });

    later(function () {
      write(totalJobs + ' of ' + totalJobs + ' steps (100%) done');
      running = false;
      K.setToggle(toggle, false, 'Replay run', 'Stop');
    }, t);
  }

  function stop() {
    timers.forEach(clearTimeout);
    timers = [];
    running = false;
    RULES.forEach(function (r) { r.el.classList.remove('is-running'); });
    K.setToggle(toggle, false, 'Replay run', 'Stop');
  }

  if (toggle) toggle.addEventListener('click', function () { if (running) stop(); else replay(); });

  if (!K.reduced) K.onceVisible(fig, replay);
})();
