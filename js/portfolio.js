(() => {
  const root = document.documentElement;
  const body = document.body;
  const themeToggle = document.querySelector('[data-theme-toggle]');
  const menuButton = document.querySelector('[data-menu-toggle]');
  const navLinks = document.querySelector('.nav-links');
  const nav = document.querySelector('.site-nav');

  const storedTheme = localStorage.getItem('portfolio-theme');
  if (storedTheme === 'light') root.dataset.theme = 'light';

  function updateThemeButton() {
    if (!themeToggle) return;
    const isLight = root.dataset.theme === 'light';
    themeToggle.setAttribute('aria-label', isLight ? 'Switch to dark theme' : 'Switch to light theme');
    themeToggle.textContent = isLight ? '☾' : '☼';
  }
  updateThemeButton();

  themeToggle?.addEventListener('click', () => {
    const next = root.dataset.theme === 'light' ? 'dark' : 'light';
    root.dataset.theme = next;
    localStorage.setItem('portfolio-theme', next);
    updateThemeButton();
    renderVolcano();
    loadMolecule(currentPdb);
  });

  menuButton?.addEventListener('click', () => {
    const open = navLinks.classList.toggle('open');
    menuButton.setAttribute('aria-expanded', String(open));
  });

  navLinks?.querySelectorAll('a').forEach(link => link.addEventListener('click', () => {
    navLinks.classList.remove('open');
    menuButton?.setAttribute('aria-expanded', 'false');
  }));

  const sections = [...document.querySelectorAll('main section[id]')];
  const navAnchors = [...document.querySelectorAll('.nav-links a[href^="#"]')];
  function onScroll() {
    nav?.classList.toggle('scrolled', window.scrollY > 24);
    const y = window.scrollY + 130;
    let active = sections[0]?.id;
    sections.forEach(section => { if (section.offsetTop <= y) active = section.id; });
    navAnchors.forEach(a => a.classList.toggle('active', a.getAttribute('href') === `#${active}`));
  }
  window.addEventListener('scroll', onScroll, { passive: true });
  onScroll();

  const reveals = document.querySelectorAll('.reveal');
  if ('IntersectionObserver' in window) {
    const observer = new IntersectionObserver(entries => {
      entries.forEach(entry => {
        if (entry.isIntersecting) {
          entry.target.classList.add('visible');
          observer.unobserve(entry.target);
        }
      });
    }, { threshold: 0.12 });
    reveals.forEach(el => observer.observe(el));
  } else {
    reveals.forEach(el => el.classList.add('visible'));
  }

  // Ambient scientific field.
  const canvas = document.getElementById('field-canvas');
  const ctx = canvas?.getContext('2d');
  let particles = [];
  let raf = 0;
  const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;

  function resizeField() {
    if (!canvas || !ctx) return;
    const rect = canvas.getBoundingClientRect();
    const dpr = Math.min(devicePixelRatio || 1, 2);
    canvas.width = Math.round(rect.width * dpr);
    canvas.height = Math.round(rect.height * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const count = Math.max(34, Math.min(90, Math.round(rect.width / 18)));
    particles = Array.from({ length: count }, () => ({
      x: Math.random() * rect.width,
      y: Math.random() * rect.height,
      vx: (Math.random() - .5) * .16,
      vy: (Math.random() - .5) * .16,
      r: Math.random() * 1.6 + .6
    }));
  }

  function drawField() {
    if (!canvas || !ctx) return;
    const rect = canvas.getBoundingClientRect();
    ctx.clearRect(0, 0, rect.width, rect.height);
    const light = root.dataset.theme === 'light';
    ctx.fillStyle = light ? 'rgba(8,127,98,.38)' : 'rgba(98,230,189,.42)';
    ctx.strokeStyle = light ? 'rgba(8,127,98,.10)' : 'rgba(98,230,189,.10)';
    for (let i = 0; i < particles.length; i++) {
      const p = particles[i];
      if (!reduceMotion) {
        p.x += p.vx; p.y += p.vy;
        if (p.x < 0 || p.x > rect.width) p.vx *= -1;
        if (p.y < 0 || p.y > rect.height) p.vy *= -1;
      }
      ctx.beginPath(); ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2); ctx.fill();
      for (let j = i + 1; j < particles.length; j++) {
        const q = particles[j];
        const dx = p.x - q.x, dy = p.y - q.y;
        const d = Math.hypot(dx, dy);
        if (d < 92) {
          ctx.globalAlpha = 1 - d / 92;
          ctx.beginPath(); ctx.moveTo(p.x, p.y); ctx.lineTo(q.x, q.y); ctx.stroke();
          ctx.globalAlpha = 1;
        }
      }
    }
    if (!reduceMotion) raf = requestAnimationFrame(drawField);
  }
  if (canvas) {
    resizeField(); drawField();
    addEventListener('resize', () => { cancelAnimationFrame(raf); resizeField(); drawField(); }, { passive: true });
  }

  const atlasData = {
    immune: {
      kicker: 'Immune profiling',
      title: 'PhIP-Seq antibody landscapes',
      text: 'Large-scale serological profiling where sparse peptide reactivity is converted into reproducible biological signals across patient cohorts.',
      bullets: ['Cohort-aware QC and statistical modeling', 'Peptide and organism-level reactivity', 'Stability-selection framework for robust signatures'],
      tags: ['PhIP-Seq', 'R', 'edgeR', 'limma', 'GLM', 'ML']
    },
    influenza: {
      kicker: 'Virology + NGS',
      title: 'Influenza neutralization at scale',
      text: 'Computational workflows for a highly multiplexed sequencing-based influenza A neutralization assay, from barcode-aware read processing to titration models.',
      bullets: ['Snakemake + SLURM workflow design', 'Cross-run QC and plate-aware validation', 'Percent-neutralization and titration modeling'],
      tags: ['Influenza', 'NGS', 'Snakemake', 'SLURM', 'Python']
    },
    structure: {
      kicker: 'Structural bioinformatics',
      title: 'From sequence to molecular mechanism',
      text: 'Protein modeling, docking, molecular dynamics, epitope localization, and structure-aware interpretation across viral and therapeutic systems.',
      bullets: ['Molecular dynamics and free-energy analysis', 'Protein/ligand and protein/protein docking', 'Structure-aware epitope mapping'],
      tags: ['GROMACS', 'PyMOL', 'ChimeraX', 'AutoDock', 'FoldDisco']
    },
    software: {
      kicker: 'Research software',
      title: 'Reproducible tools built for real analyses',
      text: 'Small utilities, research applications, and production pipelines designed to make computational biology easier to rerun, inspect, and extend.',
      bullets: ['Python and R tooling', 'HPC and cloud execution', 'Reproducible environments and manifests'],
      tags: ['Python', 'R', 'Bash', 'Git', 'AWS', 'HPC']
    },
    evolution: {
      kicker: 'Viral evolution',
      title: 'Computational virology and adaptation',
      text: 'Sequence and structure-informed studies of influenza adaptation, reassortment, and host-range signals with a focus on biologically interpretable evidence.',
      bullets: ['H5N1 mammalian adaptation signatures', 'HA/NA sequence and structural context', 'Comparative and evolutionary analysis'],
      tags: ['H5N1', 'Influenza', 'Genomics', 'Phylogenetics']
    }
  };

  const atlasPanel = document.querySelector('[data-atlas-panel]');
  const atlasNodes = document.querySelectorAll('[data-atlas]');
  function setAtlas(key) {
    const item = atlasData[key];
    if (!item || !atlasPanel) return;
    atlasNodes.forEach(n => n.classList.toggle('active', n.dataset.atlas === key));
    atlasPanel.innerHTML = `<span class="kicker">${item.kicker}</span><h3>${item.title}</h3><p>${item.text}</p><ul>${item.bullets.map(x => `<li>${x}</li>`).join('')}</ul><div class="atlas-tags">${item.tags.map(x => `<span class="tag">${x}</span>`).join('')}</div>`;
  }
  atlasNodes.forEach(node => node.addEventListener('click', () => setAtlas(node.dataset.atlas)));
  setAtlas('immune');

  // Representative, deterministic scientific plot. Not study data.
  const plotButtons = document.querySelectorAll('[data-plot-view]');
  let plotView = 'volcano';
  const seed = 1337;
  function rngFactory(s) {
    let t = s;
    return () => {
      t += 0x6D2B79F5;
      let r = Math.imul(t ^ t >>> 15, 1 | t);
      r ^= r + Math.imul(r ^ r >>> 7, 61 | r);
      return ((r ^ r >>> 14) >>> 0) / 4294967296;
    };
  }
  function representativeData() {
    const rnd = rngFactory(seed);
    return Array.from({ length: 180 }, (_, i) => {
      const x = (rnd() - .5) * 8;
      const base = Math.abs(x) * .55 + rnd() * 2.7;
      const y = Math.min(7.5, Math.max(.08, base + (rnd() > .92 ? rnd() * 3 : 0)));
      return { x, y, prevalence: 3 + rnd() * 42, label: `peptide_${String(i + 1).padStart(3, '0')}` };
    });
  }
  const demoData = representativeData();

  function renderVolcano() {
    const el = document.getElementById('volcano-plot');
    if (!el || typeof Plotly === 'undefined') return;
    const light = root.dataset.theme === 'light';
    const text = light ? '#38554c' : '#b2c7c0';
    const grid = light ? 'rgba(22,78,64,.11)' : 'rgba(174,231,212,.11)';
    const bg = 'rgba(0,0,0,0)';
    let data;
    let layout;
    if (plotView === 'volcano') {
      const sig = demoData.map(d => d.y > 3 && Math.abs(d.x) > 1.4);
      data = [{
        x: demoData.map(d => d.x), y: demoData.map(d => d.y), text: demoData.map(d => d.label),
        mode: 'markers', type: 'scattergl', hovertemplate: '%{text}<br>log2 FC: %{x:.2f}<br>-log10 p: %{y:.2f}<extra></extra>',
        marker: { size: 7, color: sig.map(v => v ? (light ? '#087f62' : '#62e6bd') : (light ? '#9ab0aa' : '#567069')), opacity: .86 }
      }];
      layout = { xaxis: { title: 'log2 fold change' }, yaxis: { title: '-log10 p-value' } };
    } else {
      const ordered = [...demoData].sort((a, b) => b.prevalence - a.prevalence).slice(0, 36);
      data = [{ x: ordered.map(d => d.prevalence), y: ordered.map((_, i) => `feature ${i + 1}`), type: 'bar', orientation: 'h', hovertemplate: '%{y}<br>prevalence: %{x:.1f}%<extra></extra>', marker: { color: light ? '#087f62' : '#62e6bd', opacity: .72 } }];
      layout = { xaxis: { title: 'Representative prevalence (%)' }, yaxis: { autorange: 'reversed', showticklabels: false } };
    }
    Plotly.react(el, data, {
      ...layout, paper_bgcolor: bg, plot_bgcolor: bg, font: { family: 'Inter, sans-serif', color: text, size: 11 },
      margin: { l: 58, r: 20, t: 24, b: 52 }, showlegend: false,
      xaxis: { ...layout.xaxis, gridcolor: grid, zerolinecolor: grid, tickcolor: grid, linecolor: grid },
      yaxis: { ...layout.yaxis, gridcolor: grid, zerolinecolor: grid, tickcolor: grid, linecolor: grid }
    }, { responsive: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d'] });
  }
  plotButtons.forEach(btn => btn.addEventListener('click', () => {
    plotView = btn.dataset.plotView;
    plotButtons.forEach(b => b.classList.toggle('active', b === btn));
    renderVolcano();
  }));
  renderVolcano();

  // 3D molecular viewer using public PDB structures.
  const pdbButtons = document.querySelectorAll('[data-pdb]');
  let currentPdb = '4K62';
  function loadMolecule(pdb) {
    currentPdb = pdb;
    const el = document.getElementById('molecule-viewer');
    const status = document.querySelector('[data-viewer-status]');
    if (!el || typeof $3Dmol === 'undefined') return;
    pdbButtons.forEach(btn => btn.classList.toggle('active', btn.dataset.pdb === pdb));
    status.textContent = `Loading PDB ${pdb}…`;
    const light = root.dataset.theme === 'light';
    const viewer = $3Dmol.createViewer(el, { backgroundColor: light ? '#edf5f2' : '#050b0a' });
    fetch(`https://files.rcsb.org/download/${pdb}.pdb`)
      .then(r => { if (!r.ok) throw new Error('PDB unavailable'); return r.text(); })
      .then(pdbText => {
        viewer.addModel(pdbText, 'pdb');
        viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
        viewer.zoomTo();
        viewer.render();
        viewer.zoom(0.86, 600);
        status.textContent = `PDB ${pdb} · drag to rotate · scroll to zoom`;
      })
      .catch(() => { status.textContent = `Could not load PDB ${pdb}. The viewer requires network access to RCSB PDB.`; });
  }
  pdbButtons.forEach(btn => btn.addEventListener('click', () => loadMolecule(btn.dataset.pdb)));
  loadMolecule(currentPdb);

  const skillData = {
    omics: ['PhIP-Seq', 'RNA-seq', 'scRNA-seq', 'ChIP-seq', 'WGS', 'edgeR', 'DESeq2', 'limma', 'Seurat', 'Scanpy', 'STAR', 'BWA', 'Minimap2', 'SAMtools', 'BEDTools'],
    programming: ['Python', 'R', 'Bash', 'SQL', 'Git', 'Linux', 'SLURM', 'AWS', 'EC2', 'S3', 'Snakemake', 'Docker', 'Singularity', 'Jupyter'],
    structural: ['GROMACS', 'AMBER', 'PyMOL', 'ChimeraX', 'VMD', 'MODELLER', 'AlphaFold', 'Rosetta', 'HADDOCK', 'AutoDock', 'RDKit', 'Open Babel', 'MM/PBSA'],
    stats: ['scikit-learn', 'XGBoost', 'LASSO', 'SVM', 'randomForest', 'PCA', 'UMAP', 't-SNE', 'Clustering', 'Differential analysis', 'Enrichment analysis', 'Network analysis'],
    reporting: ['ggplot2', 'Plotly', 'R Markdown', 'Quarto', 'Shiny', 'Dash', 'LaTeX', 'Cytoscape', 'BioRender']
  };
  const skillMeta = {
    omics: ['Bioinformatics & omics', 'Sequencing, immune profiling, and statistical analysis.'],
    programming: ['Programming & platforms', 'Reproducible research software from workstation to HPC and cloud.'],
    structural: ['Structural biology', 'Molecular modeling, docking, dynamics, and structure interpretation.'],
    stats: ['Statistics & machine learning', 'Interpretable models and exploratory methods for biological data.'],
    reporting: ['Visualization & reporting', 'Interactive and publication-ready scientific communication.']
  };
  const skillPanel = document.querySelector('[data-skill-panel]');
  const skillTabs = document.querySelectorAll('[data-skill-tab]');
  function setSkill(key) {
    if (!skillPanel) return;
    skillTabs.forEach(t => t.classList.toggle('active', t.dataset.skillTab === key));
    const [title, desc] = skillMeta[key];
    skillPanel.innerHTML = `<h3>${title}</h3><p>${desc}</p><div class="skill-cloud">${skillData[key].map(x => `<span class="skill-chip">${x}</span>`).join('')}</div>`;
  }
  skillTabs.forEach(tab => tab.addEventListener('click', () => setSkill(tab.dataset.skillTab)));
  setSkill('omics');
})();
