(function () {
  "use strict";

  const TYPE_COLORS = ["#078f82", "#347db8", "#d7782f", "#7758b6", "#24966e", "#c94f70", "#9a741c", "#377e8b", "#674d98", "#435fa5", "#ad4e87", "#c85236", "#60718f"];
  const STATUS_COLORS = ["#8b9491", "#0b8b6f", "#359a84", "#c98019", "#4d74b8", "#c84961", "#a34286"];
  const TAXON_COLORS = { 1: "#208aaa", 2: "#c97828", 4: "#7558b5", 8: "#c84955", mixed: "#2f9661", other: "#818b89" };
  const STATUS_SHORT = ["Unresolved", "Direct · strict", "Direct · permissive", "RiboGRAM-88", "Cluster consensus", "Novel candidate", "Evidence conflict"];
  const STATUS_GROUPS = {
    direct: new Set([1, 2, 6]), model: new Set([3]), propagated: new Set([4]),
    novel: new Set([5]), unresolved: new Set([0])
  };
  const clamp = (value, low, high) => Math.max(low, Math.min(high, value));
  const escapeHTML = (value) => String(value).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;").replaceAll('"', "&quot;").replaceAll("'", "&#39;");
  const compact = (value) => new Intl.NumberFormat("en", { notation: "compact", maximumFractionDigits: 1 }).format(value || 0);

  function hexRGB(hex) {
    const value = Number.parseInt(hex.slice(1), 16);
    return [(value >> 16) & 255, (value >> 8) & 255, value & 255];
  }

  function hashColor(value) {
    let hash = 2166136261;
    const text = String(value);
    for (let i = 0; i < text.length; i += 1) { hash ^= text.charCodeAt(i); hash = Math.imul(hash, 16777619); }
    const hue = ((hash >>> 0) % 300) + 20;
    const saturation = 45 + ((hash >>> 9) % 24);
    const lightness = 53 + ((hash >>> 17) % 14);
    const c = (1 - Math.abs(2 * lightness / 100 - 1)) * saturation / 100;
    const x = c * (1 - Math.abs((hue / 60) % 2 - 1));
    const m = lightness / 100 - c / 2;
    let rgb;
    if (hue < 60) rgb = [c, x, 0]; else if (hue < 120) rgb = [x, c, 0];
    else if (hue < 180) rgb = [0, c, x]; else if (hue < 240) rgb = [0, x, c];
    else if (hue < 300) rgb = [x, 0, c]; else rgb = [c, 0, x];
    return rgb.map((part) => Math.round((part + m) * 255));
  }

  async function loadCompressedJSON(url) {
    const response = await fetch(url, { cache: "force-cache" });
    if (!response.ok) throw new Error(`Could not load ${url} (${response.status})`);
    const buffer = await response.arrayBuffer();
    let text;
    if (new URL(url, location.href).pathname.endsWith(".gz")) {
      if (!("DecompressionStream" in window)) throw new Error("This browser cannot decompress the atlas data package.");
      const stream = new Blob([buffer]).stream().pipeThrough(new DecompressionStream("gzip"));
      text = await new Response(stream).text();
    } else {
      text = new TextDecoder().decode(buffer);
    }
    return JSON.parse(text);
  }

  function createShader(gl, type, source) {
    const shader = gl.createShader(type); gl.shaderSource(shader, source); gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(shader));
    return shader;
  }

  class UniverseAtlas {
    constructor(root) {
      this.root = root;
      this.stage = root.querySelector("[data-map-stage]");
      this.canvas = root.querySelector("[data-map-canvas]");
      this.tooltip = root.querySelector("[data-tooltip]");
      this.loading = root.querySelector("[data-loading]");
      this.controls = root.querySelector("[data-controls]");
      this.intro = root.querySelector("[data-intro]");
      this.search = root.querySelector("[data-search]");
      this.colorMode = root.querySelector("[data-color-mode]");
      this.statusFilter = root.querySelector("[data-status-filter]");
      this.legend = root.querySelector("[data-legend]");
      this.labelsLayer = root.querySelector("[data-cluster-labels]");
      this.inspector = root.querySelector("[data-inspector]");
      this.backdrop = root.querySelector("[data-inspector-backdrop]");
      this.fatal = root.querySelector("[data-fatal]");
      this.zoom = .88; this.panX = 0; this.panY = 0;
      this.drag = null; this.hovered = -1; this.selected = -1; this.query = "";
      this.filter = "all"; this.mode = "type13"; this.detailCache = new Map();
      this.gridSize = 220; this.grid = new Map(); this.queryMatches = null; this.labels = [];
      this.entered = false; this.resizeFrame = null; this.pointerFrame = null;
      this.bindUI();
    }

    async start() {
      try {
        this.data = await loadCompressedJSON(this.root.dataset.source);
        this.points = this.data.points; this.families = this.data.dictionaries.families;
        this.count = this.points.id.length;
        this.positions = new Float32Array(this.count * 2);
        for (let i = 0; i < this.count; i += 1) { this.positions[i * 2] = this.points.x[i]; this.positions[i * 2 + 1] = this.points.y[i]; }
        this.bounds = this.computeBounds();
        this.initGL(); this.buildGrid(); this.buildLabels(); this.updateStats(); this.updateColors();
        this.resize(); this.draw();
        this.loading.classList.add("is-hidden");
        this.restoreHash();
      } catch (error) {
        console.error(error); this.fatal.hidden = false; this.fatal.textContent = `Atlas unavailable: ${error.message}`;
        this.loading.classList.add("is-hidden");
      }
    }

    bindUI() {
      rootEvent(this.root, "[data-enter-map]", "click", () => this.enter());
      rootEvent(this.root, "[data-reset]", "click", () => { this.zoom = .88; this.panX = 0; this.panY = 0; this.draw(); });
      rootEvent(this.root, "[data-inspector-close]", "click", () => this.closeInspector());
      this.backdrop.addEventListener("click", () => this.closeInspector());
      const dialog = this.root.querySelector("[data-method-dialog]");
      this.root.querySelectorAll("[data-method-open]").forEach((element) => element.addEventListener("click", () => dialog.showModal()));
      rootEvent(this.root, "[data-method-close]", "click", () => dialog.close());
      this.colorMode.addEventListener("change", () => { this.mode = this.colorMode.value; this.updateColors(); });
      this.statusFilter.addEventListener("change", () => { this.filter = this.statusFilter.value; this.updateColors(); });
      let searchTimer;
      this.search.addEventListener("input", () => { clearTimeout(searchTimer); searchTimer = setTimeout(() => { this.query = this.search.value.trim().toLowerCase(); this.updateQuery(); }, 120); });
      this.search.addEventListener("keydown", (event) => {
        if (event.key !== "Enter") return; event.preventDefault(); clearTimeout(searchTimer);
        this.query = this.search.value.trim().toLowerCase(); this.updateQuery();
        if (this.queryMatches?.length) this.select(this.queryMatches[0]);
      });
      window.addEventListener("keydown", (event) => {
        if (event.key === "/" && document.activeElement !== this.search) { event.preventDefault(); this.enter(); this.search.focus(); }
        if (event.key === "Escape") this.closeInspector();
      });
      window.addEventListener("resize", () => { cancelAnimationFrame(this.resizeFrame); this.resizeFrame = requestAnimationFrame(() => { this.resize(); this.draw(); }); });
      this.canvas.addEventListener("wheel", (event) => this.wheel(event), { passive: false });
      this.canvas.addEventListener("pointerdown", (event) => this.pointerDown(event));
      window.addEventListener("pointermove", (event) => this.pointerMove(event));
      window.addEventListener("pointerup", (event) => this.pointerUp(event));
      this.canvas.addEventListener("click", (event) => { if (!this.drag?.moved) { const hit = this.nearest(event.clientX, event.clientY, 16); if (hit >= 0) this.select(hit); } });
    }

    enter() { this.entered = true; this.intro.classList.add("is-dismissed"); this.controls.classList.add("is-visible"); }

    computeBounds() {
      let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
      for (let i = 0; i < this.count; i += 1) { const x = this.positions[i * 2], y = this.positions[i * 2 + 1]; minX = Math.min(minX, x); maxX = Math.max(maxX, x); minY = Math.min(minY, y); maxY = Math.max(maxY, y); }
      return { minX, maxX, minY, maxY };
    }

    initGL() {
      const gl = this.canvas.getContext("webgl", { alpha: true, antialias: true, premultipliedAlpha: false });
      if (!gl) throw new Error("WebGL is required to draw this atlas."); this.gl = gl;
      const vertex = createShader(gl, gl.VERTEX_SHADER, `
        attribute vec2 a_position; attribute vec4 a_color; uniform vec2 u_pan;
        uniform float u_zoom; uniform float u_aspect; uniform float u_size;
        varying vec4 v_color; void main(){ vec2 p=(a_position+u_pan)*u_zoom;
        gl_Position=vec4(p.x/u_aspect,p.y,0.,1.); gl_PointSize=u_size; v_color=a_color; }`);
      const fragment = createShader(gl, gl.FRAGMENT_SHADER, `
        precision mediump float; varying vec4 v_color; void main(){ vec2 q=gl_PointCoord-.5;
        float d=dot(q,q); if(d>.25) discard; float edge=smoothstep(.25,.12,d);
        gl_FragColor=vec4(v_color.rgb,v_color.a*edge); }`);
      const program = gl.createProgram(); gl.attachShader(program, vertex); gl.attachShader(program, fragment); gl.linkProgram(program);
      if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
      this.program = program; gl.useProgram(program);
      this.positionBuffer = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, this.positionBuffer); gl.bufferData(gl.ARRAY_BUFFER, this.positions, gl.STATIC_DRAW);
      const position = gl.getAttribLocation(program, "a_position"); gl.enableVertexAttribArray(position); gl.vertexAttribPointer(position, 2, gl.FLOAT, false, 0, 0);
      this.colorBuffer = gl.createBuffer(); this.colorAttribute = gl.getAttribLocation(program, "a_color"); gl.enableVertexAttribArray(this.colorAttribute);
      this.uniforms = { pan: gl.getUniformLocation(program, "u_pan"), zoom: gl.getUniformLocation(program, "u_zoom"), aspect: gl.getUniformLocation(program, "u_aspect"), size: gl.getUniformLocation(program, "u_size") };
      gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); gl.disable(gl.DEPTH_TEST);
    }

    resize() {
      const ratio = Math.min(devicePixelRatio || 1, 2), rect = this.stage.getBoundingClientRect();
      this.width = rect.width; this.height = rect.height; this.aspect = Math.max(.62, rect.width / rect.height * .86);
      this.canvas.width = Math.round(rect.width * ratio); this.canvas.height = Math.round(rect.height * ratio);
      this.gl?.viewport(0, 0, this.canvas.width, this.canvas.height);
    }

    draw() {
      if (!this.gl) return; const gl = this.gl; gl.useProgram(this.program);
      gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
      gl.uniform2f(this.uniforms.pan, this.panX, this.panY); gl.uniform1f(this.uniforms.zoom, this.zoom);
      gl.uniform1f(this.uniforms.aspect, this.aspect); gl.uniform1f(this.uniforms.size, clamp(1.7 + Math.log2(this.zoom + 1) * .58, 1.8, 5.2) * Math.min(devicePixelRatio || 1, 2));
      gl.bindBuffer(gl.ARRAY_BUFFER, this.positionBuffer); const p = gl.getAttribLocation(this.program, "a_position"); gl.vertexAttribPointer(p, 2, gl.FLOAT, false, 0, 0);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.colorBuffer); gl.vertexAttribPointer(this.colorAttribute, 4, gl.UNSIGNED_BYTE, true, 0, 0);
      gl.drawArrays(gl.POINTS, 0, this.count); this.positionLabels();
    }

    pointColor(index) {
      if (this.mode === "status") return hexRGB(STATUS_COLORS[this.points.status[index]] || STATUS_COLORS[0]);
      if (this.mode === "taxonomy") { const mask = this.points.taxon_mask[index]; return hexRGB(TAXON_COLORS[mask] || (mask && (mask & (mask - 1)) ? TAXON_COLORS.mixed : TAXON_COLORS.other)); }
      if (this.mode === "family") { const family = this.points.family[index]; return family >= 0 ? hashColor(`f${family}`) : [103, 117, 123]; }
      if (this.mode === "cluster") return hashColor(`c${this.points.cluster[index]}`);
      return hexRGB(TYPE_COLORS[this.points.type13[index] % TYPE_COLORS.length]);
    }

    matches(index) {
      if (this.filter !== "all" && !STATUS_GROUPS[this.filter]?.has(this.points.status[index])) return false;
      if (this.queryMatches && !this.querySet.has(index)) return false;
      return true;
    }

    updateColors() {
      if (!this.data) return; const colors = new Uint8Array(this.count * 4); let visible = 0;
      for (let i = 0; i < this.count; i += 1) {
        const rgb = this.pointColor(i), show = this.matches(i); if (show) visible += 1;
        colors[i * 4] = rgb[0]; colors[i * 4 + 1] = rgb[1]; colors[i * 4 + 2] = rgb[2]; colors[i * 4 + 3] = show ? 205 : 3;
      }
      this.colors = colors; const gl = this.gl; gl.bindBuffer(gl.ARRAY_BUFFER, this.colorBuffer); gl.bufferData(gl.ARRAY_BUFFER, colors, gl.DYNAMIC_DRAW);
      this.root.querySelector("[data-visible-count]").textContent = `${compact(visible)} RNAs visible`;
      this.renderLegend(); this.draw();
    }

    updateQuery() {
      if (!this.query) { this.queryMatches = null; this.querySet = null; this.updateColors(); return; }
      const familyMatches = new Set();
      this.families.forEach((family, index) => { if (`${family.accession} ${family.description || ""}`.toLowerCase().includes(this.query)) familyMatches.add(index); });
      let requestedCluster = null; const match = this.query.match(/^(?:cluster\s*|c)(\d+)$/); if (match) requestedCluster = Number(match[1]);
      const results = [];
      for (let i = 0; i < this.count; i += 1) {
        if (this.points.id[i].toLowerCase().includes(this.query) || familyMatches.has(this.points.family[i]) || this.points.cluster[i] === requestedCluster) results.push(i);
      }
      this.queryMatches = results; this.querySet = new Set(results); this.updateColors();
    }

    renderLegend() {
      let items;
      if (this.mode === "status") items = STATUS_SHORT.map((label, index) => [label, STATUS_COLORS[index]]);
      else if (this.mode === "taxonomy") items = [["Bacteria", TAXON_COLORS[1]], ["Archaea", TAXON_COLORS[2]], ["Eukaryota", TAXON_COLORS[4]], ["Viruses", TAXON_COLORS[8]], ["Mixed", TAXON_COLORS.mixed]];
      else if (this.mode === "type13") items = this.data.dictionaries.ncrna13.slice(0, 8).map((label, index) => [label, TYPE_COLORS[index]]);
      else items = [[this.mode === "cluster" ? "Stable color per Leiden cluster" : "Stable color per RNA family", "#7cc7c2"], ["Unannotated", "#67757b"]];
      this.legend.innerHTML = items.map(([label, color]) => `<span class="legend-item" style="--legend-color:${color}"><i></i>${escapeHTML(label)}</span>`).join("");
    }

    buildGrid() {
      const { minX, maxX, minY, maxY } = this.bounds;
      for (let i = 0; i < this.count; i += 1) {
        const gx = clamp(Math.floor((this.positions[i * 2] - minX) / (maxX - minX || 1) * (this.gridSize - 1)), 0, this.gridSize - 1);
        const gy = clamp(Math.floor((this.positions[i * 2 + 1] - minY) / (maxY - minY || 1) * (this.gridSize - 1)), 0, this.gridSize - 1);
        const key = gy * this.gridSize + gx; let bucket = this.grid.get(key); if (!bucket) { bucket = []; this.grid.set(key, bucket); } bucket.push(i);
      }
    }

    screenToData(clientX, clientY) {
      const rect = this.canvas.getBoundingClientRect(), clipX = (clientX - rect.left) / rect.width * 2 - 1, clipY = 1 - (clientY - rect.top) / rect.height * 2;
      return { x: clipX * this.aspect / this.zoom - this.panX, y: clipY / this.zoom - this.panY };
    }

    dataToScreen(x, y) { return { x: ((x + this.panX) * this.zoom / this.aspect + 1) * this.width / 2, y: (1 - (y + this.panY) * this.zoom) * this.height / 2 }; }

    nearest(clientX, clientY, pixels) {
      if (!this.data) return -1; const target = this.screenToData(clientX, clientY), b = this.bounds;
      const gx = clamp(Math.floor((target.x - b.minX) / (b.maxX - b.minX || 1) * (this.gridSize - 1)), 0, this.gridSize - 1);
      const gy = clamp(Math.floor((target.y - b.minY) / (b.maxY - b.minY || 1) * (this.gridSize - 1)), 0, this.gridSize - 1);
      const radius = Math.max(1, Math.ceil(pixels / this.width * this.gridSize * this.aspect / this.zoom * 2)); let best = -1, bestD = pixels * pixels;
      for (let y = Math.max(0, gy - radius); y <= Math.min(this.gridSize - 1, gy + radius); y += 1) for (let x = Math.max(0, gx - radius); x <= Math.min(this.gridSize - 1, gx + radius); x += 1) {
        const bucket = this.grid.get(y * this.gridSize + x); if (!bucket) continue;
        for (const index of bucket) { if (!this.matches(index)) continue; const s = this.dataToScreen(this.positions[index * 2], this.positions[index * 2 + 1]); const dx = s.x - (clientX - this.canvas.getBoundingClientRect().left), dy = s.y - (clientY - this.canvas.getBoundingClientRect().top), d = dx * dx + dy * dy; if (d < bestD) { bestD = d; best = index; } }
      }
      return best;
    }

    wheel(event) {
      event.preventDefault(); this.enter(); const before = this.screenToData(event.clientX, event.clientY);
      this.zoom = clamp(this.zoom * Math.exp(-event.deltaY * .00115), .45, 36);
      const rect = this.canvas.getBoundingClientRect(), clipX = (event.clientX - rect.left) / rect.width * 2 - 1, clipY = 1 - (event.clientY - rect.top) / rect.height * 2;
      this.panX = clipX * this.aspect / this.zoom - before.x; this.panY = clipY / this.zoom - before.y; this.draw();
    }

    pointerDown(event) { if (event.button !== 0) return; this.enter(); this.drag = { id: event.pointerId, x: event.clientX, y: event.clientY, panX: this.panX, panY: this.panY, moved: false }; this.canvas.setPointerCapture?.(event.pointerId); this.stage.classList.add("is-dragging"); }
    pointerMove(event) {
      if (this.drag && this.drag.id === event.pointerId) { const dx = event.clientX - this.drag.x, dy = event.clientY - this.drag.y; if (Math.abs(dx) + Math.abs(dy) > 3) this.drag.moved = true; this.panX = this.drag.panX + dx * 2 * this.aspect / (this.width * this.zoom); this.panY = this.drag.panY - dy * 2 / (this.height * this.zoom); this.draw(); return; }
      cancelAnimationFrame(this.pointerFrame); this.pointerFrame = requestAnimationFrame(() => this.hover(event));
    }
    pointerUp(event) { if (this.drag?.id === event.pointerId) { setTimeout(() => { this.drag = null; }, 0); this.stage.classList.remove("is-dragging"); } }

    hover(event) {
      const index = this.nearest(event.clientX, event.clientY, 11); this.hovered = index;
      if (index < 0) { this.tooltip.hidden = true; return; }
      const family = this.points.family[index] >= 0 ? this.families[this.points.family[index]] : null;
      this.tooltip.innerHTML = `<strong>${escapeHTML(this.points.id[index])}</strong><span>${family ? `${escapeHTML(family.accession)} · ${escapeHTML(family.description || "RNA family")}` : "No family label"}</span><span><b>${escapeHTML(STATUS_SHORT[this.points.status[index]])}</b> · cluster ${this.points.cluster[index].toLocaleString()} · ${this.points.length[index]} nt</span>`;
      this.tooltip.hidden = false; const pad = 14, box = this.tooltip.getBoundingClientRect(); let x = event.clientX + 14, y = event.clientY + 14; if (x + box.width > innerWidth - pad) x = event.clientX - box.width - 14; if (y + box.height > innerHeight - pad) y = event.clientY - box.height - 14; this.tooltip.style.left = `${x}px`; this.tooltip.style.top = `${y}px`;
    }

    buildLabels() {
      const candidates = this.data.clusters.map((row, id) => ({ id, row })).filter((item) => item.row[3] >= 0 || item.row[6]);
      candidates.sort((a, b) => b.row[0] - a.row[0]); this.labelsLayer.innerHTML = "";
      const chosen = [], seenFamilies = new Set();
      for (const item of candidates) {
        if (item.row[3] >= 0) {
          if (seenFamilies.has(item.row[3])) continue;
          seenFamilies.add(item.row[3]);
        }
        chosen.push(item); if (chosen.length >= 28) break;
      }
      for (const item of chosen) {
        const family = item.row[3] >= 0 ? this.families[item.row[3]] : null, element = document.createElement("span"); element.className = "cluster-label";
        element.innerHTML = `<b>${escapeHTML(family?.accession || `Cluster ${item.id}`)}</b>${escapeHTML(family?.description || "novel-candidate neighborhood")}`;
        this.labelsLayer.appendChild(element); this.labels.push({ x: item.row[7], y: item.row[8], element });
      }
    }
    positionLabels() {
      if (!this.labels.length) return; const occupied = [];
      for (const label of this.labels) {
        const p = this.dataToScreen(label.x, label.y), width = Math.max(100, label.element.offsetWidth), height = Math.max(29, label.element.offsetHeight);
        const box = { left: p.x - width / 2 - 10, right: p.x + width / 2 + 10, top: p.y - height / 2 - 7, bottom: p.y + height / 2 + 7 };
        const collision = occupied.some((other) => !(box.right < other.left || box.left > other.right || box.bottom < other.top || box.top > other.bottom));
        const underIntro = !this.entered && box.left < Math.min(500, this.width * .44) && box.top < this.height * .74;
        const on = !this.query && this.filter === "all" && !collision && !underIntro && p.x > 80 && p.x < this.width - 80 && p.y > 80 && p.y < this.height - 60 && this.zoom > .55;
        label.element.style.left = `${p.x}px`; label.element.style.top = `${p.y}px`; label.element.style.opacity = on ? String(clamp(.72 - Math.log2(this.zoom) * .12, .25, .72)) : "0";
        if (on) occupied.push(box);
      }
    }

    updateStats() {
      this.root.querySelector("[data-stat-pilot]").textContent = compact(this.data.summary.pilot_sequences);
      this.root.querySelector("[data-stat-clusters]").textContent = compact(this.data.summary.leiden_clusters);
      this.root.querySelector("[data-stat-display]").textContent = compact(this.data.summary.displayed_sequences);
    }

    taxonomy(mask) { const values = []; if (mask & 1) values.push("Bacteria"); if (mask & 2) values.push("Archaea"); if (mask & 4) values.push("Eukaryota"); if (mask & 8) values.push("Viruses"); return values.join(" · ") || "Unresolved"; }

    async select(index) {
      this.enter(); this.selected = index; const familyIndex = this.points.family[index], family = familyIndex >= 0 ? this.families[familyIndex] : null;
      const cluster = this.data.clusters[this.points.cluster[index]], status = this.points.status[index], type = this.data.dictionaries.ncrna13[this.points.type13[index]];
      this.inspector.innerHTML = `<button type="button" class="inspector-close" data-inspector-close aria-label="Close inspector">×</button>
        <p class="eyebrow">RNAcentral · cluster ${this.points.cluster[index].toLocaleString()}</p><h2>${escapeHTML(this.points.id[index])}</h2>
        <div class="family-title">${family ? `${escapeHTML(family.accession)} · ${escapeHTML(family.description || "")}` : "No resolved RNA family"}</div>
        <span class="evidence-chip">${escapeHTML(STATUS_SHORT[status])}</span>
        <dl><div><dt>Length</dt><dd>${this.points.length[index].toLocaleString()} nt</dd></div><div><dt>RiboGRAM class</dt><dd>${escapeHTML(type)} · ${this.points.type13_conf[index]}%</dd></div>
        <div><dt>Cluster size</dt><dd>${cluster[0].toLocaleString()} RNAs</dd></div><div><dt>Cluster evidence</dt><dd>${cluster[2]} anchors · ${(cluster[4] * 100).toFixed(1)}% purity</dd></div>
        <div><dt>Taxonomy</dt><dd>${escapeHTML(this.taxonomy(this.points.taxon_mask[index]))}</dd></div><div><dt>Observed taxa</dt><dd>${this.points.taxa_count[index].toLocaleString()}</dd></div></dl>
        <div data-detail-target>${this.points.detail_chunk[index] >= 0 ? "<div class='detail-unavailable'>Loading nucleotide-resolution RiboGRAM signals…</div>" : "<div class='detail-unavailable'>Per-nucleotide tracks are published for representative RNAs. This point remains explorable through its cluster and evidence metadata.</div>"}</div>`;
      this.inspector.querySelector("[data-inspector-close]").addEventListener("click", () => this.closeInspector());
      this.inspector.classList.add("is-open"); this.backdrop.classList.add("is-visible"); history.replaceState(null, "", `#u=${encodeURIComponent(this.points.id[index])}`);
      if (this.points.detail_chunk[index] >= 0) {
        try { const record = await this.loadDetail(this.points.detail_chunk[index], this.points.id[index]); if (this.selected === index) this.renderDetail(record); }
        catch (error) { this.inspector.querySelector("[data-detail-target]").innerHTML = `<div class="detail-unavailable">Nucleotide tracks could not be loaded: ${escapeHTML(error.message)}</div>`; }
      }
    }

    async loadDetail(chunk, id) {
      if (!this.detailCache.has(chunk)) {
        const pattern = this.data.details.url_pattern.replace("{chunk:04d}", String(chunk).padStart(4, "0"));
        const base = new URL(this.root.dataset.source, location.href); const url = new URL(pattern, base).href;
        this.detailCache.set(chunk, loadCompressedJSON(url));
      }
      const payload = await this.detailCache.get(chunk), record = payload.records.find((item) => item.id === id);
      if (!record) throw new Error("representative record is absent from its detail chunk"); return record;
    }

    renderDetail(record) {
      const target = this.inspector.querySelector("[data-detail-target]");
      target.innerHTML = `<section class="signal-section"><h3>Context and pairing grammar</h3><p>Layer-7 (teal) and Layer-12 (blue) received attention; Head-20 zero-shot pairing (amber). The dashed line marks the manuscript pairing threshold, 0.4.</p><canvas class="signal-canvas"></canvas><pre class="sequence-block">${escapeHTML(this.formatSequence(record.sequence))}</pre></section>`;
      const canvas = target.querySelector("canvas"); this.drawSignals(canvas, record);
    }

    drawSignals(canvas, record) {
      const ratio = Math.min(devicePixelRatio || 1, 2), width = canvas.clientWidth, height = canvas.clientHeight; canvas.width = width * ratio; canvas.height = height * ratio;
      const ctx = canvas.getContext("2d"); ctx.scale(ratio, ratio); ctx.fillStyle = "#faf9f5"; ctx.fillRect(0, 0, width, height);
      const left = 30, right = 10, top = 14, bottom = 22, plotW = width - left - right, plotH = height - top - bottom;
      ctx.strokeStyle = "#d8dad5"; ctx.lineWidth = 1; for (let i = 0; i <= 3; i += 1) { const y = top + plotH * i / 3; ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(width - right, y); ctx.stroke(); }
      const peak = Math.max(2, ...record.attention7, ...record.attention12), tracks = [[record.attention7, "#2e9b83", peak], [record.attention12, "#467fb6", peak], [record.pair_score, "#d28c2f", .65]];
      for (const [values, color, maximum] of tracks) { ctx.beginPath(); for (let i = 0; i < values.length; i += 1) { const x = left + i / Math.max(1, values.length - 1) * plotW, y = top + plotH - clamp(values[i] / maximum, 0, 1) * plotH; if (!i) ctx.moveTo(x, y); else ctx.lineTo(x, y); } ctx.strokeStyle = color; ctx.lineWidth = 1.25; ctx.globalAlpha = .88; ctx.stroke(); }
      ctx.globalAlpha = 1; const thresholdY = top + plotH - .4 / .65 * plotH; ctx.setLineDash([4, 3]); ctx.strokeStyle = "#bd7b2a"; ctx.beginPath(); ctx.moveTo(left, thresholdY); ctx.lineTo(width - right, thresholdY); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = "#68777a"; ctx.font = "8px Plex, sans-serif"; ctx.fillText("1", left, height - 7); ctx.fillText(String(record.sequence.length), width - right - ctx.measureText(String(record.sequence.length)).width, height - 7); ctx.save(); ctx.translate(10, top + plotH / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("relative signal", -25, 0); ctx.restore();
    }

    formatSequence(sequence) { const parts = []; for (let i = 0; i < sequence.length; i += 10) parts.push(sequence.slice(i, i + 10)); const lines = []; for (let i = 0; i < parts.length; i += 5) lines.push(parts.slice(i, i + 5).join(" ")); return lines.join("\n"); }
    closeInspector() { this.inspector.classList.remove("is-open"); this.backdrop.classList.remove("is-visible"); this.selected = -1; if (location.hash.startsWith("#u=")) history.replaceState(null, "", location.pathname + location.search); }
    restoreHash() { const match = location.hash.match(/^#u=(.+)$/); if (!match) return; const id = decodeURIComponent(match[1]).toLowerCase(), index = this.points.id.findIndex((value) => value.toLowerCase() === id); if (index >= 0) this.select(index); }
  }

  function rootEvent(root, selector, event, handler) { const element = root.querySelector(selector); if (element) element.addEventListener(event, handler); }
  document.querySelectorAll("[data-universe]").forEach((root) => new UniverseAtlas(root).start());
})();
