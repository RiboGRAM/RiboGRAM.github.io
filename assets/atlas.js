(function () {
    "use strict";

    // Fill the desktop canvas with the learned landscape while keeping the
    // outer islands visible.  The previous conservative fit left a large,
    // empty halo around the data on widescreen displays.
    const DEFAULT_ZOOM = 1.22;

    // COMPASS-inspired scientific palette.  The map reads first as three
    // biological classes: cool tones for gene-derived ncRNAs, warm tones for
    // cis-regulators, and violet for introns.  Each subtype occupies a distinct
    // light-to-dark step within its class ramp; deliberately varied alpha values
    // preserve the finer annotation without producing a rainbow effect.
    const TYPE_STYLES = {
        "miRNA": { color: "#59638f", alpha: .62 },
        "sRNA": { color: "#2d8092", alpha: .60 },
        "C/D-box snoRNA": { color: "#1f5c91", alpha: .68 },
        "H/ACA-box snoRNA": { color: "#4c86bf", alpha: .64 },
        "lncRNA": { color: "#3f8fa1", alpha: .57 },
        "Gene": { color: "#74838f", alpha: .52 },
        "CRISPR": { color: "#465077", alpha: .68 },
        "Antisense RNA": { color: "#466e85", alpha: .59 },
        "Ribozyme": { color: "#237169", alpha: .63 },
        "snRNA": { color: "#745b90", alpha: .64 },
        "Spliceosomal RNA": { color: "#55456f", alpha: .67 },
        "rRNA": { color: "#173a5e", alpha: .68 },
        "Antitoxin RNA": { color: "#3d817b", alpha: .59 },
        "tRNA": { color: "#2e6d9b", alpha: .68 },
        "Cis-regulatory RNA": { color: "#9b6120", alpha: .68 },
        "Leader": { color: "#7c4c18", alpha: .74 },
        "Riboswitch": { color: "#e48925", alpha: .80 },
        "Thermoregulator": { color: "#c6761d", alpha: .75 },
        "IRES": { color: "#a9783d", alpha: .68 },
        "Frameshift element": { color: "#c25735", alpha: .75 },
        "Intron": { color: "#624572", alpha: .84 },
        "Unclassified": { color: "#8c969e", alpha: .45 }
    };
    const TYPE_COLORS = Object.fromEntries(Object.entries(TYPE_STYLES).map(([key, value]) => [key, value.color]));
    const CLASS_STYLES = {
        "Gene": { color: "#2d5f83", alpha: .50 },
        "Cis-reg.": { color: "#b97020", alpha: .60 },
        "Intron": { color: "#694b7c", alpha: .68 }
    };
    const SELECTED_COLOR = "#0073e6";
    const NEUTRAL_COLOR = "#9aa7b2";
    const FAMILY_PALETTE = ["#0073e6", "#2d4773", "#e0a800", "#d9623b", "#1f9e89", "#8e5bc9",
        "#16a3c2", "#c44e8a", "#6f8f3a", "#b0413e", "#5b7fbf", "#a67c52",
        "#3a8f5c", "#9c6ade", "#d98e04", "#4a6378", "#e07a9a", "#2a9d8f"];
    // Matplotlib-like magma palette used by the released attention notebook.
    const HEAT_STOPS = [[0, 0, 4], [28, 16, 68], [79, 18, 123], [129, 37, 129],
        [181, 54, 122], [229, 80, 100], [251, 135, 97], [254, 194, 135], [252, 253, 191]];
    const BASE_COLORS = { A: "#1a9a7a", C: "#2979bd", G: "#d99a24", U: "#d75d75", T: "#d75d75", N: "#8b95a0" };
    const clamp = (value, minimum, maximum) => Math.max(minimum, Math.min(maximum, value));

    function escapeHTML(value) {
        return String(value).replaceAll("&", "&amp;").replaceAll("<", "&lt;")
            .replaceAll(">", "&gt;").replaceAll('"', "&quot;").replaceAll("'", "&#039;");
    }

    function stableIndex(value, length) {
        let hash = 0;
        const text = String(value);
        for (let index = 0; index < text.length; index += 1) hash = ((hash << 5) - hash + text.charCodeAt(index)) | 0;
        return Math.abs(hash) % length;
    }

    function familyColor(value) { return FAMILY_PALETTE[stableIndex(value, FAMILY_PALETTE.length)]; }

    function rfamClass(type) {
        const value = String(type || "");
        if (value.startsWith("Cis-reg")) return "Cis-reg.";
        if (value.startsWith("Intron")) return "Intron";
        return "Gene";
    }

    function heatColor(t) {
        const value = clamp(t, 0, 1) * (HEAT_STOPS.length - 1);
        const left = Math.floor(value), right = Math.min(HEAT_STOPS.length - 1, left + 1), mix = value - left;
        return HEAT_STOPS[left].map((channel, index) => Math.round(channel + (HEAT_STOPS[right][index] - channel) * mix));
    }

    function oneBasedPositionTicks(length, maxTicks = 7) {
        if (length <= 1) return [{ position: 0, label: 1 }];
        const rawStep = Math.max(1, length / maxTicks);
        const magnitude = 10 ** Math.floor(Math.log10(rawStep));
        const multiplier = [1, 2, 5, 10].find((value) => value * magnitude >= rawStep) || 10;
        const step = Math.max(1, Math.round(multiplier * magnitude));
        const labels = [1];
        for (let label = step; label <= length; label += step) labels.push(label);
        if (labels[labels.length - 1] !== length) {
            if (length - labels[labels.length - 1] < step * .6 && labels.length > 1) {
                labels[labels.length - 1] = length;
            } else {
                labels.push(length);
            }
        }
        return [...new Set(labels)].map((label) => ({ position: label - 1, label }));
    }

    function scoreLabel(value) {
        const number = Number(value) || 0;
        if (Math.abs(number) >= 100) return number.toFixed(0);
        if (Math.abs(number) >= 10) return number.toFixed(1);
        return number.toFixed(2);
    }

    function formatSequence(sequence) {
        const groups = [];
        for (let index = 0; index < sequence.length; index += 10) groups.push(sequence.slice(index, index + 10));
        const lines = [];
        for (let index = 0; index < groups.length; index += 5) lines.push(groups.slice(index, index + 5).join(" "));
        return lines.join("\n");
    }

    function decodeUint16(values, scaleMaximum) {
        const binary = atob(values), buffer = new ArrayBuffer(binary.length), bytes = new Uint8Array(buffer);
        for (let index = 0; index < binary.length; index += 1) bytes[index] = binary.charCodeAt(index);
        const view = new DataView(buffer), decoded = new Float32Array(binary.length / 2), scale = Number(scaleMaximum) || 0;
        for (let index = 0; index < decoded.length; index += 1) decoded[index] = view.getUint16(index * 2, true) / 65535 * scale;
        return decoded;
    }

    class RfamExplorer {
        constructor(root) {
            this.root = root;
            this.canvas = root.querySelector("[data-rfx-canvas]");
            this.context = this.canvas.getContext("2d", { alpha: true });
            this.map = root.querySelector("[data-rfx-map]");
            this.tooltip = root.querySelector("[data-rfx-tooltip]");
            this.inspector = root.querySelector("[data-rfx-inspector]");
            this.search = root.querySelector("[data-rfx-search]");
            this.familySelect = root.querySelector("[data-rfx-family]");
            this.colorSelect = root.querySelector("[data-rfx-color]");
            this.legend = root.querySelector("[data-rfx-legend]");
            this.classKey = root.querySelector("[data-rfx-class-key]");
            this.loading = root.querySelector("[data-rfx-loading]");
            this.visibleLabel = root.querySelector("[data-rfx-visible]");
            this.fatal = root.querySelector("[data-rfx-fatal]");
            this.legendShell = root.querySelector(".atlas-legend-shell");
            this.aboutDialog = root.querySelector("[data-about-dialog]");
            this.points = [];
            this.families = new Map();
            this.detailCache = new Map();
            this.dependencyCache = new Map();
            this.filteredPoints = [];
            this.spatialBins = new Map();
            this.binCount = 180;
            this.selected = null;
            this.hovered = null;
            this.query = "";
            this.familyFilter = "all";
            this.colorMode = "subtype";
            this.zoom = DEFAULT_ZOOM;
            this.panX = 0;
            this.panY = 0;
            this.cssWidth = 0;
            this.cssHeight = 0;
            this.drag = null;
            this.selectionRequest = 0;
            this.pointerFrame = null;
            this.subtypeAnchors = [];
            this.familyAnchors = [];
            this.selectedFamilyPoints = [];
            this.resizeObserver = new ResizeObserver(() => this.resize());
        }

        async init() {
            this.bindControls();
            this.resizeObserver.observe(this.map);
            try {
                this.manifestURL = new URL(this.root.dataset.manifest, document.baseURI);
                this.data = await this.loadManifest();
                if (Number(this.data.schema_version) !== 3) throw new Error("the complete Rfam 15.0 data bundle has not been installed");
                this.points = this.data.points;
                this.families = new Map(this.data.families.map((item) => [item.accession, item]));
                this.decoratePoints();
                this.buildLabelAnchors();
                this.setBounds();
                this.populateFamilies();
                this.updateSummary();
                this.applyFilter();
                this.renderLegend();
                this.loading.hidden = true;
                this.root.classList.add("atlas-loaded");
                this.restoreURLState();
            } catch (error) {
                this.loading.hidden = true;
                this.fatal.hidden = false;
                this.fatal.textContent = `The explorer data could not be loaded (${error.message}).`;
            }
        }

        async loadManifest() {
            if ("DecompressionStream" in window) {
                try {
                    const compressedURL = new URL(`${this.manifestURL.href}.gz`);
                    const compressed = await fetch(compressedURL, { cache: "no-store" });
                    if (compressed.ok && compressed.body) {
                        const stream = compressed.body.pipeThrough(new DecompressionStream("gzip"));
                        return await new Response(stream).json();
                    }
                } catch (_error) {
                    // Fall through to the uncompressed release asset.
                }
            }
            const response = await fetch(this.manifestURL, { cache: "no-store" });
            if (!response.ok) throw new Error(`HTTP ${response.status}`);
            return response.json();
        }

        decoratePoints() {
            for (const point of this.points) {
                const family = this.families.get(point.family);
                point._family = family;
                point._type = family ? family.subtype : "Unclassified";
                point._class = family ? rfamClass(family.type) : "Gene";
                const subtypeStyle = TYPE_STYLES[point._type] || TYPE_STYLES.Unclassified;
                const classStyle = CLASS_STYLES[point._class] || CLASS_STYLES.Gene;
                point._typeColor = subtypeStyle.color;
                point._typeAlpha = subtypeStyle.alpha;
                point._classColor = classStyle.color;
                point._classAlpha = classStyle.alpha;
                point._familyColor = familyColor(point.family);
                point._search = [point.id, point.accession, point.family, family ? family.name : "",
                    family ? family.description : "", family ? family.type : "", family ? family.subtype : "",
                    point.sequence].join(" ").toLowerCase();
            }
        }

        buildLabelAnchors() {
            const subtypeMap = new Map();
            for (const point of this.points) {
                const key = point._type;
                const anchor = subtypeMap.get(key) || { label: key, x: 0, y: 0, count: 0, color: point._typeColor };
                anchor.x += point.x;
                anchor.y += point.y;
                anchor.count += 1;
                subtypeMap.set(key, anchor);
            }
            this.subtypeAnchors = [...subtypeMap.values()]
                .map((anchor) => ({ ...anchor, x: anchor.x / anchor.count, y: anchor.y / anchor.count }))
                .sort((left, right) => right.count - left.count);
            this.familyAnchors = this.data.families
                .filter((family) => Array.isArray(family.centroid_2d))
                .sort((left, right) => right.seed_count - left.seed_count)
                .slice(0, 90)
                .map((family) => ({
                    label: family.name,
                    accession: family.accession,
                    x: family.centroid_2d[0],
                    y: family.centroid_2d[1],
                    count: family.seed_count,
                    color: TYPE_COLORS[family.subtype] || familyColor(family.accession)
                }));
        }

        restoreURLState() {
            const params = new URLSearchParams(window.location.hash.replace(/^#/, ""));
            const family = params.get("family");
            const sequence = params.get("sequence");
            this.search.value = "";
            this.query = "";
            this.familyFilter = family && this.families.has(family) ? family : "all";
            this.familySelect.value = this.familyFilter;
            this.applyFilter();
            if (sequence) {
                const point = this.points.find((item) => item.id === sequence || item.accession === sequence);
                if (point) {
                    // Keep enough of the atlas visible to preserve geographic context.  Tiny
                    // families previously zoomed all the way to individual pixels, which made
                    // deep links feel empty rather than like a place in a larger universe.
                    if (this.familyFilter === "all") this.fitPoints(this.points.filter((item) => item.family === point.family), 2.65);
                    else this.fitPoints(this.filteredPoints);
                    this.selectPoint(point, false);
                    return;
                }
            }
            this.selectionRequest += 1;
            this.selected = null;
            this.selectedFamilyPoints = [];
            this.root.classList.remove("is-inspector-open");
            if (this.familyFilter === "all") this.resetView(); else this.fitPoints(this.filteredPoints);
        }

        bindControls() {
            let searchTimer = null;
            this.search.addEventListener("input", () => {
                window.clearTimeout(searchTimer);
                searchTimer = window.setTimeout(() => {
                    this.query = this.search.value.trim().toLowerCase();
                    if (this.query) {
                        this.familySelect.value = "all";
                        this.familyFilter = "all";
                    }
                    this.applyFilter();
                    if (!this.query) {
                        this.resetView();
                    } else if (this.filteredPoints.length) {
                        this.fitPoints(this.filteredPoints.slice(0, 5000));
                        if (this.filteredPoints.length <= 25) this.selectPoint(this.filteredPoints[0]);
                    }
                }, 150);
            });
            this.search.addEventListener("keydown", (event) => {
                if (event.key !== "Enter") return;
                event.preventDefault();
                window.clearTimeout(searchTimer);
                this.query = this.search.value.trim().toLowerCase();
                if (this.query) {
                    this.familySelect.value = "all";
                    this.familyFilter = "all";
                }
                this.applyFilter();
                if (this.filteredPoints.length) {
                    this.selectPoint(this.filteredPoints[0]);
                    this.fitPoints(this.filteredPoints.slice(0, 4000));
                }
            });
            this.familySelect.addEventListener("change", () => {
                window.clearTimeout(searchTimer);
                this.search.value = "";
                this.query = "";
                this.familyFilter = this.familySelect.value;
                this.applyFilter();
                if (this.familyFilter === "all") {
                    this.resetView();
                } else {
                    this.fitPoints(this.filteredPoints);
                    if (this.filteredPoints.length) this.selectPoint(this.filteredPoints[0]);
                }
            });
            this.colorSelect.addEventListener("change", () => {
                this.colorMode = this.colorSelect.value;
                this.renderLegend();
                this.draw();
            });
            this.root.querySelector("[data-rfx-reset]").addEventListener("click", () => {
                this.search.value = "";
                this.familySelect.value = "all";
                this.query = "";
                this.familyFilter = "all";
                this.selected = null;
                this.selectedFamilyPoints = [];
                this.applyFilter();
                this.resetView();
                this.closeInspector(true);
            });
            this.canvas.addEventListener("pointerdown", (event) => this.pointerDown(event));
            this.canvas.addEventListener("pointermove", (event) => this.queuePointerMove(event));
            this.canvas.addEventListener("pointerup", (event) => this.pointerUp(event));
            this.canvas.addEventListener("pointercancel", () => this.cancelDrag());
            this.canvas.addEventListener("pointerleave", () => { if (!this.drag) this.hideTooltip(); });
            this.canvas.addEventListener("wheel", (event) => this.onWheel(event), { passive: false });
            this.canvas.addEventListener("dblclick", () => this.resetView());
            this.canvas.addEventListener("keydown", (event) => {
                if (event.key === "Enter" && this.hovered) this.selectPoint(this.hovered);
                if (event.key === "Escape") this.closeInspector();
            });
            document.addEventListener("keydown", (event) => {
                if (event.key === "/" && document.activeElement !== this.search) {
                    event.preventDefault();
                    this.search.focus();
                }
                if (event.key === "Escape" && this.root.classList.contains("is-inspector-open")) this.closeInspector();
            });
            window.addEventListener("hashchange", () => this.restoreURLState());
            this.root.querySelector("[data-inspector-backdrop]")?.addEventListener("click", () => this.closeInspector());
            this.root.querySelector("[data-legend-toggle]")?.addEventListener("click", (event) => {
                const expanded = this.legendShell.classList.toggle("is-expanded");
                event.currentTarget.setAttribute("aria-expanded", String(expanded));
                event.currentTarget.textContent = expanded ? "Hide legend" : "Legend";
            });
            this.root.querySelector("[data-atlas-about]")?.addEventListener("click", () => this.aboutDialog?.showModal());
            this.root.querySelector("[data-about-close]")?.addEventListener("click", () => this.aboutDialog?.close());
            this.aboutDialog?.addEventListener("click", (event) => {
                if (event.target === this.aboutDialog) this.aboutDialog.close();
            });
        }

        updateSummary() {
            const points = this.root.querySelector("[data-rfx-stat='points']");
            const families = this.root.querySelector("[data-rfx-stat='families']");
            if (points) points.textContent = this.data.summary.points.toLocaleString();
            if (families) families.textContent = this.data.summary.families.toLocaleString();
        }

        populateFamilies() {
            const fragment = document.createDocumentFragment();
            for (const family of this.data.families) {
                const option = document.createElement("option");
                option.value = family.accession;
                option.textContent = `${family.accession} · ${family.name} (${family.seed_count})`;
                fragment.appendChild(option);
            }
            this.familySelect.appendChild(fragment);
        }

        setBounds() {
            let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
            for (const point of this.points) {
                minX = Math.min(minX, point.x); maxX = Math.max(maxX, point.x);
                minY = Math.min(minY, point.y); maxY = Math.max(maxY, point.y);
            }
            this.bounds = { minX, maxX, minY, maxY };
        }

        resize() {
            const rect = this.map.getBoundingClientRect();
            const width = Math.max(1, Math.round(rect.width)), height = Math.max(1, Math.round(rect.height));
            const dpr = Math.min(window.devicePixelRatio || 1, 2);
            if (width === this.cssWidth && height === this.cssHeight && this.dpr === dpr) return;
            this.cssWidth = width; this.cssHeight = height; this.dpr = dpr;
            this.canvas.width = Math.round(width * dpr); this.canvas.height = Math.round(height * dpr);
            this.context.setTransform(dpr, 0, 0, dpr, 0, 0);
            this.draw();
        }

        plotGeometry() {
            const padding = 42, rangeX = Math.max(1e-6, this.bounds.maxX - this.bounds.minX),
                rangeY = Math.max(1e-6, this.bounds.maxY - this.bounds.minY);
            return {
                scale: Math.min(Math.max(1, this.cssWidth - padding * 2) / rangeX,
                    Math.max(1, this.cssHeight - padding * 2) / rangeY),
                midX: (this.bounds.minX + this.bounds.maxX) / 2,
                midY: (this.bounds.minY + this.bounds.maxY) / 2
            };
        }

        screenPosition(point) {
            const geometry = this.plotGeometry();
            return {
                x: this.cssWidth / 2 + (point.x - geometry.midX) * geometry.scale * this.zoom + this.panX,
                y: this.cssHeight / 2 - (point.y - geometry.midY) * geometry.scale * this.zoom + this.panY
            };
        }

        pointColor(point) {
            if (this.colorMode === "neutral") return NEUTRAL_COLOR;
            if (this.colorMode === "family") return point._familyColor;
            if (this.colorMode === "class") return point._classColor;
            return point._typeColor;
        }

        pointAlpha(point, selectedFamily = false) {
            let alpha = .58;
            if (this.colorMode === "neutral") alpha = .46;
            else if (this.colorMode === "family") alpha = .42;
            else if (this.colorMode === "class") alpha = point._classAlpha;
            else alpha = point._typeAlpha;
            if (this.filteredPoints.length < 30000) alpha = Math.min(.88, alpha * 1.16);
            else alpha = Math.min(.95, alpha * 1.15);
            return selectedFamily ? alpha * .30 : alpha;
        }

        draw() {
            if (!this.cssWidth || !this.cssHeight || !this.context) return;
            const ctx = this.context;
            ctx.clearRect(0, 0, this.cssWidth, this.cssHeight);
            this.drawGrid(ctx);
            if (!this.filteredPoints.length) {
                ctx.fillStyle = "#77838f";
                ctx.font = '600 13px "IBM Plex Sans", "Segoe UI", sans-serif';
                ctx.textAlign = "center"; ctx.textBaseline = "middle";
                ctx.fillText("No matching sequences", this.cssWidth / 2, this.cssHeight / 2);
                return;
            }
            const radius = this.filteredPoints.length <= 10 ? 4.2
                : this.filteredPoints.length <= 200 ? 3
                    : this.zoom > 3.5 ? 2.35 : this.zoom > 1.8 ? 1.55 : this.filteredPoints.length > 30000 ? .70 : 1.35;
            const selectedFamily = this.selected && this.filteredPoints.includes(this.selected) ? this.selected.family : null;
            for (const point of this.filteredPoints) {
                const position = this.screenPosition(point);
                if (position.x < -5 || position.y < -5 || position.x > this.cssWidth + 5 || position.y > this.cssHeight + 5) continue;
                ctx.globalAlpha = this.pointAlpha(point, Boolean(selectedFamily));
                ctx.fillStyle = this.pointColor(point);
                ctx.fillRect(position.x - radius, position.y - radius, radius * 2, radius * 2);
            }
            if (selectedFamily) {
                this.drawNeighborhood(ctx);
                ctx.globalAlpha = .94;
                for (const point of this.selectedFamilyPoints) {
                    const position = this.screenPosition(point);
                    if (position.x < -5 || position.y < -5 || position.x > this.cssWidth + 5 || position.y > this.cssHeight + 5) continue;
                    ctx.fillStyle = SELECTED_COLOR;
                    const familyRadius = Math.max(1.75, radius * 1.65);
                    ctx.fillRect(position.x - familyRadius, position.y - familyRadius, familyRadius * 2, familyRadius * 2);
                }
            }
            ctx.globalAlpha = 1;
            this.drawLabels(ctx);
            if (this.hovered) this.drawMarker(ctx, this.hovered, 5.5, "#15324b");
            if (this.selected && this.filteredPoints.includes(this.selected)) {
                this.drawMarker(ctx, this.selected, 6.5, "#ffffff", 3.6);
                this.drawMarker(ctx, this.selected, 8.5, "#16384e", 2.1);
                this.drawMarker(ctx, this.selected, 13, SELECTED_COLOR, 1.6);
            }
        }

        drawNeighborhood(ctx) {
            if (!this.selected || !this.selectedFamilyPoints.length) return;
            const source = this.screenPosition(this.selected);
            const neighbors = this.selectedFamilyPoints
                .filter((point) => point !== this.selected)
                .map((point) => ({ point, distance: (point.x - this.selected.x) ** 2 + (point.y - this.selected.y) ** 2 }))
                .sort((left, right) => left.distance - right.distance)
                .slice(0, 10);
            ctx.save();
            ctx.globalAlpha = .1;
            ctx.lineWidth = .75;
            ctx.strokeStyle = SELECTED_COLOR;
            for (const neighbor of neighbors) {
                const target = this.screenPosition(neighbor.point);
                ctx.beginPath();
                ctx.moveTo(source.x, source.y);
                ctx.lineTo(target.x, target.y);
                ctx.stroke();
            }
            ctx.restore();
        }

        drawLabels(ctx) {
            if (this.filteredPoints.length < 3) return;
            let anchors;
            if (this.selected) {
                const family = this.selected._family;
                anchors = [{
                    label: family.name,
                    accession: family.accession,
                    x: family.centroid_2d?.[0] ?? this.selected.x,
                    y: family.centroid_2d?.[1] ?? this.selected.y,
                    count: family.seed_count,
                    color: SELECTED_COLOR,
                    selected: true
                }];
            } else if (this.familyFilter !== "all") {
                const family = this.families.get(this.familyFilter);
                anchors = family && Array.isArray(family.centroid_2d) ? [{
                    label: family.name,
                    accession: family.accession,
                    x: family.centroid_2d[0], y: family.centroid_2d[1], count: family.seed_count,
                    color: TYPE_COLORS[family.subtype] || familyColor(family.accession)
                }] : [];
            } else if (this.zoom < 2.35) {
                // The overview is intentionally label-free, like an atlas rather than a
                // dashboard.  Labels appear as the visitor zooms into actual geography.
                anchors = [];
            } else {
                anchors = this.familyAnchors.slice(0, this.zoom > 4 ? 42 : 22);
            }
            const occupied = [];
            ctx.save();
            ctx.textBaseline = "middle";
            for (const anchor of anchors) {
                const position = this.screenPosition(anchor);
                if (position.x < 80 || position.y < 95 || position.x > this.cssWidth - 80 || position.y > this.cssHeight - 54) continue;
                const isSubtype = !anchor.accession;
                ctx.font = `${isSubtype ? 650 : 580} ${isSubtype ? 9 : 8}px "IBM Plex Sans", "Segoe UI", sans-serif`;
                const label = isSubtype ? `${anchor.label} · ${anchor.count.toLocaleString()}` : `${anchor.label} · ${anchor.accession}`;
                const width = ctx.measureText(label).width + (anchor.selected ? 22 : (isSubtype ? 16 : 13));
                const height = anchor.selected ? 24 : (isSubtype ? 22 : 19);
                const centerX = anchor.selected ? position.x + width / 2 + 17 : position.x;
                const box = { left: centerX - width / 2, right: centerX + width / 2, top: position.y - height / 2, bottom: position.y + height / 2 };
                if (occupied.some((item) => !(box.right < item.left || box.left > item.right || box.bottom < item.top || box.top > item.bottom))) continue;
                occupied.push(box);
                ctx.globalAlpha = anchor.selected ? .97 : (isSubtype ? .92 : .78);
                ctx.fillStyle = anchor.selected ? "rgba(255, 255, 255, .97)" : "rgba(255, 255, 255, .86)";
                ctx.beginPath();
                if (ctx.roundRect) ctx.roundRect(box.left, box.top, width, height, height / 2);
                else ctx.rect(box.left, box.top, width, height);
                ctx.fill();
                ctx.globalAlpha = 1;
                ctx.fillStyle = anchor.color;
                ctx.beginPath();
                ctx.arc(box.left + (anchor.selected ? 9 : 7), position.y, anchor.selected ? 3 : 2.2, 0, Math.PI * 2);
                ctx.fill();
                ctx.fillStyle = anchor.selected ? "rgba(18, 45, 63, .96)" : (isSubtype ? "rgba(24, 43, 59, .9)" : "rgba(62, 79, 94, .82)");
                ctx.textAlign = "left";
                ctx.fillText(label, box.left + (anchor.selected ? 16 : 12), position.y + .2);
            }
            ctx.restore();
        }

        drawMarker(ctx, point, radius, color, width = 2.2) {
            const position = this.screenPosition(point);
            ctx.beginPath(); ctx.arc(position.x, position.y, radius, 0, Math.PI * 2);
            ctx.lineWidth = width; ctx.strokeStyle = color; ctx.stroke();
        }

        drawGrid(ctx) {
            ctx.save();
            const gradient = ctx.createRadialGradient(this.cssWidth * .5, this.cssHeight * .48, 0,
                this.cssWidth * .5, this.cssHeight * .48, Math.max(this.cssWidth, this.cssHeight) * .62);
            gradient.addColorStop(0, "rgba(255, 255, 255, .46)");
            gradient.addColorStop(1, "rgba(244, 239, 228, 0)");
            ctx.fillStyle = gradient;
            ctx.fillRect(0, 0, this.cssWidth, this.cssHeight);
            ctx.restore();
        }

        applyFilter() {
            const query = this.query;
            this.filteredPoints = this.points.filter((point) => {
                if (this.familyFilter !== "all" && point.family !== this.familyFilter) return false;
                return !query || point._search.includes(query);
            });
            this.buildSpatialIndex();
            this.visibleLabel.textContent = `${this.filteredPoints.length.toLocaleString()} sequence${this.filteredPoints.length === 1 ? "" : "s"}`;
            this.draw();
        }

        buildSpatialIndex() {
            this.spatialBins = new Map();
            const rangeX = Math.max(1e-9, this.bounds.maxX - this.bounds.minX), rangeY = Math.max(1e-9, this.bounds.maxY - this.bounds.minY);
            for (const point of this.filteredPoints) {
                const bx = clamp(Math.floor((point.x - this.bounds.minX) / rangeX * this.binCount), 0, this.binCount - 1);
                const by = clamp(Math.floor((point.y - this.bounds.minY) / rangeY * this.binCount), 0, this.binCount - 1);
                const key = bx * this.binCount + by, bucket = this.spatialBins.get(key);
                if (bucket) bucket.push(point); else this.spatialBins.set(key, [point]);
            }
        }

        renderLegend() {
            if (!this.data) return;
            if (this.classKey) {
                this.classKey.hidden = !["subtype", "class"].includes(this.colorMode);
                this.classKey.innerHTML = Object.entries(CLASS_STYLES).map(([label, style]) =>
                    `<span><i style="background:${style.color}"></i>${escapeHTML(label === "Gene" ? "Gene / ncRNA" : label)}</span>`
                ).join("");
            }
            if (this.colorMode === "neutral") {
                this.legend.innerHTML = `<span class="rfx-legend-item"><i class="rfx-legend-dot" style="background:${NEUTRAL_COLOR}"></i>All Rfam Seed RNAs</span><span class="rfx-legend-item"><i class="rfx-legend-dot" style="background:${SELECTED_COLOR}"></i>Selected family</span>`;
                return;
            }
            if (this.colorMode === "family") {
                this.legend.innerHTML = `<span class="rfx-legend-item"><i class="rfx-legend-dot" style="background:#2d4773"></i>Stable color per Rfam accession</span><span class="rfx-legend-item">Official Rfam 15.0 family annotation</span>`;
                return;
            }
            if (this.colorMode === "class") {
                this.legend.innerHTML = Object.entries(CLASS_STYLES).map(([label, style]) =>
                    `<span class="rfx-legend-item"><i class="rfx-legend-dot" style="background:${style.color};opacity:.9"></i>${escapeHTML(label === "Gene" ? "Gene · ncRNA" : label)}</span>`
                ).join("");
                return;
            }
            const subtypeInfo = new Map();
            for (const family of this.data.families) {
                const current = subtypeInfo.get(family.subtype) || { label: family.subtype, className: rfamClass(family.type), count: 0 };
                current.count += Number(family.seed_count) || 0;
                subtypeInfo.set(family.subtype, current);
            }
            const groups = { "Gene": [], "Cis-reg.": [], "Intron": [] };
            for (const item of subtypeInfo.values()) (groups[item.className] || groups.Gene).push(item);
            for (const items of Object.values(groups)) items.sort((left, right) => right.count - left.count);
            this.legend.innerHTML = Object.entries(groups).map(([className, items]) => {
                const heading = className === "Gene" ? "Gene · ncRNA subtypes" : className;
                const entries = items.map((item) => {
                    const style = TYPE_STYLES[item.label] || TYPE_STYLES.Unclassified;
                    return `<span class="rfx-legend-item"><i class="rfx-legend-dot" style="background:${style.color};opacity:.9"></i>${escapeHTML(item.label)}</span>`;
                }).join("");
                return `<section class="rfx-legend-group"><strong>${escapeHTML(heading)}</strong><div>${entries}</div></section>`;
            }).join("");
        }

        nearestPoint(x, y, radius = 12) {
            const geometry = this.plotGeometry();
            const dataX = geometry.midX + (x - this.cssWidth / 2 - this.panX) / (geometry.scale * this.zoom);
            const dataY = geometry.midY - (y - this.cssHeight / 2 - this.panY) / (geometry.scale * this.zoom);
            const rangeX = Math.max(1e-9, this.bounds.maxX - this.bounds.minX), rangeY = Math.max(1e-9, this.bounds.maxY - this.bounds.minY);
            const radiusData = radius / (geometry.scale * this.zoom);
            const bx0 = clamp(Math.floor((dataX - radiusData - this.bounds.minX) / rangeX * this.binCount), 0, this.binCount - 1);
            const bx1 = clamp(Math.floor((dataX + radiusData - this.bounds.minX) / rangeX * this.binCount), 0, this.binCount - 1);
            const by0 = clamp(Math.floor((dataY - radiusData - this.bounds.minY) / rangeY * this.binCount), 0, this.binCount - 1);
            const by1 = clamp(Math.floor((dataY + radiusData - this.bounds.minY) / rangeY * this.binCount), 0, this.binCount - 1);
            let nearest = null, best = radius * radius;
            for (let bx = bx0; bx <= bx1; bx += 1) for (let by = by0; by <= by1; by += 1) {
                const bucket = this.spatialBins.get(bx * this.binCount + by) || [];
                for (const point of bucket) {
                    const position = this.screenPosition(point), distance = (position.x - x) ** 2 + (position.y - y) ** 2;
                    if (distance <= best) { best = distance; nearest = point; }
                }
            }
            return nearest;
        }

        pointerDown(event) {
            if (event.button !== 0) return;
            this.canvas.setPointerCapture(event.pointerId);
            this.drag = { id: event.pointerId, startX: event.offsetX, startY: event.offsetY, panX: this.panX, panY: this.panY, moved: false };
        }

        queuePointerMove(event) {
            const state = { pointerId: event.pointerId, offsetX: event.offsetX, offsetY: event.offsetY };
            if (this.drag) { this.pointerMove(state); return; }
            if (this.pointerFrame) cancelAnimationFrame(this.pointerFrame);
            this.pointerFrame = requestAnimationFrame(() => { this.pointerFrame = null; this.pointerMove(state); });
        }

        pointerMove(event) {
            if (this.drag && this.drag.id === event.pointerId) {
                const dx = event.offsetX - this.drag.startX, dy = event.offsetY - this.drag.startY;
                if (Math.abs(dx) + Math.abs(dy) > 3) this.drag.moved = true;
                if (this.drag.moved) {
                    this.canvas.classList.add("rfx-dragging"); this.panX = this.drag.panX + dx; this.panY = this.drag.panY + dy;
                    this.hideTooltip(); this.draw();
                }
                return;
            }
            const point = this.nearestPoint(event.offsetX, event.offsetY);
            if (point !== this.hovered) { this.hovered = point; this.draw(); }
            if (point) this.showTooltip(point, event.offsetX, event.offsetY); else this.hideTooltip();
        }

        pointerUp(event) {
            if (!this.drag || this.drag.id !== event.pointerId) return;
            const moved = this.drag.moved; this.cancelDrag();
            if (!moved) { const point = this.nearestPoint(event.offsetX, event.offsetY, 14); if (point) this.selectPoint(point); }
        }

        cancelDrag() { this.drag = null; this.canvas.classList.remove("rfx-dragging"); }

        onWheel(event) {
            event.preventDefault();
            const rect = this.canvas.getBoundingClientRect(), mouseX = event.clientX - rect.left, mouseY = event.clientY - rect.top;
            const oldZoom = this.zoom, nextZoom = clamp(oldZoom * Math.exp(-event.deltaY * .0012), .7, 18);
            const centerX = this.cssWidth / 2, centerY = this.cssHeight / 2;
            const localX = (mouseX - centerX - this.panX) / oldZoom, localY = (mouseY - centerY - this.panY) / oldZoom;
            this.panX = mouseX - centerX - localX * nextZoom; this.panY = mouseY - centerY - localY * nextZoom;
            this.zoom = nextZoom; this.draw();
        }

        resetView() { this.zoom = DEFAULT_ZOOM; this.panX = 0; this.panY = 0; this.draw(); }

        fitPoints(points, maximumZoom = 12) {
            if (!points.length || !this.cssWidth || !this.cssHeight) return;
            let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
            for (const point of points) {
                minX = Math.min(minX, point.x); maxX = Math.max(maxX, point.x);
                minY = Math.min(minY, point.y); maxY = Math.max(maxY, point.y);
            }
            const geometry = this.plotGeometry();
            const width = Math.max(1, (maxX - minX) * geometry.scale);
            const height = Math.max(1, (maxY - minY) * geometry.scale);
            this.zoom = clamp(Math.min((this.cssWidth * .62) / width, (this.cssHeight * .62) / height), 1, maximumZoom);
            const centerX = (minX + maxX) / 2, centerY = (minY + maxY) / 2;
            this.panX = -(centerX - geometry.midX) * geometry.scale * this.zoom;
            this.panY = (centerY - geometry.midY) * geometry.scale * this.zoom;
            this.draw();
        }

        showTooltip(point, x, y) {
            const family = point._family;
            this.tooltip.innerHTML = `<strong>${escapeHTML(family.name)} · ${escapeHTML(point.family)}</strong><span>${escapeHTML(point.accession)} · ${point.length.toLocaleString()} nt</span><span>${escapeHTML(family.subtype)} · ${escapeHTML(family.type)}</span>`;
            this.tooltip.hidden = false;
            const width = 240, left = x + width + 22 > this.cssWidth ? x - width - 12 : x + 12,
                top = y + 94 > this.cssHeight ? y - 86 : y + 12;
            this.tooltip.style.left = `${Math.max(8, left)}px`; this.tooltip.style.top = `${Math.max(8, top)}px`;
        }

        hideTooltip() { this.tooltip.hidden = true; }

        async selectPoint(point, updateURL = true) {
            this.selected = point;
            this.selectedFamilyPoints = this.points.filter((item) => item.family === point.family);
            this.root.classList.add("is-inspector-open");
            this.hideTooltip();
            if (updateURL) {
                const params = new URLSearchParams();
                if (this.familyFilter !== "all") params.set("family", this.familyFilter);
                params.set("sequence", point.id);
                window.history.replaceState(null, "", `#${params.toString()}`);
            }
            this.draw();
            const requestId = ++this.selectionRequest;
            this.renderInspectorShell(point);
            if (this.root.dataset.details !== "true") {
                this.inspector.querySelector("[data-rfx-analysis]")?.remove();
                return;
            }
            try {
                const familyDetails = await this.loadFamilyDetails(point._family.detail_url);
                if (requestId !== this.selectionRequest) return;
                const detail = familyDetails.points[point.id];
                if (!detail) throw new Error("detail record is missing");
                this.renderAnalysis(point, detail);
                if (this.root.dataset.dependencies) {
                    try {
                        const dependency = await this.loadDependency(point.family);
                        if (requestId !== this.selectionRequest) return;
                        this.renderDependency(point, dependency, detail);
                    } catch (_dependencyError) {
                        // Some Rfam 15.0 families do not have an eligible Rfam 15.1
                        // representative.  The two released analyses remain available.
                    }
                }
            } catch (_error) {
                if (requestId !== this.selectionRequest) return;
                const target = this.inspector.querySelector("[data-rfx-analysis]");
                if (target) target.innerHTML = `<div class="rfx-analysis-unavailable"><strong>Sequence analyses unavailable</strong><span>The released per-sequence analysis record could not be loaded.</span></div>`;
            }
        }

        closeInspector(clearURL = false) {
            this.selectionRequest += 1;
            this.root.classList.remove("is-inspector-open");
            this.selected = null;
            this.selectedFamilyPoints = [];
            this.draw();
            const params = new URLSearchParams(window.location.hash.replace(/^#/, ""));
            params.delete("sequence");
            if (clearURL || this.familyFilter === "all") params.delete("family");
            const next = params.toString();
            window.history.replaceState(null, "", next ? `#${next}` : `${window.location.pathname}${window.location.search}`);
        }

        async loadFamilyDetails(url) {
            const resolvedURL = new URL(url, this.manifestURL).href;
            if (!this.detailCache.has(resolvedURL)) this.detailCache.set(resolvedURL, fetch(resolvedURL).then((response) => {
                if (!response.ok) throw new Error(`HTTP ${response.status}`);
                if (resolvedURL.endsWith(".gz")) {
                    if (!("DecompressionStream" in window) || !response.body) {
                        throw new Error("this browser cannot decompress the detail bundle");
                    }
                    const stream = response.body.pipeThrough(new DecompressionStream("gzip"));
                    return new Response(stream).json();
                }
                return response.json();
            }).catch((error) => { this.detailCache.delete(resolvedURL); throw error; }));
            return this.detailCache.get(resolvedURL);
        }

        async loadDependency(family) {
            const base = new URL(this.root.dataset.dependencies.replace(/\/?$/, "/"), document.baseURI);
            const resolvedURL = new URL(`${encodeURIComponent(family)}.json.gz`, base).href;
            if (!this.dependencyCache.has(resolvedURL)) this.dependencyCache.set(resolvedURL, fetch(resolvedURL).then((response) => {
                if (!response.ok) throw new Error(`HTTP ${response.status}`);
                if (!("DecompressionStream" in window) || !response.body) {
                    throw new Error("this browser cannot decompress the dependency bundle");
                }
                return new Response(response.body.pipeThrough(new DecompressionStream("gzip"))).json();
            }).catch((error) => { this.dependencyCache.delete(resolvedURL); throw error; }));
            return this.dependencyCache.get(resolvedURL);
        }

        renderInspectorShell(point) {
            const family = point._family;
            const windowNote = point.truncated
                ? `<span class="rfx-window-note">Model plots show nt 1–${point.analysis_length.toLocaleString()} of ${point.length.toLocaleString()}.</span>` : "";
            const clan = family.clan_accession || "Not assigned";
            const sequence = point.sequence.toUpperCase();
            const bases = ["A", "C", "G", "U"].map((base) => ({
                base,
                count: [...sequence].filter((letter) => letter === base || (base === "U" && letter === "T")).length
            }));
            const centroid = family.centroid_2d || [0, 0];
            const centroidDistance = Math.hypot(point.x - centroid[0], point.y - centroid[1]);
            const compositionBars = bases.map((item) => `<i style="width:${(item.count / Math.max(1, sequence.length) * 100).toFixed(2)}%;background:${BASE_COLORS[item.base]}"></i>`).join("");
            const compositionLabels = bases.map((item) => `<span style="color:${BASE_COLORS[item.base]}">${item.base} ${(item.count / Math.max(1, sequence.length) * 100).toFixed(1)}%</span>`).join("");
            this.inspector.innerHTML = `<div class="rfx-inspector-head"><p class="rfx-kicker">Rfam Seed sequence · ${escapeHTML(point.id)}</p><div class="rfx-inspector-title-row"><div><h3>${escapeHTML(family.name)}</h3><div class="rfx-accession">${escapeHTML(point.family)} · ${escapeHTML(point.accession)}</div></div><div class="rfx-head-actions"><button type="button" class="rfx-copy" data-rfx-copy>Copy sequence</button><button type="button" class="rfx-drawer-close" data-rfx-close aria-label="Close sequence details">×</button></div></div><div class="rfx-tag-row"><span class="rfx-tag" style="--tag-color:${point._typeColor}">${escapeHTML(family.subtype)}</span><span class="rfx-tag" style="--tag-color:${point._classColor}">${escapeHTML(point._class)} class</span><span class="rfx-tag" style="--tag-color:#81919c">${point.length.toLocaleString()} nt</span><span class="rfx-tag" style="--tag-color:#6e9a86">${family.seed_count.toLocaleString()} Seed sequences</span></div></div><div class="rfx-inspector-body"><div data-rfx-analysis><div class="rfx-analysis-loading"><i></i><span>Loading sequence analyses…</span></div></div><details class="rfx-disclosure" open><summary>Sequence & annotation <span>5′ → 3′</span></summary><div class="rfx-disclosure-body"><div class="rfx-annotation-grid"><div class="rfx-annotation"><span>Rfam family</span><strong title="${escapeHTML(family.description)}">${escapeHTML(point.family)} · ${escapeHTML(family.name)}</strong><small>${escapeHTML(family.description)}</small></div><div class="rfx-annotation"><span>ncRNA category</span><strong>${escapeHTML(family.subtype)}</strong><small>${escapeHTML(family.type)}</small></div><div class="rfx-annotation"><span>Source sequence</span><strong title="${escapeHTML(point.accession)}">${escapeHTML(point.accession)}</strong><small>Official Rfam 15.0 Seed accession</small></div><div class="rfx-annotation"><span>Rfam clan</span><strong>${escapeHTML(clan)}</strong><small>${escapeHTML(point.id)}</small></div></div>${windowNote}<pre class="rfx-sequence">${escapeHTML(formatSequence(point.sequence))}</pre><div class="rfx-composition"><div class="rfx-composition-bar">${compositionBars}</div><div class="rfx-composition-labels">${compositionLabels}</div></div></div></details><details class="rfx-disclosure"><summary>Representation coordinates <span>Layer 12 · UMAP</span></summary><div class="rfx-disclosure-body"><div class="rfx-space-metrics"><div><span>UMAP 1</span><strong>${point.x.toFixed(3)}</strong></div><div><span>UMAP 2</span><strong>${point.y.toFixed(3)}</strong></div><div><span>To centroid</span><strong>${centroidDistance.toFixed(3)}</strong></div></div></div></details><a class="rfx-rfam-link" href="https://rfam.org/family/${encodeURIComponent(point.family)}" target="_blank" rel="noreferrer"><span>Open ${escapeHTML(point.family)} in Rfam</span><span>↗</span></a></div>`;
            this.inspector.scrollTop = 0;
            this.inspector.querySelector("[data-rfx-copy]").addEventListener("click", async (event) => {
                await navigator.clipboard.writeText(point.sequence);
                const button = event.currentTarget; button.textContent = "Copied";
                window.setTimeout(() => { button.textContent = "Copy sequence"; }, 1200);
            });
            this.inspector.querySelector("[data-rfx-close]").addEventListener("click", () => this.closeInspector());
        }

        renderAnalysis(point, detail) {
            const target = this.inspector.querySelector("[data-rfx-analysis]");
            const contextValues = Array.isArray(detail.contextual_dependency.values)
                ? detail.contextual_dependency.values
                : decodeUint16(detail.contextual_dependency.values, detail.contextual_dependency.scale_max);
            const analysisSequence = point.sequence.slice(0, detail.analysis_length || point.analysis_length);
            const context = detail.contextual_dependency || {};
            const attention = detail.attention || {};
            const pairs = Array.isArray(attention.pairs) ? attention.pairs : [];
            const strongest = pairs.reduce((best, pair) => Number(pair[2]) > Number(best?.[2] || 0) ? pair : best, null);
            const pairList = pairs.slice(0, 8).map((pair) => {
                const left = Number(pair[0]), right = Number(pair[1]);
                const leftBase = analysisSequence[left - 1] || "N", rightBase = analysisSequence[right - 1] || "N";
                return `<li><span>${leftBase}${left}—${rightBase}${right}</span><strong>${Number(pair[2]).toFixed(3)}</strong></li>`;
            }).join("");
            target.innerHTML = `<section class="rfx-interpretation-intro"><p class="rfx-kicker">RiboGRAM-derived analyses</p><h4>Contextual constraint and zero-shot pairing</h4><p>Two complementary readouts from the frozen RiboGRAM backbone: masked-nucleotide contextual recovery and an APC-corrected attention pairing signal.</p><div class="rfx-model-metrics"><div><span>Mean contextual score</span><strong>${Number(context.mean || 0).toFixed(3)}</strong></div><div><span>Maximum contextual score</span><strong>${Number(context.maximum || context.scale_max || 0).toFixed(3)}</strong></div><div><span>Candidate pairs</span><strong>${Number(attention.predicted_pair_count || pairs.length).toLocaleString()}</strong></div><div><span>Maximum APC score</span><strong>${strongest ? Number(strongest[2]).toFixed(3) : "—"}</strong></div></div></section><section class="rfx-analysis-section"><div class="rfx-analysis-heading"><div><p>01 · Masked-nucleotide context</p><h4>Contextual recovery along the sequence</h4></div><span>Position-resolved</span></div><p class="rfx-analysis-note">Confidence with which the frozen model recovers each native nucleotide when that position is masked. This is not a mutation-effect or base-pairing probability.</p><div class="rfx-chart-frame rfx-context-frame"><canvas data-rfx-context aria-label="Position-resolved sequence contextual dependency scores"></canvas><span class="rfx-chart-tip" data-rfx-context-tip hidden></span></div></section><section class="rfx-analysis-section"><div class="rfx-analysis-heading"><div><p>02 · Zero-shot base pairing</p><h4>APC-corrected attention map</h4></div><span>L12 · H20 · APC ≥ ${Number(attention.pair_threshold || .4).toFixed(1)}</span></div><p class="rfx-analysis-note">Symmetrized, APC-corrected Layer-12 Head-20 attention. Cells above the indicated threshold denote candidate base-pairing signals obtained without structure supervision.</p><div class="rfx-zero-benchmark" aria-label="Rfam-wide zero-shot reference recovery"><p>Rfam reference recovery <span>92,781 Seed sequences · not held-out</span></p><div><span><small>Precision</small><strong>67.1%</strong></span><span><small>Recall</small><strong>41.7%</strong></span><span><small>F1</small><strong>51.4%</strong></span></div></div><div class="rfx-chart-frame rfx-attention-frame"><canvas data-rfx-attention aria-label="Zero-shot base-pair signal map"></canvas><span class="rfx-chart-tip" data-rfx-attention-tip hidden></span></div>${pairList ? `<ol class="rfx-pair-list">${pairList}</ol>` : `<p class="rfx-no-pairs">No candidate pair exceeded the specified APC threshold.</p>`}</section><div data-rfx-dependency-slot></div>`;
            requestAnimationFrame(() => {
                this.drawContext(target.querySelector("[data-rfx-context]"), analysisSequence, contextValues);
                this.drawAttention(target.querySelector("[data-rfx-attention]"), attention);
            });
        }

        renderDependency(point, dependency, detail) {
            const slot = this.inspector.querySelector("[data-rfx-dependency-slot]");
            if (!slot || !dependency) return;
            if (String(dependency.sequence_accession) !== String(point.accession)) {
                const representative = this.points.find((item) => item.family === point.family
                    && String(item.accession) === String(dependency.sequence_accession));
                if (!representative) return;
                slot.innerHTML = `<section class="rfx-dependency-available"><div><p>03 · Saturation-mutagenesis map</p><strong>Available for the family representative</strong><span>${escapeHTML(dependency.sequence_accession)} · ${Number(dependency.length).toLocaleString()} nt</span></div><button type="button" data-rfx-open-dependency>Open map</button></section>`;
                slot.querySelector("[data-rfx-open-dependency]").addEventListener("click", () => this.selectPoint(representative));
                return;
            }

            const summary = dependency.summary || {};
            const arcs = Array.isArray(dependency.arcs) ? dependency.arcs : [];
            const pairKey = (left, right) => `${Math.min(Number(left), Number(right))}:${Math.max(Number(left), Number(right))}`;
            const comparisonLength = Math.min(
                Number(detail?.attention?.shape?.[0] || detail?.analysis_length || dependency.length),
                Number(dependency.length)
            );
            const zeroShotPairs = new Set((detail?.attention?.pairs || []).map((pair) => pairKey(Number(pair[0]) - 1, Number(pair[1]) - 1)));
            const q999Pairs = new Set(arcs
                .filter((arc) => Number(arc[4]) === 2 && Number(arc[1]) < comparisonLength)
                .map((arc) => pairKey(arc[0], arc[1])));
            const referencePairs = new Set((dependency.reference_pairs || [])
                .filter((pair) => Number(pair[1]) < comparisonLength)
                .map((pair) => pairKey(pair[0], pair[1])));
            const zeroDependency = [...zeroShotPairs].filter((key) => q999Pairs.has(key));
            const referenceBoth = [...referencePairs].filter((key) => zeroShotPairs.has(key) && q999Pairs.has(key));
            const referenceEither = [...referencePairs].filter((key) => zeroShotPairs.has(key) || q999Pairs.has(key));
            const zeroDependencyRate = zeroShotPairs.size ? zeroDependency.length / zeroShotPairs.size * 100 : 0;
            const referenceBothRate = referencePairs.size ? referenceBoth.length / referencePairs.size * 100 : 0;
            const referenceEitherRate = referencePairs.size ? referenceEither.length / referencePairs.size * 100 : 0;
            const strongest = [...arcs].sort((left, right) => Math.max(Number(right[2]), Number(right[3]))
                - Math.max(Number(left[2]), Number(left[3]))).slice(0, 8);
            const sequence = String(dependency.sequence || point.sequence).toUpperCase().replaceAll("T", "U");
            const arcList = strongest.map((arc) => {
                const left = Number(arc[0]), right = Number(arc[1]);
                const forward = Number(arc[2]), reverse = Number(arc[3]);
                const forwardDominant = forward >= reverse;
                const source = forwardDominant ? left : right, target = forwardDominant ? right : left;
                const label = `${sequence[source] || "N"}${source + 1} → ${sequence[target] || "N"}${target + 1}`;
                return `<li><span>${label}</span><strong>${Math.max(forward, reverse).toFixed(3)} · q${Number(arc[4]) === 2 ? "99.9" : "99"}</strong></li>`;
            }).join("");
            const recall = Number(summary.reference_recall_q999 || 0) * 100;
            slot.innerHTML = `<section class="rfx-analysis-section rfx-dependency-section"><div class="rfx-analysis-heading"><div><p>03 · Saturation mutagenesis</p><h4>Directed nucleotide-dependency architecture</h4></div><span>3 substitutions per source · all targets</span></div><p class="rfx-analysis-note">Each source nucleotide was replaced by all three alternatives. Edge weight is the maximum absolute change in target log₂ odds over the three substitutions and four target bases. Thresholds are matched to dinucleotide-preserving null RNAs by length, GC content and sequence separation.</p><div class="rfx-dependency-metrics"><div><span>q99 edges</span><strong>${Number(summary.q99_pairs || 0).toLocaleString()}</strong></div><div><span>q99.9 edges</span><strong>${Number(summary.q999_pairs || 0).toLocaleString()}</strong></div><div><span>Rfam-pair recall</span><strong>${recall.toFixed(1)}%</strong></div><div><span>Displayed arcs</span><strong>${arcs.length.toLocaleString()}</strong></div></div><section class="rfx-convergence" aria-label="Evidence convergence"><div class="rfx-convergence-head"><div><p>Evidence convergence</p><h5>Pairing signal × perturbation response</h5></div><span>${comparisonLength.toLocaleString()}-nt shared window</span></div><p class="rfx-convergence-note">Position-matched comparison of symmetric zero-shot pairs, q99.9 substitution dependencies and Rfam consensus pairs. Direction is ignored only for this comparison.</p><div class="rfx-convergence-metrics"><div><span>Zero-shot ∩ q99.9</span><strong>${zeroDependencyRate.toFixed(1)}%</strong><small>${zeroDependency.length.toLocaleString()} of ${zeroShotPairs.size.toLocaleString()} pairs</small></div><div><span>Rfam recovered by both</span><strong>${referenceBothRate.toFixed(1)}%</strong><small>${referenceBoth.length.toLocaleString()} of ${referencePairs.size.toLocaleString()} pairs</small></div><div><span>Rfam recovered by either</span><strong>${referenceEitherRate.toFixed(1)}%</strong><small>${referenceEither.length.toLocaleString()} of ${referencePairs.size.toLocaleString()} pairs</small></div></div><p class="rfx-convergence-boundary">Dependency arcs are capped at 250 per family; overlap values are conservative lower bounds.</p></section><div class="rfx-dependency-key"><span><i class="is-forward"></i>5′ → 3′ dependency</span><span><i class="is-reverse"></i>3′ → 5′ dependency</span><span><i class="is-reference"></i>Rfam consensus pair</span></div><div class="rfx-chart-frame rfx-dependency-frame"><canvas data-rfx-dependency aria-label="Directed substitution-dependency arc map"></canvas></div>${arcList ? `<ol class="rfx-pair-list rfx-dependency-list">${arcList}</ol>` : ""}<p class="rfx-method-boundary">These arcs report dependencies inside the frozen model. They are hypotheses about learned sequence organization, not demonstrated molecular contacts.</p></section>`;
            requestAnimationFrame(() => this.drawDependencyArcs(
                slot.querySelector("[data-rfx-dependency]"), dependency
            ));
        }

        drawDependencyArcs(canvas, dependency) {
            if (!canvas) return;
            const frame = canvas.parentElement;
            const width = Math.max(280, Math.round(frame.getBoundingClientRect().width));
            const height = 360, dpr = Math.min(window.devicePixelRatio || 1, 2);
            frame.style.height = `${height}px`;
            canvas.style.width = `${width}px`; canvas.style.height = `${height}px`;
            canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
            const ctx = canvas.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            ctx.clearRect(0, 0, width, height); ctx.fillStyle = "#fffdf9"; ctx.fillRect(0, 0, width, height);
            const length = Math.max(1, Number(dependency.length) || String(dependency.sequence || "").length);
            const leftMargin = 24, rightMargin = 16, baseline = 178;
            const plotWidth = width - leftMargin - rightMargin;
            const positionX = (position) => leftMargin + Number(position) / Math.max(1, length - 1) * plotWidth;
            const references = Array.isArray(dependency.reference_pairs) ? dependency.reference_pairs : [];
            ctx.save(); ctx.strokeStyle = "rgba(116, 91, 61, .24)"; ctx.lineWidth = .8; ctx.setLineDash([2, 3]);
            for (const pair of references) {
                const x1 = positionX(pair[0]), x2 = positionX(pair[1]);
                const amplitude = Math.min(118, 12 + Math.sqrt(Math.abs(x2 - x1) / Math.max(1, plotWidth)) * 125);
                ctx.beginPath(); ctx.moveTo(x1, baseline); ctx.quadraticCurveTo((x1 + x2) / 2, baseline - amplitude, x2, baseline); ctx.stroke();
            }
            ctx.restore();

            const distanceLabel = (distance) => {
                if (distance <= 3) return "2-3";
                if (distance <= 7) return "4-7";
                if (distance <= 15) return "8-15";
                if (distance <= 31) return "16-31";
                if (distance <= 63) return "32-63";
                if (distance <= 127) return "64-127";
                if (distance <= 255) return "128-255";
                if (distance <= 512) return "256-512";
                return "513-1024";
            };
            const directed = [];
            for (const arc of (dependency.arcs || [])) {
                const lower = Number(arc[0]), upper = Number(arc[1]);
                const thresholds = dependency.thresholds_by_separation?.[distanceLabel(upper - lower)] || {};
                const q99 = Number(thresholds.q99 || 0), q999 = Number(thresholds.q999 || Infinity);
                const forward = Number(arc[2]), reverse = Number(arc[3]);
                if (forward >= q99) directed.push({ source: lower, target: upper, score: forward, tier: forward >= q999 ? 2 : 1 });
                if (reverse >= q99) directed.push({ source: upper, target: lower, score: reverse, tier: reverse >= q999 ? 2 : 1 });
            }
            directed.sort((left, right) => left.score - right.score);
            const maximum = Math.max(1e-8, ...directed.map((edge) => edge.score));
            for (const edge of directed) {
                const x1 = positionX(edge.source), x2 = positionX(edge.target), forward = edge.source < edge.target;
                const span = Math.abs(x2 - x1);
                const amplitude = Math.min(128, 13 + Math.sqrt(span / Math.max(1, plotWidth)) * 132);
                const controlY = baseline + (forward ? -amplitude : amplitude);
                const strength = Math.sqrt(edge.score / maximum);
                ctx.strokeStyle = forward
                    ? `rgba(42, 104, 151, ${(.15 + .68 * strength).toFixed(3)})`
                    : `rgba(36, 126, 104, ${(.15 + .68 * strength).toFixed(3)})`;
                ctx.lineWidth = edge.tier === 2 ? 1.35 : .75;
                ctx.beginPath(); ctx.moveTo(x1, baseline); ctx.quadraticCurveTo((x1 + x2) / 2, controlY, x2, baseline); ctx.stroke();
            }

            ctx.strokeStyle = "#526575"; ctx.lineWidth = 1; ctx.beginPath();
            ctx.moveTo(leftMargin, baseline + .5); ctx.lineTo(width - rightMargin, baseline + .5); ctx.stroke();
            const sequence = String(dependency.sequence || "").toUpperCase().replaceAll("T", "U");
            if (length <= 80 && plotWidth / length >= 5.4) {
                ctx.font = '600 8px "IBM Plex Sans", "Segoe UI", sans-serif'; ctx.textAlign = "center"; ctx.textBaseline = "middle";
                for (let index = 0; index < length; index += 1) {
                    ctx.fillStyle = BASE_COLORS[sequence[index]] || BASE_COLORS.N;
                    ctx.fillText(sequence[index] || "N", positionX(index), baseline + 12);
                }
            }
            const ticks = oneBasedPositionTicks(length, Math.max(5, Math.floor(width / 70)));
            ctx.fillStyle = "#516579"; ctx.strokeStyle = "#718294"; ctx.font = '10px "IBM Plex Sans", "Segoe UI", sans-serif';
            ctx.textAlign = "center"; ctx.textBaseline = "top";
            for (const tick of ticks) {
                const x = positionX(tick.position);
                ctx.beginPath(); ctx.moveTo(x, baseline - 3); ctx.lineTo(x, baseline + 4); ctx.stroke();
                ctx.fillText(String(tick.label), x, baseline + (length <= 80 ? 22 : 7));
            }
            ctx.fillStyle = "#0a2540"; ctx.font = '11px "IBM Plex Sans", "Segoe UI", sans-serif';
            ctx.fillText("RNA position (5′ → 3′)", leftMargin + plotWidth / 2, height - 17);
        }

        drawAttention(canvas, attention) {
            const rows = Number(attention.shape[0]), columns = Number(attention.shape[1]);
            const frame = canvas.parentElement, frameRect = frame.getBoundingClientRect();
            const width = Math.max(280, Math.round(frameRect.width));
            const height = Math.max(300, Math.min(410, Math.round(width * .94)));
            const dpr = Math.min(window.devicePixelRatio || 1, 2);
            frame.style.height = `${height}px`;
            canvas.style.width = `${width}px`; canvas.style.height = `${height}px`;
            canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
            const ctx = canvas.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            ctx.clearRect(0, 0, width, height);
            const margins = { left: 43, right: 51, top: 13, bottom: 47 };
            const plotSize = Math.max(160, Math.min(width - margins.left - margins.right, height - margins.top - margins.bottom));
            const plotX = margins.left + Math.max(0, (width - margins.left - margins.right - plotSize) / 2);
            const plotY = margins.top;
            ctx.fillStyle = "#000004"; ctx.fillRect(plotX, plotY, plotSize, plotSize);
            const lookup = new Map();
            let maximum = Number(attention.scale_max) || 0;
            if (attention.encoding === "sparse-upper-threshold") {
                maximum = (attention.pairs || []).reduce((value, pair) => Math.max(value, Number(pair[2])), maximum);
            }
            maximum = Math.max(maximum, 1e-8);
            const cellWidth = plotSize / columns, cellHeight = plotSize / rows;
            const paint = (row, column, score) => {
                const color = heatColor(Number(score) / maximum);
                ctx.fillStyle = `rgb(${color[0]} ${color[1]} ${color[2]})`;
                ctx.fillRect(plotX + column * cellWidth, plotY + row * cellHeight,
                    Math.max(.8, cellWidth + .15), Math.max(.8, cellHeight + .15));
                lookup.set(`${row}:${column}`, Number(score));
            };
            if (attention.encoding === "sparse-upper-threshold") {
                for (const pair of (attention.pairs || [])) {
                    const row = Number(pair[0]) - 1, column = Number(pair[1]) - 1, score = Number(pair[2]);
                    paint(row, column, score); paint(column, row, score);
                }
            } else {
                const binary = atob(attention.values);
                for (let index = 0; index < binary.length; index += 1) {
                    const score = binary.charCodeAt(index) / 255 * maximum;
                    if (score > 0) paint(Math.floor(index / columns), index % columns, score);
                }
            }
            ctx.strokeStyle = "#26313a"; ctx.lineWidth = 1; ctx.strokeRect(plotX - .5, plotY - .5, plotSize + 1, plotSize + 1);
            const ticks = oneBasedPositionTicks(Math.max(rows, columns));
            ctx.fillStyle = "#516579"; ctx.strokeStyle = "#718294"; ctx.lineWidth = 1;
            ctx.font = '10px "IBM Plex Sans", "Segoe UI", sans-serif';
            for (const tick of ticks) {
                const x = plotX + tick.position / Math.max(1, columns - 1) * plotSize;
                const y = plotY + tick.position / Math.max(1, rows - 1) * plotSize;
                ctx.beginPath(); ctx.moveTo(x, plotY + plotSize); ctx.lineTo(x, plotY + plotSize + 4); ctx.stroke();
                ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillText(String(tick.label), x, plotY + plotSize + 6);
                ctx.beginPath(); ctx.moveTo(plotX - 4, y); ctx.lineTo(plotX, y); ctx.stroke();
                ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.fillText(String(tick.label), plotX - 7, y);
            }
            ctx.fillStyle = "#0a2540"; ctx.font = '11px "IBM Plex Sans", "Segoe UI", sans-serif';
            ctx.textAlign = "center"; ctx.textBaseline = "bottom";
            ctx.fillText("RNA position (5′ → 3′)", plotX + plotSize / 2, height - 2);
            ctx.save(); ctx.translate(10, plotY + plotSize / 2); ctx.rotate(-Math.PI / 2);
            ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillText("RNA position (5′ → 3′)", 0, 0); ctx.restore();
            const barX = plotX + plotSize + 13, barY = plotY, barWidth = 10;
            const gradient = ctx.createLinearGradient(0, barY + plotSize, 0, barY);
            HEAT_STOPS.forEach((color, index) => gradient.addColorStop(index / (HEAT_STOPS.length - 1), `rgb(${color[0]} ${color[1]} ${color[2]})`));
            ctx.fillStyle = gradient; ctx.fillRect(barX, barY, barWidth, plotSize);
            ctx.strokeStyle = "#26313a"; ctx.strokeRect(barX - .5, barY - .5, barWidth + 1, plotSize + 1);
            ctx.fillStyle = "#516579"; ctx.font = '10px "IBM Plex Sans", "Segoe UI", sans-serif'; ctx.textAlign = "left";
            ctx.textBaseline = "top"; ctx.fillText(scoreLabel(maximum), barX + barWidth + 4, barY - 2);
            ctx.textBaseline = "bottom"; ctx.fillText("0", barX + barWidth + 4, barY + plotSize + 2);
            const tip = canvas.parentElement.querySelector("[data-rfx-attention-tip]");
            canvas.addEventListener("pointermove", (event) => {
                const rect = canvas.getBoundingClientRect();
                const localX = (event.clientX - rect.left) / rect.width * width;
                const localY = (event.clientY - rect.top) / rect.height * height;
                if (localX < plotX || localX > plotX + plotSize || localY < plotY || localY > plotY + plotSize) {
                    tip.hidden = true; return;
                }
                const column = clamp(Math.floor((localX - plotX) / plotSize * columns), 0, columns - 1);
                const row = clamp(Math.floor((localY - plotY) / plotSize * rows), 0, rows - 1);
                const score = lookup.get(`${row}:${column}`);
                tip.hidden = false;
                tip.textContent = `${row + 1} × ${column + 1} · ${score === undefined ? `< ${Number(attention.pair_threshold || 0).toFixed(1)}` : score.toFixed(3)}`;
            });
            canvas.addEventListener("pointerleave", () => { tip.hidden = true; });
        }

        drawContext(canvas, sequence, sourceValues) {
            const values = Array.from(sourceValues, Number), frame = canvas.parentElement;
            const viewportWidth = Math.max(260, Math.round(frame.getBoundingClientRect().width));
            const left = 39, right = 14, top = 13, bottom = 28, perBase = 9;
            const width = Math.max(viewportWidth, left + right + values.length * perBase), height = 190;
            const dpr = Math.min(window.devicePixelRatio || 1, 2);
            frame.style.height = `${height}px`;
            canvas.style.width = `${width}px`; canvas.style.height = `${height}px`;
            canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
            const ctx = canvas.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            const plotWidth = width - left - right, plotHeight = height - top - bottom, baseline = top + plotHeight;
            const maxValue = Math.max(...values, 1e-8);
            ctx.clearRect(0, 0, width, height); ctx.fillStyle = "#fffdf9"; ctx.fillRect(0, 0, width, height);
            ctx.strokeStyle = "rgba(10, 37, 64, .10)"; ctx.lineWidth = 1;
            for (let line = 0; line <= 2; line += 1) {
                const y = top + plotHeight * line / 2;
                ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(width - right, y); ctx.stroke();
                const value = maxValue * (1 - line / 2);
                ctx.fillStyle = "#516579"; ctx.font = '10px "IBM Plex Sans", "Segoe UI", sans-serif';
                ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.fillText(scoreLabel(value), left - 6, y);
            }
            const cellWidth = plotWidth / Math.max(1, values.length);
            const drawLetter = (base, centerX, targetHeight, color, alpha = 1) => {
                if (targetHeight <= .4) return;
                const fontSize = 100;
                ctx.save(); ctx.font = `900 ${fontSize}px "Arial Black", Arial, sans-serif`; ctx.textAlign = "left"; ctx.textBaseline = "alphabetic";
                const glyphWidth = Math.max(1, ctx.measureText(base).width);
                const targetWidth = Math.max(3, cellWidth * .96);
                ctx.translate(centerX - targetWidth / 2, baseline);
                ctx.scale(targetWidth / glyphWidth, targetHeight / fontSize);
                ctx.fillStyle = color; ctx.globalAlpha = alpha; ctx.fillText(base, 0, 0); ctx.restore();
            };
            values.forEach((value, index) => {
                const base = (sequence[index] || "N").toUpperCase().replace("T", "U");
                const centerX = left + (index + .5) * cellWidth;
                drawLetter(base, centerX, 6, "#cfd5da", .72);
                const normalized = clamp(Number(value) / maxValue, 0, 1);
                drawLetter(base, centerX, normalized * plotHeight, BASE_COLORS[base] || BASE_COLORS.N, .38 + .62 * Math.sqrt(normalized));
            });
            ctx.strokeStyle = "#5d6977"; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(left, baseline + .5); ctx.lineTo(width - right, baseline + .5); ctx.stroke();
            const ticks = oneBasedPositionTicks(values.length, Math.max(6, Math.floor(viewportWidth / 65)));
            ctx.fillStyle = "#516579"; ctx.font = '10px "IBM Plex Sans", "Segoe UI", sans-serif'; ctx.textAlign = "center"; ctx.textBaseline = "top";
            for (const tick of ticks) {
                const x = left + (tick.position + .5) * cellWidth;
                ctx.beginPath(); ctx.moveTo(x, baseline); ctx.lineTo(x, baseline + 4); ctx.stroke();
                ctx.fillText(String(tick.label), x, baseline + 6);
            }
            const tip = canvas.parentElement.querySelector("[data-rfx-context-tip]");
            canvas.addEventListener("pointermove", (event) => {
                const current = canvas.getBoundingClientRect();
                const logicalX = (event.clientX - current.left) / current.width * width;
                const localX = clamp(logicalX - left, 0, plotWidth - .001);
                const index = clamp(Math.floor(localX / plotWidth * values.length), 0, values.length - 1);
                tip.hidden = false;
                tip.textContent = `${sequence[index] || "N"}${index + 1} · ${Number(values[index]).toFixed(3)}`;
            });
            canvas.addEventListener("pointerleave", () => { tip.hidden = true; });
        }
    }

    document.addEventListener("DOMContentLoaded", () => {
        document.querySelectorAll("[data-rfam-explorer]").forEach((root) => new RfamExplorer(root).init());
    });
}());
