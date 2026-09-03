const DATA_ROOT = "data/final";
const IMAGE_ROOT = "images";
const REACTION_TYPES = ["C-C", "C-O", "C-S", "HAT", "Si-X", "C-Hal", "N-X"];
const SAMPLE_REACTIONS = [1, 13, 20, 22, 31, 38, 42];
let categoryComparisonRows = [];

// Shared handler for the collapsible dataset and model cards in index.html.
// It is declared at top level because the cards call it from inline onclick handlers.
function toggleCard(contentId, arrowId) {
    const content = document.getElementById(contentId);
    const arrow = document.getElementById(arrowId);

    if (!content || !arrow) {
        console.error(`Dropdown target not found: ${contentId} / ${arrowId}`);
        return;
    }

    const isOpen = content.classList.toggle("open");
    arrow.classList.toggle("rotated", isOpen);
    arrow.setAttribute("aria-expanded", String(isOpen));
}

function numeric(value) {
    const parsed = Number.parseFloat(value);
    return Number.isFinite(parsed) ? parsed : null;
}

function formatNumber(value, digits = 3) {
    const parsed = numeric(value);
    return parsed === null ? "—" : parsed.toFixed(digits);
}

function formatTime(value) {
    const parsed = numeric(value);
    if (parsed === null) return "—";
    if (parsed < 0.01) return parsed.toExponential(2);
    if (parsed < 10) return parsed.toFixed(3);
    return parsed.toFixed(2);
}

function maeColor(value, minimum = 0, maximum = 10) {
    const parsed = numeric(value);
    if (parsed === null) return "transparent";
    const ratio = Math.max(0, Math.min(1, (parsed - minimum) / (maximum - minimum)));
    const hue = 120 * (1 - ratio);
    return `hsl(${hue}, 65%, 90%)`;
}

function maeColumn(field) {
    return {
        data: field,
        render: value => formatNumber(value),
        createdCell: (cell, value) => {
            cell.style.backgroundColor = maeColor(value);
        },
    };
}

function modelSlug(model) {
    return model.toLowerCase().replace("ω", "w").replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
}

function showParityPlot(plotType, model, button) {
    const image = document.getElementById(`${plotType}-parity-image`);
    const viewer = document.querySelector(`[data-parity-type="${plotType}"]`);
    image.src = `${IMAGE_ROOT}/parity_plots/${modelSlug(model)}_${plotType}_parity.png`;
    image.alt = `${model} ${plotType === "thermo" ? "thermodynamic" : "kinetic"} parity plot`;
    viewer.querySelectorAll(".parity-tab").forEach(tab => {
        tab.classList.toggle("active", tab === button);
    });
}

function showRadarPlot(model) {
    const image = document.getElementById("radar-plot-image");
    image.src = `${IMAGE_ROOT}/radar_plots/radar_${modelSlug(model)}.png`;
    image.alt = `${model} reaction-category MAE radar plot`;
}

function initializeModelPlotControls(rows) {
    const models = rows.map(row => row.model);
    const initialModel = models.includes("B2GP-PLYP-D4") ? "B2GP-PLYP-D4" : models[0];

    document.querySelectorAll(".parity-viewer").forEach(viewer => {
        const plotType = viewer.dataset.parityType;
        const tabs = viewer.querySelector(".parity-tabs");
        tabs.replaceChildren();
        models.forEach(model => {
            const button = document.createElement("button");
            button.type = "button";
            button.className = `parity-tab${model === initialModel ? " active" : ""}`;
            button.textContent = model;
            button.addEventListener("click", () => showParityPlot(plotType, model, button));
            tabs.append(button);
        });
    });

    const select = document.getElementById("radar-model-select");
    select.replaceChildren();
    models.forEach(model => {
        const option = document.createElement("option");
        option.value = model;
        option.textContent = model;
        select.append(option);
    });
    select.value = initialModel;
    select.addEventListener("change", () => showRadarPlot(select.value));
    showRadarPlot(initialModel);
}

function initializeCoordinateViewer() {
    const tabs = document.querySelector(".coordinate-tabs");
    SAMPLE_REACTIONS.forEach((reaction, index) => {
        const button = document.createElement("button");
        button.type = "button";
        button.className = `coordinate-tab${index === 0 ? " active" : ""}`;
        button.textContent = `Rxn ${reaction}`;
        button.addEventListener("click", () => {
            document.getElementById("reaction-coordinate-image").src =
                `${IMAGE_ROOT}/reaction_coordinates/rxn${reaction}/method_comparison.png`;
            document.getElementById("reaction-coordinate-image").alt =
                `Reaction ${reaction} coordinate diagram comparing methods`;
            tabs.querySelectorAll("button").forEach(tab => tab.classList.toggle("active", tab === button));
        });
        tabs.append(button);
    });
}

function svgElement(name, attributes = {}, text = "") {
    const element = document.createElementNS("http://www.w3.org/2000/svg", name);
    Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
    if (text) element.textContent = text;
    return element;
}

function radarPoint(index, value, maximum) {
    const angle = -Math.PI / 2 + 2 * Math.PI * index / REACTION_TYPES.length;
    const radius = 175 * value / maximum;
    return [310 + radius * Math.cos(angle), 320 + radius * Math.sin(angle)];
}

function updateCategoryComparison() {
    if (!categoryComparisonRows.length) return;
    const names = [
        document.getElementById("comparison-model-one").value,
        document.getElementById("comparison-model-two").value,
    ];
    const selection = document.getElementById("comparison-property").value;
    const properties = selection === "both" ? ["thermo", "kinetic"] : [selection];
    const colors = ["#1e4d2b", "#2e86de", "#b45f06", "#7d3c98"];
    const series = [];
    names.forEach((name, modelIndex) => {
        const row = categoryComparisonRows.find(item => item.model === name);
        properties.forEach((property, propertyIndex) => series.push({
            label: `${name} — ${property === "thermo" ? "Thermodynamic" : "Kinetic"}`,
            values: REACTION_TYPES.map(category => Number(row[property][category])),
            color: colors[modelIndex * properties.length + propertyIndex],
        }));
    });
    const maximum = Math.max(5, Math.ceil(Math.max(...series.flatMap(item => item.values)) / 5) * 5);
    const svg = document.getElementById("model-category-radar");
    svg.replaceChildren();
    for (let ring = 1; ring <= 5; ring += 1) {
        const value = maximum * ring / 5;
        const points = REACTION_TYPES.map((_, index) => radarPoint(index, value, maximum));
        svg.append(svgElement("polygon", {
            points: points.map(point => point.join(",")).join(" "), fill: "none",
            stroke: "#d4d4ca", "stroke-width": 1,
        }));
        svg.append(svgElement("text", {x: 318, y: 324 - 175 * ring / 5, fill: "#555", "font-size": 12}, value.toFixed(0)));
    }
    REACTION_TYPES.forEach((category, index) => {
        const [x, y] = radarPoint(index, maximum, maximum);
        const [lx, ly] = radarPoint(index, maximum * 1.14, maximum);
        svg.append(svgElement("line", {x1: 310, y1: 320, x2: x, y2: y, stroke: "#d4d4ca"}));
        svg.append(svgElement("text", {x: lx, y: ly + 5, "text-anchor": "middle", "font-size": 16}, category));
    });
    series.forEach(item => {
        const points = item.values.map((value, index) => radarPoint(index, value, maximum));
        svg.append(svgElement("polygon", {
            points: points.map(point => point.join(",")).join(" "), fill: item.color,
            "fill-opacity": 0.08, stroke: item.color, "stroke-width": 3,
        }));
    });
    series.forEach((item, index) => {
        const y = 72 + index * 27;
        svg.append(svgElement("line", {x1: 550, y1: y, x2: 578, y2: y, stroke: item.color, "stroke-width": 4}));
        svg.append(svgElement("text", {x: 589, y: y + 5, "font-size": 14}, item.label));
    });
    document.getElementById("model-category-status").textContent =
        `Category MAE comparison; radial maximum ${maximum} kcal/mol.`;
}

function initializeCategoryComparison(rows) {
    categoryComparisonRows = rows;
    const models = rows.map(row => row.model).sort((a, b) => a.localeCompare(b));
    const selects = [document.getElementById("comparison-model-one"), document.getElementById("comparison-model-two")];
    selects.forEach(select => models.forEach(model => {
        const option = document.createElement("option"); option.value = model; option.textContent = model; select.append(option);
    }));
    selects[0].value = models.includes("B2GP-PLYP-D4") ? "B2GP-PLYP-D4" : models[0];
    selects[1].value = models.includes("esen-md-direct-all-omol") ? "esen-md-direct-all-omol" : models[1];
    [...selects, document.getElementById("comparison-property")].forEach(control => control.addEventListener("change", updateCategoryComparison));
    updateCategoryComparison();
}

function parseCsv(path) {
    return new Promise((resolve, reject) => {
        Papa.parse(path, {
            download: true,
            header: true,
            skipEmptyLines: true,
            complete: results => {
                if (results.errors.length) {
                    reject(new Error(results.errors[0].message));
                    return;
                }
                resolve(results.data);
            },
            error: reject,
        });
    });
}

function showLoadError(statusId, path, error) {
    const status = document.getElementById(statusId);
    status.classList.add("error");
    status.textContent = `Could not load ${path}. Serve TRIP50 from its repository root (for example: python3 -m http.server). ${error.message}`;
}

async function loadFinalResults() {
    const path = `${DATA_ROOT}/final_data_table.csv`;
    try {
        const rows = await parseCsv(path);
        new DataTable("#results-table", {
            data: rows,
            pageLength: 25,
            order: [[5, "asc"]],
            columns: [
                { data: "model" },
                { data: "model_category" },
                { data: "average_run_time_seconds", render: formatTime },
                maeColumn("mae_thermo_kcal_mol"),
                maeColumn("mae_kinetic_kcal_mol"),
                maeColumn("combined_mae_kcal_mol"),
            ],
        });
        initializeModelPlotControls(rows);
        document.getElementById("results-status").textContent = `${rows.length} models loaded from the final benchmark table.`;
    } catch (error) {
        showLoadError("results-status", path, error);
    }
}

function pivotCategoryRows(rows, finalRows) {
    const byModel = new Map();
    rows.forEach(row => {
        if (!byModel.has(row.model)) {
            byModel.set(row.model, {
                model: row.model,
                model_category: row.model_category,
                thermo: {},
                kinetic: {},
            });
        }
        const model = byModel.get(row.model);
        model.thermo[row.reaction_type] = row.mae_thermo_kcal_mol;
        model.kinetic[row.reaction_type] = row.mae_kinetic_kcal_mol;
    });

    const overall = new Map(finalRows.map(row => [row.model, row]));
    return [...byModel.values()].map(row => ({
        ...row,
        overall_thermo: overall.get(row.model)?.mae_thermo_kcal_mol,
        overall_kinetic: overall.get(row.model)?.mae_kinetic_kcal_mol,
    }));
}

function categoryColumns(property) {
    return [
        { data: "model" },
        { data: "model_category" },
        maeColumn(`overall_${property}`),
        ...REACTION_TYPES.map(category => ({
            data: `${property}.${category}`,
            render: (_value, _type, row) => formatNumber(row[property][category]),
            createdCell: (cell, _value, row) => {
                cell.style.backgroundColor = maeColor(row[property][category]);
            },
        })),
    ];
}

async function loadCategoryTables() {
    const categoryPath = `${DATA_ROOT}/model_results_by_reaction_type.csv`;
    const finalPath = `${DATA_ROOT}/final_data_table.csv`;
    const summaryPath = `${DATA_ROOT}/reaction_type_summary.csv`;
    try {
        const [categoryRows, finalRows, summaryRows] = await Promise.all([
            parseCsv(categoryPath),
            parseCsv(finalPath),
            parseCsv(summaryPath),
        ]);
        const pivoted = pivotCategoryRows(categoryRows, finalRows);
        initializeCategoryComparison(pivoted);

        new DataTable("#thermo-table", {
            data: pivoted,
            pageLength: 25,
            order: [[2, "asc"]],
            columns: categoryColumns("thermo"),
        });
        new DataTable("#kinetic-table", {
            data: pivoted,
            pageLength: 25,
            order: [[2, "asc"]],
            columns: categoryColumns("kinetic"),
        });
        new DataTable("#category-summary-table", {
            data: summaryRows,
            paging: false,
            searching: false,
            info: false,
            order: [[0, "asc"]],
            columns: [
                { data: "reaction_type" },
                maeColumn("mae_thermo_kcal_mol"),
                maeColumn("mae_kinetic_kcal_mol"),
                { data: "average_run_time_seconds", render: formatTime },
            ],
        });
        document.getElementById("category-status").textContent =
            `${categoryRows.length} model/category combinations loaded across ${summaryRows.length} reaction types.`;
    } catch (error) {
        showLoadError("category-status", "the final reaction-category tables", error);
    }
}

function categorySlug(category) {
    return category.toLowerCase().replaceAll("-", "_");
}

function currentParetoMode() {
    return document.querySelector('input[name="pareto-mode"]:checked').value;
}

function updateParetoImage(measure, category) {
    const mode = currentParetoMode();
    const directory = mode === "absolute" ? "Paretofrontsbycat_absolute" : "Paretofrontsbycat";
    const image = document.getElementById(`${measure}-pareto-image`);
    image.src = `${IMAGE_ROOT}/${directory}/pareto_${measure}_${categorySlug(category)}.png`;
    image.alt = `${measure === "thermo" ? "Thermodynamic" : "Kinetic"} ${mode} MAE Pareto front for ${category} reactions`;
    image.dataset.category = category;
}

function initializeParetoViewers() {
    document.querySelectorAll(".pareto-viewer").forEach(viewer => {
        const measure = viewer.dataset.paretoType;
        const tabs = viewer.querySelector(".pareto-tabs");
        REACTION_TYPES.forEach((category, index) => {
            const button = document.createElement("button");
            button.type = "button";
            button.className = `pareto-tab${index === 0 ? " active" : ""}`;
            button.textContent = category;
            button.addEventListener("click", () => {
                tabs.querySelectorAll("button").forEach(tab => tab.classList.remove("active"));
                button.classList.add("active");
                updateParetoImage(measure, category);
            });
            tabs.append(button);
        });
    });

    document.querySelectorAll('input[name="pareto-mode"]').forEach(input => {
        input.addEventListener("change", () => {
            ["thermo", "kinetic"].forEach(measure => {
                const image = document.getElementById(`${measure}-pareto-image`);
                updateParetoImage(measure, image.dataset.category || "C-C");
            });
        });
    });
}

document.addEventListener("DOMContentLoaded", () => {
    initializeParetoViewers();
    initializeCoordinateViewer();
    loadFinalResults();
    loadCategoryTables();
});
