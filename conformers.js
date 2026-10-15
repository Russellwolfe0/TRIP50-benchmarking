const STUDY_ROOT = "calculations/conformer_study/conformer_analysis";

function numeric(value) {
    const parsed = Number.parseFloat(value);
    return Number.isFinite(parsed) ? parsed : null;
}

function formatNumber(value, digits = 2) {
    const parsed = numeric(value);
    return parsed === null ? "—" : parsed.toFixed(digits);
}

// Green at 100% match, red at 40% or below.
function matchColor(value) {
    const parsed = numeric(value);
    if (parsed === null) return "transparent";
    const ratio = Math.max(0, Math.min(1, (parsed - 40) / 60));
    return `hsl(${120 * ratio}, 65%, 90%)`;
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

async function loadModelSummary() {
    const path = `${STUDY_ROOT}/model_summary.csv`;
    try {
        const rows = await parseCsv(path);
        new DataTable("#conformer-model-table", {
            data: rows,
            paging: false,
            order: [[1, "desc"]],
            columns: [
                { data: "model" },
                {
                    data: "pct_matches_reference",
                    render: value => formatNumber(value, 1),
                    createdCell: (cell, value) => { cell.style.backgroundColor = matchColor(value); },
                },
                { data: "n_lower_than_reference" },
                { data: "n_higher_than_reference" },
                { data: "n_linear_artifact" },
                { data: "n_not_converged" },
                { data: "n_rearranged" },
                { data: "n_conformers_converged", render: (value, _type, row) => `${value} / ${row.n_conformers}` },
            ],
        });
        document.getElementById("conformer-model-status").textContent =
            `${rows.length} models loaded from ${path}.`;
    } catch (error) {
        showLoadError("conformer-model-status", path, error);
    }
}

async function loadSpeciesConsensus() {
    const path = `${STUDY_ROOT}/species_consensus.csv`;
    try {
        const rows = await parseCsv(path);
        new DataTable("#conformer-species-table", {
            data: rows,
            pageLength: 25,
            order: [[3, "asc"]],
            columns: [
                { data: "species" },
                { data: "multiplicity" },
                { data: "n_conformers" },
                { data: "n_matches_reference" },
                { data: "n_lower_than_reference" },
                { data: "n_higher_than_reference" },
                { data: "min_dE_vs_ref_kcal", render: value => formatNumber(value) },
                { data: "mean_dE_vs_ref_kcal", render: value => formatNumber(value) },
            ],
        });
        document.getElementById("conformer-species-status").textContent =
            `${rows.length} species loaded from ${path}.`;
    } catch (error) {
        showLoadError("conformer-species-status", path, error);
    }
}

document.addEventListener("DOMContentLoaded", () => {
    loadModelSummary();
    loadSpeciesConsensus();
});
