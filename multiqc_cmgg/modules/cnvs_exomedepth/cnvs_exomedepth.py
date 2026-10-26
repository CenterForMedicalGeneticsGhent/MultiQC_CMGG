import html
import json
import logging
import re
from collections import defaultdict
from typing import Any, Dict, List, Tuple, Union

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table

log = logging.getLogger(__name__)

Value = Union[int, float, str]


def parse_tsv(content: str) -> Tuple[List[str], List[Dict[str, str]]]:
    """Parse a tab-separated file with a header row into column names and rows."""
    lines = content.splitlines()
    if not lines:
        return [], []
    header_cols = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        values = line.split("\t")
        rows.append(dict(zip(header_cols, values)))
    return header_cols, rows


def convert_value(value: str) -> Value:
    """Convert a raw string cell value to int/float where possible."""
    if value in ("NA", ""):
        return value
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def build_data(
    rows: List[Dict[str, str]], key_col: str
) -> Dict[str, Dict[str, Value]]:
    """Turn parsed rows into a dict keyed by a (de-duplicated) row identifier."""
    data: Dict[str, Dict[str, Value]] = {}
    counts: Dict[str, int] = defaultdict(int)
    for row in rows:
        raw_key = row.get(key_col, "unknown")
        counts[raw_key] += 1
        key = raw_key if counts[raw_key] == 1 else f"{raw_key} ({counts[raw_key]})"
        data[key] = {
            col: convert_value(val) for col, val in row.items() if col != key_col
        }
    return data


def natural_sort_key(s: str) -> List[Union[int, str]]:
    """Key for natural sorting of alphanumeric strings (e.g. patient2 before patient10)."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r"(\d+)", s)]


def build_patient_summary(
    panels_rows: List[Dict[str, str]], designs_rows: List[Dict[str, str]]
) -> Dict[str, Dict[str, Any]]:
    """Aggregate CNV counts per patient across panel and design summary files."""
    summary: Dict[str, Dict[str, Any]] = {}

    for r in panels_rows:
        p = r.get("Patient", "").strip()
        if not p:
            continue
        if p not in summary:
            summary[p] = {"panel_count": 0, "design_count": 0}
        summary[p]["panel_count"] += 1

    for r in designs_rows:
        p = r.get("Patient", "").strip()
        if not p:
            continue
        if p not in summary:
            summary[p] = {"panel_count": 0, "design_count": 0}
        summary[p]["design_count"] += 1

    return dict(sorted(summary.items(), key=lambda item: natural_sort_key(item[0])))


def build_patient_filter_html(patient_summary: Dict[str, Dict[str, Any]]) -> str:
    """Build the interactive Patient Filter UI card, styling, and JavaScript logic."""
    if not patient_summary:
        return ""

    total_patients = len(patient_summary)
    meta_json = json.dumps(
        {
            p: {"p": s["panel_count"], "d": s["design_count"]}
            for p, s in patient_summary.items()
        }
    )

    options = [f'<option value="">-- Toon alle stalen ({total_patients}) --</option>']
    pills = []

    for p, s in patient_summary.items():
        p_cnt = s["panel_count"]
        d_cnt = s["design_count"]

        parts = []
        if p_cnt > 0:
            parts.append(f"{p_cnt} panel")
        if d_cnt > 0:
            parts.append(f"{d_cnt} design")
        desc = ", ".join(parts)
        options.append(
            f'<option value="{html.escape(p, quote=True)}">{html.escape(p)} ({desc})</option>'
        )

        pills.append(
            f'<button type="button" class="btn btn-sm btn-outline-primary ed-patient-pill" '
            f'data-patient="{html.escape(p, quote=True)}" title="{html.escape(p)}: {desc}">'
            f'{html.escape(p)}</button>'
        )

    options_html = "\n        ".join(options)
    pills_html = "\n      ".join(pills)

    return f"""
<div class="card ed-patient-filter-card mb-3 shadow-sm">
  <div class="card-header bg-body-tertiary d-flex justify-content-between align-items-center py-2 flex-wrap gap-2">
    <div>
      <strong class="me-2">Patiënten Filter</strong>
      <span class="text-muted small">({total_patients} stalen met CNVs)</span>
    </div>
    <div class="d-flex align-items-center gap-2">
      <select id="ed-patient-select" class="form-select form-select-sm" style="width: auto; min-width: 200px;">
        {options_html}
      </select>
      <button type="button" class="btn btn-outline-secondary btn-sm" id="ed-reset-filter-btn" style="display: none;">
        Toon alle stalen
      </button>
    </div>
  </div>
  <div class="card-body py-2">
    <div id="ed-active-filter-alert" class="alert alert-primary py-2 px-3 mb-2 d-flex justify-content-between align-items-center" style="display: none;">
      <div>
        <span class="badge bg-primary me-2">Gefilterd staal</span>
        <strong id="ed-active-patient-name" class="fs-6"></strong>
        <span id="ed-active-patient-counts" class="ms-2 text-muted small"></span>
      </div>
      <button type="button" class="btn btn-sm btn-outline-primary" id="ed-clear-filter-btn">
        ✕ Toon alle stalen
      </button>
    </div>
    <div class="small text-muted mb-1">Klik op een patiëntnaam om de tabellen te filteren:</div>
    <div class="ed-patient-pills-container d-flex flex-wrap gap-1 p-1" style="max-height: 180px; overflow-y: auto;">
      {pills_html}
    </div>
  </div>
</div>

<style>
.ed-filtered-out {{
    display: none !important;
}}
.ed-patient-pill {{
    font-size: 0.82rem;
    padding: 0.2rem 0.5rem;
    transition: all 0.15s ease-in-out;
}}
.ed-patient-pill:hover {{
    transform: translateY(-1px);
}}
.ed-patient-filter-card {{
    border: 1px solid var(--bs-border-color, #dee2e6);
    border-radius: 0.375rem;
}}
#cnvs_exomedepth_panels_table_table tbody th.rowheader,
#cnvs_exomedepth_designs_table_table tbody th.rowheader {{
    cursor: pointer;
    text-decoration: underline dotted;
    text-underline-offset: 3px;
}}
#cnvs_exomedepth_panels_table_table tbody th.rowheader:hover,
#cnvs_exomedepth_designs_table_table tbody th.rowheader:hover {{
    color: #0d6efd;
    background-color: rgba(13, 110, 253, 0.08);
}}
</style>

<script>
window.edPatientMeta = {meta_json};

function getEdRowPatient(tr) {{
    var key = tr.getAttribute('data-sample-group');
    if (!key) {{
        var th = tr.querySelector('th.rowheader');
        key = th ? (th.getAttribute('data-sorting-val') || th.textContent.trim()) : '';
    }}
    return key.replace(/ \\(\\d+\\)$/, '').trim();
}}

function setExomeDepthPatientFilter(selectedPatient) {{
    selectedPatient = (selectedPatient || '').trim();

    // Dropdown update
    $('#ed-patient-select').val(selectedPatient);

    // Pill active states
    $('.ed-patient-pill').each(function() {{
        var p = $(this).attr('data-patient');
        if (selectedPatient && p === selectedPatient) {{
            $(this).removeClass('btn-outline-primary').addClass('btn-primary active');
        }} else {{
            $(this).removeClass('btn-primary active').addClass('btn-outline-primary');
        }}
    }});

    // Banner and reset button
    if (selectedPatient) {{
        $('#ed-active-patient-name').text(selectedPatient);
        var meta = window.edPatientMeta && window.edPatientMeta[selectedPatient];
        var countStr = '';
        if (meta) {{
            countStr = '(' + meta.p + ' panel, ' + meta.d + ' design)';
        }}
        $('#ed-active-patient-counts').text(countStr);
        $('#ed-active-filter-alert').show();
        $('#ed-reset-filter-btn').show();
    }} else {{
        $('#ed-active-filter-alert').hide();
        $('#ed-reset-filter-btn').hide();
    }}

    // Filter Panel CNVs table
    var visiblePanels = 0;
    var totalPanels = 0;
    $('#cnvs_exomedepth_panels_table_table tbody tr:not(.ed-no-data-row)').each(function() {{
        totalPanels++;
        var rowPatient = getEdRowPatient(this);
        if (!selectedPatient || rowPatient === selectedPatient) {{
            $(this).removeClass('ed-filtered-out');
            visiblePanels++;
        }} else {{
            $(this).addClass('ed-filtered-out');
        }}
    }});

    $('#cnvs_exomedepth_panels_no_data').remove();
    if (selectedPatient && visiblePanels === 0 && totalPanels > 0) {{
        var cols = $('#cnvs_exomedepth_panels_table_table thead th').length || 15;
        $('#cnvs_exomedepth_panels_table_table tbody').append(
            '<tr id="cnvs_exomedepth_panels_no_data" class="ed-no-data-row"><td colspan="' + cols + '" class="text-center text-muted fst-italic py-3">' +
            'Geen Panel CNVs gevonden voor staal <strong>' + selectedPatient + '</strong>.</td></tr>'
        );
    }}

    if ($('#cnvs_exomedepth_panels_table_table_numrows').length) {{
        $('#cnvs_exomedepth_panels_table_table_numrows').text(selectedPatient ? visiblePanels : totalPanels);
    }}

    // Filter Design CNVs table
    var visibleDesigns = 0;
    var totalDesigns = 0;
    $('#cnvs_exomedepth_designs_table_table tbody tr:not(.ed-no-data-row)').each(function() {{
        totalDesigns++;
        var rowPatient = getEdRowPatient(this);
        if (!selectedPatient || rowPatient === selectedPatient) {{
            $(this).removeClass('ed-filtered-out');
            visibleDesigns++;
        }} else {{
            $(this).addClass('ed-filtered-out');
        }}
    }});

    $('#cnvs_exomedepth_designs_no_data').remove();
    if (selectedPatient && visibleDesigns === 0 && totalDesigns > 0) {{
        var cols = $('#cnvs_exomedepth_designs_table_table thead th').length || 15;
        $('#cnvs_exomedepth_designs_table_table tbody').append(
            '<tr id="cnvs_exomedepth_designs_no_data" class="ed-no-data-row"><td colspan="' + cols + '" class="text-center text-muted fst-italic py-3">' +
            'Geen Design CNVs gevonden voor staal <strong>' + selectedPatient + '</strong>.</td></tr>'
        );
    }}

    if ($('#cnvs_exomedepth_designs_table_table_numrows').length) {{
        $('#cnvs_exomedepth_designs_table_table_numrows').text(selectedPatient ? visibleDesigns : totalDesigns);
    }}

    // Update Design CNVs summary text
    var $summary = $('#cnvs_exomedepth_designs_summary');
    if ($summary.length) {{
        if (selectedPatient) {{
            $summary.html('Click to show the Design CNVs table <span class="badge bg-secondary ms-2">' + visibleDesigns + ' CNV' + (visibleDesigns === 1 ? '' : 's') + ' voor ' + selectedPatient + '</span>');
        }} else {{
            $summary.text('Click to show the Design CNVs table');
        }}
    }}
}}

$(function() {{
    // Attach tooltips to table headers
    setTimeout(function() {{
        $('#cnvs_exomedepth_panels_table_table tbody th.rowheader, #cnvs_exomedepth_designs_table_table tbody th.rowheader')
            .attr('title', 'Klik om te filteren op dit staal');
    }}, 500);

    // Pill click handler (toggle if already active)
    $(document).on('click', '.ed-patient-pill', function(e) {{
        e.preventDefault();
        var p = $(this).attr('data-patient');
        var current = $('#ed-patient-select').val();
        if (current === p) {{
            setExomeDepthPatientFilter('');
        }} else {{
            setExomeDepthPatientFilter(p);
        }}
    }});

    // Dropdown change handler
    $(document).on('change', '#ed-patient-select', function() {{
        setExomeDepthPatientFilter($(this).val());
    }});

    // Reset buttons
    $(document).on('click', '#ed-clear-filter-btn, #ed-reset-filter-btn', function(e) {{
        e.preventDefault();
        setExomeDepthPatientFilter('');
    }});

    // Rowheader click in tables
    $(document).on('click', '#cnvs_exomedepth_panels_table_table tbody th.rowheader, #cnvs_exomedepth_designs_table_table tbody th.rowheader', function() {{
        var rawText = $(this).attr('data-sorting-val') || $(this).text().trim();
        var p = rawText.replace(/ \\(\\d+\\)$/, '').trim();
        if (p) {{
            var current = $('#ed-patient-select').val();
            if (current === p) {{
                setExomeDepthPatientFilter('');
            }} else {{
                setExomeDepthPatientFilter(p);
            }}
        }}
    }});

    // Ensure filter persists after table column sorting
    $(document).on('click', '#cnvs_exomedepth_panels_table_table thead th, #cnvs_exomedepth_designs_table_table thead th', function() {{
        setTimeout(function() {{
            var active = $('#ed-patient-select').val();
            if (active) {{
                setExomeDepthPatientFilter(active);
            }}
        }}, 50);
    }});
}});
</script>
"""


# Color legend shown above the Reads.ratio tables (hover tooltips don't work here: MultiQC's
# table CSS renders cell values on a negative z-index layer behind the <td>, so the mouse
# always hits the <td> and never the inner text where a title attribute would live).
READS_RATIO_LEGEND = (
    "Reads.ratio colour legend - predicted cnv  type: "
    '<span class="badge" style="background-color:rgb(128,0,0)">&lt; 0.2 homozygote deletie</span> '
    '<span class="badge" style="background-color:rgb(255,0,0)">0.2\u20131 heterozygote deletie</span> '
    '<span class="badge" style="background-color:rgb(0,0,255)">1\u20131.83 duplicatie</span> '
    '<span class="badge" style="background-color:rgb(0,0,128)">&gt;= 1.83 triplicatie</span>'
)


def build_headers(data: Dict[str, Dict[str, Value]], columns: List[str]) -> Dict:
    """Build MultiQC table headers, formatting columns containing floats with 2 decimals."""
    headers = {}
    for col in columns:
        col_values = [row[col] for row in data.values() if col in row]
        is_float = any(isinstance(v, float) for v in col_values)
        header = {"title": col}
        if is_float:
            header["format"] = "{:,.2f}"
        if col == "Reads.ratio":
            header["cond_formatting_rules"] = {
                "homozygote_deletie": [{"lt": 0.2}],
                "heterozygote_deletie": [{"ge": 0.2}],
                "duplicatie": [{"ge": 1}],
                "triplicatie": [{"ge": 1.83}],
            }
            # Order matters: later entries override earlier ones for values matching multiple rules
            header["cond_formatting_colours"] = [
                {"homozygote_deletie": "rgb(128,0,0)"},
                {"heterozygote_deletie": "rgb(255,0,0)"},
                {"duplicatie": "rgb(0,0,255)"},
                {"triplicatie": "rgb(0,0,128)"},
            ]
        headers[col] = header
    return headers


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ExomeDepth CNVs",
            info="Summary tables of CNV calls detected by ExomeDepth, per SeqCap panel and per design.",
        )

        panels_data, panels_columns, panels_rows = self.parse_cnv_files(
            "cnvs_exomedepth/panels"
        )
        designs_data, designs_columns, designs_rows = self.parse_cnv_files(
            "cnvs_exomedepth/designs"
        )

        if not panels_data and not designs_data:
            log.debug("No ExomeDepth CNV summary files found")
            return

        patient_summary = build_patient_summary(panels_rows, designs_rows)
        patient_filter_html = build_patient_filter_html(patient_summary)

        # Panel CNVs table is added first, Design CNVs second
        if panels_data:
            self.write_data_file(panels_data, "cnvs_exomedepth_panels")
            self.add_section(
                name="Panel CNVs",
                anchor="cnvs_exomedepth_panels",
                description="CNVs called by ExomeDepth, summarised per panel.",
                content_before_plot=f"{patient_filter_html}<p>{READS_RATIO_LEGEND}</p>",
                plot=table.plot(
                    data=panels_data,
                    headers=build_headers(panels_data, panels_columns),
                    pconfig={
                        "id": "cnvs_exomedepth_panels_table",
                        "title": "ExomeDepth: Panel CNVs",
                        "col1_header": "Patient",
                        "sort_rows": True,
                        "no_violin": True,
                    },
                ),
            )

        if designs_data:
            self.write_data_file(designs_data, "cnvs_exomedepth_designs")
            # If panels_data wasn't present, place the patient filter at the top of Design CNVs
            design_before_plot = (
                f"{patient_filter_html}<p>{READS_RATIO_LEGEND}</p>"
                if not panels_data
                else f"<p>{READS_RATIO_LEGEND}</p>"
            )
            self.add_section(
                name="Design CNVs",
                anchor="cnvs_exomedepth_designs",
                description="CNVs called by ExomeDepth, summarised per design.",
                # Collapsed by default: wrap the plot in a native <details> disclosure
                content_before_plot=(
                    f'{design_before_plot}<details id="cnvs_exomedepth_designs_details">'
                    f'<summary id="cnvs_exomedepth_designs_summary">Click to show the Design CNVs table</summary>'
                ),
                plot=table.plot(
                    data=designs_data,
                    headers=build_headers(designs_data, designs_columns),
                    pconfig={
                        "id": "cnvs_exomedepth_designs_table",
                        "title": "ExomeDepth: Design CNVs",
                        "col1_header": "Patient",
                        "sort_rows": True,
                        "no_violin": True,
                    },
                ),
                content="</details>",
            )

    def parse_cnv_files(
        self, sp_key: str, key_col: str = "Patient"
    ) -> Tuple[Dict[str, Dict[str, Value]], List[str], List[Dict[str, str]]]:
        rows: List[Dict[str, str]] = []
        header_cols: List[str] = []
        for f in self.find_log_files(sp_key, filecontents=True, filehandles=False):
            self.add_data_source(f)
            file_header_cols, file_rows = parse_tsv(f["f"])
            if not file_header_cols:
                continue
            header_cols = file_header_cols
            rows.extend(file_rows)

        if not rows:
            return {}, [], []

        data = build_data(rows, key_col)
        columns = [c for c in header_cols if c != key_col]
        return data, columns, rows
