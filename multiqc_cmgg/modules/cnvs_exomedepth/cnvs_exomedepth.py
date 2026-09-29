import logging
from collections import defaultdict
from typing import Dict, List, Tuple, Union

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

        panels_data, panels_columns = self.parse_cnv_files("cnvs_exomedepth/panels")
        designs_data, designs_columns = self.parse_cnv_files("cnvs_exomedepth/designs")

        if not panels_data and not designs_data:
            log.debug("No ExomeDepth CNV summary files found")
            return

        # Panel CNVs table is added first, Design CNVs second
        if panels_data:
            self.write_data_file(panels_data, "cnvs_exomedepth_panels")
            self.add_section(
                name="Panel CNVs",
                anchor="cnvs_exomedepth_panels",
                description="CNVs called by ExomeDepth, summarised per panel.",
                content_before_plot=f"<p>{READS_RATIO_LEGEND}</p>",
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
            self.add_section(
                name="Design CNVs",
                anchor="cnvs_exomedepth_designs",
                description="CNVs called by ExomeDepth, summarised per design.",
                # Collapsed by default: wrap the plot in a native <details> disclosure
                content_before_plot=f"<p>{READS_RATIO_LEGEND}</p><details><summary>Click to show the Design CNVs table</summary>",
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
    ) -> Tuple[Dict[str, Dict[str, Value]], List[str]]:
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
            return {}, []

        data = build_data(rows, key_col)
        columns = [c for c in header_cols if c != key_col]
        return data, columns
