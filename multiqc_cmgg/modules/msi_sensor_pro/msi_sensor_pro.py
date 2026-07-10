import logging
import re
from collections import defaultdict
from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc import report
from multiqc.utils.util_functions import update_dict
from multiqc.plots import table, bargraph
from typing import Dict, Union, List, Optional
from collections import OrderedDict

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="msisensor-pro",
            info="This table show a summary of the output of the msisensor-pro pro command for MSI detection.",
        )

        # Load configurable thresholds (default values if not specified)
        self.coverage_threshold = getattr(
            config, "msi_sensor_pro_coverage_threshold", 30
        )
        self.min_sites_threshold = getattr(config, "msi_sensor_pro_min_sites", 30)
        self.msi_high_threshold = getattr(config, "msi_high_threshold", 30.0)
        self.low_coverage_sites_threshold = getattr(
            config, "msi_sensor_pro_low_coverage_sites_threshold", 10
        )

        # Parsing and loading data from msiSensorPro summary and all files
        data_dicts_summary = self.parse_summary()
        log.info(f"Summary samples: {list(data_dicts_summary.keys())}")
        data_dicts_all = self.parse_all()
        log.info(f"All-loci samples: {list(data_dicts_all.keys())}")
        self.annotate_summary_low_coverage(data_dicts_summary, data_dicts_all)
        msisensorpro_data, all_zero = self.prepare_msisensorpro_data(data_dicts_summary)

        # Table configuration
        config_table = {
            "id": "msi_summary",
            "title": "msi_summary",
        }
        headers = {
            "num_sites": {
                "title": "Number of sites",
            },
            "num_unstable_sites": {
                "title": "Number of unstable sites",
            },
            "low_coverage_sites": {
                "title": "Number of low-coverage sites",
            },
            "perc": {
                "title": "Percentage of unstable sites",
                "format": "{:.2f}",
                "suffix": "%",
            },
        }

        # headers for all table
        all_loci = set()
        for sample in data_dicts_all.values():
            all_loci.update(sample.keys())

        headers2 = {}
        for locus in sorted(all_loci):
            headers2[locus] = {"title": locus, "description": f"MSI status at {locus}"}

        # summary table - also convert to HTML to avoid violin warnings
        try:
            # Sort summary data: primary by perc (descending), secondary by low_coverage_sites (ascending)
            sorted_samples = sorted(
                data_dicts_summary.items(),
                key=lambda x: (-x[1].get('perc', 0), x[1].get('low_coverage_sites', 0))
            )
            
            # Build table with sortable columns
            table_id = "msi-summary-table"
            html_parts = [
                f'<style>',
                f'  .bg-light-success {{ background-color: #d4edda !important; }}',
                f'  .sortable {{ cursor: pointer; user-select: none; }}',
                f'  .sortable::after {{ content: " ⇅"; font-size: 0.8em; opacity: 0.5; }}',
                f'</style>',
                f'<table id="{table_id}" class="table table-striped table-hover">',
                f'<thead><tr>',
                f'  <th class="sortable" onclick="sortTable(\'{table_id}\', 0)">Sample</th>'
            ]
            
            # Header row with sortable columns
            for idx, (col_key, col_info) in enumerate(headers.items(), start=1):
                title = col_info.get('title', col_key)
                html_parts.append(f'  <th class="sortable" onclick="sortTable(\'{table_id}\', {idx})">{title}</th>')
            html_parts.append('</tr></thead><tbody>')
            
            # Data rows (already sorted)
            for sample_name, sample_data in sorted_samples:
                html_parts.append('<tr>')
                html_parts.append(f'<td>{sample_name}</td>')
                for col_key in headers.keys():
                    value = sample_data.get(col_key, '')
                    # Apply formatting and conditional colors
                    if col_key == 'perc':
                        try:
                            val = float(value)
                            css_class = ''
                            if val >= self.msi_high_threshold:
                                css_class = 'class="bg-danger"'  # Red for high percentage
                            else:
                                css_class = 'class="bg-light-success"'  # Light green for normal
                            html_parts.append(f'<td {css_class}>{value:.2f}%</td>')
                        except:
                            html_parts.append(f'<td>{value}</td>')
                    elif col_key == 'low_coverage_sites':
                        try:
                            val = int(value)
                            css_class = ''
                            if val > self.low_coverage_sites_threshold:
                                css_class = 'class="bg-warning"'  # Orange for high
                            html_parts.append(f'<td {css_class}>{value}</td>')
                        except:
                            html_parts.append(f'<td>{value}</td>')
                    elif col_key == 'num_sites':
                        try:
                            val = int(value)
                            css_class = ''
                            if val <= self.min_sites_threshold:
                                css_class = 'class="bg-warning"'  # Orange for low
                            html_parts.append(f'<td {css_class}>{value}</td>')
                        except:
                            html_parts.append(f'<td>{value}</td>')
                    else:
                        html_parts.append(f'<td>{value}</td>')
                html_parts.append('</tr>')
            html_parts.append('</tbody></table>')
            
            # Add JavaScript for sortable columns
            html_parts.append('''
            <script>
            function sortTable(tableId, columnIdx) {
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                
                const isNumeric = (str) => !isNaN(parseFloat(str)) && isFinite(str);
                
                rows.sort((a, b) => {
                    const aCell = a.cells[columnIdx].textContent.trim();
                    const bCell = b.cells[columnIdx].textContent.trim();
                    
                    if (isNumeric(aCell) && isNumeric(bCell)) {
                        return parseFloat(bCell) - parseFloat(aCell);  // Numeric descending
                    }
                    return aCell.localeCompare(bCell);  // String ascending
                });
                
                tbody.innerHTML = '';
                rows.forEach(row => tbody.appendChild(row));
            }
            </script>
            ''')
            
            summary_html = '\n'.join(html_parts)
            self.add_section(
                content=summary_html
            )
        except Exception as e:
            log.debug(f"Failed to create summary HTML: {e}", exc_info=True)
            # Fallback to original table.plot
            self.add_section(
                plot=table.plot(
                    data=data_dicts_summary, headers=headers, pconfig=config_table
                ),
            )

        # all table
        try:
            try:
                num_samples = len(data_dicts_all)
                num_loci = len(headers2)
                total_cells = sum(len(v) for v in data_dicts_all.values())
                log.info(f"Writing all-loci data: samples={num_samples}, loci={num_loci}, cells={total_cells}")
            except Exception:
                log.debug("Failed to compute all-loci stats", exc_info=True)

            # Write raw TSV using BaseModule method so MultiQC tracks it
            self.write_data_file(data_dicts_all, "msiSensorPro_all_table")

            # Also write JSON as backup
            try:
                report.write_data_file(data_dicts_all, "msiSensorPro_all_table_json", data_format="json")
            except Exception:
                log.debug("Failed to write all-loci JSON data file", exc_info=True)
        except Exception:
            log.debug("Failed to write explicit all-loci data file", exc_info=True)

        # Create HTML table for all loci
        try:
            table_id = "msi-all-loci-table"
            html_parts = [
                '<style>',
                '  .bg-light-success { background-color: #d4edda !important; }',
                '</style>',
                f'<table id="{table_id}" class="table table-striped table-hover">',
                '<thead><tr><th>Sample</th>'
            ]
            
            for locus in sorted(headers2.keys()):
                html_parts.append(f'<th>{locus}</th>')
            html_parts.append('</tr></thead>')
            
            # Data rows
            html_parts.append('<tbody>')
            for sample_name in sorted(data_dicts_all.keys()):
                html_parts.append('<tr>')
                html_parts.append(f'<td>{sample_name}</td>')
                for locus in sorted(headers2.keys()):
                    status = data_dicts_all[sample_name].get(locus, 'N/A')
                    css_class = ''
                    if isinstance(status, str):
                        if 'Stable' in status:
                            css_class = 'class="bg-light-success"'  # Light green
                        elif 'Unstable' in status:
                            css_class = 'class="bg-danger"'   # Red
                        elif 'Low-coverage' in status:
                            css_class = 'class="bg-warning"'  # Orange
                    html_parts.append(f'<td {css_class}>{status}</td>')
                html_parts.append('</tr>')
            html_parts.append('</tbody>')
            html_parts.append('</table>')
            
            html_str = '\n'.join(html_parts)
            
            self.add_section(
                name="msisensor-pro - All Loci",
                anchor="msisensorpro_all_loci",
                description="Detailed MSI status for every loci with the coverage in brackets.",
                content=html_str,
            )
        except Exception as e:
            log.debug(f"Failed to create custom HTML table: {e}", exc_info=True)
            self.add_section(
                name="msisensor-pro - All Loci",
                anchor="msisensorpro_all_loci",
                description="Detailed MSI status for every loci with the coverage in brackets.",
                plot=table.plot(
                    data=data_dicts_all,
                    headers=headers2,
                    pconfig={
                        "id": "msiSensorPro_all_table",
                        "title": "msiSensorPro - All Site Metrics",
                        "no_violin": True,
                        "parse_numeric": False,
                        "save_file": True,
                        "save_data_file": True,
                    },
                ),
            )
        if not all_zero:
            categories = {
                "MSS": {"name": "MSS", "color": "#2ecc71"},
                "MSI-high": {"name": "MSI-high", "color": "#e74c3c"},
                "Low-coverage": {"name": "Low-coverage", "color": "#f39c12"},
            }

            # Bargraph configuration
            self.add_section(
                name="msisensor-pro - Bargraph",
                anchor="msisensorpro_bargraph",
                description="This graph visualizes the MSI status per sample in a bargraph.",
                plot=bargraph.plot(
                    data=msisensorpro_data,
                    cats=categories,
                    pconfig={
                        "id": "msiSensorPro_bargraph_v2",
                        "title": "MSI Sensor Pro Summary",
                        "ylab": "Percentage of unstable sites",
                        "ymin": 0,
                        "ymax": 100,
                        "cpswitch_counts_label": "Number of Sites",
                    },
                ),
            )
        else:
            log.info("Skipping bargraph: All samples have 0% unstable sites")

    def prepare_msisensorpro_data(self, data_summary):
        """
        Transform summary data to msisensorpro score with MSI classification
        """
        msisensorpro_data = {}
        all_zero = True
        min_bar = 0.9  # Minimum bar height for visibility in graph

        for sample_name, sample_data in data_summary.items():
            msisensorpro_score = sample_data["perc"]
            low_cov_sites = sample_data.get("low_coverage_sites", 0)

            if low_cov_sites >= self.low_coverage_sites_threshold or sample_data["num_sites"] <= self.min_sites_threshold:
                msi_status = "Low-coverage"
            elif msisensorpro_score >= self.msi_high_threshold:
                msi_status = "MSI-high"
            else:
                msi_status = "MSS"

            display_score = msisensorpro_score  
            if display_score == 0.0:
                display_score = min_bar

            sample_entry = {
                msi_status: display_score,
            }
            
            sample_label = f"{sample_name} ({msi_status})"
            msisensorpro_data[sample_label] = sample_entry
            if any(value != 0.0 for value in sample_entry.values()):
                all_zero = False

        return msisensorpro_data, all_zero

    def normalize_sample_name(self, s_name: str) -> str:
        """
        Normalize sample names for both summary and all files by removing file-specific suffixes.
        """
        normalized = s_name.replace(".txt", "")
        normalized = re.sub(r"_(summary|all)_msi$", "", normalized)
        return normalized

    def parse_summary(self):
        """
        Parse the msiSensorPro summary file.
        """
        data_summary: Dict[str, Dict[str, float]] = {}
        for f in self.find_log_files(
            "msi_sensor_pro/summary", filecontents=True, filehandles=False
        ):
            raw_name = self.clean_s_name(f["fn"], f)
            s_name = self.normalize_sample_name(raw_name)
            
            log.debug(f"parse_summary: raw_fn='{f['fn']}', raw_name='{raw_name}', s_name='{s_name}'")
            lines = f["f"].splitlines()
            header = lines[0]
            for line in lines:
                if line != header:
                    num_sites, num_unstable_sites, perc = line.split("\t")
                    data_summary[s_name] = {
                        "num_sites": int(num_sites),
                        "num_unstable_sites": int(num_unstable_sites),
                        "perc": float(perc),
                    }
            log.info(data_summary)
        return data_summary

    def annotate_summary_low_coverage(
        self, data_summary: Dict[str, Dict[str, Union[int, float]]], data_all: Dict[str, Dict]
    ) -> None:
        for sample_name, summary in data_summary.items():
            low_coverage_count = sum(
                1
                for status in data_all.get(sample_name, {}).values()
                if isinstance(status, str) and status.startswith("Low-coverage")
            )
            summary["low_coverage_sites"] = low_coverage_count

    def parse_all(self) -> Dict[str, Dict]:
        """
        Parse the msiSensorPro all file into a loci-centric structure.
        """
        sample_data = defaultdict(dict)
        all_loci = set()

        for f in self.find_log_files(
            "msi_sensor_pro/all", filecontents=True, filehandles=False
        ):
            raw_name = self.clean_s_name(f["fn"], f)
            s_name = self.normalize_sample_name(raw_name)
            
            log.debug(f"parse_all: raw_fn='{f['fn']}', raw_name='{raw_name}', s_name='{s_name}'")
            lines = f["f"].splitlines()
            
            sample_data.setdefault(s_name, {})
            
            for line in lines[1:]:
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    log.warning(
                        f"Skipping line in {s_name} due to insufficient columns"
                    )
                    continue

                chrom, pos = parts[0], parts[1]
                locus_id = f"{chrom}:{pos}"
                all_loci.add(locus_id)
                pro_p = float(parts[6])
                coverage = int(parts[8])
                threshold = float(parts[9])

                if pro_p > threshold:
                    status = f"Unstable ({coverage})"
                elif coverage < self.coverage_threshold:
                    status = f"Low-coverage ({coverage})"
                else:
                    status = f"Stable ({coverage})"

                sample_data[s_name][locus_id] = status

        for sample in sample_data.values():
            missing = all_loci - set(sample.keys())
            for locus in missing:
                sample[locus] = "Low-coverage (0)"

        return sample_data