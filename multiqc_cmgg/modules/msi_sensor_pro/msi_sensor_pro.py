import logging
import re
from collections import defaultdict
from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc import report
from multiqc.plots import table, bargraph
from typing import Dict, Union

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
        log.info("Parsed %d summary samples", len(data_dicts_summary))
        data_dicts_all = self.parse_all()
        log.info("Parsed %d all-loci samples", len(data_dicts_all))
        self.annotate_summary_low_coverage(data_dicts_summary, data_dicts_all)
        msisensorpro_data, all_zero = self.prepare_msisensorpro_data(data_dicts_summary)

        sorted_summary = dict(
            sorted(
                data_dicts_summary.items(),
                key=lambda item: (
                    -item[1].get("perc", 0),
                    item[1].get("low_coverage_sites", 0),
                ),
            )
        )

        headers = {
            "num_sites": {
                "title": "Number of sites",
                "cond_formatting_rules": {
                    "low": [{"n_le": self.min_sites_threshold}]
                },
                "cond_formatting_colours": [{"low": "#f39c12"}],
            },
            "num_unstable_sites": {
                "title": "Number of unstable sites",
            },
            "low_coverage_sites": {
                "title": "Number of low-coverage sites",
                "cond_formatting_rules": {
                    "high": [{"n_gt": self.low_coverage_sites_threshold}]
                },
                "cond_formatting_colours": [{"high": "#f39c12"}],
            },
            "perc": {
                "title": "Percentage of unstable sites",
                "format": "{:.2f}",
                "suffix": "%",
                "cond_formatting_rules": {
                    "high": [{"n_ge": self.msi_high_threshold}],
                    "normal": [{"n_lt": self.msi_high_threshold}],
                },
                "cond_formatting_colours": [
                    {"high": "#e74c3c"},
                    {"normal": "#d4edda"},
                ],
            },
        }

        # headers for all table
        all_loci = set()
        for sample in data_dicts_all.values():
            all_loci.update(sample.keys())

        headers2 = {}
        for locus in sorted(all_loci):
            headers2[locus] = {"title": locus, "description": f"MSI status at {locus}"}

        self.add_section(
            plot=table.plot(
                data=sorted_summary,
                headers=headers,
                pconfig={
                    "id": "msi_summary",
                    "title": "msi_summary",
                    "no_violin": True,
                },
            ),
        )

        # all table
        self.log_all_loci_stats(data_dicts_all, headers2)
        self.write_all_loci_data_files(data_dicts_all)

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

    def log_all_loci_stats(self, data_dicts_all, headers2):
        num_samples = len(data_dicts_all)
        num_loci = len(headers2)
        total_cells = sum(len(v) for v in data_dicts_all.values())
        log.info(
            f"Writing all-loci data: samples={num_samples}, loci={num_loci}, cells={total_cells}"
        )

    def write_all_loci_data_files(self, data_dicts_all):
        self.write_data_file(data_dicts_all, "msiSensorPro_all_table")
        report.write_data_file(
            data_dicts_all,
            "msiSensorPro_all_table_json",
            data_format="json",
        )

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

                if coverage < self.coverage_threshold:
                    status = f"Low-coverage ({coverage})"
                elif pro_p > threshold:
                    status = f"Unstable ({coverage})"
                else:
                    status = f"Stable ({coverage})"

                sample_data[s_name][locus_id] = status

        for sample in sample_data.values():
            missing = all_loci - set(sample.keys())
            for locus in missing:
                sample[locus] = "Low-coverage (0)"

        return sample_data