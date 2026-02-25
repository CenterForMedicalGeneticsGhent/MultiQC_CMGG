import logging
from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table
from typing import Dict, Union

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="targeted_MSH2",
            info="Module to report MSH2 hotspot variants NM_000251.3:c.942+3A>T,c.942+2T>A,c.942+2T>C,c.942+2T>G",
        )

        # Get config with defaults
        msh2_config = getattr(config, "targeted_MSH2_config", {"sanger_threshold": 28})
        log.info(f"threshold Sanger has been set to {msh2_config.get('sanger_threshold', 28)}")

        # Find and load any input files for this module
        MSH2_varcount_data: Dict[str, Dict[str, Union[float, str]]] = dict()

        log.info("Searching for targeted_MSH2 files...")
        for f in self.find_log_files("targeted_MSH2"):
            log.info(f"Found file: {f['fn']} with s_name: {f['s_name']}")
            self.add_data_source(f)
            s_name = f["s_name"]
            parsed = parse_file(f["f"], msh2_config)
            if s_name not in MSH2_varcount_data:
                MSH2_varcount_data[s_name] = parsed

            if s_name in MSH2_varcount_data:
                MSH2_varcount_data[s_name].update(parsed)
                # Filter to strip out ignored sample names
                MSH2_varcount_data = self.ignore_samples(MSH2_varcount_data)

        # Debug for amount of reports found
        n_reports_found = len(MSH2_varcount_data)
        if n_reports_found > 0:
            log.debug(f"Found {len(MSH2_varcount_data)} targeted_MSH2 reports")

        if n_reports_found == 0:
            log.debug("No targeted_MSH2 reports found")

        # Write parsed report data to a file
        self.write_data_file(MSH2_varcount_data, "multiqc_targeted_MSH2")
        self.add_software_version(None)

        # Add targeted_MSH2 Table
        config_table = {
            "id": "targeted_MSH2",
            "title": "Targeted: MSH2",
            "sort_rows": True,
            "no_violin": True,
        }

        headers = {
            "MSH2_c.942+3_wt": {
                "title": "WT readcount",
                "description": "wild type readcount",
                "scale": "PuBu",
            },
            "MSH2_c.942+3A>T": {
                "title": "c.942+3A>T frequency (readcount)",
                "description": "frequency c.942+3A>T (readcount)",
                "cond_formatting_rules": {"sanger": [{"s_contains": " "}]},
                "cond_formatting_colours": [{"sanger": "#EE4B2B"}],
                "scale": False,
            },
            "MSH2_c.942+2T>A": {
                "title": "c.942+2T>A frequency (readcount)",
                "description": "frequency c.942+2T>A (readcount)",
                "cond_formatting_rules": {"sanger": [{"s_contains": " "}]},
                "cond_formatting_colours": [{"sanger": "#EE4B2B"}],
                "scale": False,
            },
            "MSH2_c.942+2T>C": {
                "title": "c.942+2T>C frequency (readcount)",
                "description": "frequency c.942+2T>C (readcount)",
                "cond_formatting_rules": {"sanger": [{"s_contains": " "}]},
                "cond_formatting_colours": [{"sanger": "#EE4B2B"}],
                "scale": False,
            },
            "MSH2_c.942+2T>G": {
                "title": "c.942+2T>G frequency (readcount)",
                "description": "frequency c.942+2T>G (readcount)",
                "cond_formatting_rules": {"sanger": [{"s_contains": " "}]},
                "cond_formatting_colours": [{"sanger": "#EE4B2B"}],
                "scale": False,
            },
        }
        self.add_section(
            plot=table.plot(
                data=MSH2_varcount_data, headers=headers, pconfig=config_table
            ),
        )


def parse_file(f: str, config: Dict[str, int]) -> Dict[str, Union[float, str]]:
    """
    Parses a single samplegender TSV file content and returns a dictionary
    with the relevant data from columns 2-6.
    """
    parsed_data: Dict[str, Union[float, str]] = {}
    lines = f.splitlines()

    if len(lines) < 3:
        # Not enough data, return an empty dictionary
        log.warning(
            "Not enough lines in the file to parse MSH2 hotspot variant counts."
        )
        return parsed_data
    headers = lines[2].strip().split(" ")[1:6]
    values = lines[3].strip().split(" ")[1:6]

    for key, value in zip(headers, values):
        parsed_data[key] = value

    # Calculating frequency of variants and determining need for sanger sequencing:
    for variant, counts in parsed_data.items():
        if variant != "MSH2_c.942+3_wt":
            freq = round(
                (int(counts))
                / (int(parsed_data["MSH2_c.942+3_wt"]) + int(counts))
                * 100,
                2,
            )

            if freq >= float(config["sanger_threshold"]):
                parsed_data[variant] = f"{freq}% ({counts})"
            else:
                parsed_data[variant] = f"{freq}% ({counts})"
        else:
            parsed_data[variant] = int(parsed_data[variant])

    return parsed_data
