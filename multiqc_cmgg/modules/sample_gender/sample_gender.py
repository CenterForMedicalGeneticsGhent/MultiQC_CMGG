import logging
from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table
from typing import Dict, Union, Any
import math

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="Sex prediction",
            anchor="sample_gender",
            info="This table show the results of the sex prediction with the expected sex and some metrics.",
        )

        # Find and load any input files for this module
        self.samplegender_data: Dict[str, Dict[str, Union[float, str]]] = dict()

        # Define search keys and their extensions
        method_dict = {
            "sample_gender/xy": "_xy",
            "sample_gender/hetx": "_hetx",
            "sample_gender/sry": "_sry",
        }

        for method, extension in method_dict.items():
            for f in self.find_log_files(method):
                self.add_data_source(f)
                # Clean sample name
                s_name = f["s_name"]
                if s_name.endswith(extension):
                    s_name = s_name[: -len(extension)]

                parsed = parse_file(f["f"])

                if not parsed:
                    continue

                if s_name not in self.samplegender_data:
                    self.samplegender_data[s_name] = parsed
                else:
                    self.samplegender_data[s_name].update(parsed)

        # Filter to strip out ignored sample names
        self.samplegender_data = self.ignore_samples(self.samplegender_data)

        if len(self.samplegender_data) == 0:
            raise UserWarning

        log.info(f"Found {len(self.samplegender_data)} SampleGender reports")

        # Calculate certainty and predicted sex
        self._calculate_gender_consensus()

        # Write parsed report data to a file
        self.write_data_file(self.samplegender_data, "multiqc_ngsbits_samplegender")

        # Table configuration
        self._create_table()

    def _calculate_gender_consensus(self):
        list_gender_methods = ["gender_xy", "gender_hetx", "gender_sry"]

        for sample, data in self.samplegender_data.items():
            count_M = 0
            count_F = 0

            for method in list_gender_methods:
                gender = data.get(method)
                if gender == "M":
                    count_M += 1
                elif gender == "F":
                    count_F += 1

            total = 3.0  # number of methods
            percentage = 0.0
            calc_sex = "Unknown"

            if count_M >= count_F:
                percentage = count_M / total
                if count_M > count_F:
                    calc_sex = "M"
            elif count_F >= count_M:
                percentage = count_F / total
                if count_F > count_M:
                    calc_sex = "F"
            self.samplegender_data[sample]["certainty"] = percentage
            self.samplegender_data[sample]["calc_gender"] = calc_sex

    def _create_table(self):
        config_table = {
            "id": "sex_prediction",
            "title": "Sex prediction",
            "namespace": "SampleGender",
        }

        headers = {
            "certainty": {
                "title": "Certainty",
                "description": "Certainty of sex match",
                "format": "{:.0%}",
                "scale": False,
                "max": 1,
                "cond_formatting_rules": {
                    "pass": [{"eq": 1}],
                    "warn": [{"lt": 1}],
                    "fail": [{"lt": 0.4}],
                },
                "cond_formatting_colours": [
                    {"pass": "#5cb85c"},
                    {"warn": "#f0ad4e"},
                    {"fail": "#d9534f"},
                ],
            },
            "calc_gender": {
                "title": "Calculated Sex",
                "description": "Consensus predicted sex",
                "scale": False,
            },
            "gender_xy": {
                "title": "Sex (XY)",
                "description": "Predicted gender based on chromosome read ratios",
                "scale": False,
            },
            "gender_sry": {
                "title": "Sex (SRY)",
                "description": "Predicted gender based on SRY gene coverage",
                "scale": False,
            },
            "gender_hetx": {
                "title": "Sex (HETX)",
                "description": "Predicted gender based on heterozygous variants on X chromosome",
                "scale": False,
                "cond_formatting_rules": {
                    "unknown": [{"s_eq": "unknown (too few SNPs)"}]
                },
                "cond_formatting_colours": [{"unknown": "#808080"}],
            },
            "ratio_chry_chrx": {
                "title": "ChrY/ChrX Ratio",
                "description": "Ratio of reads mapped to ChrY vs ChrX",
                "min": 0,
                "format": "{:.4f}",
                "scale": "Purples",
            },
            "coverage_sry": {
                "title": "Coverage SRY",
                "description": "Coverage of SRY in chrY",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.2f}",
            },
            "het_fraction": {
                "title": "Fraction HETX",
                "description": "Fraction of heterozygous SNPs in chrX",
                "min": 0,
                "scale": "Reds",
                "format": "{:,.4f}",
            },
        }

        self.add_section(
            plot=table.plot(
                data=self.samplegender_data, headers=headers, pconfig=config_table
            ),
        )


def parse_file(f: str) -> Dict[str, Union[float, str]]:
    """
    Parses a single samplegender TSV file content and returns a dictionary
    with the relevant data from columns 2-6.
    """
    parsed_data: Dict[str, Union[float, str]] = {}
    lines = f.splitlines()

    if len(lines) < 2:
        return parsed_data

    headers = lines[0].strip().split("\t")[1:6]
    values = lines[1].strip().split("\t")[1:6]

    # Map full names to abbreviations
    if values[0].lower() == "male":
        values[0] = "M"
    elif values[0].lower() == "female":
        values[0] = "F"

    # Map original header names to internal keys
    param_map = {
        "reads_chry": "gender_xy",
        "het_fraction": "gender_hetx",
        "coverage_sry": "gender_sry",
    }

    # Update headers based on map
    # Note: The original code logic was a bit fragile: "if param in headers: headers[0]=test"
    # It assumes gender_xy/hetx/sry is always at index 0?
    # Let's preserve original intent but make it safer if possible.
    # The original code iterated through paramdict and if any key was in headers, it replaced headers[0] with the value.
    # This implies the file type is identified by the presence of a specific column.

    # Reconstructing logic:
    for col_name, new_key in param_map.items():
        if col_name in headers:
            # The first column is the gender result for that method
            headers[0] = new_key

    for key, value in zip(headers, values):
        try:
            val_float = float(value)
            if math.isnan(val_float):
                parsed_data[key] = "N/A"
            else:
                parsed_data[key] = val_float
        except ValueError:
            parsed_data[key] = value

    return parsed_data
