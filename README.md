<p float="left">
<img src="docs/images/CMGG_logo.png" width="250" title="CMMG">
<img src="docs/images/MultiQC_logo.png" width="250" title="MultiQC">
</p>

**MultiQC_CMGG is a plugin for MultiQC containing customized modules and templates**

For more information about MultiQC, see [http://multiqc.info](http://multiqc.info)

## Description

The MultiQC_CMGG plugin add a custom template, used for internal QC reports and module customizations that can't be merged in core MultiQC.

## Installation

This plugin can be installed using the following methods

- using `pip`:

```bash
pip install --upgrade --force-reinstall git+https://github.com/CenterForMedicalGeneticsGhent/MultiQC_CMGG.git
```

- using `setuptools`:

```bash
git clone https://github.com/CenterForMedicalGeneticsGhent/MultiQC_CMGG
cd MultiQC_CMGG
python setup.py install
```

## Usage

### Modules

---

Name: Sampletracking
Description: >
Parse Picard CrosscheckFingerprints output and format sensible tables and plots

---

### Templates

---

Name: cmgg
Description: >
CMGG specfic template with custom logo's and affiliate links. To enable this template, add the `-t/--template cmgg` option to the command line

---
