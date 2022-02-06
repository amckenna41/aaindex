# aaindex
Python package for working with AAIndex database (https://www.genome.jp/aaindex/)

![alt text](https://raw.githubusercontent.com/amckenna41/aaindex/main/images/aaindex_logo.png)

# pySAR <a name="TOP"></a>
[![AAIndex](https://img.shields.io/pypi/v/aaindex)](https://pypi.org/project/aaindex/)
[![pytest](https://github.com/amckenna41/pySAR/workflows/Building%20and%20Testing%20%F0%9F%90%8D/badge.svg)](https://github.com/amckenna41/pySAR/actions?query=workflowBuilding%20and%20Testing%20%F0%9F%90%8D)
[![Platforms](https://img.shields.io/badge/platforms-linux%2C%20macOS%2C%20Windows-green)](https://pypi.org/project/aaindex/)
[![PythonV](https://img.shields.io/pypi/pyversions/pySAR?logo=2)](https://pypi.org/project/aaindex/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)
[![Build](https://img.shields.io/github/workflow/status/amckenna41/pySAR/Deploy%20to%20PyPI%20%F0%9F%93%A6)](https://github.com/amckenna41/pySAR/actions)
[![Build Status](https://travis-ci.com/amckenna41/aaindex.svg?branch=main)](https://travis-ci.com/amckenna41/aaindex)
[![Issues](https://img.shields.io/github/issues/amckenna41/aaindex)](https://github.com/amckenna41/aaindex/issues)
[![Size](https://img.shields.io/github/repo-size/amckenna41/aaindex)](https://github.com/amckenna41/aaindex)
[![Commits](https://img.shields.io/github/commit-activity/w/amckenna41/aaindex)](https://github.com/amckenna41/aaindex)

Table of Contents
-----------------

  * [Introduction](#introduction)
  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Tests](#tests)
  * [Contact](#contact)

Introduction
------------

Requirements
------------
* [Python][python] >= 3.6
* [numpy][numpy] >= 1.16.0
* [pandas][pandas] >= 1.1.0
* [sklearn][sklearn] >= 0.24
* [requests][requests] >= 2.24.0


Installation
-----------------
Install the latest version of aaindex using pip:

```bash
pip3 install aaindex
```

Install by cloning repository:
```bash
git clone https://github.com/amckenna41/aaindex.git
python3 setup.py install
```


Usage
-----

```python
import aaindex as aai

```

Directories
-----------
* `/pySAR/PyBioMed` - package partially forked from https://github.com/gadsbyfly/PyBioMed, used in
the calculation of the protein descriptors.
* `/Results` - stores the associated research paper created alongside the software. Also includes all results generated from the research as well as the article's supplementary materials.
* `/tests` - unit and integration tests for pySAR.
* `/data` - all required data and datasets are stored in this folder.


Tests
-----
To run all tests, from the main aaindex folder run:
```
python3 -m unittest discover
```

To run main test module, from the main aaindex folder run:
```
python -m unittest tests.test_aaindex -v
```

Contact
-------
If you have any questions or comments, please contact amckenna41@qub.ac.uk or raise an issue on the [Issues][Issues] tab.

[Back to top](#TOP)


[python]: https://www.python.org/downloads/release/python-360/
[numpy]: https://numpy.org/
[pandas]: https://pandas.pydata.org/
[sklearn]: https://scikit-learn.org/stable/
[requests]: https://docs.python-requests.org/en/latest/
[Issues]: https://github.com/amckenna41/pySAR/issues

