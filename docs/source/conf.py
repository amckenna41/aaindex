#!/usr/bin/env python3
# Configuration file for the Sphinx documentation builder.
#
# For a full list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys

# Add the project root so autodoc can import the aaindex package.
sys.path.insert(0, os.path.abspath("../.."))

import aaindex  # noqa: E402

# -- Project information -----------------------------------------------------

project = "aaindex"
copyright = "2021-2026, AJ McKenna"
author = "AJ McKenna"
version = aaindex.__version__
release = aaindex.__version__


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = []

# Napoleon settings — project uses Google-style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = False


# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
master_doc = "index"

# -- Intersphinx mapping -----------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}