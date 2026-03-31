aaindex Documentation
=====================

.. image:: ../../images/aaindex_logo.png
   :align: center
   :alt: aaindex logo
   :height: 400px

**aaindex** is a lightweight Python package for accessing the data in the
`AAindex <https://www.genome.jp/aaindex/>`_ databases, which represent the
physicochemical, biochemical, and structural properties of amino acids as
numerical indices. The AAindex database, maintained by the
`Kyoto Encyclopedia of Genes and Genomes (KEGG) <https://www.kegg.jp/>`_,
is a widely-used reference resource in bioinformatics and computational biology,
providing curated amino acid property data derived from published literature.

This package offers a simple, dependency-free interface for querying any record
across all three AAindex databases. Records can be retrieved by their accession
number or searched by keyword, and individual fields — including numerical values,
descriptions, references, and correlation coefficients — are all directly accessible.
It is well-suited for use in protein sequence analysis, feature engineering for
machine learning models, and any application that requires standardised amino acid
physicochemical data.

The package is used other key custom-built projects: `pySAR <https://github.com/amckenna41/pySAR/>`_ and `protpy <https://github.com/amckenna41/protpy/>`_.

The package supports all three AAindex databases:

- **AAindex1** — 566 amino acid indices (single numerical values per amino acid)
- **AAindex2** — 94 amino acid substitution matrices (20 x 20)
- **AAindex3** — 47 amino acid contact potential matrices (20 x 20)

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   userGuide

.. toctree::
   :maxdepth: 2
   :caption: Database Reference

   api/aaindex1
   api/aaindex2
   api/aaindex3
   api/aaindex_matrix

.. toctree::
   :maxdepth: 1
   :caption: Project

   changelog


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`