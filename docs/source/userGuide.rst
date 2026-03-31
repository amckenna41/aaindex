User Guide
==========

Installation
------------

Install the latest release from PyPI:

.. code-block:: bash

   pip install aaindex --upgrade

Or install from source:

.. code-block:: bash

   git clone https://github.com/amckenna41/aaindex.git
   cd aaindex
   pip install .


Quick Start
-----------

Import any of the three pre-built singleton instances:

.. code-block:: python

   from aaindex import aaindex1, aaindex2, aaindex3


AAindex1 — Amino Acid Indices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

AAindex1 contains 566 numerical indices representing various physicochemical and
biochemical properties of amino acids.

.. code-block:: python

   from aaindex import aaindex1

   # Number of records in the database
   aaindex1.num_records()   # 566

   # Access a record by its accession number
   record = aaindex1["ANDN920101"]
   record.description       # "alpha-CH chemical shifts (Andersen et al., 1992)"
   record.values            # {'A': 4.35, 'R': 4.38, ...}
   record.references
   record.pmid
   record.category
   record.correlation_coefficients

   # Get only the amino acid values for a record
   aaindex1.values("ANDN920101")  # {'A': 4.35, 'R': 4.38, ...}

   # Search records by keyword
   results = aaindex1.search("hydrophobicity")
   len(results)  # number of matching records

   # List all accession numbers and descriptions
   aaindex1.record_codes()
   aaindex1.record_names()

   # List valid amino acid codes
   aaindex1.amino_acids()  # ['-', 'A', 'C', 'D', ...]

   # Filter records by category
   aaindex1.get_record_by_category("sec_struct")

   # Visualise a record as a bar chart (requires matplotlib)
   aaindex1.plot("ANDN920101")


AAindex2 — Substitution Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

AAindex2 contains 94 amino acid substitution matrices, each providing a 20 x 20
matrix of pairwise scores.

.. code-block:: python

   from aaindex import aaindex2

   # Number of records
   aaindex2.num_records()   # 94

   # Get the full 20x20 matrix for a record
   matrix = aaindex2.values("ALTS910101")  # {'A': {'A': 0.0, 'R': ...}, ...}

   # Get a single pairwise score
   aaindex2.get("ALTS910101", "A", "R")  # -2.0

   # Access full record metadata
   record = aaindex2["ALTS910101"]
   record.description
   record.matrix

   # Search by keyword
   results = aaindex2.search("mutation")

   # Visualise as a heatmap (requires matplotlib)
   aaindex2.plot("ALTS910101")


AAindex3 — Contact Potential Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

AAindex3 contains 47 amino acid pairwise contact potential matrices, with the same
interface as AAindex2.

.. code-block:: python

   from aaindex import aaindex3

   # Number of records
   aaindex3.num_records()   # 47

   # Get a single pairwise score
   aaindex3.get("TANS760101", "A", "R")

   # Get the full 20x20 matrix
   aaindex3.values("TANS760101")

   # Visualise as a heatmap (requires matplotlib)
   aaindex3.plot("TANS760101")


Common Operations
-----------------

All three classes support Python container protocols:

.. code-block:: python

   # Length
   len(aaindex1)            # 566

   # Membership testing
   "ANDN920101" in aaindex1  # True

   # Iteration over record codes
   for code in aaindex1:
       print(code)

Records returned by ``__getitem__`` are :class:`~aaindex._aaindex_matrix.Map`
objects, which support both dict-style and attribute-style access:

.. code-block:: python

   record = aaindex1["ANDN920101"]
   record["description"]    # dict-style
   record.description       # attribute-style (identical result)


Visualisation
-------------

Each database provides a ``plot()`` method. AAindex1 produces a bar chart of the
20 amino acid values; AAindex2 and AAindex3 produce heatmaps of the 20 x 20
matrix. Requires ``matplotlib``:

.. code-block:: bash

   pip install matplotlib

.. code-block:: python

   aaindex1.plot("ANDN920101")
   aaindex2.plot("ALTS910101")
   aaindex3.plot("TANS760101")