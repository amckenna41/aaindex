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

   # Encode a sequence as a numeric vector using one index
   aaindex1.encode("ACDEFGHIKLM", "KYTJ820101")
   # [1.8, 2.5, -3.5, -3.5, 2.8, ...]

   # Specify a fill value for non-canonical amino acids (default: None)
   aaindex1.encode("ACXZ", "KYTJ820101", unknown_value=0.0)
   # [1.8, 2.5, 0.0, 0.0]

   # Encode against multiple indices in one call (multi-index feature extraction)
   aaindex1.multi_encode("ACDE", ["KYTJ820101", "EISD840101"])
   # {'KYTJ820101': [1.8, 2.5, -3.5, -3.5], 'EISD840101': [...]}

   # Encode a batch of sequences
   aaindex1.batch_encode(["ACDE", "FGHI"], "KYTJ820101")
   # [[1.8, 2.5, -3.5, -3.5], [-3.2, -0.4, ...]]

   # Export as DataFrame (set include_gap=False to drop the '-' column)
   aaindex1.to_dataframe("KYTJ820101", include_gap=False)

   # Plot amino acid value distribution (requires matplotlib)
   ax = aaindex1.plot_distribution("KYTJ820101")
   import matplotlib.pyplot as plt; plt.show()


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

   # Search by keyword (results sorted alphabetically)
   results = aaindex2.search("mutation")

   # Plot the 20×20 matrix as a heatmap (requires matplotlib)
   ax = aaindex2.plot_heatmap("ALTS910101")
   import matplotlib.pyplot as plt; plt.show()


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

   # Plot the 20×20 contact potential matrix as a heatmap (requires matplotlib)
   ax = aaindex3.plot_heatmap("TANS760101")
   import matplotlib.pyplot as plt; plt.show()


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


