AAindex1
========

The AAindex1 section currently contains 566 amino acid indices representing the various
physicochemical, structural and biochemical properties of amino acids. Each entry consists
of an accession number, a short description on the index, the reference information, notes,
PMID (pubmed ID) and the numerical values for the property of 20 amino acids. In addition,
it contains neighbour information; namely, the cross-links to other entries with an absolute
value for the correlation coefficient of 0.8 or larger, allowing users to identify entries
describing similar properties.

.. rubric:: Record format

.. code-block:: text

   ************************************************************************
   *                                                                      *
   * H Accession number                                                   *
   * D Data description                                                   *
   * R Pub med article ID (PMID)                                          *
   * A Author(s)                                                          *
   * T Title of the article                                               *
   * J Journal reference                                                  *
   * * Comment or missing                                                 *
   * C Accession numbers of similar entries with the correlation          *
   *   coefficients of 0.8 (-0.8) or more (less).                         *
   *   Notice: The correlation coefficient is calculated with zeros       *
   *   filled for missing values.                                         *
   * I Amino acid index data in the following order                       *
   *   Ala    Arg    Asn    Asp    Cys    Gln    Glu    Gly    His    Ile *
   *   Leu    Lys    Met    Phe    Pro    Ser    Thr    Trp    Tyr    Val *
   * //                                                                   *
   ************************************************************************

.. automodule:: aaindex.aaindex1
   :no-members:

.. autoclass:: aaindex.aaindex1.AAIndex1
   :members: num_records, record_codes, record_names, values, search,
             amino_acids, get_record_by_category, get_all_categories,
             parse_categories, parse_aaindex
   :undoc-members:
   :show-inheritance:
   :special-members: __getitem__, __len__, __contains__, __iter__, __repr__, __sizeof__

.. rubric:: References

[1] Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database.
*Nucleic Acids Res.* 28, 374 (2000).