AAindex2
========

The AAindex2 section currently contains 94 amino acid mutation matrices: 47 symmetric
matrices and 19 non-symmetric matrices. The format of the entry is almost the same as
that of AAindex1 except that it contains 210 numerical values (20 diagonal and
20 × 19/2 off-diagonal elements) for a symmetric matrix and 400 or more numerical values
for a non-symmetric matrix (some matrices include a gap or distinguish two states of
cysteine).

.. rubric:: Record format

.. code-block:: text

   ************************************************************************
   *                                                                      *
   * Each entry has the following format.                                 *
   *                                                                      *
   * H Accession number                                                   *
   * D Data description                                                   *
   * R PMID                                                               *
   * A Author(s)                                                          *
   * T Title of the article                                               *
   * J Journal reference                                                  *
   * * Comment or missing                                                 *
   * M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV           *
   *   AA                                                                 *
   *   AR RR                                                              *
   *   AN RN NN                                                           *
   *   AD RD ND DD                                                        *
   *   AC RC NC DC CC                                                     *
   *   AQ RQ NQ DQ CQ QQ                                                  *
   *   AE RE NE DE CE QE EE                                               *
   *   AG RG NG DG CG QG EG GG                                            *
   *   AH RH NH DH CH QH EH GH HH                                         *
   *   AI RI NI DI CI QI EI GI HI II                                      *
   *   AL RL NL DL CL QL EL GL HL IL LL                                   *
   *   AK RK NK DK CK QK EK GK HK IK LK KK                                *
   *   AM RM NM DM CM QM EM GM HM IM LM KM MM                             *
   *   AF RF NF DF CF QF EF GF HF IF LF KF MF FF                          *
   *   AP RP NP DP CP QP EP GP HP IP LP KP MP FP PP                       *
   *   AS RS NS DS CS QS ES GS HS IS LS KS MS FS PS SS                    *
   *   AT RT NT DT CT QT ET GT HT IT LT KT MT FT PT ST TT                 *
   *   AW RW NW DW CW QW EW GW HW IW LW KW MW FW PW SW TW WW              *
   *   AY RY NY DY CY QY EY GY HY IY LY KY MY FY PY SY TY WY YY           *
   *   AV RV NV DV CV QV EV GV HV IV LV KV MV FV PV SV TV WV YV VV        *
   * //                                                                   *
   ************************************************************************

.. automodule:: aaindex.aaindex2
   :no-members:

.. autoclass:: aaindex.aaindex2.AAIndex2
   :members: num_records, record_codes, record_names, values, get, search,
             amino_acids, parse_aaindex
   :inherited-members:
   :undoc-members:
   :show-inheritance:
   :special-members: __getitem__, __len__, __contains__, __iter__, __repr__, __sizeof__

.. rubric:: References

[1] Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database.
*Nucleic Acids Res.* 28, 374 (2000).
