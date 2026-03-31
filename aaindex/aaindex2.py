################################################################################
################                    AAindex2                   #################
################################################################################

#importing required modules and dependencies
from typing import List
from ._aaindex_matrix import _AAIndexMatrix

__all__: List[str] = ['AAIndex2', 'aaindex2']

class AAIndex2(_AAIndexMatrix):
    """Python parser for AAindex2: Amino Acid Substitution Matrix Database.

    Inherits all parsing, search, and lookup functionality from _AAIndexMatrix.
    Stores the 94 known 20x20 substitution matrices from AAindex2
    (http://www.genome.jp/aaindex/).

    References:
        [1]: Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database.
             Nucleic Acids Res. 28, 374 (2000).
    """

    def __init__(self) -> None:
        super().__init__("aaindex2")


#create instance of AAIndex2 class
aaindex2 = AAIndex2()

