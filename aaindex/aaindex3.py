################################################################################
################                    AAindex3                   #################
################################################################################

#importing required modules and dependencies
from typing import List
from ._aaindex_matrix import _AAIndexMatrix

__all__: List[str] = ['AAIndex3', 'aaindex3']

class AAIndex3(_AAIndexMatrix):
    """Python parser for AAindex3: Amino Acid Contact Potential Matrix Database.

    Inherits all parsing, search, and lookup functionality from _AAIndexMatrix.
    Stores the 47 known 20x20 contact potential matrices from AAindex3
    (http://www.genome.jp/aaindex/).

    References:
        [1]: Kawashima, S. and Kanehisa, M.; AAindex: amino acid index database.
             Nucleic Acids Res. 28, 374 (2000).
    """

    def __init__(self) -> None:
        super().__init__("aaindex3")


#create instance of AAIndex3 class
aaindex3 = AAIndex3()