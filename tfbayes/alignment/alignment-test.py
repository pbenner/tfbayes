
# test native alignment class and fasta parser
# ------------------------------------------------------------------------------

from tfbayes.alignment import *
from tfbayes.interface import *
from tfbayes.phylotree import *

tree      = pt_root_t  ("alignment-test.nh")
alignment = alignment_t("alignment-test.fa", tree)

print alignment["mm9"]

for column in alignment:
    print column

# test converter from biopython
# ------------------------------------------------------------------------------

from tfbayes.alignment import *
from tfbayes.interface import *
from tfbayes.phylotree import *
from tfbayes           import alignio

tree           = pt_root_t  ("alignment-test.nh")
alignment_list = list(alignio.parse("alignment-test.maf", "maf"))
alignment      = alignment_list[0]
alignment      = alignment_t(alignment, tree)
