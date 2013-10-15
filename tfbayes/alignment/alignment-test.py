
from tfbayes.alignment import *
from tfbayes.interface import *
from tfbayes.phylotree import *

tree      = pt_root_t  ("alignment-test.nh")
alignment = alignment_t("alignment-test.fa", tree)

print str(alignment["mm9"])

for column in alignment:
    print column

