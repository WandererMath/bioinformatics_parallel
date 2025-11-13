from func import *
from mylib.utils import get_paths_ends_with_something

# sam2bedgraph("test/A1.txt", "test/cov", "test/norm")

merge_bedgraphs(get_paths_ends_with_something("test/norm", ".bedgraph"), "test/norm_merged/test.bedgrph")