import re
import os
from mylib.utils import get_paths_ends_with_something

BASE_DIR='Ribo-seq/intermediate'

if __name__ == '__main__':
    RAW_DIR=os.path.join(BASE_DIR, 'raw')
    raw_names= get_paths_ends_with_something(RAW_DIR, '.gz')
    raw_names=[os.path.basename(name) for name in raw_names]
    os.chdir(RAW_DIR)
    for filename in raw_names:
        # match = re.search(r'KFredrick016_([A-Z]R\d)_.*_(R\d)_', filename)
        # if match:
        #     sample= match.group(1)
        #     r= match.group(2)
        #     new_name = f"{sample}_{r}.fastq.gz"
        #     os.rename(filename, new_name)
        # else:
        #     raise Exception()
        new_name=filename[0]+filename[2:]
        os.rename(filename, new_name)
