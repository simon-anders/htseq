import sys
import subprocess as sp

call = [
        'htseq-count',
        '-m', 'intersection-nonempty',
        '--nonunique', 'none',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ]
print(' '.join(call))
output = sp.check_output(call)

with open('example_data/yeast_RNASeq_excerpt_withNH_counts.tsv', 'r') as f:
    expected = f.read()

assert output == expected
