import os
import subprocess as sp

tests = [
    {'call': [
        'htseq-count',
        'example_data/bamfile_no_qualities.sam',
        'example_data/bamfile_no_qualities.gtf',
        ],
     'expected_fn': 'example_data/bamfile_no_qualities.tsv'},
    {'call': [
        'htseq-count',
        '-c', 'test_output.tsv',
        'example_data/bamfile_no_qualities.sam',
        'example_data/bamfile_no_qualities.gtf',
        ],
     'expected_fn': 'example_data/bamfile_no_qualities.tsv'},
    # Testing multiple cores on travis makes a mess
    #{'call': [
    #    'htseq-count',
    #    '-n', '2',
    #    'example_data/bamfile_no_qualities.sam',
    #    'example_data/bamfile_no_qualities.gtf',
    #    ],
    # 'expected_fn': 'example_data/bamfile_no_qualities.tsv'},
    {'call': [
        'htseq-count',
        'example_data/bamfile_no_qualities.bam',
        'example_data/bamfile_no_qualities.gtf',
        ],
     'expected_fn': 'example_data/bamfile_no_qualities.tsv'},
    {'call': [
        'htseq-count',
        '-m', 'intersection-nonempty',
        '--nonunique', 'none',
        '--secondary-alignments', 'score',
        '--supplementary-alignments', 'score',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ],
     'expected_fn': 'example_data/yeast_RNASeq_excerpt_withNH_counts.tsv'},
    {'call': [
        'htseq-count',
        '-m', 'intersection-nonempty',
        '--nonunique', 'none',
        '--secondary-alignments', 'score',
        '--supplementary-alignments', 'score',
        '--feature-query', 'gene_id == "YPR036W-A"',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ],
     'expected_fn': 'example_data/yeast_RNASeq_excerpt_withNH_counts_YPR036W-A.tsv'},
    {'call': [
        'htseq-count-barcodes',
        '-m', 'intersection-nonempty',
        '--nonunique', 'none',
        '--secondary-alignments', 'score',
        '--supplementary-alignments', 'score',
        'example_data/yeast_RNASeq_excerpt_withbarcodes.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ],
     'expected_fn': 'example_data/yeast_RNASeq_excerpt_withbarcodes.tsv'},
    {'call': [
        'htseq-count',
        '-m', 'intersection-nonempty',
        '--nonunique', 'none',
        '--secondary-alignments', 'score',
        '--supplementary-alignments', 'score',
        '--additional-attr', 'gene_name',
        '--additional-attr', 'exon_number',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ],
     'expected_fn': 'example_data/yeast_RNASeq_excerpt_withNH_counts_additional_attributes.tsv'},
    {'call': [
        'htseq-count',
        '-m', 'intersection-nonempty',
        '--nonunique', 'fraction',
        '--secondary-alignments', 'score',
        '--supplementary-alignments', 'score',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ],
     'expected_fn': 'example_data/yeast_RNASeq_excerpt_withNH_counts_nonunique_fraction.tsv'},
    {'call': [
        'htseq-count',
        '-m', 'intersection-nonempty',
        '--nonunique', 'all',
        '--secondary-alignments', 'score',
        '--supplementary-alignments', 'score',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ],
     'expected_fn': 'example_data/yeast_RNASeq_excerpt_withNH_counts_nonunique.tsv'},
    {'call': [
        'htseq-count',
        '-m', 'intersection-nonempty',
        '-i', 'gene_id',
        '--additional-attr', 'gene_name',
        '--nonunique', 'none',
        '--secondary-alignments', 'score',
        '--supplementary-alignments', 'score',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ],
     'expected_fn': 'example_data/yeast_RNASeq_excerpt_withNH_counts_twocolumns.tsv'},
    {'call': [
        'htseq-count',
        '-m', 'intersection-nonempty',
        '--nonunique', 'none',
        '--secondary-alignments', 'ignore',
        '--supplementary-alignments', 'score',
        'example_data/yeast_RNASeq_excerpt_withNH.sam',
        'example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
        ],
     'expected_fn': 'example_data/yeast_RNASeq_excerpt_withNH_counts_ignore_secondary.tsv'},
    ]


# Run the tests
for t in tests:
    expected_fn = t['expected_fn']
    call = t['call']
    # local testing
    #call = ['python', 'python3/HTSeq/scripts/count.py'] + call[1:]

    print(' '.join(call))
    output = sp.check_output(call).decode()

    if '-c' in call:
        output_fn = call[call.index('-c') + 1]
        with open(output_fn, 'r') as f:
            output = f.read()
    else:
        output_fn = None

    with open(expected_fn, 'r') as f:
        expected = f.read()

    try:
        assert output == expected
    except AssertionError:
        for out, exp in zip(output.split('\n'), expected.split('\n')):
            print(out, exp)
            if out != exp:
                break

        raise
    finally:
        if output_fn is not None:
            os.remove(output_fn)

