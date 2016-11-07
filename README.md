qfilt
=====

This simple program is meant to filter sequencing data,
optionally removing or splitting reads with poor quality scores
and to optionally _only_ retain fragments from reads that are tagged with a given 5' sequence.

BUILD & INSTALL
---------------

The build process requires [CMake](http://www.cmake.org/). To build, type:

    cmake [-DINSTALL_PREFIX=/install/path (default=/usr/local)] .
    make install

USAGE
-----

    qfilt [-h] [-v] [-o OUTPUT] [-q QSCORE] [-l LENGTH] [-m MODE] [-T PREFIX] [-t MISMATCH] (-F FASTQ | -Q FASTA QUAL)

To try it using the example data provided:

    qfilt -F data/test.fna data/test.qual -q 15 -l 30 -T ATATCGCGAGGA

### OUTPUT: ####

#### stderr: ####

	run settings:
		input fasta:         data/test.fna
		input qual:          data/test.qual
		min q-score:         15
		min fragment length: 30
		run mode:            0 (truncate/don't tolerate homopolymers/don't tolerate ambigs)
		5' tag:              ATATCGCGAGGA
		max tag mismatches:  0

	run summary:
		total bases       :  29570
		original reads    :  305
		q10               :  0.995502
		q20               :  0.860298
		q30               :  0.728745
		mean q-score      :  33.2273
		contributing reads:  5
		retained fragments:  5

	original read length distribution:
		mean:                96.9508
		median:              77
		variance             3743.03
		standard deviation:  61.1803
		min:                 49
		2.5%:                54
		97.5%:               332
		max:                 497

	retained fragment length distribution:
		mean:                41
		median:              37
		variance             72.5
		standard deviation:  8.51469
		min:                 33
		2.5%:                33
		97.5%:               54
		max:                 54

#### stdout: ####

    >GM98SRO01B77KU rank=0000671 x=796.0 y=1996.0 length=58
    CCACGCGTATCGATGTCGACTTTTTTTTCTTTTCTTACATAGTAG
    >GM98SRO01BA3RP rank=0000953 x=419.5 y=1603.5 length=87
    CTGATGCTGCACCAACTGTACTCCCTCGCGATA
    >GM98SRO01E1BKW rank=0001233 x=1948.0 y=846.0 length=66
    TACAGTTGGTGCAGCATCAGAAAAGTACGACATCGATACGCGTGGTCCTCGCGA
    >GM98SRO01DVVNY rank=0001304 x=1476.0 y=636.5 length=84
    ACGGCTGATGCTGCACCAACTGTACTCCCTCGCGATA
    >GM98SRO01D6FIX rank=0001416 x=1596.0 y=1415.0 length=91
    CGGCTGATGCTGCACCAACTGTACTCCCTCGCGATA

ARGUMENTS
---------

    optional arguments:
    -h, --help               show this help message and exit
    -v, --version            show version and exit
    -o OUTPUT                direct retained fragments to a file named OUTPUT (default=stdout)
    -q QSCORE                minimum per-base quality score below which a read will be split
                             or truncated (default=20)
    -l LENGTH                minimum retained fragment LENGTH (default=50)
    -m MODE                  MODE is a 3-bitmask (an integer in [0-7], default=0):
                             if the lowest bit is set, it is like passing -s;
                             if the middle bit is set, it is like passing -p;
                             and if the highest bit is set, it is like passing -a
    -s                       when encountering a low q-score, split instead of truncate
    -p                       tolerate low q-score homopolymeric regions
    -a                       tolerate low q-score ambiguous nucleotides
    -T PREFIX                if supplied, only reads with this PREFIX are retained,
                             and the PREFIX is stripped from each contributing read
    -t MISMATCH              if PREFIX is supplied, prefix matching tolerates at most
                             MISMATCH mismatches (default=0)
    -f FORMAT                output in FASTA or FASTQ format (default=FASTA)
    -j                       output run diagnostics to stderr as JSON (default is to write ASCII text)
