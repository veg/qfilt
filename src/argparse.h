
#ifndef ARGPARSE_H
#define ARGPARSE_H

// program name
#define QFILT "qfilt"

// argument defaults
#define DEFAULT_MIN_LENGTH 50
#define DEFAULT_MIN_QSCORE 20
#define DEFAULT_MODE 0
#define DEFAULT_TAG_MISMATCH 0

typedef struct {
    const char * fastq;
    const char * fasta;
    const char * qual;
    long min_length;
    long min_qscore;
    bool split; // split not truncate
    bool hpoly; // tolerate homopolymers
    bool ambig; // tolerate ambigs ('N')
    char tag[256];
    long tag_length;
    long tag_mismatch;
} args_t;

void parse_args(args_t &, int, const char **);

#endif // ARGPARSE_H
