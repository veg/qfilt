
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "argparse.h"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

const char usage[] =
    "usage: " QFILT " [-h] "
    "[-q MIN_QSCORE] "
    "[-l MIN_LENGTH] "
    "[-m MODE] "
    "[-T PREFIX] "
    "[-t MAX_MISMATCH] "
    "(-F FASTQ | -Q FASTA QUAL)\n";

const char help_msg[] =
    "filter sequencing data using some simple heuristics\n"
    "\n"
    "required arguments:\n"
    "  -F FASTQ                 FASTQ file\n"
    "  -Q FASTA QUAL            FASTA and QUAL files\n"
    "\n"
    "optional arguments:\n"
    "  -h, --help               show this help message and exit\n"
    "  -q MIN_QSCORE            minimum q-score (default=" TO_STR(DEFAULT_MIN_QSCORE) ")\n"
    "  -l MIN_LENGTH            minimum retained fragment length (default=" TO_STR(DEFAULT_MIN_LENGTH) ")\n"
    "  -m MODE                  MODE is a 3-bitmask (an integer in [0-7], default=" TO_STR(DEFAULT_MODE) "):\n"
    "                           if the lowest bit is set, a low q-score causes reads to be split,\n"
    "                           otherwise they are truncated;\n"
    "                           if the second bit is set, low q-score homopolymers are tolerated;\n"
    "                           and if the third bit is set, low q-score 'N's are tolerated\n"
    "  -T PREFIX                if PREFIX is supplied, only reads with this prefix are retained\n"
    "  -t MAX_MISMATCH          if PREFIX is supplied, prefix matching tolerates at most\n"
    "                           MAX_MISMATCH mismatches (default=" TO_STR(DEFAULT_TAG_MISMATCH) ")\n";

inline
void help()
{
    fprintf(stderr, "%s\n%s", usage, help_msg);
    exit(1);
}

#define PARSE_ERROR(msg, args...) \
{ \
    fprintf(stderr, "%s" QFILT ": error: " msg "\n", usage , ##args); \
    exit(1); \
}

void parse_qual(args_t & args, const char * fasta, const char * qual)
{
    if (args.fastq)
        PARSE_ERROR("-Q and -F are mutually exclusive");

    args.fasta = fasta;
    args.qual = qual;
}

void parse_fastq(args_t & args, const char * fastq)
{
    if (args.fasta || args.qual)
        PARSE_ERROR("-F and -Q are mutually exclusive");

    args.fastq = fastq;
}

void parse_minlength(args_t & args, const char * min_length)
{
    args.min_length = atoi(min_length);
    if (args.min_length <= 1)
        PARSE_ERROR("minimum length expected a positive integer, had: %s", min_length);
}

void parse_minqscore(args_t & args, const char * min_qscore)
{
    args.min_qscore = atoi(min_qscore);
    if (args.min_qscore < 0)
        PARSE_ERROR("min q-score expected a non-negative integer, had: %s", min_qscore);
}

void parse_mode(args_t & args, const char * mode)
{
    int m = atoi(mode);

    if (m < 0 || m > 7)
        PARSE_ERROR("mode must be an integer in [0, 7], had: %s", mode);

    args.split = (m & 1);
    args.hpoly = (m & 2);
    args.ambig = (m & 4);
}

void parse_tag(args_t & args, const char * tag)
{
    int nvar = sscanf(tag, "%256s", args.tag);
    if (nvar != 1)
        PARSE_ERROR("failed to process tag argument %s", tag);
    args.tag_length = strlen(tag);
}

void parse_tagmismatch(args_t & args, const char * tag_mismatch)
{
    args.tag_mismatch = atoi(tag_mismatch);
    if (args.tag_mismatch < 0)
        PARSE_ERROR("maximum tag mismatch expected non-negative integer, had: %s", tag_mismatch);
}

void parse_args(args_t & args, int argc, const char * argv[])
{
    int i;

    // defaults
    args.fastq = NULL;
    args.fasta = NULL;
    args.qual = NULL;
    args.min_length = DEFAULT_MIN_LENGTH;
    args.min_qscore = DEFAULT_MIN_QSCORE;
    args.tag[0] = '\0';
    args.tag_length = 0;
    args.tag_mismatch = DEFAULT_TAG_MISMATCH;

    // handle the mode separately
    parse_mode(args, TO_STR(DEFAULT_MODE));

    // skip arg[0], it's just the program name
    for (i = 1; i < argc; ++i) {
        const char * arg = argv[i];

        if (arg[0] == '-' && arg[1] == '-') {
            if (!strcmp(&arg[2], "help"))
                help();
#if 0
                 if (!strcmp(&arg[2], "fastq")) parse_fastq(args, argv[++i]);
            else if (!strcmp(&arg[2], "qual")) { parse_qual(args, argv[i+1], argv[i+2]); i += 2; }
            else if (!strcmp(&arg[2], "minlength")) parse_minlength(args, argv[++i]);
            else if (!strcmp(&arg[2], "minqscore")) parse_minqscore(args, argv[++i]);
            else if (!strcmp(&arg[2], "mode")) parse_mode(args, argv[++i]);
            else if (!strcmp(&arg[2], "tag")) parse_tag(args, argv[++i]);
            else if (!strcmp(&arg[2], "tagmismatch")) parse_tagmismatch(args, argv[++i]);
#endif
            else
                PARSE_ERROR("unknown argument: %s", arg);
        }
        else if (arg[0] == '-') {
                 if (!strcmp(&arg[1], "F")) parse_fastq(args, argv[++i]);
            else if (!strcmp(&arg[1], "Q")) { parse_qual(args, argv[i+1], argv[i+2]); i += 2; }
            else if (!strcmp(&arg[1], "l")) parse_minlength(args, argv[++i]);
            else if (!strcmp(&arg[1], "q")) parse_minqscore(args, argv[++i]);
            else if (!strcmp(&arg[1], "m")) parse_mode(args, argv[++i]);
            else if (!strcmp(&arg[1], "T")) parse_tag(args, argv[++i]);
            else if (!strcmp(&arg[1], "t")) parse_tagmismatch(args, argv[++i]);
            else if (!strcmp(&arg[1], "h")) help();
            else
                PARSE_ERROR("unknown argument: %s", arg);
        }
        else
            PARSE_ERROR("unknown argument: %s", arg);
    }

    if (!args.fastq && (!args.fasta || !args.qual))
        PARSE_ERROR("missing required argument -F FASTQ or -Q FASTA QUAL");
}
