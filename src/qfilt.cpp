#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "vector.h"

void init_genrand(unsigned long);
unsigned long genrand_int32(void);

static const char * const usage = "qfilt <fasta (.fna) file> <quality scores (.qual) file> <min phred score> <min run length> <filtering mode> [5' tag] [tag mismatch: default 0]\n"
                                  "filtering mode can be 0 through 7 (see README for details)\n"
                                  "Viewed as a bit mask, if the highest bit is set, all 'N' characters as skipped\n"
                                  "If the second bit is set, low q-scores in homopolymers are tolerated.\n"
                                  "If the lowest bit is set, then reads are SPLIT at low score points, otherwise they are truncated.\n";

/*
static const char * const valid_chars = "ACGTNacgtn";

static long char_lookup[256];
*/

typedef struct {
    const char * fastq;
    const char * fasta;
    const char * qual;
    long min_length;
    long min_qscore;
    bool split;
    bool homo;
    bool skipN;
    char tag[256];
    long tag_length;
    long tag_mismatch;
} args_t;

// this has to be a class because C++ is stupid
class pos_t {
  public:
    const char * file;
    unsigned long line;
    unsigned long col;
    pos_t(const char * file, unsigned long line, unsigned long col);
};

// idiotic constructor because C++ continues to be dumb
pos_t::pos_t(const char * file, unsigned long line, unsigned long col) : file(file), line(line), col(col) {}

void fprint_vector_stats(FILE * file, vector<long> & vec, const char * hdr)
{
    double sum = 0.,
           var = 0.;
    long i;

    for (i = 0; i < vec.length(); ++i)
    {
        sum += vec[i];
        var += vec[i] * vec[i];
    }

    // remember, i == vec.length()

    var = (var - (sum * sum) / vec.length()) / (vec.length() - 1);

    fprintf(file, "%s\n"
                  "    mean:   %g\n"
                  "    median: %g\n"
                  "    var:    %g\n"
                  "    stdev:  %g\n"
                  "    min:    %ld\n"
                  "    2.5%%:   %ld\n"
                  "    97.5%%:  %ld\n"
                  "    max:    %ld\n"
                  "\n",
        hdr,
        sum / i,
        (i % 2) ? 1.0 * vec[i / 2] : 0.5 * (vec[i / 2] + vec[i / 2 - 1]),
        var,
        sqrt(var),
        vec[0],
        vec[(long) (0.025 * i)],
        vec[(long) (0.975 * i)],
        vec[i - 1]
    );
}


inline void parse_error(const char * msg, pos_t & pos) 
{
    fprintf(stderr, "\nERROR (file %s, line %ld, column %ld): %s\n", pos.file, pos.line, pos.col, msg); 
    exit(1);
}

#define MALFUNCTION() \
{ \
    fprintf(stderr, "\nERROR (file %s, line %d): state machine malfunction\n", __FILE__, __LINE__); \
    exit(1); \
}

void parse_qual(args_t & args, const char * fasta, const char * qual)
{
    if (args.fastq)
    {
        fprintf(stderr, "--qual incompatible with --fastq\n");
        exit(1);
    }

    args.fasta = fasta;
    args.qual = qual;
}

void parse_fastq(args_t & args, const char * fastq)
{
    if (args.fasta || args.qual)
    {
        fprintf(stderr, "--fastq incompatible with --qual\n");
        exit(1);
    }

    args.fastq = fastq; 
}

void parse_minlength(args_t & args, const char * min_length)
{
    args.min_length = (long) atoi(min_length);
    if (args.min_length <= 1)
    {
        fprintf(stderr, "Expected a positive integer to specify the minimum high quality run length (50 is a good default): had %s\n", min_length);
        exit(1);
    }
}

void parse_minqscore(args_t & args, const char * min_qscore)
{
    args.min_qscore = (long) atoi(min_qscore);
    if (args.min_qscore < 0)
    {
        fprintf(stderr, "Expected a non-negative integer to specify the minimum q-score (20 is a good default): had %s\n", min_qscore);
        exit(1);
    }
}

void parse_mode(args_t & args, const char * mode)
{
    int m = atoi(mode);
    args.split = (m & 1);
    args.homo = (m & 2);
    args.skipN = (m & 4);
}

void parse_tag(args_t & args, const char * tag)
{
    int nvar = sscanf(tag, "%256s", args.tag);
    if (nvar != 1)
    {
        fprintf(stderr, "Failed to process the tag argument %s\n", tag);
        exit(1);
    }
    args.tag_length = strlen(tag);
}

void parse_tagmismatch(args_t & args, const char * tag_mismatch)
{
    args.tag_mismatch = atoi(tag_mismatch);
}

void parse_args(int argc, const char * argv[], args_t & args)
{
    int i;

    // defaults
    args.fastq = NULL;
    args.fasta = NULL;
    args.qual = NULL;
    args.min_length = 50;
    args.min_qscore = 20;
    args.split = false;
    args.homo = false;
    args.skipN = false;
    args.tag[0] = '\n';
    args.tag_length = 0;
    args.tag_mismatch = 0;

    // skip arg[0], it's just the program name
    for (i = 1; i < argc; ++i)
    {
        const char * arg = argv[i];

        if (arg[0] == '-' && arg[1] == '-')
        {
                 if (!strcmp(&arg[2], "fastq")) parse_fastq(args, argv[++i]);
            else if (!strcmp(&arg[2], "qual")) { parse_qual(args, argv[i+1], argv[i+2]); i += 2; }
            else if (!strcmp(&arg[2], "minlength")) parse_minlength(args, argv[++i]);
            else if (!strcmp(&arg[2], "minqscore")) parse_minqscore(args, argv[++i]);
            else if (!strcmp(&arg[2], "mode")) parse_mode(args, argv[++i]);
            else if (!strcmp(&arg[2], "tag")) parse_tag(args, argv[++i]); 
            else if (!strcmp(&arg[2], "tagmismatch")) parse_tagmismatch(args, argv[++i]);
            else { fprintf(stderr, "Unknown argument %s\n", arg); exit(1); }
        }
        else if (arg[0] == '-')
        {
                 if (!strcmp(&arg[1], "F")) parse_fastq(args, argv[++i]); 
            else if (!strcmp(&arg[1], "Q")) { parse_qual(args, argv[i+1], argv[i+2]); i += 2; }
            else if (!strcmp(&arg[1], "l")) parse_minlength(args, argv[++i]); 
            else if (!strcmp(&arg[1], "q")) parse_minqscore(args, argv[++i]); 
            else if (!strcmp(&arg[1], "m")) parse_mode(args, argv[++i]); 
            else if (!strcmp(&arg[1], "T")) parse_tag(args, argv[++i]); 
            else if (!strcmp(&arg[1], "t")) parse_tagmismatch(args, argv[++i]); 
            else { fprintf(stderr, "Unknown argument %s\n", arg); exit(1); }
        }
        else
        {
            fprintf(stderr, "Unknown argument %s\n", arg);
            exit(1);
        }
    }

    if (!args.fastq && (!args.fasta || !args.qual))
    {
        fprintf(stderr, "Neither --fastq FASTQ file nor --qual FASTA QUAL files provided, ergo nothing to filter\n");
        exit(1);
    }
}

typedef struct {
    string id;
    string seq;
    vector<long> quals;
    long length;
} seq_t;

void clear_seq(seq_t & seq)
{
    seq.id.clear();
    seq.seq.clear();
    seq.quals.clear();
    seq.length = 0;
}

enum state_t
{
    UNKNOWN,
    ID,
    SEQUENCE,
    QUALITY
};

enum filetype_t
{
    FASTA,
    QUAL,
    FASTQ
};

class state_machine
{
  private:
    FILE * fastq;
    FILE * fasta;
    FILE * qual;
    pos_t fpos;
    pos_t qpos;
    state_t fstate;
    state_t qstate;
  public:
    state_machine(const char *);
    state_machine(const char *, const char *);
    bool next(seq_t &);
};

state_machine::state_machine(const char * fastq_file):
    fasta(NULL),
    qual(NULL),
    fpos(pos_t(fastq_file, 0, 0)),
    qpos(pos_t(NULL, 0, 0)),
    fstate(UNKNOWN),
    qstate(UNKNOWN)
{
    fastq = fopen(fastq_file, "rb");
    if (!fastq)
    {
        fprintf(stderr, "\nERROR: failed to open the FASTQ file %s\n", fastq_file);
        exit(1);
    }
}

state_machine::state_machine(const char * fasta_file, const char * qual_file):
    fasta(NULL),
    qual(NULL),
    fpos(pos_t(fasta_file, 0, 0)),
    qpos(pos_t(qual_file, 0, 0)),
    fstate(UNKNOWN),
    qstate(UNKNOWN)
{
    fasta = fopen(fasta_file, "rb");
    if (!fasta)
    {
        fprintf(stderr, "\nERROR: failed to open the FASTA file %s\n", fasta_file);
        exit(1);
    }

    qual = fopen(qual_file, "rb");
    if (!qual)
    {
        fprintf(stderr, "\nERROR: failed to open the QUAL file %s\n", qual_file);
        exit(1);
    }
}

inline char skip_whitespace(FILE * file, pos_t & pos)
{
    char c;
    do {
        c = fgetc(file);
        if (c == '\n') {
            pos.line += 1;
            pos.col = 0;
        }
        else
            pos.col += 1;
    } while (c == ' ' || c == '\t' || c == '\n' || c == '\r');
    return c;
}

int extend_until(string & str, const char until, FILE * file, pos_t & pos)
{
    char buf[1024],
         c = '\0';
    long len = 0, total = 0;

    while (fgets(buf, 1024, file)) {
        len = strlen(buf);
        if (buf[len - 1] == '\n') {
            pos.line += 1;
            pos.col = 0;
            // extend the string, but without the newline
            str.extend((const char *) buf, len - 1);
            total += len - 1;
            // if newline is our until ...
            if (buf[len -1 ] == c)
                break;
            // keep excising whitespace
            c = skip_whitespace(file, pos);
            // check the until character
            if (c == until)
                break;
            // if not, append and try again
            else {
                str.append(c);
                total += 1;
            }
        }
        else {
            str.extend(buf);
            total += len;
        }
    }
    return total;
}

bool state_machine::next(seq_t & seq)
{
    const char hdr = fastq ? '@' : '>',
               sep = fastq ? '+' : '>';
    FILE * & file = fastq ? fastq : fasta;
    filetype_t filetype = fastq ? FASTQ : FASTA;
    pos_t & pos = fpos;
    state_t & state = fstate;

begin:
    do {
        switch (state) {
            case UNKNOWN: {
                const char c = skip_whitespace(file, pos);
                if (c == hdr)
                    state = ID;
                else if (c != EOF)
                    parse_error("malformed file\n", pos);
                else
                    return false;
                break;
            }
            case ID: {
                int nelem;
                nelem = extend_until(seq.id, '\n', file, pos);
                if (nelem < 1)
                    parse_error("malformed file: missing ID\n", pos);
                state = SEQUENCE;
                break;
            }
            case SEQUENCE: {
                switch (filetype) {
                    case FASTA:
                    case FASTQ: {
                        int nelem;
                        nelem = extend_until(seq.seq, sep, file, pos);
                        if (nelem < 1)
                            parse_error("malformed file: missing ID\n", pos);
                        if (filetype == FASTA) {
                            pos.col -= 1;
                            fseek(file, -1, SEEK_CUR);
                            fstate = UNKNOWN;
                        }
                        else // FASTQ
                            fstate = QUALITY;
                        break;
                    }
                    case QUAL: {
                        break;
                    }
                    default:
                        MALFUNCTION();
                }
                break;
            }
            case QUALITY: {
                int nelem;
                string qs;

                nelem = extend_until(qs, hdr, file, pos);
                if (nelem < 1)
                    parse_error("malformed file: missing quality scores\n", pos);

                if (filetype == QUAL) {
                    char * buf = strtok(qs.c_str(), " \t\n\r");
                    while (buf != NULL) {
                        seq.quals.append(atoi(buf));
                        buf = strtok(NULL, " \t\n\r");
                    }
                }
                else { // FASTQ
                    int i;
                    for (i = 0; i < qs.length(); ++i) {
                        // encoding: chr(phred+33)
                        seq.quals.append(((long) qs[i]) - 33);
                    }
                }
                // reset the state to unkonw and just prior to the header
                pos.col -= 1;
                fseek(file, -1, SEEK_CUR);
                state = UNKNOWN;
                break;
            }
            default:
                MALFUNCTION();
        }
    } while (state != UNKNOWN);

    if (file == fasta) {
        file = qual;
        filetype = QUAL;
        pos = qpos;
        state = qstate;
        goto begin;
    }

    if (file == qual && seq.seq.length() != seq.quals.length())
        parse_error("malformed file: sequence length does not match the number of quality scores\n", pos);
    else
        seq.length = seq.seq.length();

    return true;
}

            
/*---------------------------------------------------------------------------------------------------- */

int main(int argc, const char * argv[])
{
    args_t args;

    state_machine * machine = NULL;

    seq_t * seq = NULL;

    /*
    for (int i = 0; i < 256; ++i)
        char_lookup[i] = -1;

    for (int i = 0; i < valid_char_count; ++i)
        char_lookup[valid_chars[i]] = i;
    */

    parse_args(argc, argv, args);
    
    if (args.fastq)
        machine = new state_machine(args.fastq);
    else
        machine = new state_machine(args.fasta, args.qual);

    vector<long> read_lengths;
    vector<long> contig_lengths;

    while (machine->next(*seq)) {
        char buf[256] = "";
        long to = 0, contig = 0;

        read_lengths.append(seq->length);

        // compare the sequence prefix to the tag,
        // if it matches by at least tag_mismatch,
        // keep the sequence, otherwise discard
        if (args.tag_length) {
            long mismatch = 0;

            for (to = 0; to < args.tag_length; ++to) {
                if (seq->seq[to] != args.tag[to])
                    mismatch += 1;
            }

            if (mismatch > args.tag_mismatch)
                continue;
        }

        // if we're splitting,
        // continue the following process until we reach the end of the sequence,
        // but only continue if there's enough left to produce a minimum-sized contig
        do {
            bool abort = false;
            char curr = '\0',
                 last = '\0';
            long from = 0,
                 i = 0;

            // reset the read contig identifier
            memset(buf, '\0', 256);

            // push through the sequence until the quality score meets the minimum
            for (; (to < (seq->length - args.min_length)) && (seq->quals[to] < args.min_qscore); ++to)

            // begin with positive quality score
            from = to;

            // build a read until we hit a low quality score,
            // that is, unless we're skipping Ns or retaining homopolymers
            for (; (to < seq->length) && (!abort); ++to) {
                curr = seq->seq[to];

                if (seq->quals[to] < args.min_qscore) {
                    // if homopolymer, continue
                    if (!args.homo || last != curr)
                        goto next;
                    // if skipping Ns, continue
                    else if (args.skipN && (curr == 'N' || curr == 'n'))
                        goto next;
                    // otherwise, ABORT!!!
                    else
                        abort = true;
                }
next:
                last = curr;
            }

            if (from - to)
            // give the read a CONTIG identifier (if there's more than 1) 
            if (contig > 0)
                sprintf(buf, " CONTIG=%ld", contig);

            // print the read ID (with CONTIG identifier)
            fprintf(stdout, ">%s%s\n", seq->id.c_str(), buf);

            // print the read sequence
            for (i = from; i < to; i += 80) {
                strncpy(buf, seq->seq.c_str() + i, 80);
                buf[80] = '\0';
                fprintf(stdout, "%s", buf);
            }

            contig_lengths.append(to - from);

            contig += 1;
        } while (args.split && to < (seq->length - args.min_length));
      
        clear_seq(*seq);
    }

    read_lengths.sort();
    contig_lengths.sort();

    fprintf(stderr, "run settings:\n");
    if (args.fasta)
        fprintf(stderr,
                    "    input fasta: %s\n"
                    "    input qual: %s\n",
            args.fasta,
            args.qual
        );
    else
        fprintf(stderr,
                    "    input fastq: %s\n",
            args.fastq
        );
    fprintf(stderr, "    q-score min: %ld\n"
                    "    min contig length: %ld\n"
                    "    run mode: %d\n"
                    "    5' tag: %s\n"
                    "    max tag mismatches: %ld\n"
                    "\n",
        args.min_qscore,
        args.min_length,
        ((args.split ? 1 : 0) | (args.homo ? 2 : 0) | (args.skipN ? 4 : 0)),
        args.tag,
        args.tag_mismatch
    );

    fprint_vector_stats(stderr, read_lengths, "original read length statistics:");
    fprint_vector_stats(stderr, contig_lengths, "retained contig length statistics:");

    delete seq;
    delete machine;

    return 0;
}
