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

// vec must be sorted
void fprint_vector_stats(FILE * file, vec_t<long> & vec, const char * hdr)
{
    double sum = 0.,
           var = 0.,
           mean = 0.,
           median = 0.;
    long i,
         min = 0,
         two5 = 0,
         ninetyseven5 = 0,
         max = 0;

    for (i = 0; i < vec.length(); ++i)
    {
        sum += vec[i];
        var += vec[i] * vec[i];
    }

    // remember, i == vec.length()

    if (vec.length()) {
        var = (var - (sum * sum) / vec.length()) / (vec.length() - 1);
        mean = sum / i;
        median = (i % 2) ? 1.0 * vec[i / 2] : 0.5 * (vec[i / 2] + vec[i / 2 - 1]);
        min = vec[0];
        two5 = vec[long(0.025 * i)];
        ninetyseven5 = vec[long(0.975 * i)];
        max = vec[i - 1];
    }

    fprintf(file, "%s\n"
                  "    mean:                %g\n"
                  "    median:              %g\n"
                  "    variance             %g\n"
                  "    standard deviation:  %g\n"
                  "    min:                 %ld\n"
                  "    2.5%%:                %ld\n"
                  "    97.5%%:               %ld\n"
                  "    max:                 %ld\n",
        hdr,
        mean,
        median,
        var,
        sqrt(var),
        max,
        two5,
        ninetyseven5,
        max
    );
}

inline
void parse_error(const char * msg, pos_t & pos) {
    fprintf(stderr, "\nERROR (file: %s, line: %ld, column: %ld): %s\n", pos.file, pos.line, pos.col, msg);
    exit(1);
}

#define MALFUNCTION() \
{ \
    fprintf(stderr, "\nERROR (file: %s, line: %d): state machine malfunction\n", __FILE__, __LINE__); \
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
    args.min_length = atoi(min_length);
    if (args.min_length <= 1)
    {
        fprintf(stderr, "Expected a positive integer to specify the minimum high quality run length (50 is a good default): had %s\n", min_length);
        exit(1);
    }
}

void parse_minqscore(args_t & args, const char * min_qscore)
{
    args.min_qscore = atoi(min_qscore);
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
    args.tag[0] = '\0';
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

class seq_t {
  public:
    str_t * id;
    str_t * seq;
    vec_t<long> * quals;
    long length;
    seq_t() : length(0) {
        id = new str_t();
        seq = new str_t();
        quals = new vec_t<long>();
    }
    void clear() {
        id->clear();
        seq->clear();
        quals->clear();
        length = 0;
    }
    ~seq_t() {
        if (id) delete id;
        if (seq) delete seq;
        if (quals) delete quals;
    }
};

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

class machine_t
{
  private:
    FILE * fasta;
    FILE * fastq;
    FILE * qual;
    pos_t fpos;
    pos_t qpos;
    state_t fstate;
    state_t qstate;

    // these are permanent static buffer
    str_t * qid;
    str_t * qs;

    void init() {
        qid = new str_t();
        qs = new str_t();
    }
  public:
    machine_t(const char *);
    machine_t(const char *, const char *);
    bool next(seq_t &);
    ~machine_t() {
        if (fasta) fclose(fasta);
        if (fastq) fclose(fastq);
        if (qual) fclose(qual);
        if (qid) delete qid;
        if (qs) delete qs;
    }
};

machine_t::machine_t(const char * fastq_file):
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

    init();
}

machine_t::machine_t(const char * fasta_file, const char * qual_file):
    fastq(NULL),
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

    init();
}

inline
char fgetc_skip_ws(FILE * file, pos_t & pos)
{
    char c;
    do {
        c = fgetc(file);
        if (c == '\n') {
            pos.line += 1;
            pos.col = 0;
        }
        else if (c != EOF) {
            pos.col += 1;
        }
    } while (c == ' ' || c == '\t' || c == '\n' || c == '\r');
    return c;
}

int extend_until(str_t & str, const char until, FILE * file, pos_t & pos, bool trim)
{
    char buf[256];
    long total = 0;

    while (fgets(buf, 256, file)) {
        long len = strlen(buf),
             i;

        for (i = 0; i < len; ++i) {
            const char c = buf[i];
            // if we hit either a newline or until
            // extend str up until this point
            if (c == '\n' || c == until) {
                str.extend(buf, i);
                total += i;
            }
            // if we hit until,
            // rewind by what's left over
            if (c == until) {
                fseek(file, i - len + 1, SEEK_CUR);
                return total;
            }
            // if we hit newline,
            // break out for EOL handling below
            else if (c == '\n')
                break;
            // otherwise increment the position column
            else
                pos.col += 1;
        }

        // handle EOL
        if (buf[len - 1] == '\n') {
            pos.line += 1;
            pos.col = 0;

            // if we're not trimming, append the newline
            if (!trim)
                str.append('\n');
            // otherwise, strip all whitespace,
            // collect the next character,
            // and if it's until or EOF, return;
            // if not, append it and continue
            else {
                const char c = fgetc_skip_ws(file, pos);
                if (c == until || c == EOF)
                    return total;
                else {
                    str.append(c);
                    total += 1;
                }
            }
        }
    }
    return total;
}

// by default trim whitespace from ends
int extend_until(str_t & str, const char until, FILE * file, pos_t & pos)
{
    return extend_until(str, until, file, pos, true);
}

bool machine_t::next(seq_t & seq)
{
    const char hdr = fastq ? '@' : '>',
               sep = fastq ? '+' : '>';
    FILE * file = fastq ? fastq : fasta;
    filetype_t filetype = fastq ? FASTQ : FASTA;
    pos_t * pos = &fpos;
    state_t * state = &fstate;

begin:
    do {
        switch (*state) {
            case UNKNOWN: {
                const char c = fgetc_skip_ws(file, *pos);
                if (c == hdr)
                    *state = ID;
                else if (c == EOF)
                    return false;
                else
                    parse_error("malformed file", *pos);
                break;
            }
            case ID: {
                switch (filetype) {
                    case FASTA:
                    case FASTQ: {
                        int nelem = extend_until(*seq.id, '\n', file, *pos);
                        if (nelem < 1)
                            parse_error("malformed file: missing ID", *pos);
                        *state = SEQUENCE;
                        break;
                    }
                    case QUAL: {
                        int nelem = extend_until(*qid, '\n', file, *pos);
                        if (nelem < 1)
                            parse_error("malformed file: missing ID", *pos);
                        // clear the qid after use
                        qid->clear();
                        *state = QUALITY;
                        break;
                    }
                    default:
                        MALFUNCTION();
                }
                break;
            }
            case SEQUENCE: {
                switch (filetype) {
                    case FASTA:
                    case FASTQ: {
                        int nelem;
                        nelem = extend_until(*seq.seq, sep, file, *pos);
                        if (nelem < 1)
                            parse_error("malformed file: missing sequence", *pos);
                        if (filetype == FASTA) {
                            if (!feof(file))
                                fseek(file, -1, SEEK_CUR);
                            pos->col -= 1;
                            *state = UNKNOWN;
                        }
                        else // FASTQ
                            *state = QUALITY;
                        break;
                    }
                    default:
                        MALFUNCTION();
                }
                break;
            }
            case QUALITY: {
                int nelem = extend_until(*qs, hdr, file, *pos, false);
                if (nelem < 1)
                    parse_error("malformed file: missing quality scores", *pos);

                if (filetype == QUAL) {
                    char * buf = strtok(qs->c_str(), " \t\n\r");
                    while (buf != NULL) {
                        seq.quals->append(atoi(buf));
                        buf = strtok(NULL, " \t\n\r");
                    }
                }
                else { // FASTQ
                    int i;
                    for (i = 0; i < qs->length(); ++i) {
                        // encoding: chr(phred+33)
                        seq.quals->append(long((*qs)[i]) - 33);
                    }
                }
                // clear the qual data after use
                qs->clear();
                // reset the state to UNKNOWN and just prior to the header
                if (!feof(file))
                    fseek(file, -1, SEEK_CUR);
                pos->col -= 1;
                *state = UNKNOWN;
                break;
            }
            default:
                MALFUNCTION();
        }
    } while (*state != UNKNOWN);

    if (qual && file == fasta) {
        file = qual;
        filetype = QUAL;
        pos = &qpos;
        state = &qstate;
        goto begin;
    }

    if (file == qual && seq.seq->length() != seq.quals->length()) {
        char buf[512];
        sprintf(
            buf,
            "malformed file: sequence length (%ld) does not match the number of quality scores (%ld)",
            seq.seq->length(),
            seq.quals->length()
        );
        parse_error(buf, *pos);
    }
    else
        seq.length = seq.seq->length();

    return true;
}


/*---------------------------------------------------------------------------------------------------- */

int main(int argc, const char * argv[])
{
    args_t args;

    machine_t * machine = NULL;

    seq_t seq = seq_t();

    long ncontrib = 0;

    /*
    for (int i = 0; i < 256; ++i)
        char_lookup[i] = -1;

    for (int i = 0; i < valid_char_count; ++i)
        char_lookup[valid_chars[i]] = i;
    */

    parse_args(argc, argv, args);

    if (args.fastq)
        machine = new machine_t(args.fastq);
    else
        machine = new machine_t(args.fasta, args.qual);

    vec_t<long> read_lengths = vec_t<long>();
    vec_t<long> fragment_lengths = vec_t<long>();

    for (; machine->next(seq); seq.clear()) {
        // maxto is the maximum value of "to",
        // NOT THE UPPER BOUND
        const long maxto = seq.length - args.min_length;
        long nfragment = 0,
             to = 0;

        read_lengths.append(seq.length);

        // compare the sequence prefix to the tag,
        // if it matches by at least tag_mismatch,
        // keep the sequence, otherwise discard
        if (args.tag_length) {
            long mismatch = 0;

            for (to = 0; to < args.tag_length; ++to) {
                if ((*seq.seq)[to] != args.tag[to])
                    mismatch += 1;
            }

            if (mismatch > args.tag_mismatch)
                continue;
        }

        // if we're splitting,
        // continue the following process until we reach the end of the sequence,
        // but only continue if there's enough left to produce a minimum-sized fragment
        while (true) {
            char buf[256] = "\n";
            long from = 0,
                 i = 0,
                 nambigs = 0;

            // push through the sequence until the quality score meets the minimum
            while ((to <= maxto) && ((*seq.quals)[to] < args.min_qscore)) {
                to += 1;
            }

            // if we don't have enough length left,
            // skip to the next sequence
            if (to > maxto) {
                fprintf(stderr, "skipping %s because insufficient sequence length remains (to: %ld, maxto: %ld)\n", seq.id->c_str(), to, maxto);
                break;
            }

            // begin with positive quality score
            from = to;

            // build a read until we hit a low quality score,
            // that is, unless we're skipping Ns or retaining homopolymers
            for (; to < seq.length; ++to) {
                char curr = (*seq.seq)[to],
                     last = -1;

                if ((*seq.quals)[to] < args.min_qscore) {
                    // if homopolymer (toupper => equiv), continue
                    if (args.homo && toupper(last) == toupper(curr))
                        goto next;
                    // if skipping Ns, continue
                    else if (args.skipN && (curr == 'N' || curr == 'n')) {
                        nambigs += 1;
                        goto next;
                    }
                    // otherwise, ABORT!!!
                    else
                        break;
                }
next:
                last = curr;
            }

            // "to" is now the upper bound

            // if our fragment isn't long enough,
            // skip to the next fragment
            if (to - from - nambigs < args.min_length) {
                fprintf(stderr, "skipping %s because of insufficient fragment length (from: %ld, to: %ld, nambigs: %ld)\n", seq.id->c_str(), from, to, nambigs);
                continue;
            }

            // give the read a fragment identifier (if there's more than 1)
            if (nfragment > 0)
                sprintf(buf, " fragment=%ld\n", nfragment);
            // if its the first fragment, count the contributing read
            else
                ncontrib += 1;

            // print the read ID (with fragment identifier)
            fprintf(stdout, ">%s%s", seq.id->c_str(), buf);

            // print the read sequence
            for (i = from; i < to; i += 60) {
                const int nitem = (to - i < 60) ? to - i : 60;
                strncpy(buf, seq.seq->c_str() + i, nitem);
                buf[nitem] = '\0';
                fprintf(stdout, "%s\n", buf);
            }
#if 0
            // for printing quality scores
            fprintf(stdout, "+\n");
            for (i = from; i < to; ++i) {
                char s[] = " ";
                if (i == from)
                    s[0] = '\0';
                fprintf(stdout, "%s%ld", s, (*seq.quals)[i]);
            }
            fprintf(stdout, "\n");
#endif
            fragment_lengths.append(to - from - nambigs);

            if (!args.split)
                break;

            // only increment fragment identifier after printing
            nfragment += 1;
        }
    }

    read_lengths.sort();
    fragment_lengths.sort();

    fprintf(stderr, "run settings:\n");
    if (args.fasta)
        fprintf(stderr,
                    "    input fasta:         %s\n"
                    "    input qual:          %s\n",
            args.fasta,
            args.qual
        );
    else
        fprintf(stderr,
                    "    input fastq:         %s\n",
            args.fastq
        );
    fprintf(stderr, "    min q-score:         %ld\n"
                    "    min fragment length: %ld\n"
                    "    run mode:            %d (%s/%s/%s)\n",
        args.min_qscore,
        args.min_length,
        ((args.split ? 1 : 0) | (args.homo ? 2 : 0) | (args.skipN ? 4 : 0)),
        args.split ? "split" : "truncate",
        args.homo ? "retain homopolymers" : "don't retain homopolymers",
        args.skipN ? "skip ambigs" : "don't skip ambigs"
    );
    if (args.tag_length)
        fprintf(stderr,
                    "    5' tag:              %s\n"
                    "    max tag mismatches:  %ld\n",
            args.tag,
            args.tag_mismatch
        );

    fprintf(stderr, "\nrun summary:\n"
                    "    original reads:      %ld\n"
                    "    contributing reads:  %ld\n"
                    "    retained fragments:  %ld\n",
        read_lengths.length(),
        ncontrib,
        fragment_lengths.length()
    );

    fprint_vector_stats(stderr, read_lengths, "\noriginal read length distribution:");
    fprint_vector_stats(stderr, fragment_lengths, "\nretained fragment length distribution:");

    delete machine;

    return 0;
}
