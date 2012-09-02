
#include "seq.h"

// helper classes and functions

#define MALFUNCTION() \
{ \
    fprintf(stderr, "\nERROR (file: %s, line: %d): state machine malfunction\n", __FILE__, __LINE__); \
    exit(1); \
}

// pos_t methods

pos_t::pos_t(const char * file) : file(file), line(1), col(1)
{
    cols = new vec_t<unsigned long>; 
}

pos_t::~pos_t()
{
    delete cols;
}

void pos_t::get(const char * & f, long & l, long & c)
{
    f = file;
    l = line;
    c = col;
}

// seq_t methods

seq_t::seq_t() : length(0)
{
    id = new str_t();
    seq = new str_t();
    quals = new vec_t<long>();
}

seq_t::~seq_t()
{
    if (id) delete id;
    if (seq) delete seq;
    if (quals) delete quals;
}

void seq_t::clear()
{
    id->clear();
    seq->clear();
    quals->clear();
    length = 0;
}

// helper functions for parser_t

inline
void parse_error(const char * msg, pos_t & pos)
{
    const char * file;
    long line, col;
    pos.get(file, line, col);
    fprintf(stderr, "\nERROR (file: %s, line: %ld, column: %ld): %s\n", file, line, col, msg);
    exit(1);
}

inline
char fgetc_skip_ws(FILE * file, pos_t & pos)
{
    char c;
    do {
        c = fgetc(file);
        if (c == '\n') {
            pos.next_line();
        }
        else if (c != EOF) {
            pos.next_col();
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
                pos.next_col();
        }

        // handle EOL
        if (buf[len - 1] == '\n') {
            pos.next_line();

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
inline
int extend_until(str_t & str, const char until, FILE * file, pos_t & pos)
{
    return extend_until(str, until, file, pos, true);
}

// parser methods

// protected parser methods

void parser_t::init()
{
    qid = new str_t();
    qs = new str_t();
}

// public parser methods

parser_t::parser_t(const char * fastq_file):
    fasta(NULL),
    qual(NULL),
    fpos(pos_t(fastq_file)),
    qpos(pos_t(NULL)),
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

parser_t::parser_t(const char * fasta_file, const char * qual_file):
    fastq(NULL),
    fpos(pos_t(fasta_file)),
    qpos(pos_t(qual_file)),
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

parser_t::~parser_t()
{
    if (fasta)
        fclose(fasta);
    if (fastq)
        fclose(fastq);
    if (qual)
        fclose(qual);
    if (qid)
        delete qid;
    if (qs)
        delete qs;
}

bool parser_t::next(seq_t & seq)
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
                            pos->prev_col();
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
                pos->prev_col();
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
