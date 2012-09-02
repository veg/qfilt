
#include <cstdio>

#include "vec.h"

#ifndef SEQ_H
#define SEQ_H

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

class pos_t {
  public:
    const char * file;
    unsigned long line;
    unsigned long col;
    pos_t(const char * file, unsigned long line, unsigned long col) :
        file(file), line(line), col(col) {}
};

class parser_t
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
    parser_t(const char *);
    parser_t(const char *, const char *);
    bool next(seq_t &);
    ~parser_t() {
        if (fasta) fclose(fasta);
        if (fastq) fclose(fastq);
        if (qual) fclose(qual);
        if (qid) delete qid;
        if (qs) delete qs;
    }
};

#endif // SEQ_H
