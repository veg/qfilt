
#ifndef SEQ_H
#define SEQ_H

#include <cstdio>
#include <string>
#include <vector>

enum state_t {
    UNKNOWN,
    ID,
    SEQUENCE,
    QUALITY
};

enum filetype_t {
    FASTA,
    QUAL,
    FASTQ
};

class pos_t
{
private:
    const char * file;
    unsigned long line;
    unsigned long col;
    std::vector<unsigned long> cols;

public:
    pos_t( const char * file );
    void get( const char *& f, long & l, long & c ) const;

    inline
    void next_col( const long ncol ) {
        col += ncol;
    }

    inline
    void next_col() {
        next_col( 1 );
    }

    inline
    void next_line() {
        line += 1;
        cols.push_back( col );
        col = 1;
    }

    inline
    void prev_col( const long ncol ) {
        col -= ncol;
    }

    inline
    void prev_col() {
        prev_col( 1 );
    }

    inline
    void prev_line() {
        line -= 1;
        col = cols[line - 1];
    }
};

class seq_t
{
public:
    std::string id;
    std::string seq;
    std::vector<long> quals;
    long length;
    seq_t();
    void clear();
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
    char fchr;
    char qchr;

    // these are for parsing
    const char * const hdr;
    const char * const sep;

    // these are permanent static buffer
    std::string qid;
    std::string qs;

public:
    parser_t( const char * );
    parser_t( const char *, const char * );
    ~parser_t();
    bool next( seq_t & );
};

#endif // SEQ_H
