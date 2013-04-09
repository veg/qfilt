
#include "common.hpp"
#include "seq.hpp"
#include "strtok.hpp"

// helper classes and functions

#define MALFUNCTION() \
{ \
    fprintf( stderr, "\nERROR (file: %s, line: %d): state machine malfunction\n", basename( __FILE__ ), __LINE__ ); \
    exit( 1 ); \
}

#define PARSE_ERROR( pos, msg, ... ) \
{ \
    const char * file; \
    long line, col; \
    (pos).get( file, line, col ); \
    fprintf( stderr, "\nERROR (file: %s, line: %ld, column: %ld): " msg "\n", file, line, col , ## __VA_ARGS__ ); \
    exit( 1 ); \
}

#define IS_WHITESPACE( chr ) ( ( chr ) == ' ' || ( chr ) == '\t' || ( chr ) == '\r' || ( chr ) == '\n' ) 

// pos_t ------------------------------------------------------------------------------------------------------------ //

pos_t::pos_t( const char * file ) : file( file ), line( 1 ), col( 1 )
{
}

void pos_t::get( const char *& f, long & l, long & c ) const
{
    f = file;
    l = line;
    c = col;
}

// seq_t ------------------------------------------------------------------------------------------------------------ //

seq_t::seq_t() : length( 0 )
{
}

void seq_t::clear()
{
    id.clear();
    seq.clear();
    quals.clear();
    length = 0;
}

// parser_t helper functions ---------------------------------------------------------------------------------------- //

inline
void skip_ws( FILE * file, pos_t & pos, char & chr )
{
    if ( chr == EOF )
        return;

    if ( !IS_WHITESPACE( chr ) )
        return;
    
    for ( chr = fgetc( file ); chr != EOF; chr = fgetc( file ) ) {
        if ( !IS_WHITESPACE( chr ) )
            break;
        else if ( chr == '\n' )
            pos.next_line();
        else if ( chr != '\r' )
            pos.next_col();
    }
}

char extend_until( std::string & str, const char * delim, FILE * file, pos_t & pos, bool trim )
{
    char chr;

    for ( chr = fgetc( file ); chr != EOF; chr = fgetc( file ) ) {
        if ( chr == '\n' )
            pos.next_line();
        else if ( chr != '\r' )
            pos.next_col();

        if ( strrchr( delim, chr ) )
            break;
        else if ( !trim || !strrchr( "\r\n", chr ) )
            str.push_back( chr );
    }
        
    if ( trim )
        skip_ws( file, pos, chr );

    return chr;
}

// by default trim whitespace from the beginning of lines
inline
char extend_until( std::string & str, const char * delim, FILE * file, pos_t & pos )
{
    return extend_until( str, delim, file, pos, true );
}

// parser_t --------------------------------------------------------------------------------------------------------- //

const char chr[] = ">\0@\0+";

parser_t::parser_t( const char * fastq_file ):
    fasta( NULL ), 
    qual( NULL ), 
    fpos( pos_t( fastq_file ) ),
    qpos( pos_t( NULL ) ),
    fstate( UNKNOWN ),
    qstate( UNKNOWN ),
    fchr( ' ' ),
    qchr( ' ' ),
    hdr( chr + 2 ),
    sep( chr + 4 )
{
    if ( !strcmp( fastq_file, "-" ) )
        fastq = stdin;
    else
        fastq = fopen( fastq_file, "rb" ); 

    if ( !fastq ) {
        fprintf( stderr, "\nERROR: failed to open the FASTQ file %s\n", fastq_file );
        exit( 1 );
    }
}

parser_t::parser_t( const char * fasta_file, const char * qual_file ):
    fastq( NULL ), 
    fpos( pos_t( fasta_file ) ),
    qpos( pos_t( qual_file ) ),
    fstate( UNKNOWN ),
    qstate( UNKNOWN ),
    fchr( ' ' ),
    qchr( ' ' ),
    hdr( chr + 0 ),
    sep( chr + 0 )
{

    if ( !strcmp( fasta_file, "-" ) )
        fasta = stdin;
    else
        fasta = fopen( fasta_file, "rb" ); 

    if ( !fasta ) {
        fprintf( stderr, "\nERROR: failed to open the FASTA file %s\n", fasta_file );
        exit( 1 );
    }

    if ( qual_file ) {
        if ( !strcmp( qual_file, "-" ) )
            qual = stdin; 
        else
            qual = fopen( qual_file, "rb" ); 
    
        if ( !qual ) {
            fprintf( stderr, "\nERROR: failed to open the QUAL file %s\n", qual_file );
            exit( 1 );
        }
    }
    else
        qual = NULL;

    if ( fasta == stdin && qual == stdin ) {
        fprintf( stderr, "\nERROR: FASTA and QUAL file may not both be STDIN" );
        exit( 1 );
    }
}

parser_t::~parser_t()
{
    if ( fasta && fasta != stdin )
        fclose(fasta);

    if ( fastq && fastq != stdin )
        fclose(fastq);

    if ( qual && qual != stdin )
        fclose(qual);
}

bool parser_t::next( seq_t & seq )
{
    FILE * file = fastq ? fastq : fasta;
    filetype_t filetype = fastq ? FASTQ : FASTA;
    pos_t * pos = &fpos;
    state_t * state = &fstate;
    char * chr = &fchr;
begin:

    do {
        switch ( *state ) {
        case UNKNOWN: {
            skip_ws( file, *pos, *chr );

            if ( strchr( hdr, *chr ) )
                *state = ID;
            else if ( *chr == EOF )
                return false;
            else
                PARSE_ERROR( *pos, "malformed file: %c", *chr )
            
            break;
        }

        case ID: {
            std::string & str = ( filetype == QUAL ) ? qid : seq.id;
            *chr = extend_until( str, "\r\n", file, *pos );

            if ( str.length() < 1 )
                PARSE_ERROR( *pos, "malformed file: missing ID" )
                
            if ( filetype == QUAL ) {
                qid.clear();
                *state = QUALITY;
            }
            else
                *state = SEQUENCE;

            break;
        }

        case SEQUENCE: {
            switch ( filetype ) {
            case FASTA:
            case FASTQ: {
                if ( *chr != EOF )
                    seq.seq.push_back( *chr );

                *chr = extend_until( seq.seq, sep, file, *pos );

                if ( !seq.seq.length() )
                    PARSE_ERROR( *pos, "malformed file: missing sequence" )
            
                if ( filetype == FASTA )
                    *state = UNKNOWN;
                else { // FASTQ
                    // skip the whitespace after the + separator
                    *chr = fgetc( file );
                    skip_ws( file, *pos, *chr );
                    *state = QUALITY;
                }

                break;
            }

            default:
                MALFUNCTION()
            }

            break;
        }

        case QUALITY: {
            // from the previous skip_ws
            if ( *chr != EOF )
                qs.push_back( *chr );

            // FASTQ files permit '@' to appear as a valid quality value (31),
            // so look for a newline instead of the hdr
            *chr = extend_until( qs, ( filetype == QUAL ) ? hdr : "\r\n", file, *pos, false );

            if ( qs.length() < 1 )
                PARSE_ERROR( *pos, "malformed file: missing quality scores" )
                
            if ( filetype == QUAL ) {
                char * buf = NULL;

                strtok_t tok( qs.c_str() );

                while ( ( buf = tok.next( " \t\r\n" ) ) )
                    seq.quals.push_back( atoi( buf ) );
            }
            else { // FASTQ
                size_t i;

                for ( i = 0; i < qs.length(); ++i ) {
                    // encoding: chr(phred+33)
                    seq.quals.push_back( long( qs[i] ) - 33 );
                }
            }

            // clear the qual data after use
            qs.clear();

            *state = UNKNOWN;
            break;
        }

        default:
            MALFUNCTION()
        }
    }
    while ( *state != UNKNOWN );

    if ( qual && file == fasta ) {
        file = qual;
        filetype = QUAL;
        pos = &qpos;
        state = &qstate;
        chr = &qchr;
        goto begin;
    }    

    if ( qual && seq.seq.length() != seq.quals.size() )
        PARSE_ERROR(
            *pos,
            "malformed file: sequence length (%ld) does not match the number of quality scores (%ld)",
            seq.seq.length(),
            seq.quals.size()
        )
    else
        seq.length = seq.seq.length();

    return true;
}
