
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <stdlib.h>

#include "argparse.hpp"
#include "seq.hpp"

#if 0
static const char * const valid_chars = "ACGTNacgtn";
static long char_lookup[256];
#endif

const size_t BUF_LEN = 60;

// vec must be sorted
void fprint_vector_stats( FILE * file, std::vector<size_t> & vec, const char * hdr, bool do_json )
{
    double sum = 0.,
           var = 0.,
           mean = 0.,
           median = 0.;
           
    long min = 0,
         two5 = 0,
         ninetyseven5 = 0,
         max = 0;
         
    size_t i;

    for ( i = 0; i < vec.size(); ++i ) {
        sum += vec[i];
        var += vec[i] * vec[i];
    }

    // remember, i == vec.size()

    if ( vec.size() ) {
        var = ( var - ( sum * sum ) / vec.size() ) / ( vec.size() - 1 );
        mean = sum / i;
        median = ( i % 2 ) ? 1.0 * vec[i / 2] : 0.5 * ( vec[i / 2] + vec[i / 2 - 1] );
        min = vec[0];
        two5 = vec[long( 0.025 * i )];
        ninetyseven5 = vec[long( 0.975 * i )];
        max = vec[i - 1];
    }

    if (do_json) {
        fprintf( file, ",\n\t\"%s\": {"
                 "\n\t\t\"mean\":                %g,"
                 "\n\t\t\"median\":              %g,"
                 "\n\t\t\"variance\":            %g,"
                 "\n\t\t\"standard deviation\":  %g,"
                 "\n\t\t\"min\":                 %ld,"
                 "\n\t\t\"2.5%%\":                %ld,"
                 "\n\t\t\"97.5%%\":               %ld,"
                 "\n\t\t\"max\":                 %ld}",
                 hdr,
                 mean,
                 median,
                 var,
                 sqrt( var ),
                 min,
                 two5,
                 ninetyseven5,
                 max
                );
   }
    else {
        fprintf( file, "\n%s\n"
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
                 sqrt( var ),
                 min,
                 two5,
                 ninetyseven5,
                 max
               );
        }
}

// main ------------------------------------------------------------------------------------------------------------- //

int main( int argc, const char * argv[] )
{
    argparse::args_t args = argparse::args_t( argc, argv );
    seq::parser_t * parser = NULL;
    seq::seq_t seq = seq::seq_t();
    long ncontrib = 0L;
    
    long      total_bases = 0L,
              q_over10    = 0L,
              q_over20    = 0L,
              q_over30    = 0L;
    
    double    q_score_sum = 0.0;
    
    // initialize the parser
    if ( args.fastq )
        parser = new seq::parser_t( args.fastq );
    else
        parser = new seq::parser_t( args.fasta, args.qual );

    if ( !parser ) {
        fprintf( stderr, "\nERROR: failed to initialize parser\n" );
        exit( 1 );
    }

#if 0

    for ( int i = 0; i < 256; ++i )
        char_lookup[i] = -1;

    for ( int i = 0; i < valid_char_count; ++i )
        char_lookup[valid_chars[i]] = i;

#endif
    std::vector<size_t> read_lengths;
    std::vector<size_t> fragment_lengths;

    for ( ; parser->next( seq ); seq.clear() ) {
        // maxto is the maximum value of "to",
        // NOT THE UPPER BOUND
        
        if (seq.length == 0) continue;
        
        const size_t maxto = seq.length - args.min_length;
        size_t nfragment = 0,
               to = 0;

        read_lengths.push_back( seq.length );
        total_bases += seq.length;
        
        for (size_t i = 0; i < seq.length; ++i ) {
            q_score_sum += seq.quals[i];
            if (seq.quals[i] >= 10L) {
                q_over10 ++;
                if (seq.quals[i] >= 20L) {
                    q_over20++;
                    if (seq.quals[i] >= 30L) {
                        q_over30++;
                    }
                }
            }
        }
        

        if ( seq.length < args.min_length )
            continue;

        // compare the sequence prefix to the tag,
        // if it matches by at least tag_mismatch,
        // keep the sequence, otherwise discard
        if ( args.tag_length ) {
            size_t mismatch = 0;

            if ( maxto < args.tag_length )
                continue;

            for ( to = 0; to < args.tag_length; ++to ) {
                // tolower -> case insensitive
                if ( toupper( seq.seq[to] ) != toupper( args.tag[to] ) )
                    mismatch += 1;
            }

            if ( mismatch > args.tag_mismatch )
                continue;
        }
        
        // gather read stats
        
        

        if ( args.punch ) {
            size_t i;
            unsigned long nambigs = 0UL;

            //char buf[BUF_LEN + 1];
            //buf[BUF_LEN] = '\0';

            // print the read ID
            
            if (seq.length >= args.min_length) {

              char * buffer = new char[seq.length+1];

            
              for ( i = 0; to < seq.length; ++i, ++to ) {
                  if ( seq.quals[to] < args.min_qscore ) {
                      buffer[i] = args.punch;
                      nambigs ++;
                      if (nambigs >= args.remove_count) {
                        break;
                      }
                  }
                  else
                      buffer[i] = seq.seq[to];
              }
              
              if (to == seq.length) {
                buffer[to] = '\0';
                
              // print the remaining portion of the sequence
                ncontrib ++;
                
                fragment_lengths.push_back( seq.length );
              
                fprintf(
                  args.output,
                  "%c%s\n",
                  ( args.format == argparse::FASTQ ) ? '@' : '>',
                  seq.id.c_str()
                  );

                 fprintf(
                    args.output,
                    ( args.format == argparse::FASTQ ) ? "%s" : "%s\n",
                    buffer
                    );
                 
                 
                if ( args.format == argparse::FASTQ ) {
                    fprintf( args.output, "\n+\n" );
                    for ( i = 0; i < to; i += BUF_LEN ) {
                        char buf[BUF_LEN + 1];
                        const int nitem = ( to - i < BUF_LEN ) ? to - i : BUF_LEN;
                        for ( int j = 0; j < nitem; ++j )
                            buf[j] = ( char ) ( seq.quals[i + j] + 33 );
                        buf[nitem] = '\0';
                        fprintf( args.output, "%s", buf );
                    }
                    fprintf( args.output, "\n" );
                }
              }
              delete [] buffer; 
            }
        }
        // if we're splitting,
        // continue the following process until we reach the end of the sequence,
        // but only continue if there's enough left to produce a minimum-sized fragment
        else while ( true ) {
            size_t from = 0,
                   i = 0,
                   nambigs = 0;

            // push through the sequence until the quality score meets the minimum
            while ( ( to <= maxto ) && ( seq.quals[to] < args.min_qscore ) ) {
                to += 1;
            }

            // if we don't have enough length left,
            // skip to the next sequence
            if ( to > maxto )
                break;

            // begin with positive quality score
            from = to;

            // build a read until we hit a low quality score,
            // that is, unless we're skipping Ns or retaining homopolymers
            for ( ; to < seq.length; ++to ) {
                char curr = seq.seq[to],
                     last = -1;

                if ( seq.quals[to] < args.min_qscore ) {
                    // if homopolymer (toupper -> case insensitive), continue (last == curr)
                    if ( args.hpoly && toupper( last ) == toupper( curr ) )
                        continue;
                    // if skipping Ns, continue (without assigning last)
                    else if ( args.ambig && ( curr == 'N' || curr == 'n' ) ) {
                        nambigs += 1;
                        continue;
                    }
                    // otherwise, ABORT!!!
                    else
                        break;
                }

                last = curr;
            }

            // "to" is now the upper bound

            // if our fragment isn't long enough,
            // skip to the next fragment
            if ( to - from - nambigs < args.min_length )
                continue;

            // print the read ID
            fprintf(
                args.output,
                "%c%s",
                ( args.format == argparse::FASTQ ) ? '@' : '>',
                seq.id.c_str()
                );

            // print the fragment identifier
            if ( nfragment > 0 )
                fprintf( args.output, " fragment=%ld\n", nfragment + 1 );
            else {
                fprintf( args.output, "\n" );
                // if it's the first fragment,
                // count the contributing read
                ncontrib += 1;
            }

            // print the read sequence
            for ( i = from; i < to; i += BUF_LEN ) {
                char buf[BUF_LEN + 1];
                const size_t nitem = ( to - i < BUF_LEN ) ? to - i : BUF_LEN;
                strncpy( buf, seq.seq.c_str() + i, nitem );
                buf[nitem] = '\0';
                fprintf(
                    args.output,
                    ( args.format == argparse::FASTQ ) ? "%s" : "%s\n",
                    buf
                    );
            }

            if ( args.format == argparse::FASTQ ) {
                fprintf( args.output, "\n+\n" );
                for ( i = from; i < to; i += BUF_LEN ) {
                    char buf[BUF_LEN + 1];
                    const int nitem = ( to - i < BUF_LEN ) ? to - i : BUF_LEN;
                    for ( int j = 0; j < nitem; ++j )
                        buf[j] = ( char ) ( seq.quals[i + j] + 33 );
                    buf[nitem] = '\0';
                    fprintf( args.output, "%s", buf );
                }
                fprintf( args.output, "\n" );
            }
#if 0
            // for printing quality scores
            fprintf( args.output, "+\n" );

            for ( i = from; i < to; ++i ) {
                char s[] = " ";

                if ( i == from )
                    s[0] = '\0';

                fprintf( args.output, "%s%ld", s, ( *seq.quals )[i] );
            }

            fprintf( args.output, "\n" );
#endif
            fragment_lengths.push_back( to - from - nambigs );

            if ( !args.split )
                break;

            // only increment fragment identifier after printing
            nfragment += 1;
        }
    }

    if ( args.json ) {
        fprintf( stderr, "{\"Settings\":{\n\t" );
        if ( args.fasta )
            fprintf( stderr,
                "\"fasta\": \"%s\",\n\t"
                "\"qual\":  \"%s\",\n\t",
                args.fasta->path,
                args.qual->path
                );
        else
            fprintf( stderr,
                "\"fastq\": \"%s\",\n\t",
                args.fastq->path
                );

        fprintf( stderr,
            "\"min q-score\": %ld,\n\t"
            "\"min fragment length\": %ld,\n\t",
            args.min_qscore,
            args.min_length
        );
            
        if (args.punch) {
           fprintf( stderr,
                 "\"punch low scores with\":    \"%c\", \n\t"
                 "\"skip sequence if more than\":  %ld \n\t",
                  args.punch,
                  args.remove_count
                 );
        
        } else {
          fprintf( stderr,
              "\"on low scores\": \"%s\",\n\t"
              "\"on homopolymers\": \"%s\",\n\t"
              "\"on ambiguities\": \"%s\"",
              args.split ? "split" : "truncate",
              args.hpoly ? "tolerate homopolymers" : "don't tolerate homopolymers",
              args.ambig ? "tolerate ambigs" : "don't tolerate ambigs"
              );

        }

        if ( args.tag_length )
            fprintf( stderr,
                ",\n\t\"tag\": \"%s\","
                ",\n\t\"max tag mismatches\":  %ld",
                args.tag,
                args.tag_mismatch
                );

       fprintf( stderr,
            "},\n"
            "\"run summary\":{"
            "\n\t\"total bases\":      %ld,"
            "\n\t\"original reads\":      %ld,"
            "\n\t\"q10\":      %g,"
            "\n\t\"q20\":      %g,"
            "\n\t\"q30\":      %g,"
            "\n\t\"mean q-score\":      %g,"
            "\n\t\"contributing reads\":  %ld,"
            "\n\t\"retained fragments\":  %ld",
            total_bases,
            read_lengths.size(),
            q_over10 / (double) total_bases,
            q_over20 / (double) total_bases,
            q_over30 / (double) total_bases,
            q_score_sum / (double) total_bases,
            ncontrib,
            fragment_lengths.size()
            );
    } else {
        fprintf( stderr, "run settings:\n" );

        if ( args.fasta )
            fprintf( stderr,
                     "    input fasta:         %s\n"
                     "    input qual:          %s\n",
                     args.fasta->path,
                     args.qual->path
                   );
        else
            fprintf( stderr,
                     "    input fastq:         %s\n",
                     args.fastq->path
                   );

        fprintf( stderr,
                 "    min q-score:         %ld\n"
                 "    min fragment length: %ld\n" ,
                  args.min_qscore,
                  args.min_length
                );
                 
        if ( args.punch ) {
           fprintf( stderr,
                 "    punch low scores with:    %c \n"
                 "    skip sequence if more than:  %ld \n",
                  args.punch,
                  args.remove_count
                 );
         
        } else {
        
          fprintf( stderr,
                 "    run mode:            %d (%s/%s/%s)\n",
                 ( ( args.split ? 1 : 0 ) | ( args.hpoly ? 2 : 0 ) | ( args.ambig ? 4 : 0 ) ),
                 args.split ? "split" : "truncate",
                 args.hpoly ? "tolerate homopolymers" : "don't tolerate homopolymers",
                 args.ambig ? "tolerate ambigs" : "don't tolerate ambigs"
               );
        }

        if ( args.tag_length )
            fprintf( stderr,
                     "    5' tag:              %s\n"
                     "    max tag mismatches:  %ld\n",
                     args.tag,
                     args.tag_mismatch
                   );

        fprintf( stderr,
                 "\n"
                 "run summary:\n"
                 "    total bases       :  %ld\n"
                 "    original reads    :  %ld\n"
                 "    q10               :  %g\n"
                 "    q20               :  %g\n"
                 "    q30               :  %g\n"
                 "    mean q-score      :  %g\n"
                 "    contributing reads:  %ld\n"
                 "    retained fragments:  %ld\n",
                 total_bases,
                 read_lengths.size(),
                 q_over10 / (double) total_bases,
                 q_over20 / (double) total_bases,
                 q_over30 / (double) total_bases,
                 q_score_sum / (double) total_bases,
                 ncontrib,
                 fragment_lengths.size()
               );
        // print original read length and retained fragment length statistics
    }

    std::sort( read_lengths.begin(), read_lengths.end() );
    std::sort( fragment_lengths.begin(), fragment_lengths.end() );

    fprint_vector_stats( stderr, read_lengths, "original read length distribution:" , args.json);
    fprint_vector_stats( stderr, fragment_lengths, "retained fragment length distribution:", args.json );

    if ( args.json )
        fprintf( stderr, "\n\t}\n}\n");

    delete parser;

    return 0;
}
