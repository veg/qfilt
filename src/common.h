
#ifndef COMMON_H
#define COMMON_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

inline
const char * basename(const char * path)
{
    long i;
    
    if (path == NULL)
        return NULL;

    for (i = strlen(path) - 1; i >= 0; --i) {
#ifdef _WIN32
        if (path[i] == '\\')
#else
        if (path[i] == '/')
#endif
            return path + i + 1; 
    }

    return path;
}

inline
void __check_ptr(const void * ptr, const char * file, const unsigned line)
{
    if (ptr == NULL) {
        fprintf(stderr, "\nERROR (file: %s, line: %d): memory allocation failure\n", basename(file), line);
        exit(1);
    }
}

#endif // COMMON_H
