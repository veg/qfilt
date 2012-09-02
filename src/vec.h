
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef VEC_H
#define VEC_H

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
            return &path[i + 1]; 
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

inline
int __elem_cmp(const void * a, const void * b)
{
    if (*(long *) a < *(long *) b)
        return -1;
    else if (* (long *) a > *(long *) b)
        return 1;
    return 0;
}

template<class T>
class vec_t {
  private:
    void __resize(long);

  protected:
    T * data;
    long len, capacity;
  
    void init();

  public:
    vec_t();
    ~vec_t();
    long length() const;
    void append(const T);
    void clear();
    int compare(vec_t<T> &) const;
    void extend(const T *, long);
    void extend(vec_t<T> &, long, long);
    void extend(vec_t<T> &);
    void sort();
    T operator[](long i) const;
};

#include "vec.tpp"

#endif // VEC_H
