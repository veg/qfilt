
#ifndef VEC_H
#define VEC_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "common.hpp"

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
