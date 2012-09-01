
// #include "vector.h"

// #include <cstdlib> // included by parent: vector.h

#define DEFAULT_SZ 1024L

// private methods

template<class T>
inline
void vec_t<T>::__resize(long rlen)
{
    // don't forget the null terminating byte
    rlen += 1;

    if (rlen > capacity) {
        const long incr = (capacity / 8 < DEFAULT_SZ) ? DEFAULT_SZ : capacity / 8,
                   rcap = incr * (rlen / incr + 1);
        void * ptr = realloc(data, sizeof(T) * rcap);
        __check_ptr(ptr, __FILE__, __LINE__);
        data = reinterpret_cast<T *>(ptr);
        capacity = rcap;
    }
}

// protected methods

template<class T>
void vec_t<T>::init()
{
    data = reinterpret_cast<T *>(malloc(sizeof(T) * DEFAULT_SZ));
    __check_ptr(data, __FILE__, __LINE__);
    capacity = DEFAULT_SZ;
    len = 0;
    data[len] = NULL;
}

// public methods

template<class T>
void vec_t<T>::append(const T item)
{
    __resize(len + 1);
    data[len] = item;
    len += 1;
    if (len >= capacity)
        fprintf(stderr, "len >= capacity");
    data[len] = NULL;
}

template<class T>
void vec_t<T>::clear()
{
    if (capacity > DEFAULT_SZ) {
        void * ptr = realloc(data, sizeof(T) * DEFAULT_SZ);
        __check_ptr(ptr, __FILE__, __LINE__);
        data = reinterpret_cast<T *>(ptr);
        capacity = DEFAULT_SZ;
    }
    len = 0;
    data[len] = NULL;
}

template<class T>
int vec_t<T>::compare(vec_t<T> & vec) const
{
    long minlen, i;

    if (len <= vec.length())
        minlen = len;
    else
        minlen = vec.length();

    for (i = 0; i < minlen; ++i) {
        if (data[i] < vec[i])
            return -1;
        else if (data[i] > vec[i])
            return 1;
    }

    if (len == vec.length())
        return 0;

    return (len < vec.length) ? -1 : 1;
}

template<class T>
void vec_t<T>::extend(const T * vec, long nitem)
{
    if (nitem <= 0)
        return;
    else {
        __resize(len + nitem);
        memcpy(data + len, vec, sizeof(T) * nitem);
        len += nitem;
        data[len] = NULL;
    }
}

template<class T>
void vec_t<T>::extend(vec_t<T> & vec, long from, long to)
{
    long nitem = to - from;
    if (vec.length() <= 0 || nitem <= 0)
        return;
    else
        extend(vec.data + from, nitem);
}

template<class T>
void vec_t<T>::extend(vec_t<T> & vec)
{
    extend(vec.data, vec.length());
}
