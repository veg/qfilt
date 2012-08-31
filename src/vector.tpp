
// #include "vector.h"

// #include <cstdlib> // included by parent: vector.h

// private methods

template<class T>
inline
void vector<T>::__resize(long rlen)
{
    if (rlen >= capacity) {
        const long incr = (capacity / 8 < DEFAULT_SZ) ? DEFAULT_SZ : capacity / 8,
                   rcap = incr * ((rlen + 1) / incr + 1);
        fprintf(stderr, "resizing capacity from: %ld to: %ld\n", capacity, rcap);
        void * ptr = realloc(data, sizeof(T) * rcap);
        __check_ptr(ptr, __FILE__, __LINE__);
        data = reinterpret_cast<T *>(ptr);
        capacity = rcap;
    }
}

// protected methods

template<class T>
void vector<T>::init()
{
    data = reinterpret_cast<T *>(calloc(DEFAULT_SZ, sizeof(T)));
    __check_ptr(data, __FILE__, __LINE__);
    capacity = DEFAULT_SZ;
    len = 0;
    data[len] = NULL;
}

// public methods

template<class T>
void vector<T>::append(const T item)
{
    __resize(len + 1);
    data[len] = item;
    len += 1;
    data[len] = NULL;
}

template<class T>
void vector<T>::clear()
{
    if (capacity > DEFAULT_SZ) {
        void * ptr = realloc(data, sizeof(T) * DEFAULT_SZ);
        __check_ptr(ptr, __FILE__, __LINE__);
        data = reinterpret_cast<T *>(ptr);
    }
    len = 0;
    memset(data, NULL, sizeof(T) * capacity);
}

template<class T>
int vector<T>::compare(vector & vec) const
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
void vector<T>::extend(const T * vec, long nitem)
{
    if (nitem <= 0)
        return;
    else {
        __resize(len + nitem);
        fprintf(stderr, "memcpying %ld items into vector with %ld remaining slots\n", nitem, capacity - len);
        memcpy(data + len, vec, sizeof(T) * nitem);
        len += nitem;
        fprintf(stderr, "len: %ld, cap: %ld\n", len, capacity);
        data[len] = NULL;
    }
}

template<class T>
void vector<T>::extend(vector & vec, long from, long to)
{
    long nitem = to - from;
    if (vec.length() <= 0 || nitem <= 0)
        return;
    else
        extend(vec.data + from, nitem);
}

template<class T>
void vector<T>::extend(vector & vec)
{
    extend(vec.data, vec.length());
}
