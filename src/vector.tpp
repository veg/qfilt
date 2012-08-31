
// #include "vector.h"

// #include <cstdlib> // included by parent: vector.h

// protected methods

template<class T>
void vector<T>::init()
{
    data = (T *) malloc(sizeof(T) * DEFAULT_SZ);
    __check_ptr(data, __FILE__, __LINE__);
    size = DEFAULT_SZ;
    len = 0;
    data[len] = NULL;
}

// public methods

template<class T>
void vector<T>::append(const T item)
{
    if (size == len) {
        long increase = size / 8;
        if (increase < DEFAULT_SZ)
            increase = DEFAULT_SZ;
        data = (T *) realloc(data, sizeof(T) * (size + increase));
        __check_ptr(data, __FILE__, __LINE__);
        size += increase;
    }
    data[len] = item;
    len += 1;
    data[len] = NULL;
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
        // don't forget the extra null byte (1 + ...)
        const long diff = 1 + len + nitem - size;
        if (diff < 0) {
            long inc = size / 8,
                 increase;
            if (inc < DEFAULT_SZ)
                inc = DEFAULT_SZ;
            increase = inc;
            while (increase < diff)
                increase += inc;
            data = (char *) realloc(data, sizeof(T) * (size + increase));
            __check_ptr(data, __FILE__, __LINE__);
            size += increase;
        }
        memcpy(&data[len], vec, sizeof(T) * nitem);
        len += nitem;
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
        extend(&vec.data[from], nitem);
}

template<class T>
void vector<T>::extend(vector & vec)
{
    extend(vec.data, vec.length());
}
