
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define DEFAULT_SZ 1024L

inline void __check_ptr(const void * ptr, const char * file, const unsigned line)
{
    if (ptr == NULL) {
        fprintf(stderr, "\nERROR (file: %s, line %d): memory allocation failure\n", file, line);
        exit(1);
    }
}

int __elem_cmp(const void * a, const void * b)
{
    if (*(long *) a < *(long *) b)
        return -1;
    else if (* (long *) a > *(long *) b)
        return 1;
    return 0;
}

template<class T>
class vector {
  protected:
    T * data;
    long len, size;
  
    void init();
    void destroy() {
        free(data);
    }

  public:
    vector() {
        init();
    }
    long length() const {
        return len;
    }
    void append(const T);
    void clear() {
        len = 0;
        data[0] = NULL;
    }
    int compare(vector<T> &) const;
    void extend(const T *, long);
    void extend(vector<T> &, long, long);
    void extend(vector<T> &);
    void sort() {
        qsort(data, len, sizeof(T), __elem_cmp); 
    }
    T operator[](long i) const {
        return data[i];
    };
    ~vector() {
        destroy();
    }
};

class string : public vector<char> {
  public:
    using vector<char>::extend;
    string() {
        init();
    }
    string(const char * str) {
        init();
        extend(str);
    }
    char * c_str() const {
        return data;
    }
    void extend(const char * str) {
        extend(str, (long) strlen(str));
    }
    void lower() {
        long i;
        for (i = 0; i < len; ++i)
            data[i] = tolower(data[i]);
    }
    void upper() {
        long i;
        for (i = 0; i < len; ++i)
            data[i] = toupper(data[i]);
    }
    ~string() {
        destroy();
    }
};

#include "vector.tpp"
