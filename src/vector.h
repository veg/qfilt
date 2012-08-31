
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define DEFAULT_SZ 1024L

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
  private:
    void __resize(long);

  protected:
    T * data;
    long len, capacity;
  
    void init();

  public:
    vector() {
        init();
    }
    long length() const {
        return len;
    }
    void append(const T);
    void clear();
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
        if (data != NULL)
            free(data);
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
};

#include "vector.tpp"
