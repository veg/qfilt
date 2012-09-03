
#ifndef STRTOK_H
#define STRTOK_H

class strtok_t {
  private:
    char * str,
         * ptr;
    const long len;

  public:
    strtok_t(const char *);
    ~strtok_t();
    char * next(const char *);
};

#endif // STRTOK_H
