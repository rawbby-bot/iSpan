#ifndef __UTIL_H__
#define __UTIL_H__
#include <algorithm>
#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <cstdint>

#define LOCK(vert, lock) while (!__sync_bool_compare_and_swap(lock + vert, 0, -1))
#define UNLOCK(vert, lock) lock[vert] = 0

typedef std::int64_t vertex_t;
typedef std::int64_t index_t;
typedef double path_t;
typedef std::int64_t depth_t;
typedef std::int64_t color_t;

#define INFTY (float)10000000
#define NEGATIVE (int)-1
#define ORPHAN (unsigned char)254
#define UNVIS (long)-1
#define DEBUG 0
#define VERBOSE 0
#define OUT 1
#define OUTPUT_TIME 1

#define TRIM_TIMES 3
inline off_t
fsize(const char* filename)
{
  struct stat st;
  if (stat(filename, &st) == 0)
    return st.st_size;
  return -1;
}

#endif
