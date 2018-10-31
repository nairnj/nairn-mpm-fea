#include "Read_XML/HashString.hpp"

hash_t hash_(char const* str)
{
  hash_t ret{ basis };

  while (*str) {
    ret ^= *str;
    ret *= prime;
    str++;
  }

  return ret;
}


