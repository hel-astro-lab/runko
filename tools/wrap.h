#pragma once

#include <cmath>

/* wrap x -> [0,max) */
double wrap_max(double x, double max)
{
    return fmod(max + fmod(x, max), max);
}
/* wrap x -> [min,max) */
double wrap(double x, double min, double max)
{
    return min + wrap_max(x - min, max - min);
}

/* integer math: `(max + x % max) % max` */
int wrap_max(int x, int max)
{
  return (max + x % max) % max;
}

/* wrap x -> [min,max) */
int wrap(int x, int min, int max)
{
  return min + wrap_max(x - min, max - min);
}


