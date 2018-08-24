#pragma once

#include <cmath>

/* wrap x -> [0,max) */
inline double wrap_max(double x, double max)
{
    return fmod(max + fmod(x, max), max);
}
/* wrap x -> [min,max) */
inline double wrap(double x, double min, double max)
{
    return min + wrap_max(x - min, max - min);
}

/* wrap x -> [0,max) */
inline float wrap_max(float x, float max)
{
    return fmod(max + fmod(x, max), max);
}
/* wrap x -> [min,max) */
inline float wrap(float x, float min, float max)
{
    return min + wrap_max(x - min, max - min);
}



/* integer math: `(max + x % max) % max` */
inline int wrap_max(int x, int max)
{
  return (max + x % max) % max;
}

/* wrap x -> [min,max) */
inline int wrap(int x, int min, int max)
{
  return min + wrap_max(x - min, max - min);
}


