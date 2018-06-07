#pragma once

#include "math.h"

/* wrap x -> [0,max) */
double wrap_max(double x, double max)
{
    /* integer math: `(max + x % max) % max` */
    return fmod(max + fmod(x, max), max);
}
/* wrap x -> [min,max) */
double wrap(double x, double min, double max)
{
    return min + wrap_max(x - min, max - min);
}
