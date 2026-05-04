// Copyright 2016 - 2026, Joonas Nättilä and the hel-astro-lab contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cmath>

/* limit x -> [min,max] */
inline int limit(int x, int min, int max)
{
  x = x < min ? min : x;
  x = x > max ? max : x;
  return x;
}

/* limit x -> [min,max] */
inline float limit(float x, float min, float max)
{
  x = x < min ? min : x;
  x = x > max ? max : x;
  return x;
}

/* limit x -> [min,max] */
inline double limit(double x, double min, double max)
{
  x = x < min ? min : x;
  x = x > max ? max : x;
  return x;
}

