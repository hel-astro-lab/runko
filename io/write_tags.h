#pragma once

#include "../corgi/tags.h"

/// I/O Write modes, 
//
// Note: 
//    Naming should be Name{} name;
//
// Usage: 
//    Define functions that depend on tags as
//    func(xx, yy, ..., WriteMode::TagName)
//    and call the function as
//    yy = func(aa, bb,..., WriteMode::tagName)

namespace h5io {


class WriteMode :
  public corgi::tags::WriteMode
{
  public:

    /// write only analysis part of fields::Tile
    static struct Analysis{} analysis;
};



} // ns h5io
