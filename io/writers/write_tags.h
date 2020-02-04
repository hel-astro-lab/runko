#pragma once

#include "../corgi/tags.h"

/// I/O Write modes, 
//
// Note: 
//    Naming should be Name{} name;
//
// Usage: 
//    Define functions that depend on tags as
//    func(xx, yy, ..., Write_mode::Tag_name)
//    and call the function as
//    yy = func(aa, bb,..., Write_mode::tag_name)
//class Write_mode :
//  public corgi::tags::Write_mode

namespace h5io {


class Write_mode 
{
  public:

    /// Standard/default writing mode
    static struct Standard{} standard;

    /// write only analysis part of fields::Tile
    static struct Analysis{} analysis;
};



} // ns h5io
