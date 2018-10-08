#pragma once

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
//
class WriteMode
{
  public:

    /// write only analysis part of fields::Tile
    static struct Analysis{} analysis;
};






