#pragma once

void interpolateDevEntry( 
        toolbox::Mesh<real_short,3>& f,
        toolbox::Mesh<real_short,0>& fi,
        const std::array<int,3>& in,
        const std::array<int,3>& out
      );

void push_ebDevEntry(ffe::Tile<3>& tile);

//void add_jperpXDevEntry(fields::YeeLattice& m, ffe::SkinnyYeeLattice& dm);