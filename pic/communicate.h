#pragma once

#include <map>
#include <algorithm>

#include "cell.h"
#include "../tools/cppitertools/reversed.hpp"
#include "../tools/wrap.h"

using std::min;
using std::max;
using iter::reversed;

namespace pic {

// expose particle memory of external neighbors
pic::ParticleBlock& get_external_data(
  int i, int j,
  pic::PicCell& cell,
  corgi::Node& grid)
{ 
  auto ind = cell.neighs(i, j); 
  uint64_t cid = grid.cellId( std::get<0>(ind), std::get<1>(ind) );
  pic::PicCell& external_cell = dynamic_cast<pic::PicCell&>( grid.getCell(cid) );

  return external_cell.container;
}

//! General communicator dealing with inter-tile particle transfer
class Communicator {
  public:




  void check_outgoing_particles( pic::PicCell& cell)
  {

    // initialize pointers to particle arrays
    int nparts = cell.container.size();

    // block limits
    auto mins = cell.mins;
    auto maxs = cell.maxs;

    double xmin = mins[0];
    double xmax = maxs[0];

    double ymin = mins[1];
    double ymax = maxs[1];

    double zmin = mins[2];
    double zmax = maxs[2];


    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( cell.container.loc(i,0) );

    double x0, y0, z0;

    // map of particles traveling beyond tile boundaries
    auto& to_others = cell.container.to_other_tiles;
    to_others.clear();


    // loop and check particles
    int n1 = 0;
    int n2 = nparts;

    int i=0,j=0,k=0; // relative indices

    for(int n=n1; n<n2; n++) {
      i = 0;
      j = 0;
      k = 0;

      x0 = loc[0][n];
      y0 = loc[1][n];
      z0 = loc[2][n];

      if(x0 <  xmin) i--; // left wrap
      if(x0 >= xmax) i++; // right wrap

      if(y0 <  ymin) j--; // bottom wrap
      if(y0 >= ymax) j++; // top wrap

      if(z0 <  zmin) k--; // back
      if(z0 >= zmax) k++; // front

      if ((i == 0) && (j == 0)) continue; //TODO: hack to make this work with 2D corgi tiles

      if (i != 0 | j != 0 | k != 0) {
        to_others.insert( std::make_pair( std::make_tuple(i,j,k), n) );
      }
    }

    // next let's see if there are any particles outflowing diagonally
    //for(auto p : to_others) {
    //  i = std::get<0>(p.first); 
    //  j = std::get<1>(p.first); 
    //  k = std::get<2>(p.first); 
    //  std::cout << "to: (" << i << "," << j << "," << k << ") n: " << p.second << '\n';
    //}

  }


  //! get incoming particles from neighboring tiles
  void get_incoming_particles( 
      pic::PicCell& cell, 
      corgi::Node& grid)
  {

    // local tile limits
    //auto mins = cell.mins;
    //auto maxs = cell.maxs;

    std::vector<double> mins = {
      grid.getXmin(),
      grid.getYmin(),
      0.0
    };

    std::vector<double> maxs = {
      grid.getXmax(),
      grid.getYmax(),
      1.0
    };



    // fetch incoming particles from neighbors around me
    int k = 0;
    for (int i=-1; i<=1; i++)
    for (int j=-1; j<=1; j++) {
    //for (int k=-1; k<=1; k++) { // TODO: hack to get 2d tiles working
      //std::cout << "from: (" << i << "," << j << "," << k << ")" << '\n';
      pic::ParticleBlock& neigh  = get_external_data(i, j, cell, grid);

      // indices as seen by sender
      std::tuple<int,int,int> nindx(-i, -j, -k);

      auto& to_others = neigh.to_other_tiles;
      for (auto&& elem : to_others) {
        
        // TODO: easy way to do this once corgi supports 3D tiles
        //if( elem.first == nindx ) {

        //}
        //TODO: collapsed z-dimension due to 2D corgi tiles
        if ( std::get<0>(elem.first) == 0 &&
             std::get<1>(elem.first) == 0 ) {
          continue;
        }
        if ( std::get<0>(elem.first) == std::get<0>(nindx) &&
             std::get<1>(elem.first) == std::get<1>(nindx) ) {

          //std::cout << " incoming particle in loop indices of: ";
          //std::cout << "(" << i << "," << j << "," << k << ")";
          //std::cout << " that is sending to: ";
          //std::cout << "(" << std::get<0>(elem.first) << "," << std::get<1>(elem.first) << "," << std::get<2>(elem.first) << ")";
          //std::cout << " and compared to mirror indices of: ";
          //std::cout << "(" << std::get<0>(nindx) << "," << std::get<1>(nindx) << "," << std::get<2>(nindx) << ")";
          //std::cout << " particle #" << elem.second << '\n';

          // get particle info and wrap location to tile coordinates
          //std::cout << "wrapping " <<
            
          //std::cout << " incoming particle in loop indices of: ";
          std::vector<double> loc = {
            wrap( neigh.loc(0, elem.second), mins[0], maxs[0] ),
            wrap( neigh.loc(1, elem.second), mins[1], maxs[1] ),
            wrap( neigh.loc(2, elem.second), mins[2], maxs[2] ),
          };
          //std::cout << "x_old: " << neigh.loc(0, elem.second) << "x_new: " << loc[0] << "\n";
          //std::cout << "y_old: " << neigh.loc(1, elem.second) << "y_new: " << loc[1] << "\n";
          //std::cout << "z_old: " << neigh.loc(2, elem.second) << "z_new: " << loc[2] << "\n";

          std::vector<double> vel = {
            neigh.vel(0, elem.second),
            neigh.vel(1, elem.second),
            neigh.vel(2, elem.second),
          };

          cell.container.add_particle(loc, vel);
        }
      }
    }
  }


  //!  Finalize communication by deleting particles that went beyond tile boundaries
  void delete_transferred_particles( pic::PicCell& cell)
  {

    std::vector<int> to_be_deleted;

    //for(auto& elem : cell.container.to_other_tiles)  {
    //  std::cout << "to be removed particle # " << elem.second << " ";
    //  std::cout << "(" << std::get<0>(elem.first) << "," << std::get<1>(elem.first) << "," << std::get<2>(elem.first) << ")";
    //  std::cout << '\n';
    //}

    // get and reverse sort deleted particles
    for(auto& elem : cell.container.to_other_tiles) to_be_deleted.push_back( elem.second );
    std::sort(to_be_deleted.begin(), to_be_deleted.end(), std::greater<int>() );

    // quick delete algorithm for indx
    //vec[idx] = vec.back();
    //vec.pop_back();

    double* loc[3];
    for( int i=0; i<3; i++)
      loc[i] = &( cell.container.loc(i,0) );

    double* vel[3];
    for( int i=0; i<3; i++)
      vel[i] = &( cell.container.vel(i,0) );

    // overwrite particles with the last one on the array and then resize the array
    int last = cell.container.size();
    for(int indx : to_be_deleted) {
      last--;

      //std::cout << "deleting " << indx << " by putting it to " << last << '\n';
      for(int i=0; i<3; i++) {
        loc[i][indx] = loc[i][last];
        vel[i][indx] = vel[i][last];
      }

    }
    cell.container.resize(last);

  }




};




} // end of namespace pic
