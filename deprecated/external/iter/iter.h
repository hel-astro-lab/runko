#ifndef ITER_H
#define ITER_H
#pragma once

#include <tuple>
#include <unordered_map>
#include <map>
#include <optional>

#include <execinfo.h>
#include <stdio.h>



// vectorization is broken at the point where we add values back to the grid.
// This auxiliary function tries to hide that complexity.
template<typename T, typename S>
 inline void atomic_add(T& lhs, S rhs) 
{
    //NOTE: need to use #pragma omp atomic if vectorizing these
    #pragma omp atomic update
    lhs += static_cast<S>(rhs);
}



namespace {

class ExecPlan{
public:
    std::string stream;
    std::optional<std::tuple<int, int, int>> block;

    ExecPlan(){
        stream = "default";
    }
    ExecPlan(std::string stream_)
    {
        stream = stream_;
    }
    ExecPlan(std::tuple<int, int, int> block_, std::string stream_)
    {
        stream = stream_;
        block = block_;
    }
    ExecPlan(std::tuple<int, int, int> block_)
    {
        stream = "default";
        block = block_;
    }

};

class UniIter{
    public:

    class UniIterHost{
        template<int Dims, class F, class... Args>
        static void iterateND(F fun, int xMax, int yMax, int zMax, Args& ... args)
        {
            //
        }


    public:

        template<class F, class... Args>
        static void iterate2D(F fun, int xMax, int yMax, Args& ... args)
        {
          #pragma omp parallel for
          for (int y = 0; y < yMax; y++) {
            #pragma omp simd
            for (int x = 0; x < xMax; x++) {
              fun(x, y, args...);
            }
          }
        }

        template<class F, class... Args>
        static void iterate2D_nonvec(F fun, int xMax, int yMax, Args& ... args)
        {
          for (int y = 0; y < yMax; y++) {
            for (int x = 0; x < xMax; x++) {
              fun(x, y, args...);
            }
          }
        }

        template<class F, class... Args>
        static void iterate3D(F fun, int xMax, int yMax, int zMax, Args& ... args)
        {
            #pragma omp parallel for
            for (int z = 0; z < zMax; z++) {
                for (int y = 0; y < yMax; y++) {
                    #pragma omp simd
                    for (int x = 0; x < xMax; x++) {
                        fun(x, y, z, args...);
                    }
                }
            }
        }

        template<class F, class... Args>
        static void iterate3D_nonvec(F fun, int xMax, int yMax, int zMax, Args& ... args)
        {
            for (int z = 0; z < zMax; z++) {
                for (int y = 0; y < yMax; y++) {
                    for (int x = 0; x < xMax; x++) {
                        fun(x, y, z, args...);
                    }
                }
            }
        }

        template<class F, class... Args>
        static void iterate(F fun, int max, Args& ... args)
        {
            #pragma omp parallel for simd
            for (int x = 0; x < max; x++)
            {
                fun(x, args...);
            }
        }
    };

public:

    template<class F, class... Args>
    static void iterate2D(F fun, int xMax, int yMax, Args& ... args)
    {
        //
        //std::cout << "3d iterate" << std::endl;
        UniIterHost::iterate2D(fun, xMax, yMax, args...);
    }

    template<class F, class... Args>
    static void iterate2D(F fun, ExecPlan plan, int xMax, int yMax, Args& ... args)
    {
        //std::cout << "2d iterate w plan" << std::endl;
        UniIterHost::iterate2D(fun, xMax, yMax, args...);
    }


    template<class F, class... Args>
    static void iterate3D(F fun, int xMax, int yMax, int zMax, Args& ... args)
    {
        //
        //std::cout << "3d iterate" << std::endl;
        UniIterHost::iterate3D(fun, xMax, yMax, zMax, args...);
    }

    template<class F, class... Args>
    static void iterate3D(F fun, ExecPlan plan, int xMax, int yMax, int zMax, Args& ... args)
    {
        //
        //std::cout << "3d iterate w plan" << std::endl;
        UniIterHost::iterate3D(fun, xMax, yMax, zMax, args...);
    }
    
    template<class F, class... Args>
    static void iterate(F fun, int max, Args& ... args)
    {
        //
        //std::cout << "1d iterate" << std::endl;
        UniIterHost::iterate(fun, max, args...);
    }
};
}

#endif
