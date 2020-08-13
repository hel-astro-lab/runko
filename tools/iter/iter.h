#ifndef ITER_H
#define ITER_H
#pragma once

#ifdef GPU 
#include <cuda_runtime_api.h>
#else
#endif

#include <tuple>
#include <unordered_map>
#include <map>
#include <optional>
//#define DEBUG_DEV
#define GPU

#include "devcall.h"

#include <execinfo.h>
#include <stdio.h>

namespace {

#ifdef GPU
namespace DevKernels
{
    template<class F, class... Args>
    __global__ void iterateKern(F fun, int max, Args* ... args)
    {
        //
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        if(idx >= max) return;

        fun(idx, *args...);
    }

    template<class F, class... Args>
    __global__ void iterate3dKern(F fun, int xMax, int yMax, int zMax, Args* ... args)
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        int idy = blockIdx.y * blockDim.y + threadIdx.y;
        int idz = blockIdx.z * blockDim.z + threadIdx.z;
        if(idx >= xMax) return;
        if(idy >= yMax) return;
        if(idz >= zMax) return;

        fun(idx, idy, idz, *args...);
    }

    template<class T>
    __global__ void copyItoDMulti(int max, int **indexes, T **from, T **to)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int ptrInd = threadIdx.y;
        if(i >= max) return;
        
        int indirectInd = indexes[ptrInd][i];

        to[ptrInd][i] = from[ptrInd][indirectInd];
    }

    template<class T>
    __global__ void copyDtoIMulti(int max, int **indexes, T **from, T **to)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int ptrInd = threadIdx.y;
        if(i >= max) return;
        
        int indirectInd = indexes[ptrInd][i];

        from[ptrInd][indirectInd] = to[ptrInd][i];
    }
}
#endif

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

    #ifdef GPU
    class UniIterCU{
    public:

        static std::unordered_map<std::string, cudaStream_t> streams;
        static std::map<void*, void*> registeredAddesses;



        template<class T>
        static T *getDevPtr(T &t)
        {
          //printf("getting %p \n", &t);
          //std::cout << "getting dev ptr" << std::endl;
          T *ptr;
          auto err = cudaHostGetDevicePointer(&ptr, &t, 0);
          if(err != cudaSuccess) {
            std::cout << cudaGetErrorString(cudaGetLastError()) << " in file "
                      << __FILE__ << " at line " << __LINE__ << std::endl;
            //return reg(t);

            }
          return ptr;
        }

        static auto getGrid3(std::tuple<int, int, int> max)
        {
            dim3 block = { 8, 8, 1 };
            dim3 grid  = { 1 + (std::get<0>(max) / block.x), 1 + (std::get<1>(max) / block.y), 1 + (std::get<2>(max) / block.z) };

            return std::make_tuple(block, grid);
        }

        static auto getGrid3(std::tuple<int, int, int> max, std::tuple<int, int, int> block_)
        {
            dim3 block = { std::get<0>(block_), std::get<1>(block_), std::get<2>(block_) };
            dim3 grid  = { 1 + (std::get<0>(max) / block.x), 1 + (std::get<1>(max) / block.y), 1 + (std::get<2>(max) / block.z) };

            return std::make_tuple(block, grid);
        }

        static auto getGrid(std::tuple<int> max)
        {
            dim3 block = { 128};
            dim3 grid  = { 1 + (std::get<0>(max) / block.x)};

            return std::make_tuple(block, grid);
        }

        template<class F, class... Args>
        static void iterate(F fun, int max, Args& ... args)
        {
            auto [block, grid] = getGrid({max});

            getErrorCuda(((DevKernels::iterateKern<<<grid, block>>>(fun, max, getDevPtr(args)...))));

         }

        template<class F, class... Args>
        static void iterate3D(F fun, int xMax, int yMax, int zMax, Args& ... args)
        {
            auto [block, grid] = getGrid3({xMax, yMax, zMax});

            getErrorCuda(((DevKernels::iterate3dKern<<<grid, block>>>(fun, xMax, yMax, zMax, getDevPtr(args)...))));

        }
        template<class F, class... Args>
        static void iterate3D(F fun, ExecPlan plan, int xMax, int yMax, int zMax, Args& ... args)
        {
            std::tuple<dim3, dim3> gridArgs;

            if(plan.block)
            {
                gridArgs = getGrid3({xMax, yMax, zMax}, plan.block.value());
            }
            else
            {
                gridArgs = getGrid3({xMax, yMax, zMax});
            }

            // unpack and handle exec plan here
            if (plan.stream == "default")
            {
                // default stream
                //std::cout << "Using default stream using it" << std::endl;
                getErrorCuda(((DevKernels::iterate3dKern<<<std::get<1>(gridArgs), std::get<0>(gridArgs)>>>(fun, xMax, yMax, zMax, getDevPtr(args)...))));
            }
            else
            {
                // non default stream
                // check if it exsists
                // if not create it
                auto search = streams.find(plan.stream);
                if (search != streams.end()) {
                    //
                    std::cout << "stream found using it" << std::endl;
                } else {
                    streams[plan.stream] = cudaStream_t();
                    cudaStreamCreate(&streams[plan.stream]);
                    std::cout << "Not found stream crated" << std::endl;
                }
                getErrorCuda(((DevKernels::iterate3dKern<<<std::get<1>(gridArgs), std::get<0>(gridArgs), 0, streams[plan.stream]>>>(fun, xMax, yMax, zMax, getDevPtr(args)...))));
    
            }
        }
    };
    #endif


    class UniIterHost{
        template<int Dims, class F, class... Args>
        static void iterateND(F fun, int xMax, int yMax, int zMax, Args& ... args)
        {
            //
        }


    public:

        template<class F, class... Args>
        static void iterate3D(F fun, int xMax, int yMax, int zMax, Args& ... args)
        {
            #pragma omp parallel for
            for (int x = 0; x < xMax; x++)
            {
                for (int y = 0; y < yMax; y++)
                {
                    #pragma omp simd
                    for (int z = 0; z < zMax; z++)
                    {
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
    static void iterate3D(F fun, int xMax, int yMax, int zMax, Args& ... args)
    {
        //
        //std::cout << "3d iterate" << std::endl;
    #ifdef GPU
        UniIterCU::iterate3D(fun, xMax, yMax, zMax, args...);
    #else
        UniIterHost::iterate3D(fun, xMax, yMax, zMax, args...);
    #endif
    }

    template<class F, class... Args>
    static void iterate3D(F fun, ExecPlan plan, int xMax, int yMax, int zMax, Args& ... args)
    {
        //
        //std::cout << "3d iterate w plan" << std::endl;
    #ifdef GPU
        UniIterCU::iterate3D(fun, plan, xMax, yMax, zMax, args...);
    #else
        UniIterHost::iterate3D(fun, xMax, yMax, zMax, args...);
    #endif
    }
    
    template<class F, class... Args>
    static void iterate(F fun, int max, Args& ... args)
    {
        //
        //std::cout << "1d iterate" << std::endl;
    #ifdef GPU
        UniIterCU::iterate(fun, max, args...);
    #else
        UniIterHost::iterate(fun, max, args...);
    #endif
    }

    static void sync(){
        #ifdef GPU
        auto err = cudaDeviceSynchronize();
        // todo: check the error code
        //std::cout << err << std::endl;
        #else
        #endif
        //
    }
};

#ifdef GPU
std::unordered_map<std::string, cudaStream_t> UniIter::UniIterCU::streams = std::unordered_map<std::string, cudaStream_t>();
std::map<void*, void*> UniIter::UniIterCU::registeredAddesses = std::map<void*, void*>();
#endif
}

#endif