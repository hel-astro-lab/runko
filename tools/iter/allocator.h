#pragma once

#include <cuda_runtime_api.h>
#ifdef GPU 
#include "devcall.h"
#endif

#include <map>
#include <cstddef>

namespace 
{
    class UniAllocator{   

    public:

        template<class T>
        static T *allocate(int count){
            T *ptr;
            #ifdef GPU
            getErrorCuda((cudaMallocManaged((void**)&ptr, count * sizeof(T))));
            // todo: check the error code
            #else
            ptr = new T[count];
            #endif
            //
            return ptr;
        }
        template<class T>
        static void deallocate(T *ptr)
        {
            #ifdef GPU
            getErrorCuda(cudaFree(ptr));
            // todo: check the error code
            #else
            delete[] ptr;
            #endif
        }

    // optional allocator that will always allocate memory on the device when running on the GPU
    // useful for MPI buffers to enable MPI to directly move data to the NIC from the GPU
        template<class T>
        static T *allocateScratch(int count){
            T *ptr;
            #ifdef GPU
            getErrorCuda((cudaMalloc((void**)&ptr, count * sizeof(T))));
            // todo: check the error code
            #else
            // todo: change to omp5 allocation when available
            ptr = new T[count];
            #endif
            //
            return ptr;
        }
        template<class T>
        static void deallocateScratch(T *ptr)
        {
            #ifdef GPU
            auto err = cudaFree(ptr);
            // todo: check the error code
            #else
            delete[] ptr;
            #endif
        }
    };


}


class ManagedParent {
    public:
    void *operator new(size_t len) {
        void *ptr = UniAllocator::allocate<unsigned char>(len);
        return ptr;
    }
    void operator delete(void *ptr) {
        UniAllocator::deallocate<unsigned char>((unsigned char*)ptr);
    }
};


class ManagedParent2 {
    public:
    void *operator new(size_t len) {
        std::cout << "allocating " << std::endl;
        void *ptr = UniAllocator::allocate<unsigned char>(len);
        return ptr;
    }
    void operator delete(void *ptr) {
        UniAllocator::deallocate<unsigned char>((unsigned char*)ptr);
    }
};