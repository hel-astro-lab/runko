#pragma once

#ifdef GPU 
#include <cuda_runtime_api.h>
#include "devcall.h"
#endif

#include <map>

#ifdef GPU
#define DEV_REGISTER UniAllocator::registerClass(*this);
#else
#define DEV_REGISTER
#endif

namespace 
{

    class UniAllocator{   
        static std::map<void*, void*> registeredAddesses;

    public:

        template<class T>
        // register values passed in
        static void registerClass(T &t)
        {
            //std::cout << "found value registering " << std::endl;
            void* addrHost = &t;
            auto search = registeredAddesses.find(addrHost);
            if (search != registeredAddesses.end()) {

            } else {
                //https://migocpp.wordpress.com/2018/06/08/cuda-memory-access-global-zero-copy-unified/
                cudaHostRegister(&t, sizeof(T), cudaHostRegisterMapped);
                auto err = cudaGetLastError();
                if (err == cudaErrorHostMemoryAlreadyRegistered)
                {
                    // figure out what part is already regisered and free it then reregister the whole thing
                    auto start = registeredAddesses.lower_bound(addrHost);
                    auto end = registeredAddesses.lower_bound(addrHost+sizeof(T));

                    int numUnregistered = 0;
                    for(auto iter = start; iter != end; iter++)
                    {
                        if(true) // if it is in within the range of hostaddr and gostaddr+sizeof(T)
                        {
                            getErrorCuda((cudaHostUnregister(iter->first)));
                            numUnregistered ++;
                        }
                    }
                    if(numUnregistered != 0)
                    {
                        getErrorCuda((cudaHostRegister(&t, sizeof(T), cudaHostRegisterMapped)));
                    }
                }
                else if(err != cudaSuccess)
                {
                    std::cout << "other error while registering memory" << std::endl;
                }
            }
        }
        
        template<class T>
        static void unRegisterClass(T &t)
        {
            // todo fix
        }

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
            auto err = cudaFree(ptr);
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

    std::map<void*, void*> UniAllocator::registeredAddesses = std::map<void*, void*>();

}