// default size in bytes
#define DEFAULTSIZE 4096

#include <cstring>
#include <iostream>

#include <cuda_runtime_api.h>
#include "devcall.h"

    template <class T>
    class DevVec{

        T *ptr;
        size_t count;
        size_t cap;

        inline void realloc(size_t newCap)
        {
            //
            //std::cout << "reallocing to " << newCap << std::endl;
            T *ptrTemp;// = new T[newCap];
            #ifdef GPU
            getErrorCuda((cudaMallocManaged((void**)&ptrTemp, newCap * sizeof(T))));
            #else
            ptrTemp = new T[newCap];
            #endif

            size_t toCopyCount = cap;

            if(newCap < cap)
                toCopyCount = newCap;
            
            #ifdef GPU
            cudaMemcpy(ptrTemp, ptr, sizeof(T)*toCopyCount, cudaMemcpyDefault);
            #else
            std::memcpy(ptrTemp, ptr, sizeof(T)*toCopyCount);
            #endif

            cap = newCap;
            
            #ifdef GPU
            cudaFree(ptr);
            #else
            delete[] ptr;
            #endif

            ptr = ptrTemp;
        }

        public:
        //
        DevVec(){
            //
            cap = DEFAULTSIZE / sizeof(T);
            #ifdef GPU
            getErrorCuda((cudaMallocManaged((void**)&ptr, cap * sizeof(T))));
            #else
            ptrTemp = new T[cap];
            #endif

            count = 0;
        }

        inline void push_back(T val)
        {
            //
            if(count+1 <= cap)
            {
                ptr[count] = val;
                count++;
                return;
            }
            else{
                realloc(cap*2);
                push_back(val);
            }
        }

        inline void reserve(size_t newCap)
        {
            //std::cout << "reserve" << std::endl;
            if(newCap > cap)
                realloc(newCap);
        }

        inline void resize(size_t newCap)
        {
            //
            if(newCap < cap)
            {
                count = newCap;
                // todo go clear unused elements
            }
            else
            {
                realloc(newCap);
            }
        }

        inline void shrink_to_fit()
        {
            // does f all since it does not have to...
        }

        inline T & operator[](size_t ind)
        {
            //
            return ptr[ind];
        }

        inline const T &operator[](size_t ind) const
        {
            //
            return ptr[ind];
        }

        inline size_t size()
        {
            return count;
        }
        inline size_t capacity()
        {
            return cap;
        }
        inline void clear()
        {
            count = 0;
        }

        // used to fix the size
        inline void set_size(size_t newSize)
        {
            count = newSize;
        }


        inline std::vector<T> toVector() const
        {
            //
            //std::cout << "getting vector" << std::endl;
            std::vector<T> ret;
            for (size_t i = 0; i < count; i++)
            {
                ret.push_back(ptr[i]);
            }
            return ret;
        }
    };