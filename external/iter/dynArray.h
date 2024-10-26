//
// Created by Fredrik Robertsen on 07/11/2020.
//

#ifndef PHELL_MANVEC_H
#define PHELL_MANVEC_H


// default size in bytes
#define DEFAULTSIZE 4096

#include <cstring>
#include <iostream>

#ifdef GPU
#include <cuda_runtime_api.h>
#endif
#include "devcall.h"

    template <class T>
    class ManVec{

        T *ptr;
        size_t count;
        size_t cap;
        float overAllocFactor = 1.2;
        bool allocated;

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

            size_t toCopyCount = count;

            if(newCap < cap)
                toCopyCount = newCap;

            #ifdef GPU
            cudaMemcpy(ptrTemp, ptr, sizeof(T)*toCopyCount, cudaMemcpyDefault);
            #else
            std::memcpy(ptrTemp, ptr, sizeof(T)*toCopyCount);
            #endif

            cap = newCap;

/*
            for (size_t i = toCopyCount; i < newCap; i++)
            {
                ptrTemp[i] = T{};
            }
            
*/
            #ifdef GPU
            cudaFree(ptr);
            #else
            delete[] ptr;
            #endif

            ptr = ptrTemp;
            //std::cout << "reallocing out " << std::endl;

        }

        public:
        //
        ManVec(): allocated(true){
            //std::cout << "ManVec create " << std::endl;

            cap = DEFAULTSIZE / sizeof(T);
            #ifdef GPU
            getErrorCuda((cudaMallocManaged((void**)&ptr, cap * sizeof(T))));
            #else
            ptr = new T[cap];
            #endif

            count = 0;
            //std::cout << "ManVec create out" << std::endl;
        }
        ~ManVec(){
            //std::cout << "ManVec deconstruct " << std::endl;

            if(allocated)
            {
                #ifdef GPU
                cudaFree(ptr);
                #else
                delete[] ptr;
                #endif
            }
            //std::cout << "ManVec deconstruct out" << std::endl;
        }

        // copy constructor
        ManVec (const ManVec &old_obj)
        {
            //
            //std::cout << "ManVec copy ctr "<< old_obj.count << " " << old_obj.cap << std::endl;

            cap = old_obj.cap;
            count = old_obj.count;
            allocated = old_obj.allocated;

            #ifdef GPU
            getErrorCuda((cudaMallocManaged((void**)&ptr, cap * sizeof(T))));
            cudaMemcpy(ptr, old_obj.ptr, sizeof(T)*count, cudaMemcpyDefault);
            #else
            ptr = new T[cap];
            std::memcpy(ptr, old_obj.ptr, sizeof(T)*count);
            #endif
            
            /*
            for (size_t i = count; i < cap; i++)
            {
                ptr[i] = T{};
            }
            */
            
        }

        // move constructor
        ManVec ( ManVec && other) :allocated(false)
        {
            //std::cout << "ManVec move ctr "<< other.count << " " << other.cap << std::endl;
            count = 0;
            cap = 0;
            std::swap(count, other.count);
            std::swap(cap, other.cap);
            std::swap(ptr, other.ptr);
            std::swap(allocated, other.allocated);
        }


        // v0
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
                realloc((cap+1)*overAllocFactor);
                push_back(val);
            }
        }

        //v1 and v2 with emplace_back as the lowest operation

        //v1
        //template <typename... Args>
        //inline void emplace_back(Args&&... args) {
        //  //if (size == capacity) grow();
        //  if(count+1 > cap) {
        //    realloc((cap+1)*overAllocFactor);
        //  }
        //  
        //  return *new (start + size++) T(std::forward<Args>(args)...);
        //}
        //inline void push_back(const T & t) { return emplace_back(t); }
        //inline void push_back(T && t) { return emplace_back(std::move(t)); }

        //v2
        //template<class... Args>
        //inline void emplace_back(Args&&... args) {

        //  if(count+1 > cap){
        //    realloc((cap+1)*overAllocFactor);
        //  }

        //  ptr[count] = T(std::forward<Args>(args)...);
        //  count++;
        //  return;
        //}

        //inline void push_back(T&& val) {
        //  emplace_back(std::move(val));
        //}

        //inline void push_back(const T& val) {
        //  emplace_back(val);
        //}



        inline void reserve(size_t newCap)
        {
            //std::cout << "reserve" << std::endl;
            if(newCap > cap) realloc(newCap);
        }

        inline void resize(size_t newCap)
        {
            if(newCap <= cap)
            {
                count = newCap;
                // todo go clear unused elements
            }
            else
            {
                //std::cout << "reallocing " << newCap << " old cap " << cap << std::endl;
                //realloc((newCap+1)*overAllocFactor); // resize with overAllocFactor // NOTE: this is a bug; not a feature
                realloc(newCap);  // resize strictly what is requested
            }
            count = newCap;
        }

        inline void shrink_to_fit()
        {
            // does f all since it does not have to...

            // TODO: check if it works
            realloc(count);
        }

        DEVCALLABLE
        inline T & operator[](const size_t &ind)
        {
            //
            //if(ind >= cap)
            //    std::cout << "error " << ind << " " << cap << std::endl;
                
            return ptr[ind];
        }

        DEVCALLABLE
        inline const T &operator[](const size_t &ind) const
        {
            //
            
            //if(ind >= cap)
            //    std::cout << "error " << ind << " " << cap << std::endl;
                
            return ptr[ind];
        }

        DEVCALLABLE
        inline size_t size() const
        {
            return count;
        }

        DEVCALLABLE
        inline size_t capacity() const
        {
            return cap;
        }

        inline void clear()
        {
            count = 0;
        }

        DEVCALLABLE
        inline T*& data()
        {
            return ptr;
        }

        // used to fix the size
        inline void set_size(size_t newSize)
        {
            count = newSize;
        }

        DEVCALLABLE
        inline bool empty()
        {
            return count == 0;
        }

        DEVCALLABLE
        inline T* begin()
        {
            return ptr;
        }

        DEVCALLABLE
        inline T* end()
        {
            return ptr+count;
        }

        DEVCALLABLE
        inline T const* cbegin() const
        {
            return ptr;
        }

        DEVCALLABLE
        inline T const* cend() const
        {
            return ptr+count;
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

#endif
