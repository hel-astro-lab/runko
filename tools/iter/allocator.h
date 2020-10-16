#pragma once

#include <cuda_runtime_api.h>
#ifdef GPU 
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
            #ifdef GPU
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
            #endif
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


    //https://howardhinnant.github.io/allocator_boilerplate.html

template <class T>
struct ManagedAlloc
{
// based on https://howardhinnant.github.io/allocator_boilerplate.html
	using value_type    = T;

//     using pointer       = value_type*;
//     using const_pointer = typename std::pointer_traits<pointer>::template
//                                                     rebind<value_type const>;
//     using void_pointer       = typename std::pointer_traits<pointer>::template
//                                                           rebind<void>;
//     using const_void_pointer = typename std::pointer_traits<pointer>::template
//                                                           rebind<const void>;

//     using difference_type = typename std::pointer_traits<pointer>::difference_type;
//     using size_type       = std::make_unsigned_t<difference_type>;

//     template <class U> struct rebind {typedef ManagedAlloc<U> other;};
	
	ManagedAlloc() noexcept {}  // not required, unless used

	template <class U> constexpr ManagedAlloc(const ManagedAlloc <U>&) noexcept {}

	[[nodiscard]] T* allocate(std::size_t n)
	{
		if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
			throw std::bad_alloc();

		T* ptr;
		auto err = cudaMallocManaged((void**)&ptr, n * sizeof(T));
		if (err == cudaSuccess)
			return ptr;

		throw std::bad_alloc();
	}
	void deallocate(T* p, std::size_t) noexcept {
		cudaFree(p);
	}

//     value_type*
//     allocate(std::size_t n, const_void_pointer)
//     {
//         return allocate(n);
//     }

//     template <class U, class ...Args>
//     void
//     construct(U* p, Args&& ...args)
//     {
//         ::new(p) U(std::forward<Args>(args)...);
//     }

//     template <class U>
//     void
//     destroy(U* p) noexcept
//     {
//         p->~U();
//     }

//     std::size_t
//     max_size() const noexcept
//     {
//         return std::numeric_limits<size_type>::max();
//     }

//     allocator
//     select_on_container_copy_construction() const
//     {
//         return *this;
//     }

//     using propagate_on_container_copy_assignment = std::false_type;
//     using propagate_on_container_move_assignment = std::false_type;
//     using propagate_on_container_swap            = std::false_type;
//     using is_always_equal                        = std::is_empty<allocator>;
};

template <class T, class U>
bool
operator==(ManagedAlloc<T> const&, ManagedAlloc<U> const&) noexcept
{
    return true;
}

template <class T, class U>
bool
operator!=(ManagedAlloc<T> const& x, ManagedAlloc<U> const& y) noexcept
{
    return !(x == y);
}

}