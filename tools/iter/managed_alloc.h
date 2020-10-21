#pragma once

#include <cuda_runtime_api.h>

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
		T* ptr;
		auto err = cudaMallocManaged((void**)&ptr, n * sizeof(T));
		if (err == cudaSuccess)
			return ptr;

		throw std::bad_alloc();
	}
	void deallocate(T* p, std::size_t) noexcept 
    {
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
