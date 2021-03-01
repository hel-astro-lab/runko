#pragma once

#ifdef GPU 
#define DEVCALLABLE __host__ __device__ 
#else
#define DEVCALLABLE
#endif

#ifdef GPU 
#define DEVFUNIFGPU __device__ 
#else
#define DEVFUNIFGPU
#endif


//#define DEBUG_DEV
#ifdef DEBUG_DEV
#define getErrorCuda(command)\
		command;\
		cudaDeviceSynchronize();\
		if (cudaPeekAtLastError() != cudaSuccess){\
			std::cout << #command << " : " << cudaGetErrorString(cudaGetLastError())\
			 << " in file " << __FILE__ << " at line " << __LINE__ << " in function " << __PRETTY_FUNCTION__ << std::endl;\
             exit((1));\
		}
#endif
#ifndef DEBUG_DEV
#define getErrorCuda(command) command;
#endif


