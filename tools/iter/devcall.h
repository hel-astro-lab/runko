#pragma once

#ifdef GPU 
#define DEVCALLABLE __host__ __device__ 
#else
#define DEVCALLABLE
#endif


#ifdef DEBUG_DEV
#define getErrorCuda(command)\
		command;\
		cudaDeviceSynchronize();\
		if (cudaPeekAtLastError() != cudaSuccess){\
			std::cout << #command << " : " << cudaGetErrorString(cudaGetLastError())\
			 << " in file " << __FILE__ << " at line " << __LINE__ << std::endl;\
             exit((1));\
		}
#endif
#ifndef DEBUG_DEV
#define getErrorCuda(command) command;
#endif


