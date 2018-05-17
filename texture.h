#ifndef TEXTURE_H
#define TEXTURE_H

#include<cuda.h>

// These textures are (optionally) used to cache the 'x' vector in y += A*x

texture<float,1>  tex_x;


inline void bind_x(const float* x)
{     
    size_t offset = size_t(-1);
    gpuErrchk(cudaBindTexture(&offset, tex_x, x));
}

inline void unbind_x(const float* x)
{
    gpuErrchk(cudaUnbindTexture(tex_x));
}


__inline__ __device__ float fetch_x(const int i, const float* x)
{
	return tex1Dfetch(tex_x, i);
    //int2 v = tex1Dfetch(tex_x, i);
    //return __hiloint2double(v.y, v.x);

}
#endif
