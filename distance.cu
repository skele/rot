#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <cstdlib>
#include <omp.h>

//#define N (2 << 20)
//#define N (1<<10)
#define L 6
#define NT 1024

// CUDA Kernel
//   (x, y, z, d) are thrust vectors of size N
//   (p) is a device pointer array size 3
__global__ void distance(float *x, float *y, float *z, float *d, int i){
    int tx = threadIdx.x;
    int j  = threadIdx.x + blockDim.x * blockIdx.x;
    
    __shared__ float s_x[NT];
    __shared__ float s_y[NT];
    __shared__ float s_z[NT];

    s_x[tx] = x[j];
    s_y[tx] = y[j];
    s_z[tx] = z[j];
    __syncthreads();
    
    //so that the particle itself will not be selected as a nearest neighbor
    d[i] = 10000000.0;

    if (j != i)
    {
        d[j] = sqrt( (s_x[tx] - x[i]) * (s_x[tx] - x[i]) +
                     (s_y[tx] - y[i]) * (s_y[tx] - y[i]) +
                     (s_z[tx] - z[i]) * (s_z[tx] - z[i]));
    }
}

void obtain_densities(float *mx, float *my, float *mz, float *mradius, float *mdensity, int istart, int icount, int ikid, int N, float mass)
{

    int i, j;
    double ini, end;
    float radius,density;

    cudaSetDevice(ikid);

    // Host vectors (x, y, z, r)
    thrust::host_vector<float> h_x(N), h_y(N), h_z(N), h_r(N);
 
    // Fillong host vectors (x, y, z) with random numbers
    /*    thrust::generate(h_x.begin(), h_x.end(), rand);
    thrust::generate(h_y.begin(), h_y.end(), rand);
    thrust::generate(h_z.begin(), h_z.end(), rand);*/
    std::fill(h_r.begin(), h_r.end(), 0);

    //copy into host vectors;
    for (i = 0; i < N; i++)
      {
	h_x[i] = mx[i];
	h_y[i] = my[i];
	h_z[i] = mz[i];
	
      }

    // Device vectors (x, y, z, r ) using the memory of the Host vectors
    thrust::device_vector<float> d_x = h_x, d_y = h_y, d_z = h_z, d_r = h_r; 
   
    // Raw pointer of the vectors to use a normal CUDA kernel (x, y, z, r)
    float *x = thrust::raw_pointer_cast( &d_x[0] );
    float *y = thrust::raw_pointer_cast( &d_y[0] );
    float *z = thrust::raw_pointer_cast( &d_z[0] );
    float *r = thrust::raw_pointer_cast( &d_r[0] );

    // CUDA kernel configuration
    int nthread = NT;
    int nblocks = ceil(N/nthread);

    // GPU procedure

    ini = omp_get_wtime();
    
    for(i = istart; i < istart+icount; i++)
    {
        distance <<< nblocks, nthread >>> (x,y,z,r,i);
        cudaThreadSynchronize();

        // Sorting d_r
        thrust::sort(d_r.begin(), d_r.end());
        
        // Copy data after sorting
        //  from d_r to h_r
        thrust::copy(d_r.begin(), d_r.begin() + L, h_r.begin());
	/*        for (j = 0; j < L; j++) {
            std::cout << h_r[j] << " ";
        }
        std::cout << std::endl;
	*/
	//need to read the L-1 element of the distances
	radius = d_r[L-1];
	density = (L-1)*mass*3.0/4.0/(radius*radius*radius);
	mradius[i] = radius;
	mdensity[i] = density;
    }

    end = omp_get_wtime();
    std::cout << "GPU time  per particle: " << (end-ini)/icount << std::endl;
    std::cout << "GPU total time: " << (end-ini) << std::endl;

    // CPU procedure
    /*
    std::vector<float> e(N);

    ini = omp_get_wtime();
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            e[j] = sqrt( (h_x[j] - h_x[i]) * (h_x[j] - h_x[i]) +
                         (h_y[j] - h_y[i]) * (h_y[j] - h_y[i]) +
                         (h_z[j] - h_z[i]) * (h_z[j] - h_z[i]));

        }
        std::sort(e.begin(), e.end());
        for (j = 0; j < L; j++) {
            std::cout << e[j] << " ";
        }
        std::cout << std::endl;
    }

    end = omp_get_wtime();
    std::cout << "CPU time  per particle: " << (end-ini)/N << std::endl;
    std::cout << "CPU total time: " << (end-ini) << std::endl;
    */

}
