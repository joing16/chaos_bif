#include "mex.h"
#include <cuda_runtime.h>
#include <math.h>

// RK4의 각 단계에서 미분값(dx/dt, dy/dt)을 계산하는 디바이스 함수
__device__ void calculate_derivatives(
    double* dx_out, double* dy_out,
    const double* x_in, const double* y_in,
    const double* alpha_x, const double* beta_x,
    const double* alpha_y, const double* beta_y,
    double ext_force, int N, int tid)
{
    double dx = 0.0;
    double dy = 0.0;
    for (int j = 0; j < N; ++j) {
        dx += alpha_x[j * N + tid] * x_in[j];
        dx += beta_x[j * N + tid]  * y_in[j] * (1.0 - fabs(y_in[j]));
        dy += alpha_y[j * N + tid] * y_in[j];
        dy += beta_y[j * N + tid]  * x_in[j] * (1.0 - fabs(x_in[j]));
    }
    dx += ext_force;
    dy += ext_force;
    dx_out[tid] = dx;
    dy_out[tid] = dy;
}


__global__ void model2_kernel(
    double* d_X, double* d_Y,
    double* d_x, double* d_y,
    const double* alpha_x, const double* beta_x,
    const double* alpha_y, const double* beta_y,
    double ext_force_amp, double ext_force_freq, double phase_shift,
    double dt, int numSteps, int N)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= N) return;

    // RK4 중간 계산을 위한 공유 메모리 할당
    extern __shared__ double temp_storage[];
    double* k1_x = &temp_storage[0 * N];
    double* k1_y = &temp_storage[1 * N];
    double* k2_x = &temp_storage[2 * N];
    double* k2_y = &temp_storage[3 * N];
    double* k3_x = &temp_storage[4 * N];
    double* k3_y = &temp_storage[5 * N];
    double* k4_x = &temp_storage[6 * N];
    double* k4_y = &temp_storage[7 * N];
    double* temp_x = &temp_storage[8 * N];
    double* temp_y = &temp_storage[9 * N];

    // 초기 상태 저장
    d_X[tid * numSteps] = d_x[tid];
    d_Y[tid * numSteps] = d_y[tid];
    __syncthreads();

    for (int step = 1; step < numSteps; ++step) {
        double t = (step - 1) * dt;

        // --- 4차 룽게-쿠타(RK4) 계산 시작 ---

        // k1 계산
        double ext_force_t0 = ext_force_amp * sin(ext_force_freq * t + phase_shift);
        calculate_derivatives(k1_x, k1_y, d_x, d_y, alpha_x, beta_x, alpha_y, beta_y, ext_force_t0, N, tid);
        __syncthreads();

        // k2 계산
        double ext_force_t1 = ext_force_amp * sin(ext_force_freq * (t + 0.5 * dt) + phase_shift);
        for (int j = 0; j < N; ++j) { temp_x[j] = d_x[j] + 0.5 * dt * k1_x[j]; temp_y[j] = d_y[j] + 0.5 * dt * k1_y[j]; }
        calculate_derivatives(k2_x, k2_y, temp_x, temp_y, alpha_x, beta_x, alpha_y, beta_y, ext_force_t1, N, tid);
        __syncthreads();

        // k3 계산
        for (int j = 0; j < N; ++j) { temp_x[j] = d_x[j] + 0.5 * dt * k2_x[j]; temp_y[j] = d_y[j] + 0.5 * dt * k2_y[j]; }
        calculate_derivatives(k3_x, k3_y, temp_x, temp_y, alpha_x, beta_x, alpha_y, beta_y, ext_force_t1, N, tid);
        __syncthreads();

        // k4 계산
        double ext_force_t2 = ext_force_amp * sin(ext_force_freq * (t + dt) + phase_shift);
        for (int j = 0; j < N; ++j) { temp_x[j] = d_x[j] + dt * k3_x[j]; temp_y[j] = d_y[j] + dt * k3_y[j]; }
        calculate_derivatives(k4_x, k4_y, temp_x, temp_y, alpha_x, beta_x, alpha_y, beta_y, ext_force_t2, N, tid);
        __syncthreads();

        // 최종 상태 업데이트
        d_x[tid] += (dt / 6.0) * (k1_x[tid] + 2.0 * k2_x[tid] + 2.0 * k3_x[tid] + k4_x[tid]);
        d_y[tid] += (dt / 6.0) * (k1_y[tid] + 2.0 * k2_y[tid] + 2.0 * k3_y[tid] + k4_y[tid]);
        __syncthreads();
        
        // 결과 저장
        d_X[tid * numSteps + step] = d_x[tid];
        d_Y[tid * numSteps + step] = d_y[tid];
    }
}


void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if (nrhs != 12) {
        mexErrMsgIdAndTxt("simulateModel2_full:invalidNumInputs",
            "Expected 12 inputs: dt, numSteps, x0, y0, alpha_x, beta_x, alpha_y, beta_y, ext_force_amp, ext_force_freq, N, phase_shift");
    }

    double dt             = mxGetScalar(prhs[0]);
    int    numSteps       = (int)mxGetScalar(prhs[1]);
    double* h_x0          = mxGetPr(prhs[2]);
    double* h_y0          = mxGetPr(prhs[3]);
    double* h_alpha_x     = mxGetPr(prhs[4]);
    double* h_beta_x      = mxGetPr(prhs[5]);
    double* h_alpha_y     = mxGetPr(prhs[6]);
    double* h_beta_y      = mxGetPr(prhs[7]);
    double ext_force_amp  = mxGetScalar(prhs[8]);
    double ext_force_freq = mxGetScalar(prhs[9]);
    int    N              = (int)mxGetScalar(prhs[10]);
    double phase_shift    = mxGetScalar(prhs[11]); 

    plhs[0] = mxCreateDoubleMatrix(N, numSteps, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N, numSteps, mxREAL);
    double* h_X = mxGetPr(plhs[0]);
    double* h_Y = mxGetPr(plhs[1]);

    double *d_X, *d_Y, *d_x, *d_y, *d_alpha_x, *d_beta_x, *d_alpha_y, *d_beta_y;
    size_t stateSize  = N * numSteps * sizeof(double);
    size_t nodeSize   = N * sizeof(double);
    size_t matrixSize = N * N * sizeof(double);

    cudaMalloc(&d_X, stateSize);
    cudaMalloc(&d_Y, stateSize);
    cudaMalloc(&d_x, nodeSize);
    cudaMalloc(&d_y, nodeSize);
    cudaMalloc(&d_alpha_x, matrixSize);
    cudaMalloc(&d_beta_x,  matrixSize);
    cudaMalloc(&d_alpha_y, matrixSize);
    cudaMalloc(&d_beta_y,  matrixSize);

    cudaMemcpy(d_x, h_x0, nodeSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, h_y0, nodeSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_alpha_x, h_alpha_x, matrixSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_beta_x,  h_beta_x,  matrixSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_alpha_y, h_alpha_y, matrixSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_beta_y,  h_beta_y,  matrixSize, cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;
    
    // --- 수정된 부분: 커널 호출 시 공유 메모리 크기 지정 ---
    // RK4 계산에 필요한 10개의 임시 배열(크기 N)을 공유 메모리에 할당
    size_t shared_mem_size = 10 * N * sizeof(double);
    model2_kernel<<<numBlocks, blockSize, shared_mem_size>>>(d_X, d_Y, d_x, d_y,
                                                            d_alpha_x, d_beta_x,
                                                            d_alpha_y, d_beta_y,
                                                            ext_force_amp, ext_force_freq, phase_shift,
                                                            dt, numSteps, N);
    // ----------------------------------------------------

    cudaMemcpy(h_X, d_X, stateSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_Y, d_Y, stateSize, cudaMemcpyDeviceToHost);

    cudaFree(d_X); cudaFree(d_Y);
    cudaFree(d_x); cudaFree(d_y);
    cudaFree(d_alpha_x); cudaFree(d_beta_x);
    cudaFree(d_alpha_y); cudaFree(d_beta_y);
}