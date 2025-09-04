// -------------------------------------------------------------
// C++ version: 28.01.2025
// Author: Maximilian Gram,
// Experimental Physics 5,
// University of Wuerzburg
// -------------------------------------------------------------

// ---------------------- use in Matlab: -----------------------
// [M, T] = mex_BLOCH_expm(M, w1x, w1y, dw0, R1, R2, dt)
// M:    4 x 1 double array
// w1x:  N x 1 double array
// w1y:  N x 1 double array
// dw0:  N x 1 double array
// R1:   1 x 1 double
// R2:   1 x 1 double
// dt:   1 x 1 double
// -------------------------------------------------------------

// -------------------------------------------------------------
// -------------------------- MEX, C++ -------------------------
// -------------------------------------------------------------

#include <mex.h>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
using namespace std;
using namespace Eigen;

// dimension of Bloch matrix: 4x4
constexpr int n = 4;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    // define pointers to the input arrays
    double* Minit = mxGetPr(prhs[0]);     // 4x1 double array
    double* w1x   = mxGetPr(prhs[1]);     // Nx1 double array
    double* w1y   = mxGetPr(prhs[2]);     // Nx1 double array
    double* dw0   = mxGetPr(prhs[3]);     // Nx1 double array
    double R1     = mxGetScalar(prhs[4]); // 1x1 double scalar
    double R2     = mxGetScalar(prhs[5]); // 1x1 double scalar
    double dt     = mxGetScalar(prhs[6]); // 1x1 double scalar

    // get number of dt steps
    int N = mxGetNumberOfElements(prhs[1]);

    // prepare Eigen arrays
    Eigen::Matrix<double, n, 1> M;
    Eigen::Matrix<double, n, n> B;
    Eigen::Matrix<double, n, n> T;
    T.setIdentity();

    // copy initial M into Eigen vector
    for (int j = 0; j < n; ++j)
    M(j) = Minit[j];

    // Compute Transition Matrix
    for (int j = 0; j < N; ++j)
    {
        B << -R2,       dw0[j],  -w1y[j],  0.0,
             -dw0[j],  -R2,       w1x[j],  0.0,
              w1y[j],  -w1x[j],  -R1,      R1,
              0.0,      0.0,      0.0,     0.0;
        B = B * dt;
        T = B.exp() * T;
    }

    // Compute the final magnetization
    M = T * M;

    // output final M
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double* output_M = mxGetPr(plhs[0]);
    for (int j = 0; j < n; ++j)
    output_M[j] = M(j);

    // output final transition matrix
    plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
    double* output_T = mxGetPr(plhs[1]);
    std::copy(T.data(), T.data() + 16, output_T);

}