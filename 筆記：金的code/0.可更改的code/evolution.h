#ifndef EVOLUTION_FILE
#define EVOLUTION_FILE
//注意:__deviece__ 在GPU調用，在GPU執行，
//注意:__host__ 在CPU調用，在CPU執行。
//注意:__global__在CPU調用，在GPU執行。
#include "interpolationHillISLBM.h"
#include "MRT_Process.h"
#include "MRT_Matrix.h"
//===========================================================================================================================================
//前三個程式碼: Fx_in, fx_old 
__device__ double ModifydRho_F378( double F3_in,  double F7_in,  double F8_in,  double f4_old,  double f9_old,  double f10_old )
{
    double drho = F3_in + F7_in + F8_in - f4_old - f9_old - f10_old;
    //正y方向 - 負y方向  (平面)
    //正y方向壁面的質量流率 
    //正y方向壁面:淨質量流率 流出  - 流入  外法線方向
    return drho;
}
    
__device__ double ModifydRho_F4910( double F4_in,  double F9_in,  double F10_in,  double f3_old,  double f7_old,  double f8_old )
{
    double drho = F4_in + F9_in + F10_in - f3_old - f7_old - f8_old;
    return drho;
}

__device__ double ModifydRho_F15( double F15_in, double f18_old )
{
    double drho = F15_in - f18_old;
    return drho;
}

__device__ double ModifydRho_F16(
    double F16_in, double f17_old)
{
    double drho = F16_in - f17_old;
    //為什麼不是 f17_old-F16_in;
    return drho;
}

__device__ double dRhoglobal( 
    double F1_in, double F2_in, double F3_in, double F4_in, double F5_in, double F6_in, double F7_in, double F8_in, double F9_in, 
    double F10_in, double F11_in, double F12_in, double F13_in, double F14_in, double F15_in, double F16_in, double F17_in, double F18_in,
    double f1_old, double f2_old, double f3_old, double f4_old, double f5_old, double f6_old, double f7_old, double f8_old, double f9_old,
    double f10_old, double f11_old, double f12_old, double f13_old, double f14_old, double f15_old, double f16_old, double f17_old, double f18_old)
{
    double globaldrho = F1_in + F2_in + F3_in + F4_in + F5_in + F6_in 
                      + F7_in + F8_in + F9_in + F10_in + F11_in + F12_in 
                      + F13_in + F14_in + F15_in + F16_in + F17_in + F18_in 
                      - f1_old - f2_old - f3_old - f4_old - f5_old - f6_old 
                      - f7_old - f8_old - f9_old - f10_old - f11_old - f12_old 
                      - f13_old - f14_old - f15_old - f16_old - f17_old - f18_old;
    return globaldrho; //控制體積單元的密度修正://少了F_0作為監測 
}
//===========================================================================================================================================
__global__ void stream_collide_Buffer(
    double *f0_old, double *f1_old, double *f2_old, double *f3_old, double *f4_old, double *f5_old, double *f6_old,
    double *f7_old, double *f8_old, double *f9_old, double *f10_old, double *f11_old, double *f12_old, double *f13_old,
    double *f14_old, double *f15_old, double *f16_old, double *f17_old, double *f18_old, 
    double *f0_new, double *f1_new, double *f2_new, double *f3_new, double *f4_new, double *f5_new, double *f6_new,
    double *f7_new, double *f8_new, double *f9_new, double *f10_new, double *f11_new, double *f12_new, double *f13_new,
    double *f14_new, double *f15_new, double *f16_new, double *f17_new, double *f18_new, 

    double *X0_0,  double *X0_1, double *X0_2,  double *X0_3,  double *X0_4,  double *X0_5,  double *X0_6,
    double *X2_0,  double *X2_1,  double *X2_2,  double *X2_3,  double *X2_4,  double *X2_5,  double *X2_6,
    double *Y0_0,  double *Y0_1, double *Y0_2,  double *Y0_3,  double *Y0_4,  double *Y0_5,  double *Y0_6,
    double *Y2_0,  double *Y2_1,  double *Y2_2,  double *Y2_3,  double *Y2_4,  double *Y2_5,  double *Y2_6,

    double *XiF3_0,  double *XiF3_1,  double *XiF3_2,  double *XiF3_3,  double *XiF3_4,  double *XiF3_5,  double *XiF3_6,
    double *XiF4_0,  double *XiF4_1,  double *XiF4_2,  double *XiF4_3,  double *XiF4_4,  double *XiF4_5,  double *XiF4_6,
    double *XiF5_0,  double *XiF5_1,  double *XiF5_2,  double *XiF5_3,  double *XiF5_4,  double *XiF5_5,  double *XiF5_6,
    double *XiF6_0,  double *XiF6_1,  double *XiF6_2,  double *XiF6_3,  double *XiF6_4,  double *XiF6_5,  double *XiF6_6,
    double *XiF15_0, double *XiF15_1, double *XiF15_2, double *XiF15_3, double *XiF15_4, double *XiF15_5, double *XiF15_6,
    double *XiF16_0, double *XiF16_1, double *XiF16_2, double *XiF16_3, double *XiF16_4, double *XiF16_5, double *XiF16_6,
    double *XiF17_0, double *XiF17_1, double *XiF17_2, double *XiF17_3, double *XiF17_4, double *XiF17_5, double *XiF17_6,
    double *XiF18_0, double *XiF18_1, double *XiF18_2, double *XiF18_3, double *XiF18_4, double *XiF18_5, double *XiF18_6,
    
    double *XBFLf37_0,   double *XBFLf37_1,   double *XBFLf37_2,   double *XBFLf37_3,   double *XBFLf37_4,   double *XBFLf37_5,   double *XBFLf37_6,
    double *XBFLf38_0,   double *XBFLf38_1,   double *XBFLf38_2,   double *XBFLf38_3,   double *XBFLf38_4,   double *XBFLf38_5,   double *XBFLf38_6,
    double *YBFLf3_0,    double *YBFLf3_1,    double *YBFLf3_2,    double *YBFLf3_3,    double *YBFLf3_4,    double *YBFLf3_5,    double *YBFLf3_6,
    double *XiBFLf3_0,   double *XiBFLf3_1,   double *XiBFLf3_2,   double *XiBFLf3_3,   double *XiBFLf3_4,   double *XiBFLf3_5,   double *XiBFLf3_6,
    double *XBFLf49_0,   double *XBFLf49_1,   double *XBFLf49_2,   double *XBFLf49_3,   double *XBFLf49_4,   double *XBFLf49_5,   double *XBFLf49_6,
    double *XBFLf410_0,  double *XBFLf410_1,  double *XBFLf410_2,  double *XBFLf410_3,  double *XBFLf410_4,  double *XBFLf410_5,  double *XBFLf410_6,
    double *YBFLf4_0,    double *YBFLf4_1,    double *YBFLf4_2,    double *YBFLf4_3,    double *YBFLf4_4,    double *YBFLf4_5,    double *YBFLf4_6,
    double *XiBFLf4_0,   double *XiBFLf4_1,   double *XiBFLf4_2,   double *XiBFLf4_3,   double *XiBFLf4_4,   double *XiBFLf4_5,   double *XiBFLf4_6,
    double *YBFLf15_0,   double *YBFLf15_1,   double *YBFLf15_2,   double *YBFLf15_3,   double *YBFLf15_4,   double *YBFLf15_5,   double *YBFLf15_6,
    double *XiBFLf15_0,  double *XiBFLf15_1,  double *XiBFLf15_2,  double *XiBFLf15_3,  double *XiBFLf15_4,  double *XiBFLf15_5,  double *XiBFLf15_6,
    double *YBFLf16_0,   double *YBFLf16_1,   double *YBFLf16_2,   double *YBFLf16_3,   double *YBFLf16_4,   double *YBFLf16_5,   double *YBFLf16_6,
    double *XiBFLf16_0,  double *XiBFLf16_1,  double *XiBFLf16_2,  double *XiBFLf16_3,  double *XiBFLf16_4,  double *XiBFLf16_5,  double *XiBFLf16_6,
    
    int *BFLReqF3_d,     int *BFLReqF4_d,     int *BFLReqF15_d,    int *BFLReqF16_d,
    double *u,           double *v,           double *w,           double *rho_d,       double *Force,       int start,     double *rho_modify, 
    double *Q3_d, double*Q4_d,  double *Q15_d, double*Q16_d)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;//定義x方向索引
    const int j = blockIdx.y*blockDim.y + threadIdx.y;//定義y方向索引
    const int k = blockIdx.z*blockDim.z + threadIdx.z;//定義z方向索引

    if( i <= 2 || i >= NX6-3 || k <= 2 || k >= NZ6-3 ) return;

    const int index = j*NX6*NZ6 + k*NX6 + i;
          int idx;
          int idx_xi = j*NZ6 + k;
    const int nface = NX6*NZ6;
    const int nline = NX6;
    double F0_in,  F1_in,  F2_in,  F3_in,  F4_in,  F5_in,  F6_in,  F7_in,  F8_in,  F9_in;
	double F10_in, F11_in, F12_in, F13_in, F14_in, F15_in, F16_in, F17_in, F18_in;

    //MRT Variable//
    double m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18;
	double meq0,meq1,meq2,meq3,meq4,meq5,meq6,meq7,meq8,meq9,meq10,meq11,meq12,meq13,meq14,meq15,meq16,meq17,meq18;
    /////////////

    int cell_z = k-3;
    if( k <= 6 ) cell_z = 3;
    if( k >= NZ6-7 ) cell_z = NZ6-10;

    // Matrix //
    Matrix;
    Inverse_Matrix;
    Relaxation;
    ////////////

    //F0_in  = f0_old[index];
	//F1_in  = f1_old[index-1];
	//F2_in  = f2_old[index+1];
	//F3_in  = f3_old[index-NX6*NZ6];
	//F4_in  = f4_old[index+NX6*NZ6];
	//F5_in  = f5_old[index-NX6];
	//F6_in  = f6_old[index+NX6];
	//F7_in  = f7_old[index-NX6*NZ6-1];
	//F8_in  = f8_old[index-NX6*NZ6+1];
	//F9_in  = f9_old[index+NX6*NZ6-1];
	//F10_in = f10_old[index+NX6*NZ6+1];
	//F11_in = f11_old[index-NX6-1];
	//F12_in = f12_old[index-NX6+1];
	//F13_in = f13_old[index+NX6-1];
	//F14_in = f14_old[index+NX6+1];
	//F15_in = f15_old[index-NX6*NZ6-NX6];
	//F16_in = f16_old[index+NX6*NZ6-NX6];
	//F17_in = f17_old[index-NX6*NZ6+NX6];
	//F18_in = f18_old[index+NX6*NZ6+NX6];
    F0_Intrpl7( f0_old,  i, j, k);
    F1_Intrpl7( f1_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6 );
    F2_Intrpl7( f2_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6 );
    F3_Intrpl7( f3_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    F4_Intrpl7( f4_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    F5_Intrpl7( f5_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    F6_Intrpl7( f6_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    //7,8,9,10需要往三個方向內插
    X_Y_XI_Intrpl7( f7_old, F7_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    X_Y_XI_Intrpl7( f8_old, F8_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    X_Y_XI_Intrpl7( f9_old, F9_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    X_Y_XI_Intrpl7(f10_old, F10_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    //x_z平面方向之分量:先對X軸做內插，再對Z軸做內插
    X_XI_Intrpl7(  f11_old, F11_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    X_XI_Intrpl7(  f12_old, F12_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    X_XI_Intrpl7(  f13_old, F13_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    X_XI_Intrpl7(  f14_old, F14_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    //y_z平面方向之分量:先對Y軸做內插，再對Z軸做內插
    Y_XI_Intrpl7(  f15_old, F15_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF15_0, XiF15_1, XiF15_2, XiF15_3, XiF15_4, XiF15_5, XiF15_6 );
    Y_XI_Intrpl7(  f16_old, F16_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF16_0, XiF16_1, XiF16_2, XiF16_3, XiF16_4, XiF16_5, XiF16_6 );
    Y_XI_Intrpl7(  f17_old, F17_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF17_0, XiF17_1, XiF17_2, XiF17_3, XiF17_4, XiF17_5, XiF17_6 );
    Y_XI_Intrpl7(  f18_old, F18_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF18_0, XiF18_1, XiF18_2, XiF18_3, XiF18_4, XiF18_5, XiF18_6 );
    

    

    if( k == 3 ){
        F5_in  = f6_old[index];
        F11_in = f14_old[index];
        F12_in = f13_old[index];
        F15_in = f18_old[index];
        F16_in = f17_old[index];
    }
    if( k == NZ6-4 ){
        F6_in  = f5_old[index];
        F13_in = f12_old[index];
        F14_in = f11_old[index];
        F17_in = f16_old[index];
        F18_in = f15_old[index];
    }
    /* if( k == 3 ){
        F5_in  = f6_old[index];
        F11_in = f14_old[index];
        F12_in = f13_old[index];
    } */

    __syncthreads();
    //BFL Linear
    if( k == 3 || k == 4 ) {
        idx_xi = (k-3)*NYD6+j;
        if ( BFLReqF3_d[(k-3)*NYD6+j] == 1 ){
            if(Q3_d[(k-3)*NYD6+j] > 0.5){
                F3_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f4_old[index]+ ((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f3_old[index];
                F7_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f10_old[index]+((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f7_old[index];
                F8_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f9_old[index]+ ((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f8_old[index];
            }
            if(Q3_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(  f4_old,  F3_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6,  XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                //X_Y_XI_Intrpl7(f10_old, F7_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf37_0, XBFLf37_1, XBFLf37_2, XBFLf37_3, XBFLf37_4, XBFLf37_5, XBFLf37_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                //X_Y_XI_Intrpl7(f9_old,  F8_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf38_0, XBFLf38_1, XBFLf38_2, XBFLf38_3, XBFLf38_4, XBFLf38_5, XBFLf38_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                Y_XI_Intrpl7(  f4_old,  F3_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6,  XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                X_Y_XI_Intrpl7(f10_old, F7_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf37_0, XBFLf37_1, XBFLf37_2, XBFLf37_3, XBFLf37_4, XBFLf37_5, XBFLf37_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                X_Y_XI_Intrpl7(f9_old,  F8_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf38_0, XBFLf38_1, XBFLf38_2, XBFLf38_3, XBFLf38_4, XBFLf38_5, XBFLf38_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
            }
            //F0_in = F0_in + ModifydRho_F378( F3_in, F7_in, F8_in, f4_old[index], f9_old[index], f10_old[index] ); 
        }
        if ( BFLReqF4_d[(k-3)*NYD6+j] == 1 ){
            if(Q4_d[(k-3)*NYD6+j] > 0.5){
                F4_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f3_old[index]+ ((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f4_old[index];
                F9_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f8_old[index]+((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f9_old[index];
                F10_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f7_old[index]+ ((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f10_old[index]; 
            }
            if(Q4_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(  f3_old, F4_in,  i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                //X_Y_XI_Intrpl7(f8_old, F9_in,  i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf49_0,  XBFLf49_1,  XBFLf49_2,  XBFLf49_3,  XBFLf49_4,  XBFLf49_5,  XBFLf49_6,  YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                //X_Y_XI_Intrpl7(f7_old, F10_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf410_0, XBFLf410_1, XBFLf410_2, XBFLf410_3, XBFLf410_4, XBFLf410_5, XBFLf410_6, YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                Y_XI_Intrpl7(  f3_old, F4_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                X_Y_XI_Intrpl7(f8_old, F9_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf49_0,  XBFLf49_1,  XBFLf49_2,  XBFLf49_3,  XBFLf49_4,  XBFLf49_5,  XBFLf49_6,  YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                X_Y_XI_Intrpl7(f7_old, F10_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf410_0, XBFLf410_1, XBFLf410_2, XBFLf410_3, XBFLf410_4, XBFLf410_5, XBFLf410_6, YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
            }
            //F0_in = F0_in + ModifydRho_F4910( F4_in, F9_in, F10_in, f3_old[index], f7_old[index], f8_old[index] );
        }
        if ( BFLReqF15_d[(k-3)*NYD6+j] == 1 ){
            if(Q15_d[(k-3)*NYD6+j] > 0.5){
                F15_in = (1/(2*Q15_d[(k-3)*NYD6+j]))*f18_old[index]+ ((2*Q15_d[(k-3)*NYD6+j]-1)/(2*Q15_d[(k-3)*NYD6+j]))*f15_old[index];
            }
            if(Q15_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(f18_old, F15_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf15_0, YBFLf15_1, YBFLf15_2, YBFLf15_3, YBFLf15_4, YBFLf15_5, YBFLf15_6, XiBFLf15_0, XiBFLf15_1, XiBFLf15_2, XiBFLf15_3, XiBFLf15_4, XiBFLf15_5, XiBFLf15_6);
                Y_XI_Intrpl7(f18_old, F15_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf15_0, YBFLf15_1, YBFLf15_2, YBFLf15_3, YBFLf15_4, YBFLf15_5, YBFLf15_6, XiBFLf15_0, XiBFLf15_1, XiBFLf15_2, XiBFLf15_3, XiBFLf15_4, XiBFLf15_5, XiBFLf15_6);
            }
            //F0_in = F0_in + ModifydRho_F15( F15_in, f18_old[index] );
        }
        if ( BFLReqF16_d[(k-3)*NYD6+j] == 1 ){
            if(Q16_d[(k-3)*NYD6+j] > 0.5){
                F16_in = (1/(2*Q16_d[(k-3)*NYD6+j]))*f17_old[index]+ ((2*Q16_d[(k-3)*NYD6+j]-1)/(2*Q16_d[(k-3)*NYD6+j]))*f16_old[index];
            }
            if(Q16_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(f17_old, F16_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf16_0, YBFLf16_1, YBFLf16_2, YBFLf16_3, YBFLf16_4, YBFLf16_5, YBFLf16_6, XiBFLf16_0, XiBFLf16_1, XiBFLf16_2, XiBFLf16_3, XiBFLf16_4, XiBFLf16_5, XiBFLf16_6);
                Y_XI_Intrpl7(f17_old, F16_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf16_0, YBFLf16_1, YBFLf16_2, YBFLf16_3, YBFLf16_4, YBFLf16_5, YBFLf16_6, XiBFLf16_0, XiBFLf16_1, XiBFLf16_2, XiBFLf16_3, XiBFLf16_4, XiBFLf16_5, XiBFLf16_6);
            }
            //F0_in = F0_in + ModifydRho_F16( F16_in, f17_old[index] );
        }
    }
    /* if( k == 3 || k == 4 ) {
        idx_xi = (k-3)*NYD6+j;
        if ( BFLReqF3_d[(k-3)*NYD6+j] == 1 ){
            F3_Intrpl7( f4_old,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf3_0,   YBFLf3_1,   YBFLf3_2,   YBFLf3_3,   YBFLf3_4,   YBFLf3_5,   YBFLf3_6,   XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6 );
            X_Y_XI_Intrpl7( f10_old, F7_in,  i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf37_0,   XBFLf37_1,   XBFLf37_2,   XBFLf37_3,   XBFLf37_4,   XBFLf37_5,   XBFLf37_6,   YBFLf3_0,   YBFLf3_1,   YBFLf3_2,   YBFLf3_3,   YBFLf3_4,   YBFLf3_5,   YBFLf3_6,   XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6 );
            X_Y_XI_Intrpl7(  f9_old, F8_in,  i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf38_0,   XBFLf38_1,   XBFLf38_2,   XBFLf38_3,   XBFLf38_4,   XBFLf38_5,   XBFLf38_6,   YBFLf3_0,   YBFLf3_1,   YBFLf3_2,   YBFLf3_3,   YBFLf3_4,   YBFLf3_5,   YBFLf3_6,   XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6 );
            //F0_in = F0_in - ModifydRho_F378( F3_in, F7_in, F8_in, f4_old[index], f9_old[index], f10_old[index] );
        }
        if ( BFLReqF4_d[(k-3)*NYD6+j] == 1 ){
            F4_Intrpl7( f3_old,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6 );
            X_Y_XI_Intrpl7( f8_old,  F9_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf49_0,   XBFLf49_1,   XBFLf49_2,   XBFLf49_3,   XBFLf49_4,   XBFLf49_5,   XBFLf49_6,   YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6 );
            X_Y_XI_Intrpl7( f7_old, F10_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf410_0,  XBFLf410_1,  XBFLf410_2,  XBFLf410_3,  XBFLf410_4,  XBFLf410_5,  XBFLf410_6,  YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6 );
            //F0_in = F0_in - ModifydRho_F4910( F4_in, F9_in, F10_in, f3_old[index], f7_old[index], f8_old[index] );
        }
        if ( BFLReqF15_d[(k-3)*NYD6+j] == 1 ){
            Y_XI_Intrpl7(  f18_old, F15_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf15_0,   YBFLf15_1,   YBFLf15_2,   YBFLf15_3,   YBFLf15_4,   YBFLf15_5,   YBFLf15_6,   XiBFLf15_0, XiBFLf15_1, XiBFLf15_2, XiBFLf15_3, XiBFLf15_4, XiBFLf15_5, XiBFLf15_6 );
            //F0_in = F0_in - ModifydRho_F15( F15_in, f18_old[index] );
        }
        if ( BFLReqF16_d[(k-3)*NYD6+j] == 1 ){
            Y_XI_Intrpl7(  f17_old, F16_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf16_0,   YBFLf16_1,   YBFLf16_2,   YBFLf16_3,   YBFLf16_4,   YBFLf16_5,   YBFLf16_6,   XiBFLf16_0, XiBFLf16_1, XiBFLf16_2, XiBFLf16_3, XiBFLf16_4, XiBFLf16_5, XiBFLf16_6 );
            //F0_in = F0_in - ModifydRho_F16( F16_in, f17_old[index] );
        }
    } */
    F0_in = F0_in + rho_modify[0];//質量修正 
    __syncthreads();//線程計算同步函數 

    double rho_s = F0_in  + F1_in  + F2_in  + F3_in  + F4_in  + F5_in  + F6_in  + F7_in  + F8_in + F9_in+
				   F10_in + F11_in + F12_in + F13_in + F14_in + F15_in + F16_in + F17_in + F18_in;
	double u1 = (F1_in+ F7_in + F9_in + F11_in+ F13_in -( F2_in+F8_in +F10_in+F12_in+F14_in ) ) / rho_s ;
	double v1 = (F3_in+ F7_in + F8_in + F15_in+ F17_in -( F4_in+F9_in +F10_in+F16_in+F18_in ) ) / rho_s ;
	double w1 = (F5_in+ F11_in+ F12_in+ F15_in+ F16_in -( F6_in+F13_in+F14_in+F17_in+F18_in ) ) / rho_s ;

    double udot = u1*u1 + v1*v1 + w1*w1;


    const double F0_eq  = (1./3.)  *rho_s*(1.0-1.5*udot);
	const double F1_eq  = (1./18.) *rho_s*(1.0+3.0*u1 +4.5*u1*u1-1.5*udot);
	const double F2_eq  = (1./18.) *rho_s*(1.0-3.0*u1 +4.5*u1*u1-1.5*udot);
	const double F3_eq  = (1./18.) *rho_s*(1.0+3.0*v1 +4.5*v1*v1-1.5*udot);
	const double F4_eq  = (1./18.) *rho_s*(1.0-3.0*v1 +4.5*v1*v1-1.5*udot);
	const double F5_eq  = (1./18.) *rho_s*(1.0+3.0*w1 +4.5*w1*w1-1.5*udot);
	const double F6_eq  = (1./18.) *rho_s*(1.0-3.0*w1 +4.5*w1*w1-1.5*udot);
	const double F7_eq  = (1./36.) *rho_s*(1.0+3.0*( u1 +v1) +4.5*( u1 +v1)*( u1 +v1)-1.5*udot);
	const double F8_eq  = (1./36.) *rho_s*(1.0+3.0*(-u1 +v1) +4.5*(-u1 +v1)*(-u1 +v1)-1.5*udot);
	const double F9_eq  = (1./36.) *rho_s*(1.0+3.0*( u1 -v1) +4.5*( u1 -v1)*( u1 -v1)-1.5*udot);
	const double F10_eq = (1./36.) *rho_s*(1.0+3.0*(-u1 -v1) +4.5*(-u1 -v1)*(-u1 -v1)-1.5*udot);
	const double F11_eq = (1./36.) *rho_s*(1.0+3.0*( u1 +w1) +4.5*( u1 +w1)*( u1 +w1)-1.5*udot);
	const double F12_eq = (1./36.) *rho_s*(1.0+3.0*(-u1 +w1) +4.5*(-u1 +w1)*(-u1 +w1)-1.5*udot);
	const double F13_eq = (1./36.) *rho_s*(1.0+3.0*( u1 -w1) +4.5*( u1 -w1)*( u1 -w1)-1.5*udot);
	const double F14_eq = (1./36.) *rho_s*(1.0+3.0*(-u1 -w1) +4.5*(-u1 -w1)*(-u1 -w1)-1.5*udot);
	const double F15_eq = (1./36.) *rho_s*(1.0+3.0*( v1 +w1) +4.5*( v1 +w1)*( v1 +w1)-1.5*udot);
	const double F16_eq = (1./36.) *rho_s*(1.0+3.0*(-v1 +w1) +4.5*(-v1 +w1)*(-v1 +w1)-1.5*udot);
	const double F17_eq = (1./36.) *rho_s*(1.0+3.0*( v1 -w1) +4.5*( v1 -w1)*( v1 -w1)-1.5*udot);
	const double F18_eq = (1./36.) *rho_s*(1.0+3.0*(-v1 -w1) +4.5*(-v1 -w1)*(-v1 -w1)-1.5*udot);

    //
    m_matrix;
    meq;
    collision;
    ////

    /* F0_in  = F0_in  + 1.0 / tau * (F0_eq  - F0_in );
	F1_in  = F1_in  + 1.0 / tau * (F1_eq  - F1_in );
	F2_in  = F2_in  + 1.0 / tau * (F2_eq  - F2_in );
	F3_in  = F3_in  + 1.0 / tau * (F3_eq  - F3_in ) + 1.0/6.0*dt*Force[0];
	F4_in  = F4_in  + 1.0 / tau * (F4_eq  - F4_in ) - 1.0/6.0*dt*Force[0];
	F5_in  = F5_in  + 1.0 / tau * (F5_eq  - F5_in );
	F6_in  = F6_in  + 1.0 / tau * (F6_eq  - F6_in );
	F7_in  = F7_in  + 1.0 / tau * (F7_eq  - F7_in ) + 1.0/12.0*dt*Force[0];
	F8_in  = F8_in  + 1.0 / tau * (F8_eq  - F8_in ) + 1.0/12.0*dt*Force[0];
	F9_in  = F9_in  + 1.0 / tau * (F9_eq  - F9_in ) - 1.0/12.0*dt*Force[0];
	F10_in = F10_in + 1.0 / tau * (F10_eq - F10_in) - 1.0/12.0*dt*Force[0];
	F11_in = F11_in + 1.0 / tau * (F11_eq - F11_in);
	F12_in = F12_in + 1.0 / tau * (F12_eq - F12_in);
	F13_in = F13_in + 1.0 / tau * (F13_eq - F13_in);
	F14_in = F14_in + 1.0 / tau * (F14_eq - F14_in);
	F15_in = F15_in + 1.0 / tau * (F15_eq - F15_in) + 1.0/12.0*dt*Force[0];
	F16_in = F16_in + 1.0 / tau * (F16_eq - F16_in) - 1.0/12.0*dt*Force[0];
	F17_in = F17_in + 1.0 / tau * (F17_eq - F17_in) + 1.0/12.0*dt*Force[0];
	F18_in = F18_in + 1.0 / tau * (F18_eq - F18_in) - 1.0/12.0*dt*Force[0];
     */
    


    __syncthreads();//記憶體同步函數 
    //Straming Step
    f0_new[index] = F0_in;
	f1_new[index] = F1_in;
	f2_new[index] = F2_in;
	f3_new[index] = F3_in;
	f4_new[index] = F4_in;
	f5_new[index] = F5_in;
	f6_new[index] = F6_in;
	f7_new[index] = F7_in;
	f8_new[index] = F8_in;
	f9_new[index] = F9_in;
	f10_new[index] = F10_in;
	f11_new[index] = F11_in;
	f12_new[index] = F12_in;
	f13_new[index] = F13_in;
	f14_new[index] = F14_in;
	f15_new[index] = F15_in;
	f16_new[index] = F16_in;
	f17_new[index] = F17_in;
	f18_new[index] = F18_in;
	u[index] = u1;
	v[index] = v1;
	w[index] = w1;
	rho_d[index] = rho_s;

}

__global__ void stream_collide(
    //一般態分佈函數_舊值
    double *f0_old, double *f1_old, double *f2_old, double *f3_old, double *f4_old, double *f5_old, double *f6_old, double *f7_old, double *f8_old, double *f9_old, double *f10_old, double *f11_old, double *f12_old, double *f13_old, double *f14_old, double *f15_old, double *f16_old, double *f17_old, double *f18_old, 
    //一般態分佈函數_新值
    double *f0_new, double *f1_new, double *f2_new, double *f3_new, double *f4_new, double *f5_new, double *f6_new, double *f7_new, double *f8_new, double *f9_new, double *f10_new, double *f11_new, double *f12_new, double *f13_new, double *f14_new, double *f15_new, double *f16_new, double *f17_new, double *f18_new, 
    double *X0_0,  double *X0_1, double *X0_2,  double *X0_3,  double *X0_4,  double *X0_5,  double *X0_6,  double *X2_0,  double *X2_1,  double *X2_2,  double *X2_3,  double *X2_4,  double *X2_5,  double *X2_6,
    double *Y0_0,  double *Y0_1, double *Y0_2,  double *Y0_3,  double *Y0_4,  double *Y0_5,  double *Y0_6,  double *Y2_0,  double *Y2_1,  double *Y2_2,  double *Y2_3,  double *Y2_4,  double *Y2_5,  double *Y2_6,
    double *XiF3_0,  double *XiF3_1,  double *XiF3_2,  double *XiF3_3,  double *XiF3_4,  double *XiF3_5,  double *XiF3_6,
    double *XiF4_0,  double *XiF4_1,  double *XiF4_2,  double *XiF4_3,  double *XiF4_4,  double *XiF4_5,  double *XiF4_6,
    double *XiF5_0,  double *XiF5_1,  double *XiF5_2,  double *XiF5_3,  double *XiF5_4,  double *XiF5_5,  double *XiF5_6,
    double *XiF6_0,  double *XiF6_1,  double *XiF6_2,  double *XiF6_3,  double *XiF6_4,  double *XiF6_5,  double *XiF6_6,
    double *XiF15_0, double *XiF15_1, double *XiF15_2, double *XiF15_3, double *XiF15_4, double *XiF15_5, double *XiF15_6,
    double *XiF16_0, double *XiF16_1, double *XiF16_2, double *XiF16_3, double *XiF16_4, double *XiF16_5, double *XiF16_6,
    double *XiF17_0, double *XiF17_1, double *XiF17_2, double *XiF17_3, double *XiF17_4, double *XiF17_5, double *XiF17_6,
    double *XiF18_0, double *XiF18_1, double *XiF18_2, double *XiF18_3, double *XiF18_4, double *XiF18_5, double *XiF18_6,
    double *XBFLf37_0,   double *XBFLf37_1,   double *XBFLf37_2,   double *XBFLf37_3,   double *XBFLf37_4,   double *XBFLf37_5,   double *XBFLf37_6,
    double *XBFLf38_0,   double *XBFLf38_1,   double *XBFLf38_2,   double *XBFLf38_3,   double *XBFLf38_4,   double *XBFLf38_5,   double *XBFLf38_6,
    double *YBFLf3_0,    double *YBFLf3_1,    double *YBFLf3_2,    double *YBFLf3_3,    double *YBFLf3_4,    double *YBFLf3_5,    double *YBFLf3_6,
    double *XiBFLf3_0,   double *XiBFLf3_1,   double *XiBFLf3_2,   double *XiBFLf3_3,   double *XiBFLf3_4,   double *XiBFLf3_5,   double *XiBFLf3_6,
    double *XBFLf49_0,   double *XBFLf49_1,   double *XBFLf49_2,   double *XBFLf49_3,   double *XBFLf49_4,   double *XBFLf49_5,   double *XBFLf49_6,
    double *XBFLf410_0,  double *XBFLf410_1,  double *XBFLf410_2,  double *XBFLf410_3,  double *XBFLf410_4,  double *XBFLf410_5,  double *XBFLf410_6,
    
    
    
    //4:YZ : BFL邊條件插值，
    double *YBFLf4_0,    double *YBFLf4_1,    double *YBFLf4_2,    double *YBFLf4_3,    double *YBFLf4_4,    double *YBFLf4_5,    double *YBFLf4_6,
    double *XiBFLf4_0,   double *XiBFLf4_1,   double *XiBFLf4_2,   double *XiBFLf4_3,   double *XiBFLf4_4,   double *XiBFLf4_5,   double *XiBFLf4_6,
    //15:YZ : BFL邊條件插值，
    double *YBFLf15_0,   double *YBFLf15_1,   double *YBFLf15_2,   double *YBFLf15_3,   double *YBFLf15_4,   double *YBFLf15_5,   double *YBFLf15_6,
    double *XiBFLf15_0,  double *XiBFLf15_1,  double *XiBFLf15_2,  double *XiBFLf15_3,  double *XiBFLf15_4,  double *XiBFLf15_5,  double *XiBFLf15_6,
    //16:YZ : BFL邊條件插值，
    double *YBFLf16_0,   double *YBFLf16_1,   double *YBFLf16_2,   double *YBFLf16_3,   double *YBFLf16_4,   double *YBFLf16_5,   double *YBFLf16_6,
    double *XiBFLf16_0,  double *XiBFLf16_1,  double *XiBFLf16_2,  double *XiBFLf16_3,  double *XiBFLf16_4,  double *XiBFLf16_5,  double *XiBFLf16_6,
    //在IsBFLBCNeeded中定義判斷式，0 or 1；在BFLReqF3_h中，填入；最後再透過cudaMemcpyHostToDevice傳入BFLReqF3_d
    int *BFLReqF3_d,     int *BFLReqF4_d,     int *BFLReqF15_d,    int *BFLReqF16_d,
    //六個宏觀參數 u , v,  x , 密度 , 外力向 ,  修正密度量
    double *u,           double *v,           double *w,           double *rho_d,       double *Force,  double *rho_modify,
    //對應e_{?}方向的Boundary node 到牆面正規化距離(?:3,4,15,16)
    double *Q3_d,        double*Q4_d,         double *Q15_d,       double*Q16_d)    
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;

    if( i <= 2 || i >= NX6-3 || j <= 6 || j >= NYD6-7 || k <= 2 || k >= NZ6-3 ) return;
    //if( i <= 2 || i >= NX6-3 || j <= 2 || j >= NYD6-3 || k <= 2 || k >= NZ6-3 ) return;

    const int index = j*NX6*NZ6 + k*NX6 + i; 
          int idx;
          int idx_xi = j*NZ6 + k;//用xi來代替變換後的格點座標 
    const int nface = NX6*NZ6;
    const int nline = NX6;

    double F0_in,  F1_in,  F2_in,  F3_in,  F4_in,  F5_in,  F6_in,  F7_in,  F8_in,  F9_in;
	double F10_in, F11_in, F12_in, F13_in, F14_in, F15_in, F16_in, F17_in, F18_in;

    //MRT Variable//
    //m為代有物理意義
    double m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18;
	double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18;
	double meq0,meq1,meq2,meq3,meq4,meq5,meq6,meq7,meq8,meq9,meq10,meq11,meq12,meq13,meq14,meq15,meq16,meq17,meq18;
    /////////////

    int cell_z = k-3;
    if( k <= 6 ) cell_z = 3;
    if( k >= NZ6-7 ) cell_z = NZ6-10;
    
    // Matrix //
    Matrix;
    Inverse_Matrix;
    Relaxation;//計算鬆弛矩陣 
    ////////////

    //F0_in  = f0_old[index];
	//F1_in  = f1_old[index-1];
	//F2_in  = f2_old[index+1];
	//F3_in  = f3_old[index-NX6*NZ6];
	//F4_in  = f4_old[index+NX6*NZ6];
	//F5_in  = f5_old[index-NX6];
	//F6_in  = f6_old[index+NX6];
	//F7_in  = f7_old[index-NX6*NZ6-1];
	//F8_in  = f8_old[index-NX6*NZ6+1];
	//F9_in  = f9_old[index+NX6*NZ6-1];
	//F10_in = f10_old[index+NX6*NZ6+1];
	//F11_in = f11_old[index-NX6-1];
	//F12_in = f12_old[index-NX6+1];
	//F13_in = f13_old[index+NX6-1];
	//F14_in = f14_old[index+NX6+1];
	//F15_in = f15_old[index-NX6*NZ6-NX6];
	//F16_in = f16_old[index+NX6*NZ6-NX6];
	//F17_in = f17_old[index-NX6*NZ6+NX6];
	//F18_in = f18_old[index+NX6*NZ6+NX6];
    //IBLBM 第一個步驟: 內插(interpolation)
    F0_Intrpl7( f0_old,  i, j, k);
    //只要相對應離散速度場有y分量者都要做z方向的插值因爲晶格座標為三維直角坐標系 但是在山丘坐標系中，同樣一個Z座標會產生偏移，不在同一個高度上
    //一維插值
    F1_Intrpl7( f1_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6 );
    F2_Intrpl7( f2_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6 );
    //二維插值
    F3_Intrpl7( f3_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    F4_Intrpl7( f4_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    //一維插值
    F5_Intrpl7( f5_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    F6_Intrpl7( f6_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    //三維插值
    X_Y_XI_Intrpl7( f7_old, F7_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    X_Y_XI_Intrpl7( f8_old, F8_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    X_Y_XI_Intrpl7( f9_old, F9_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    X_Y_XI_Intrpl7(f10_old, F10_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    
    X_XI_Intrpl7(  f11_old, F11_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    X_XI_Intrpl7(  f12_old, F12_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    X_XI_Intrpl7(  f13_old, F13_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    X_XI_Intrpl7(  f14_old, F14_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    Y_XI_Intrpl7(  f15_old, F15_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF15_0, XiF15_1, XiF15_2, XiF15_3, XiF15_4, XiF15_5, XiF15_6);
    Y_XI_Intrpl7(  f16_old, F16_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF16_0, XiF16_1, XiF16_2, XiF16_3, XiF16_4, XiF16_5, XiF16_6);
    Y_XI_Intrpl7(  f17_old, F17_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF17_0, XiF17_1, XiF17_2, XiF17_3, XiF17_4, XiF17_5, XiF17_6);
    Y_XI_Intrpl7(  f18_old, F18_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF18_0, XiF18_1, XiF18_2, XiF18_3, XiF18_4, XiF18_5, XiF18_6);
    
    
//====邊界處理============================================================================================================================================================================================
    //上壁面與下壁面邊界處理 
    //邊界條件處理錯誤，原分佈函數並沒有經過碰撞，就給予遷移(邊界條件)
    //邊界條件應該作為遷移步的一種特例，在Half-wayBounce Back 中，應該先進行碰撞，不須插值，再給予轉向(遷移步被邊界條件改)
    if( k == 3 ){//這就是遷移步，處理轉向問題。
        F5_in  = f6_old[index];
        //這個也錯，因為碰撞方程式的輸入應該是6跟6平衡態，若按此寫法，是6跟5平衡態進行碰撞
        F11_in = f14_old[index];
        F12_in = f13_old[index];
        F15_in = f18_old[index];
        F16_in = f17_old[index];
    }
    if( k == NZ6-4 ){
        F6_in  = f5_old[index];
        F13_in = f12_old[index];
        F14_in = f11_old[index];
        F17_in = f16_old[index];
        F18_in = f15_old[index];
    }
    /* if( k == 3 ){
        F5_in  = f6_old[index];
        F11_in = f14_old[index];
        F12_in = f13_old[index];
    } */
    __syncthreads();//線程計算同步函數
    //BFL 
    //邊界處理在於碰撞步之後，應視為遷移步的特例
    if( k == 3 || k == 4 ) {
        idx_xi = (k-3)*NYD6+j;
        //BFL: 插值邊界處理，Req:是否需要Requirment，F3與streaming wise的方向，_d代表device
        //(1).
        if ( BFLReqF3_d[idx_xi] == 1 ){//曲面處理判斷式
            //正規化距離判斷
            if(Q3_d[(k-3)*NYD6+j] > 0.5){//若>0.5，則做線性內插 
                //註記:(2Q-1)就是計算點到反彈後的點的距離。
                //錯誤點(1.)整條方程式應該放在碰撞步之後，也就是差值的點 第一項的timeStep = n+1(遷移後)，t+dt , 第二項:n+1 ，t+dt
                //錯誤點(2.)被更新點應該是Boundary Node且轉向，為INBB Method 的第四步驟，轉向轉錯了，4,3,4(3取為遷移前Boundary Node上的值，tinmestep = n+0.5 t = 0.5dt)
                //錯誤點(3.)(2q-1)配對的值位置應該用index-1
                F3_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f4_old[index]+ ((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f3_old[index];
                F7_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f10_old[index]+((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f7_old[index];
                F8_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f9_old[index]+ ((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f8_old[index];
            }
            if(Q3_d[(k-3)*NYD6+j] < 0.5) {//若<0.5，做高階內插
                //Y_XI_Intrpl7(  f4_old,  F3_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6,  XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                //X_Y_XI_Intrpl7(f10_old, F7_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf37_0, XBFLf37_1, XBFLf37_2, XBFLf37_3, XBFLf37_4, XBFLf37_5, XBFLf37_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                //X_Y_XI_Intrpl7(f9_old,  F8_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf38_0, XBFLf38_1, XBFLf38_2, XBFLf38_3, XBFLf38_4, XBFLf38_5, XBFLf38_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                //(1.)(2.)
                Y_XI_Intrpl7(  f4_old,  F3_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6,  XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);//插值加上轉向，寫得漂亮，但是缺少碰撞步驟
                X_Y_XI_Intrpl7(f10_old, F7_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf37_0, XBFLf37_1, XBFLf37_2, XBFLf37_3, XBFLf37_4, XBFLf37_5, XBFLf37_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                X_Y_XI_Intrpl7(f9_old,  F8_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf38_0, XBFLf38_1, XBFLf38_2, XBFLf38_3, XBFLf38_4, XBFLf38_5, XBFLf38_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
            // 這裡的轉向(遷移步寫完整)寫得很漂亮，但是因為缺少碰撞應該是失誤了。
            //此寫法為，先內插，計算內插後平衡態，碰撞，賦值(錯)
            //應該改成: 先碰撞。再內插，最後賦值
            //F0_in = F0_in + ModifydRho_F378( F3_in, F7_in, F8_in, f4_old[index], f9_old[index], f10_old[index] ); 
        }
        //(2).
       if ( BFLReqF4_d[(k-3)*NYD6+j] == 1 ){
            if(Q4_d[(k-3)*NYD6+j] > 0.5){
                F4_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f3_old[index]+ ((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f4_old[index];//首先，這裡的轉向做得很漂亮，但是後面的f4應該要用新的值，而不應該用舊的值，該有差分位置應該要用index-1，取往後走一格 
                F9_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f8_old[index]+((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f9_old[index];
                F10_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f7_old[index]+ ((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f10_old[index]; 
            }
            if(Q4_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(  f3_old, F4_in,  i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                //X_Y_XI_Intrpl7(f8_old, F9_in,  i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf49_0,  XBFLf49_1,  XBFLf49_2,  XBFLf49_3,  XBFLf49_4,  XBFLf49_5,  XBFLf49_6,  YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                //X_Y_XI_Intrpl7(f7_old, F10_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf410_0, XBFLf410_1, XBFLf410_2, XBFLf410_3, XBFLf410_4, XBFLf410_5, XBFLf410_6, YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                Y_XI_Intrpl7(  f3_old, F4_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                X_Y_XI_Intrpl7(f8_old, F9_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf49_0,  XBFLf49_1,  XBFLf49_2,  XBFLf49_3,  XBFLf49_4,  XBFLf49_5,  XBFLf49_6,  YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                X_Y_XI_Intrpl7(f7_old, F10_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf410_0, XBFLf410_1, XBFLf410_2, XBFLf410_3, XBFLf410_4, XBFLf410_5, XBFLf410_6, YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
            }
            //F0_in = F0_in + ModifydRho_F4910( F4_in, F9_in, F10_in, f3_old[index], f7_old[index], f8_old[index] );
        }
        //(3).
        if ( BFLReqF15_d[(k-3)*NYD6+j] == 1 ){
            if(Q15_d[(k-3)*NYD6+j] > 0.5){
                F15_in = (1/(2*Q15_d[(k-3)*NYD6+j]))*f18_old[index]+ ((2*Q15_d[(k-3)*NYD6+j]-1)/(2*Q15_d[(k-3)*NYD6+j]))*f15_old[index];
            }
            if(Q15_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(f18_old, F15_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf15_0, YBFLf15_1, YBFLf15_2, YBFLf15_3, YBFLf15_4, YBFLf15_5, YBFLf15_6, XiBFLf15_0, XiBFLf15_1, XiBFLf15_2, XiBFLf15_3, XiBFLf15_4, XiBFLf15_5, XiBFLf15_6);
                Y_XI_Intrpl7(f18_old, F15_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf15_0, YBFLf15_1, YBFLf15_2, YBFLf15_3, YBFLf15_4, YBFLf15_5, YBFLf15_6, XiBFLf15_0, XiBFLf15_1, XiBFLf15_2, XiBFLf15_3, XiBFLf15_4, XiBFLf15_5, XiBFLf15_6);
            }
            //F0_in = F0_in + ModifydRho_F15( F15_in, f18_old[index] );
        }
        //(4).
        if ( BFLReqF16_d[(k-3)*NYD6+j] == 1 ){
            if(Q16_d[(k-3)*NYD6+j] > 0.5){
                F16_in = (1/(2*Q16_d[(k-3)*NYD6+j]))*f17_old[index]+ ((2*Q16_d[(k-3)*NYD6+j]-1)/(2*Q16_d[(k-3)*NYD6+j]))*f16_old[index];
            }
            if(Q16_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(f17_old, F16_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf16_0, YBFLf16_1, YBFLf16_2, YBFLf16_3, YBFLf16_4, YBFLf16_5, YBFLf16_6, XiBFLf16_0, XiBFLf16_1, XiBFLf16_2, XiBFLf16_3, XiBFLf16_4, XiBFLf16_5, XiBFLf16_6);
                Y_XI_Intrpl7(f17_old, F16_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf16_0, YBFLf16_1, YBFLf16_2, YBFLf16_3, YBFLf16_4, YBFLf16_5, YBFLf16_6, XiBFLf16_0, XiBFLf16_1, XiBFLf16_2, XiBFLf16_3, XiBFLf16_4, XiBFLf16_5, XiBFLf16_6);
            }
            //F0_in = F0_in + ModifydRho_F16( F16_in, f17_old[index] );
        }
    }
    /* if( k == 3 || k == 4 ) {
        idx_xi = (k-3)*NYD6+j;
        if ( BFLReqF3_d[(k-3)*NYD6+j] == 1 ){
            F3_Intrpl7( f4_old,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf3_0,   YBFLf3_1,   YBFLf3_2,   YBFLf3_3,   YBFLf3_4,   YBFLf3_5,   YBFLf3_6,   XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6 );
            X_Y_XI_Intrpl7( f10_old, F7_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf37_0,   XBFLf37_1,   XBFLf37_2,   XBFLf37_3,   XBFLf37_4,   XBFLf37_5,   XBFLf37_6,   YBFLf3_0,   YBFLf3_1,   YBFLf3_2,   YBFLf3_3,   YBFLf3_4,   YBFLf3_5,   YBFLf3_6,   XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6 );
            X_Y_XI_Intrpl7(  f9_old, F8_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf38_0,   XBFLf38_1,   XBFLf38_2,   XBFLf38_3,   XBFLf38_4,   XBFLf38_5,   XBFLf38_6,   YBFLf3_0,   YBFLf3_1,   YBFLf3_2,   YBFLf3_3,   YBFLf3_4,   YBFLf3_5,   YBFLf3_6,   XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6 );
            //F0_in = F0_in - ModifydRho_F378( F3_in, F7_in, F8_in, f4_old[index], f9_old[index], f10_old[index] ); 
        }
        if ( BFLReqF4_d[(k-3)*NYD6+j] == 1 ){
            F4_Intrpl7( f3_old,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6 );
            X_Y_XI_Intrpl7( f8_old,  F9_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf49_0,   XBFLf49_1,   XBFLf49_2,   XBFLf49_3,   XBFLf49_4,   XBFLf49_5,   XBFLf49_6,   YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6 );
            X_Y_XI_Intrpl7( f7_old, F10_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf410_0,  XBFLf410_1,  XBFLf410_2,  XBFLf410_3,  XBFLf410_4,  XBFLf410_5,  XBFLf410_6,  YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6 );
            //F0_in = F0_in - ModifydRho_F4910( F4_in, F9_in, F10_in, f3_old[index], f7_old[index], f8_old[index] );
        }
        if ( BFLReqF15_d[(k-3)*NYD6+j] == 1 ){
            Y_XI_Intrpl7(  f18_old, F15_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf15_0,   YBFLf15_1,   YBFLf15_2,   YBFLf15_3,   YBFLf15_4,   YBFLf15_5,   YBFLf15_6,   XiBFLf15_0, XiBFLf15_1, XiBFLf15_2, XiBFLf15_3, XiBFLf15_4, XiBFLf15_5, XiBFLf15_6 );
            //F0_in = F0_in - ModifydRho_F15( F15_in, f18_old[index] );
        }
        if ( BFLReqF16_d[(k-3)*NYD6+j] == 1 ){
            Y_XI_Intrpl7(  f17_old, F16_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf16_0,   YBFLf16_1,   YBFLf16_2,   YBFLf16_3,   YBFLf16_4,   YBFLf16_5,   YBFLf16_6,   XiBFLf16_0, XiBFLf16_1, XiBFLf16_2, XiBFLf16_3, XiBFLf16_4, XiBFLf16_5, XiBFLf16_6 );
            //F0_in = F0_in - ModifydRho_F16( F16_in, f17_old[index] );
        }
    } */
//====邊界處理結束====================================================================================================================================================================================================================================================================================================================================
    F0_in = F0_in + rho_modify[0];
    __syncthreads();
//====IBLBM第二個步驟: 計算平衡態分佈函數===============================================================================================================================================================================================================================================================================================================
    double rho_s = F0_in  + F1_in  + F2_in  + F3_in  + F4_in  + F5_in  + F6_in  + F7_in  + F8_in + F9_in+
				   F10_in + F11_in + F12_in + F13_in + F14_in + F15_in + F16_in + F17_in + F18_in;
	double u1 = (F1_in+ F7_in + F9_in + F11_in+ F13_in -( F2_in+F8_in +F10_in+F12_in+F14_in ) ) / rho_s ;
	double v1 = (F3_in+ F7_in + F8_in + F15_in+ F17_in -( F4_in+F9_in +F10_in+F16_in+F18_in ) ) / rho_s ;
	double w1 = (F5_in+ F11_in+ F12_in+ F15_in+ F16_in -( F6_in+F13_in+F14_in+F17_in+F18_in ) ) / rho_s ;
    //速度場已經由插值得出//這裡的速度場要小心 
    double udot = u1*u1 + v1*v1 + w1*w1;
   
    const double F0_eq  = (1./3.)  *rho_s*(1.0-1.5*udot);
	const double F1_eq  = (1./18.) *rho_s*(1.0+3.0*u1 +4.5*u1*u1-1.5*udot);
	const double F2_eq  = (1./18.) *rho_s*(1.0-3.0*u1 +4.5*u1*u1-1.5*udot);
	const double F3_eq  = (1./18.) *rho_s*(1.0+3.0*v1 +4.5*v1*v1-1.5*udot);
	const double F4_eq  = (1./18.) *rho_s*(1.0-3.0*v1 +4.5*v1*v1-1.5*udot);
	const double F5_eq  = (1./18.) *rho_s*(1.0+3.0*w1 +4.5*w1*w1-1.5*udot);
	const double F6_eq  = (1./18.) *rho_s*(1.0-3.0*w1 +4.5*w1*w1-1.5*udot);
	const double F7_eq  = (1./36.) *rho_s*(1.0+3.0*( u1 +v1) +4.5*( u1 +v1)*( u1 +v1)-1.5*udot);
	const double F8_eq  = (1./36.) *rho_s*(1.0+3.0*(-u1 +v1) +4.5*(-u1 +v1)*(-u1 +v1)-1.5*udot);
	const double F9_eq  = (1./36.) *rho_s*(1.0+3.0*( u1 -v1) +4.5*( u1 -v1)*( u1 -v1)-1.5*udot);
	const double F10_eq = (1./36.) *rho_s*(1.0+3.0*(-u1 -v1) +4.5*(-u1 -v1)*(-u1 -v1)-1.5*udot);
	const double F11_eq = (1./36.) *rho_s*(1.0+3.0*( u1 +w1) +4.5*( u1 +w1)*( u1 +w1)-1.5*udot);
	const double F12_eq = (1./36.) *rho_s*(1.0+3.0*(-u1 +w1) +4.5*(-u1 +w1)*(-u1 +w1)-1.5*udot);
	const double F13_eq = (1./36.) *rho_s*(1.0+3.0*( u1 -w1) +4.5*( u1 -w1)*( u1 -w1)-1.5*udot);
	const double F14_eq = (1./36.) *rho_s*(1.0+3.0*(-u1 -w1) +4.5*(-u1 -w1)*(-u1 -w1)-1.5*udot);
	const double F15_eq = (1./36.) *rho_s*(1.0+3.0*( v1 +w1) +4.5*( v1 +w1)*( v1 +w1)-1.5*udot);
	const double F16_eq = (1./36.) *rho_s*(1.0+3.0*(-v1 +w1) +4.5*(-v1 +w1)*(-v1 +w1)-1.5*udot);
	const double F17_eq = (1./36.) *rho_s*(1.0+3.0*( v1 -w1) +4.5*( v1 -w1)*( v1 -w1)-1.5*udot);
	const double F18_eq = (1./36.) *rho_s*(1.0+3.0*(-v1 -w1) +4.5*(-v1 -w1)*(-v1 -w1)-1.5*udot);
    //上述為插值過後的平衡態分佈函數
//=========第二個步驟結束==============================================================================================================================================================================================================
    //
    m_matrix;
    meq;
    collision;//碰撞步 
    ////
    /* F0_in  = F0_in  + 1.0 / tau * (F0_eq  - F0_in );
	F1_in  = F1_in  + 1.0 / tau * (F1_eq  - F1_in );
	F2_in  = F2_in  + 1.0 / tau * (F2_eq  - F2_in );
	F3_in  = F3_in  + 1.0 / tau * (F3_eq  - F3_in ) + 1.0/6.0*dt*Force[0];
	F4_in  = F4_in  + 1.0 / tau * (F4_eq  - F4_in ) - 1.0/6.0*dt*Force[0];
	F5_in  = F5_in  + 1.0 / tau * (F5_eq  - F5_in );
	F6_in  = F6_in  + 1.0 / tau * (F6_eq  - F6_in );
	F7_in  = F7_in  + 1.0 / tau * (F7_eq  - F7_in ) + 1.0/12.0*dt*Force[0];
	F8_in  = F8_in  + 1.0 / tau * (F8_eq  - F8_in ) + 1.0/12.0*dt*Force[0];
	F9_in  = F9_in  + 1.0 / tau * (F9_eq  - F9_in ) - 1.0/12.0*dt*Force[0];
	F10_in = F10_in + 1.0 / tau * (F10_eq - F10_in) - 1.0/12.0*dt*Force[0];
	F11_in = F11_in + 1.0 / tau * (F11_eq - F11_in);
	F12_in = F12_in + 1.0 / tau * (F12_eq - F12_in);
	F13_in = F13_in + 1.0 / tau * (F13_eq - F13_in);
	F14_in = F14_in + 1.0 / tau * (F14_eq - F14_in);
	F15_in = F15_in + 1.0 / tau * (F15_eq - F15_in) + 1.0/12.0*dt*Force[0];
	F16_in = F16_in + 1.0 / tau * (F16_eq - F16_in) - 1.0/12.0*dt*Force[0];
	F17_in = F17_in + 1.0 / tau * (F17_eq - F17_in) + 1.0/12.0*dt*Force[0];
	F18_in = F18_in + 1.0 / tau * (F18_eq - F18_in) - 1.0/12.0*dt*Force[0]; */

    __syncthreads();
    //遷移步
    //非格點上的數值 交給 格點上的數值 
    f0_new[index] = F0_in;
	f1_new[index] = F1_in;
	f2_new[index] = F2_in;
	f3_new[index] = F3_in;
	f4_new[index] = F4_in;
	f5_new[index] = F5_in;
	f6_new[index] = F6_in;
	f7_new[index] = F7_in;
	f8_new[index] = F8_in;
	f9_new[index] = F9_in;
	f10_new[index] = F10_in;
	f11_new[index] = F11_in;
	f12_new[index] = F12_in;
	f13_new[index] = F13_in;
	f14_new[index] = F14_in;
	f15_new[index] = F15_in;
	f16_new[index] = F16_in;
	f17_new[index] = F17_in;
	f18_new[index] = F18_in;
	u[index] = u1;
	v[index] = v1;
	w[index] = w1;
	rho_d[index] = rho_s;
    /* if( i == 30 && j == 30 && k == 30){
         printf(" step = %d \t i = %d \t j = %d \t k = %d \t rho = %lf \n",step, i, j, k, rho_d[index]);
    } */
   
}

__global__ void periodicUD(//這裡應該是Stream-Wise 方向吧//應該修改為periodicS-W
    double *f0_old, double *f1_old, double *f2_old, double *f3_old, double *f4_old, double *f5_old, double *f6_old, double *f7_old, double *f8_old, double *f9_old, double *f10_old, double *f11_old, double *f12_old, double *f13_old, double *f14_old, double *f15_old, double *f16_old, double *f17_old, double *f18_old, 
    double *f0_new, double *f1_new, double *f2_new, double *f3_new, double *f4_new, double *f5_new, double *f6_new, double *f7_new, double *f8_new, double *f9_new, double *f10_new, double *f11_new, double *f12_new, double *f13_new, double *f14_new, double *f15_new, double *f16_new, double *f17_new, double *f18_new, 
    double *y_d,       double *x_d,      double *z_d,
    double *u,         double *v,        double *w,         double *rho_d)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;
          int idx, idx_buffer;
          int buffer = 3;
    
    if( i >= NX6 || k >= NZ6 ) return;//不考慮buffer之外的計算點 
    
    idx_buffer = j*NZ6*NX6 + k*NX6 + i; //j為noermal-wise之方向指向紙面內部//座標採用([y,z,x] = [NormalWise, UPandDown, StreamWsie])//這裡的idx-buffer 指的應該是:考慮外推層的座標系統
    idx = idx_buffer + (NYD6-2*buffer-1)*NZ6*NX6;
    //NYD6對於每個GPU分配到的總節點數量 
    f0_new[idx_buffer]  = f0_new[idx];    f1_new[idx_buffer]  = f1_new[idx];    f2_new[idx_buffer]  = f2_new[idx];
    f3_new[idx_buffer]  = f3_new[idx];    f4_new[idx_buffer]  = f4_new[idx];    f5_new[idx_buffer]  = f5_new[idx];
    f6_new[idx_buffer]  = f6_new[idx];    f7_new[idx_buffer]  = f7_new[idx];    f8_new[idx_buffer]  = f8_new[idx];
    f9_new[idx_buffer]  = f9_new[idx];    f10_new[idx_buffer] = f10_new[idx];   f11_new[idx_buffer] = f11_new[idx];
    f12_new[idx_buffer] = f12_new[idx];   f13_new[idx_buffer] = f13_new[idx];   f14_new[idx_buffer] = f14_new[idx];
    f15_new[idx_buffer] = f15_new[idx];   f16_new[idx_buffer] = f16_new[idx];   f17_new[idx_buffer] = f17_new[idx];
    f18_new[idx_buffer] = f18_new[idx];
    u[idx_buffer] = u[idx];
    v[idx_buffer] = v[idx];
    w[idx_buffer] = w[idx];
    rho_d[idx_buffer] = rho_d[idx];

    idx_buffer = (NYD6-1-j)*NZ6*NX6 + k*NX6 + i; //建立y方向的週期性邊界條件
    idx = idx_buffer - (NYD6-2*buffer-1)*NZ6*NX6;

    f0_new[idx_buffer]  = f0_new[idx];    f1_new[idx_buffer]  = f1_new[idx];    f2_new[idx_buffer]  = f2_new[idx];
    f3_new[idx_buffer]  = f3_new[idx];    f4_new[idx_buffer]  = f4_new[idx];    f5_new[idx_buffer]  = f5_new[idx];
    f6_new[idx_buffer]  = f6_new[idx];    f7_new[idx_buffer]  = f7_new[idx];    f8_new[idx_buffer]  = f8_new[idx];
    f9_new[idx_buffer]  = f9_new[idx];    f10_new[idx_buffer] = f10_new[idx];   f11_new[idx_buffer] = f11_new[idx];
    f12_new[idx_buffer] = f12_new[idx];   f13_new[idx_buffer] = f13_new[idx];   f14_new[idx_buffer] = f14_new[idx];
    f15_new[idx_buffer] = f15_new[idx];   f16_new[idx_buffer] = f16_new[idx];   f17_new[idx_buffer] = f17_new[idx];
    f18_new[idx_buffer] = f18_new[idx];
    u[idx_buffer] = u[idx];
    v[idx_buffer] = v[idx];
    w[idx_buffer] = w[idx];
    rho_d[idx_buffer] = rho_d[idx];

}

__global__ void periodicSW(//這裡應該是Normal-Wise方向//應該修改為periodicNW//反整此方向為展向方向
    double *f0_old, double *f1_old, double *f2_old, double *f3_old, double *f4_old, double *f5_old, double *f6_old, double *f7_old, double *f8_old, double *f9_old, double *f10_old, double *f11_old, double *f12_old, double *f13_old, double *f14_old, double *f15_old, double *f16_old, double *f17_old, double *f18_old, 
    double *f0_new, double *f1_new, double *f2_new, double *f3_new, double *f4_new, double *f5_new, double *f6_new, double *f7_new, double *f8_new, double *f9_new, double *f10_new, double *f11_new, double *f12_new, double *f13_new, double *f14_new, double *f15_new, double *f16_new, double *f17_new, double *f18_new, 
    double *y_d,       double *x_d,      double *z_d,
    double *u,         double *v,        double *w,         double *rho_d)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;
          int idx, idx_buffer;
          int buffer = 3;

    if( j >= NYD6 || k >= NZ6 ) return; //在y方向做GPU切割，換言之，NYD6為一個gpu所負責的場空間區域 (從Y方向做切片)
    //畫圖 
    idx_buffer = j*NZ6*NX6 + k*NX6 + i;//此為periodic場空間座標 [i = 0 , 2 ][3...34][35][36,37,38]
    idx = idx_buffer + (NX6-2*buffer-1);//(NX6-2*buffer-1)= NX 為物理空間總節點數 
    //idx為不考慮x方向buffer層的場空間座標;
    //0 <- 32 ; 1 <- 33 ; 3 <- 34
    f0_new[idx_buffer]  = f0_new[idx];    f1_new[idx_buffer]  = f1_new[idx];    f2_new[idx_buffer]  = f2_new[idx];
    f3_new[idx_buffer]  = f3_new[idx];    f4_new[idx_buffer]  = f4_new[idx];    f5_new[idx_buffer]  = f5_new[idx];
    f6_new[idx_buffer]  = f6_new[idx];    f7_new[idx_buffer]  = f7_new[idx];    f8_new[idx_buffer]  = f8_new[idx];
    f9_new[idx_buffer]  = f9_new[idx];    f10_new[idx_buffer] = f10_new[idx];   f11_new[idx_buffer] = f11_new[idx];
    f12_new[idx_buffer] = f12_new[idx];   f13_new[idx_buffer] = f13_new[idx];   f14_new[idx_buffer] = f14_new[idx];
    f15_new[idx_buffer] = f15_new[idx];   f16_new[idx_buffer] = f16_new[idx];   f17_new[idx_buffer] = f17_new[idx];
    f18_new[idx_buffer] = f18_new[idx];
    u[idx_buffer] = u[idx];
    v[idx_buffer] = v[idx];
    w[idx_buffer] = w[idx];
    rho_d[idx_buffer] = rho_d[idx];
   
    idx_buffer = j*NX6*NZ6 + k*NX6 + (NX6-1-i);//此idx_buffer為36,37,38 
    idx = idx_buffer - (NX6-2*buffer-1); //(NX6-2*buffer-1)=32//4,5,6
    
    f0_new[idx_buffer]  = f0_new[idx];    f1_new[idx_buffer]  = f1_new[idx];    f2_new[idx_buffer]  = f2_new[idx];
    f3_new[idx_buffer]  = f3_new[idx];    f4_new[idx_buffer]  = f4_new[idx];    f5_new[idx_buffer]  = f5_new[idx];
    f6_new[idx_buffer]  = f6_new[idx];    f7_new[idx_buffer]  = f7_new[idx];    f8_new[idx_buffer]  = f8_new[idx];
    f9_new[idx_buffer]  = f9_new[idx];    f10_new[idx_buffer] = f10_new[idx];   f11_new[idx_buffer] = f11_new[idx];
    f12_new[idx_buffer] = f12_new[idx];   f13_new[idx_buffer] = f13_new[idx];   f14_new[idx_buffer] = f14_new[idx];
    f15_new[idx_buffer] = f15_new[idx];   f16_new[idx_buffer] = f16_new[idx];   f17_new[idx_buffer] = f17_new[idx];
    f18_new[idx_buffer] = f18_new[idx];
    u[idx_buffer] = u[idx];
    v[idx_buffer] = v[idx];
    w[idx_buffer] = w[idx];
    rho_d[idx_buffer] = rho_d[idx];

}

__global__ void periodicNML(//這裡應該是Up and Down 應該修改為periodicUD 
    double *f0_old, double *f1_old, double *f2_old, double *f3_old, double *f4_old, double *f5_old, double *f6_old, double *f7_old, double *f8_old, double *f9_old, double *f10_old, double *f11_old, double *f12_old, double *f13_old, double *f14_old, double *f15_old, double *f16_old, double *f17_old, double *f18_old, 
    double *f0_new, double *f1_new, double *f2_new, double *f3_new, double *f4_new, double *f5_new, double *f6_new, double *f7_new, double *f8_new, double *f9_new, double *f10_new, double *f11_new, double *f12_new, double *f13_new, double *f14_new, double *f15_new, double *f16_new, double *f17_new, double *f18_new, 
    double *y_d,       double *x_d,      double *z_d,
    double *u,         double *v,        double *w,         double *rho_d)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;
          int idx, idx_buffer;
          int buffer = 3;
    //dim3 griddimSW(  1,      NYD6/NT+1, NZ6);
    //dim3 blockdimSW( buffer, NT,        1 );

    if( i >= NX6 || j >= NYD6 ) return;

    idx_buffer = j*NX6*NZ6 + k*NX6 + i;
    idx = idx_buffer + (NZ6-2*buffer-1)*NX6;//建立z方向的週期性邊界條件

    f0_new[idx_buffer]  = f0_new[idx];    f1_new[idx_buffer]  = f1_new[idx];    f2_new[idx_buffer]  = f2_new[idx];
    f3_new[idx_buffer]  = f3_new[idx];    f4_new[idx_buffer]  = f4_new[idx];    f5_new[idx_buffer]  = f5_new[idx];
    f6_new[idx_buffer]  = f6_new[idx];    f7_new[idx_buffer]  = f7_new[idx];    f8_new[idx_buffer]  = f8_new[idx];
    f9_new[idx_buffer]  = f9_new[idx];    f10_new[idx_buffer] = f10_new[idx];   f11_new[idx_buffer] = f11_new[idx];
    f12_new[idx_buffer] = f12_new[idx];   f13_new[idx_buffer] = f13_new[idx];   f14_new[idx_buffer] = f14_new[idx];
    f15_new[idx_buffer] = f15_new[idx];   f16_new[idx_buffer] = f16_new[idx];   f17_new[idx_buffer] = f17_new[idx];
    f18_new[idx_buffer] = f18_new[idx];
    u[idx_buffer] = u[idx];
    v[idx_buffer] = v[idx];
    w[idx_buffer] = w[idx];
    rho_d[idx_buffer] = rho_d[idx];

    idx_buffer = j*NX6*NZ6 + (NZ6-1-k)*NX6 + i;
    idx = idx_buffer - (NZ6-2*buffer-1)*NX6;

    f0_new[idx_buffer]  = f0_new[idx];    f1_new[idx_buffer]  = f1_new[idx];    f2_new[idx_buffer]  = f2_new[idx];
    f3_new[idx_buffer]  = f3_new[idx];    f4_new[idx_buffer]  = f4_new[idx];    f5_new[idx_buffer]  = f5_new[idx];
    f6_new[idx_buffer]  = f6_new[idx];    f7_new[idx_buffer]  = f7_new[idx];    f8_new[idx_buffer]  = f8_new[idx];
    f9_new[idx_buffer]  = f9_new[idx];    f10_new[idx_buffer] = f10_new[idx];   f11_new[idx_buffer] = f11_new[idx];
    f12_new[idx_buffer] = f12_new[idx];   f13_new[idx_buffer] = f13_new[idx];   f14_new[idx_buffer] = f14_new[idx];
    f15_new[idx_buffer] = f15_new[idx];   f16_new[idx_buffer] = f16_new[idx];   f17_new[idx_buffer] = f17_new[idx];
    f18_new[idx_buffer] = f18_new[idx];
    u[idx_buffer] = u[idx];
    v[idx_buffer] = v[idx];
    w[idx_buffer] = w[idx];
    rho_d[idx_buffer] = rho_d[idx];
}

__global__ void AccumulateUbulk(
    double *Ub_avg,     double *v,
    double *x,          double *z  )
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;

    if( i <= 2 || i >= NX6-3 || k <= 2 || k >= NZ6-3 ) return;

    double dx = ( x[i+1] - x[i-1] ) / 2.0;
    double dz = ( z[j*NZ6+k+1] - z[j*NZ6+k-1] ) / 2.0;

    Ub_avg[k*NX6+i] += v[j*NZ6*NX6+k*NX6+i] * dx * dz;
}

void Launch_CollisionStreaming(double *f_old[19], double *f_new[19]) {
    int buffer = 3;

    dim3 griddimSW(  1,      NYD6/NT+1, NZ6);//blockIdx.x ∈ {0}//blockIdx.y ∈ {0, 1, …, (NYD6/NT)}//blockIdx.z ∈ {0, 1, …, NZ6−1}
    dim3 blockdimSW( buffer, NT,   1);//threadIdx.x ∈ [0, buffer−1]（= [0,2] 若 buffer=3）//threadIdx.y ∈ [0, NT−1]//threadIdx.z ∈ {0}
    dim3 griddimUD(  NX6/NT+1, 1,  NZ6);
    dim3 blockdimUD( NT,       buffer, 1 );
    dim3 griddimNML( NX6/NT+1, NYD6,   1);
    dim3 blockdimNML(NT,       1,      buffer);


    dim3 griddim(  NX6/NT+1, NYD6, NZ6);
    dim3 blockdim( NT, 1, 1);

    dim3 griddimBuf(NX6/NT+1, 1, NZ6);
    dim3 blockdimBuf(NT, 4, 1);

    stream_collide_Buffer<<<griddimBuf, blockdimBuf, 0, stream1>>>(
    f_old[0],   f_old[1],   f_old[2],   f_old[3],   f_old[4],   f_old[5],   f_old[6],   
    f_old[7],   f_old[8],   f_old[9],   f_old[10],  f_old[11],  f_old[12],  f_old[13],  
    f_old[14],  f_old[15],  f_old[16],  f_old[17],  f_old[18],
    f_new[0],   f_new[1],   f_new[2],   f_new[3],   f_new[4],   f_new[5],   f_new[6],   
    f_new[7],   f_new[8],   f_new[9],   f_new[10],  f_new[11],  f_new[12],  f_new[13],  
    f_new[14],  f_new[15],  f_new[16],  f_new[17],  f_new[18],
    
    XPara0_d[0],  XPara0_d[1],  XPara0_d[2],    XPara0_d[3],    XPara0_d[4],    XPara0_d[5],    XPara0_d[6],    
    XPara2_d[0],    XPara2_d[1],    XPara2_d[2],    XPara2_d[3],    XPara2_d[4],    XPara2_d[5],    XPara2_d[6],
    YPara0_d[0],  YPara0_d[1],  YPara0_d[2],    YPara0_d[3],    YPara0_d[4],    YPara0_d[5],    YPara0_d[6],    
    YPara2_d[0],    YPara2_d[1],    YPara2_d[2],    YPara2_d[3],    YPara2_d[4],    YPara2_d[5],    YPara2_d[6],
    
    XiParaF3_d[0],  XiParaF3_d[1],  XiParaF3_d[2],  XiParaF3_d[3],  XiParaF3_d[4],  XiParaF3_d[5],  XiParaF3_d[6],
    XiParaF4_d[0],  XiParaF4_d[1],  XiParaF4_d[2],  XiParaF4_d[3],  XiParaF4_d[4],  XiParaF4_d[5],  XiParaF4_d[6],
    XiParaF5_d[0],  XiParaF5_d[1],  XiParaF5_d[2],  XiParaF5_d[3],  XiParaF5_d[4],  XiParaF5_d[5],  XiParaF5_d[6],
    XiParaF6_d[0],  XiParaF6_d[1],  XiParaF6_d[2],  XiParaF6_d[3],  XiParaF6_d[4],  XiParaF6_d[5],  XiParaF6_d[6],
    XiParaF15_d[0], XiParaF15_d[1], XiParaF15_d[2], XiParaF15_d[3], XiParaF15_d[4], XiParaF15_d[5], XiParaF15_d[6],
    XiParaF16_d[0], XiParaF16_d[1], XiParaF16_d[2], XiParaF16_d[3], XiParaF16_d[4], XiParaF16_d[5], XiParaF16_d[6],
    XiParaF17_d[0], XiParaF17_d[1], XiParaF17_d[2], XiParaF17_d[3], XiParaF17_d[4], XiParaF17_d[5], XiParaF17_d[6],
    XiParaF18_d[0], XiParaF18_d[1], XiParaF18_d[2], XiParaF18_d[3], XiParaF18_d[4], XiParaF18_d[5], XiParaF18_d[6],

    XBFLParaF37_d[0]  ,   XBFLParaF37_d[1]  ,   XBFLParaF37_d[2]  ,   XBFLParaF37_d[3]  ,   XBFLParaF37_d[4]  ,   XBFLParaF37_d[5]  ,   XBFLParaF37_d[6]  ,
    XBFLParaF38_d[0]  ,   XBFLParaF38_d[1]  ,   XBFLParaF38_d[2]  ,   XBFLParaF38_d[3]  ,   XBFLParaF38_d[4]  ,   XBFLParaF38_d[5]  ,   XBFLParaF38_d[6]  ,
    YBFLParaF378_d[0] ,   YBFLParaF378_d[1] ,   YBFLParaF378_d[2] ,   YBFLParaF378_d[3] ,   YBFLParaF378_d[4] ,   YBFLParaF378_d[5] ,   YBFLParaF378_d[6] ,
    XiBFLParaF378_d[0],   XiBFLParaF378_d[1],   XiBFLParaF378_d[2],   XiBFLParaF378_d[3],   XiBFLParaF378_d[4],   XiBFLParaF378_d[5],   XiBFLParaF378_d[6],
    
    XBFLParaF49_d[0]   ,  XBFLParaF49_d[1]   ,  XBFLParaF49_d[2]   ,  XBFLParaF49_d[3]   ,  XBFLParaF49_d[4]   ,  XBFLParaF49_d[5]   ,  XBFLParaF49_d[6]   ,
    XBFLParaF410_d[0]  ,  XBFLParaF410_d[1]  ,  XBFLParaF410_d[2]  ,  XBFLParaF410_d[3]  ,  XBFLParaF410_d[4]  ,  XBFLParaF410_d[5]  ,  XBFLParaF410_d[6]  ,
    YBFLParaF4910_d[0] ,  YBFLParaF4910_d[1] ,  YBFLParaF4910_d[2] ,  YBFLParaF4910_d[3] ,  YBFLParaF4910_d[4] ,  YBFLParaF4910_d[5] ,  YBFLParaF4910_d[6] ,
    XiBFLParaF4910_d[0],  XiBFLParaF4910_d[1],  XiBFLParaF4910_d[2],  XiBFLParaF4910_d[3],  XiBFLParaF4910_d[4],  XiBFLParaF4910_d[5],  XiBFLParaF4910_d[6],
    
    YBFLParaF15_d[0] ,    YBFLParaF15_d[1] ,    YBFLParaF15_d[2] ,    YBFLParaF15_d[3] ,    YBFLParaF15_d[4] ,    YBFLParaF15_d[5] ,    YBFLParaF15_d[6] ,
    XiBFLParaF15_d[0],    XiBFLParaF15_d[1],    XiBFLParaF15_d[2],    XiBFLParaF15_d[3],    XiBFLParaF15_d[4],    XiBFLParaF15_d[5],    XiBFLParaF15_d[6],
    
    YBFLParaF16_d[0] ,    YBFLParaF16_d[1] ,    YBFLParaF16_d[2] ,    YBFLParaF16_d[3] ,    YBFLParaF16_d[4] ,    YBFLParaF16_d[5] ,    YBFLParaF16_d[6] ,
    XiBFLParaF16_d[0],    XiBFLParaF16_d[1],    XiBFLParaF16_d[2],    XiBFLParaF16_d[3],    XiBFLParaF16_d[4],    XiBFLParaF16_d[5],    XiBFLParaF16_d[6],
    
    BFLReqF3_d,     BFLReqF4_d,     BFLReqF15_d,    BFLReqF16_d,
    u,      v,      w,      rho_d,  Force_d,    3, rho_modify_d, Q3_d, Q4_d, Q15_d, Q16_d
    );
    stream_collide_Buffer<<<griddimBuf, blockdimBuf, 0, stream1>>>(
    //一般態分佈含數據矩陣來表示各分量
    //一般態分佈函數_舊值 
    f_old[0],   f_old[1],   f_old[2],   f_old[3],   f_old[4],   f_old[5],   f_old[6],   f_old[7],   f_old[8],   f_old[9],   f_old[10],  f_old[11],  f_old[12],  f_old[13],  f_old[14],  f_old[15],  f_old[16],  f_old[17],  f_old[18],
    //一般態分佈函數_舊值 
    f_new[0],   f_new[1],   f_new[2],   f_new[3],   f_new[4],   f_new[5],   f_new[6],   f_new[7],   f_new[8],   f_new[9],   f_new[10],  f_new[11],  f_new[12],  f_new[13],  f_new[14],  f_new[15],  f_new[16],  f_new[17],  f_new[18],
    //*_h 表示host端的陣列，*_d表示device端的陣列 
    XPara0_d[0],  XPara0_d[1],  XPara0_d[2],    XPara0_d[3],    XPara0_d[4],    XPara0_d[5],    XPara0_d[6],    XPara2_d[0],    XPara2_d[1],    XPara2_d[2],    XPara2_d[3],    XPara2_d[4],    XPara2_d[5],    XPara2_d[6],
    YPara0_d[0],  YPara0_d[1],  YPara0_d[2],    YPara0_d[3],    YPara0_d[4],    YPara0_d[5],    YPara0_d[6],    YPara2_d[0],    YPara2_d[1],    YPara2_d[2],    YPara2_d[3],    YPara2_d[4],    YPara2_d[5],    YPara2_d[6],
    XiParaF3_d[0],  XiParaF3_d[1],  XiParaF3_d[2],  XiParaF3_d[3],  XiParaF3_d[4],  XiParaF3_d[5],  XiParaF3_d[6],
    XiParaF4_d[0],  XiParaF4_d[1],  XiParaF4_d[2],  XiParaF4_d[3],  XiParaF4_d[4],  XiParaF4_d[5],  XiParaF4_d[6],
    XiParaF5_d[0],  XiParaF5_d[1],  XiParaF5_d[2],  XiParaF5_d[3],  XiParaF5_d[4],  XiParaF5_d[5],  XiParaF5_d[6],
    XiParaF6_d[0],  XiParaF6_d[1],  XiParaF6_d[2],  XiParaF6_d[3],  XiParaF6_d[4],  XiParaF6_d[5],  XiParaF6_d[6],
    XiParaF15_d[0], XiParaF15_d[1], XiParaF15_d[2], XiParaF15_d[3], XiParaF15_d[4], XiParaF15_d[5], XiParaF15_d[6],
    XiParaF16_d[0], XiParaF16_d[1], XiParaF16_d[2], XiParaF16_d[3], XiParaF16_d[4], XiParaF16_d[5], XiParaF16_d[6],
    XiParaF17_d[0], XiParaF17_d[1], XiParaF17_d[2], XiParaF17_d[3], XiParaF17_d[4], XiParaF17_d[5], XiParaF17_d[6],
    XiParaF18_d[0], XiParaF18_d[1], XiParaF18_d[2], XiParaF18_d[3], XiParaF18_d[4], XiParaF18_d[5], XiParaF18_d[6],
    XBFLParaF37_d[0]  ,   XBFLParaF37_d[1]  ,   XBFLParaF37_d[2]  ,   XBFLParaF37_d[3]  ,   XBFLParaF37_d[4]  ,   XBFLParaF37_d[5]  ,   XBFLParaF37_d[6]  ,
    XBFLParaF38_d[0]  ,   XBFLParaF38_d[1]  ,   XBFLParaF38_d[2]  ,   XBFLParaF38_d[3]  ,   XBFLParaF38_d[4]  ,   XBFLParaF38_d[5]  ,   XBFLParaF38_d[6]  ,
    YBFLParaF378_d[0] ,   YBFLParaF378_d[1] ,   YBFLParaF378_d[2] ,   YBFLParaF378_d[3] ,   YBFLParaF378_d[4] ,   YBFLParaF378_d[5] ,   YBFLParaF378_d[6] ,
    XiBFLParaF378_d[0],   XiBFLParaF378_d[1],   XiBFLParaF378_d[2],   XiBFLParaF378_d[3],   XiBFLParaF378_d[4],   XiBFLParaF378_d[5],   XiBFLParaF378_d[6],
    XBFLParaF49_d[0]   ,  XBFLParaF49_d[1]   ,  XBFLParaF49_d[2]   ,  XBFLParaF49_d[3]   ,  XBFLParaF49_d[4]   ,  XBFLParaF49_d[5]   ,  XBFLParaF49_d[6]   ,
    XBFLParaF410_d[0]  ,  XBFLParaF410_d[1]  ,  XBFLParaF410_d[2]  ,  XBFLParaF410_d[3]  ,  XBFLParaF410_d[4]  ,  XBFLParaF410_d[5]  ,  XBFLParaF410_d[6]  ,
    YBFLParaF4910_d[0] ,  YBFLParaF4910_d[1] ,  YBFLParaF4910_d[2] ,  YBFLParaF4910_d[3] ,  YBFLParaF4910_d[4] ,  YBFLParaF4910_d[5] ,  YBFLParaF4910_d[6] ,
    XiBFLParaF4910_d[0],  XiBFLParaF4910_d[1],  XiBFLParaF4910_d[2],  XiBFLParaF4910_d[3],  XiBFLParaF4910_d[4],  XiBFLParaF4910_d[5],  XiBFLParaF4910_d[6],
    YBFLParaF15_d[0] ,    YBFLParaF15_d[1] ,    YBFLParaF15_d[2] ,    YBFLParaF15_d[3] ,    YBFLParaF15_d[4] ,    YBFLParaF15_d[5] ,    YBFLParaF15_d[6] ,
    XiBFLParaF15_d[0],    XiBFLParaF15_d[1],    XiBFLParaF15_d[2],    XiBFLParaF15_d[3],    XiBFLParaF15_d[4],    XiBFLParaF15_d[5],    XiBFLParaF15_d[6],
    YBFLParaF16_d[0] ,    YBFLParaF16_d[1] ,    YBFLParaF16_d[2] ,    YBFLParaF16_d[3] ,    YBFLParaF16_d[4] ,    YBFLParaF16_d[5] ,    YBFLParaF16_d[6] ,
    XiBFLParaF16_d[0],    XiBFLParaF16_d[1],    XiBFLParaF16_d[2],    XiBFLParaF16_d[3],    XiBFLParaF16_d[4],    XiBFLParaF16_d[5],    XiBFLParaF16_d[6],
    BFLReqF3_d,     BFLReqF4_d,     BFLReqF15_d,    BFLReqF16_d,
    u,      v,      w,      rho_d,  Force_d,    NYD6-7, rho_modify_d, Q3_d, Q4_d, Q15_d, Q16_d
    );

    dim3 griddim_Ubulk(  NX6/NT+1, 1, NZ6);
    dim3 blockdim_Ubulk( NT, 1, 1);

    AccumulateUbulk<<<griddim_Ubulk, blockdim_Ubulk, 0, stream1>>>(
    Ub_avg_d, v, x_d, z_d
    );

    stream_collide<<<griddim, blockdim, 0, stream0>>>(
    f_old[0],   f_old[1],   f_old[2],   f_old[3],   f_old[4],   f_old[5],   f_old[6],   f_old[7],   f_old[8],   f_old[9],   f_old[10],  f_old[11],  f_old[12],  f_old[13],  f_old[14],  f_old[15],  f_old[16],  f_old[17],  f_old[18],
    f_new[0],   f_new[1],   f_new[2],   f_new[3],   f_new[4],   f_new[5],   f_new[6],   f_new[7],   f_new[8],   f_new[9],   f_new[10],  f_new[11],  f_new[12],  f_new[13],  f_new[14],  f_new[15],  f_new[16],  f_new[17],  f_new[18],
    XPara0_d[0],  XPara0_d[1],  XPara0_d[2],    XPara0_d[3],    XPara0_d[4],    XPara0_d[5],    XPara0_d[6],    XPara2_d[0],    XPara2_d[1],    XPara2_d[2],    XPara2_d[3],    XPara2_d[4],    XPara2_d[5],    XPara2_d[6],
    YPara0_d[0],  YPara0_d[1],  YPara0_d[2],    YPara0_d[3],    YPara0_d[4],    YPara0_d[5],    YPara0_d[6],    YPara2_d[0],    YPara2_d[1],    YPara2_d[2],    YPara2_d[3],    YPara2_d[4],    YPara2_d[5],    YPara2_d[6],
    XiParaF3_d[0],  XiParaF3_d[1],  XiParaF3_d[2],  XiParaF3_d[3],  XiParaF3_d[4],  XiParaF3_d[5],  XiParaF3_d[6],
    XiParaF4_d[0],  XiParaF4_d[1],  XiParaF4_d[2],  XiParaF4_d[3],  XiParaF4_d[4],  XiParaF4_d[5],  XiParaF4_d[6],
    XiParaF5_d[0],  XiParaF5_d[1],  XiParaF5_d[2],  XiParaF5_d[3],  XiParaF5_d[4],  XiParaF5_d[5],  XiParaF5_d[6],
    XiParaF6_d[0],  XiParaF6_d[1],  XiParaF6_d[2],  XiParaF6_d[3],  XiParaF6_d[4],  XiParaF6_d[5],  XiParaF6_d[6],
    XiParaF15_d[0], XiParaF15_d[1], XiParaF15_d[2], XiParaF15_d[3], XiParaF15_d[4], XiParaF15_d[5], XiParaF15_d[6],
    XiParaF16_d[0], XiParaF16_d[1], XiParaF16_d[2], XiParaF16_d[3], XiParaF16_d[4], XiParaF16_d[5], XiParaF16_d[6],
    XiParaF17_d[0], XiParaF17_d[1], XiParaF17_d[2], XiParaF17_d[3], XiParaF17_d[4], XiParaF17_d[5], XiParaF17_d[6],
    XiParaF18_d[0], XiParaF18_d[1], XiParaF18_d[2], XiParaF18_d[3], XiParaF18_d[4], XiParaF18_d[5], XiParaF18_d[6],
    XBFLParaF37_d[0]  ,   XBFLParaF37_d[1]  ,   XBFLParaF37_d[2]  ,   XBFLParaF37_d[3]  ,   XBFLParaF37_d[4]  ,   XBFLParaF37_d[5]  ,   XBFLParaF37_d[6]  ,
    XBFLParaF38_d[0]  ,   XBFLParaF38_d[1]  ,   XBFLParaF38_d[2]  ,   XBFLParaF38_d[3]  ,   XBFLParaF38_d[4]  ,   XBFLParaF38_d[5]  ,   XBFLParaF38_d[6]  ,
    YBFLParaF378_d[0] ,   YBFLParaF378_d[1] ,   YBFLParaF378_d[2] ,   YBFLParaF378_d[3] ,   YBFLParaF378_d[4] ,   YBFLParaF378_d[5] ,   YBFLParaF378_d[6] ,
    XiBFLParaF378_d[0],   XiBFLParaF378_d[1],   XiBFLParaF378_d[2],   XiBFLParaF378_d[3],   XiBFLParaF378_d[4],   XiBFLParaF378_d[5],   XiBFLParaF378_d[6],
    XBFLParaF49_d[0]   ,  XBFLParaF49_d[1]   ,  XBFLParaF49_d[2]   ,  XBFLParaF49_d[3]   ,  XBFLParaF49_d[4]   ,  XBFLParaF49_d[5]   ,  XBFLParaF49_d[6]   ,
    XBFLParaF410_d[0]  ,  XBFLParaF410_d[1]  ,  XBFLParaF410_d[2]  ,  XBFLParaF410_d[3]  ,  XBFLParaF410_d[4]  ,  XBFLParaF410_d[5]  ,  XBFLParaF410_d[6]  ,
    YBFLParaF4910_d[0] ,  YBFLParaF4910_d[1] ,  YBFLParaF4910_d[2] ,  YBFLParaF4910_d[3] ,  YBFLParaF4910_d[4] ,  YBFLParaF4910_d[5] ,  YBFLParaF4910_d[6] ,
    XiBFLParaF4910_d[0],  XiBFLParaF4910_d[1],  XiBFLParaF4910_d[2],  XiBFLParaF4910_d[3],  XiBFLParaF4910_d[4],  XiBFLParaF4910_d[5],  XiBFLParaF4910_d[6],
    YBFLParaF15_d[0] ,    YBFLParaF15_d[1] ,    YBFLParaF15_d[2] ,    YBFLParaF15_d[3] ,    YBFLParaF15_d[4] ,    YBFLParaF15_d[5] ,    YBFLParaF15_d[6] ,
    XiBFLParaF15_d[0],    XiBFLParaF15_d[1],    XiBFLParaF15_d[2],    XiBFLParaF15_d[3],    XiBFLParaF15_d[4],    XiBFLParaF15_d[5],    XiBFLParaF15_d[6],
    YBFLParaF16_d[0] ,    YBFLParaF16_d[1] ,    YBFLParaF16_d[2] ,    YBFLParaF16_d[3] ,    YBFLParaF16_d[4] ,    YBFLParaF16_d[5] ,    YBFLParaF16_d[6] ,
    XiBFLParaF16_d[0],    XiBFLParaF16_d[1],    XiBFLParaF16_d[2],    XiBFLParaF16_d[3],    XiBFLParaF16_d[4],    XiBFLParaF16_d[5],    XiBFLParaF16_d[6],
    BFLReqF3_d,     BFLReqF4_d,     BFLReqF15_d,    BFLReqF16_d,
    u,      v,      w,      rho_d,  Force_d, rho_modify_d, Q3_d, Q4_d, Q15_d, Q16_d
    );

    /* periodicUD<<<griddimUD, blockdimUD>>>(
        f_old[0] ,f_old[1] ,f_old[2] ,f_old[3] ,f_old[4] ,f_old[5] ,f_old[6] ,f_old[7] ,f_old[8] ,f_old[9] ,f_old[10] ,f_old[11] ,f_old[12] ,f_old[13] ,f_old[14] ,f_old[15] ,f_old[16] ,f_old[17] ,f_old[18],
        f_new[0] ,f_new[1] ,f_new[2] ,f_new[3] ,f_new[4] ,f_new[5] ,f_new[6] ,f_new[7] ,f_new[8] ,f_new[9] ,f_new[10] ,f_new[11] ,f_new[12] ,f_new[13] ,f_new[14] ,f_new[15] ,f_new[16] ,f_new[17] ,f_new[18],
        y_d, x_d, z_d, u, v, w, rho_d
    ); */
    /* periodicNML<<<griddimNML, blockdimNML, 0, stream0>>>(
        f_old[0] ,f_old[1] ,f_old[2] ,f_old[3] ,f_old[4] ,f_old[5] ,f_old[6] ,f_old[7] ,f_old[8] ,f_old[9] ,f_old[10] ,f_old[11] ,f_old[12] ,f_old[13] ,f_old[14] ,f_old[15] ,f_old[16] ,f_old[17] ,f_old[18],
        f_new[0] ,f_new[1] ,f_new[2] ,f_new[3] ,f_new[4] ,f_new[5] ,f_new[6] ,f_new[7] ,f_new[8] ,f_new[9] ,f_new[10] ,f_new[11] ,f_new[12] ,f_new[13] ,f_new[14] ,f_new[15] ,f_new[16] ,f_new[17] ,f_new[18],
        y_d, x_d, z_d, u, v, w, rho_d
    ); */

    CHECK_CUDA( cudaStreamSynchronize(stream1) );

    ISend_LtRtBdry( f_new, iToLeft,    l_nbr, itag_f4, 0, 10,   4,9,10,16,18, 3,7,8,15,17  );
    IRecv_LtRtBdry( f_new, iFromRight, r_nbr, itag_f4, 1, 10,   4,9,10,16,18, 3,7,8,15,17  );
    ISend_LtRtBdry( f_new, iToRight,   r_nbr, itag_f3, 2, 10,   4,9,10,16,18, 3,7,8,15,17  );
    IRecv_LtRtBdry( f_new, iFromLeft,  l_nbr, itag_f3, 3, 10,   4,9,10,16,18, 3,7,8,15,17  );

    for( int i = 0;  i < 10; i++ ){
        CHECK_MPI( MPI_Waitall(4, request[i], status[i]) );
    }
    for( int i = 19; i < 23; i++ ){
        CHECK_MPI( MPI_Waitall(4, request[i], status[i]) );
    }

    CHECK_CUDA( cudaStreamSynchronize(stream0) );

    periodicSW<<<griddimSW, blockdimSW, 0, stream0>>>(
        f_old[0] ,f_old[1] ,f_old[2] ,f_old[3] ,f_old[4] ,f_old[5] ,f_old[6] ,f_old[7] ,f_old[8] ,f_old[9] ,f_old[10] ,f_old[11] ,f_old[12] ,f_old[13] ,f_old[14] ,f_old[15] ,f_old[16] ,f_old[17] ,f_old[18],
        f_new[0] ,f_new[1] ,f_new[2] ,f_new[3] ,f_new[4] ,f_new[5] ,f_new[6] ,f_new[7] ,f_new[8] ,f_new[9] ,f_new[10] ,f_new[11] ,f_new[12] ,f_new[13] ,f_new[14] ,f_new[15] ,f_new[16] ,f_new[17] ,f_new[18],
        y_d, x_d, z_d, u, v, w, rho_d
    );
}

void Launch_ModifyForcingTerm()
{
    const size_t nBytes = NX6 * NZ6 * sizeof(double);
    CHECK_CUDA( cudaMemcpy(Ub_avg_h, Ub_avg_d, nBytes, cudaMemcpyDeviceToHost) );
    
    double Ub_avg = 0.0;
    for( int k = 3; k < NZ6-3; k++ ){
    for( int i = 3; i < NX6-4; i++ ){
        Ub_avg = Ub_avg + Ub_avg_h[k*NX6+i];
        Ub_avg_h[k*NX6+i] = 0.0;
    }}
    Ub_avg = Ub_avg / (double)(LX*(LZ-1.0))/NDTFRC;

    CHECK_CUDA( cudaMemcpy(Ub_avg_d, Ub_avg_h, nBytes, cudaMemcpyHostToDevice) );

    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    
    double beta = max(0.001, 3.0/(double)Re);
    Force_h[0] = Force_h[0] + beta*(Uref - Ub_avg)*Uref/3.036;

    double force_avg = 0.0;
    CHECK_MPI( MPI_Reduce( (void*)Force_h, (void*)&force_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
    CHECK_MPI( MPI_Barrier( MPI_COMM_WORLD ) );

    if( myid == 0 ){
        force_avg = force_avg / (double)jp;
        Force_h[0] = force_avg;
    }

    CHECK_MPI( MPI_Bcast( (void*)Force_h, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD ) );
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );

    printf("Ub_avg = %lf\t Uref = %lf\t Force = %.5lE\n", Ub_avg, Uref , Force_h[0] );

    CHECK_CUDA( cudaMemcpy(Force_d, Force_h, sizeof(double), cudaMemcpyHostToDevice) );
    
    CHECK_CUDA( cudaDeviceSynchronize() );
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    
}
#endif
