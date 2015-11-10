
//#include <fortran.h>
#include <stdlib.h>
#include <stdio.h>

void vega_FEM_initiate_cpp(void);
void vega_interpolateIndexRatio_cpp(int *markerInterpolateIndex_,double *markerInterpolateRatio_);
void vega_deformation_cpp(double *markerPressure_,double *markerInterpolateVelocity_,double *bodyMotion_);
void vega_reNewBodyPosition_cpp(void);

void vega_fem_initiate_c_(void)
{
    printf("****** \n Initiate Vega FEM. \n");
    vega_FEM_initiate_cpp();
}

void vega_interpolateindexratio_c_(int markerInterpolateIndex_[],double markerInterpolateRatio_[])
{
    printf("****** \n Vega: get interpolate Index Ratio. \n");
    vega_interpolateIndexRatio_cpp(markerInterpolateIndex_,markerInterpolateRatio_);
}

void vega_deformation_c_(double *markerPressure_,double *markerInterpolateVelocity_,double *bodyMotion_)
{
    printf("****** \n Starting vega_deformation. \n");
    vega_deformation_cpp(markerPressure_,markerInterpolateVelocity_,bodyMotion_);
}

void vega_renewbodyposition_c_(void)
{
    printf("****** \n FSI Deformation converged. renew Body position. \n");
    vega_reNewBodyPosition_cpp();
}
