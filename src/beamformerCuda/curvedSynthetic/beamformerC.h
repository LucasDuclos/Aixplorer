#include <stdio.h>
#include <math.h>
#include "beamformer.h"

float apod_c(int n, int apertureLength)
{
    return (0.54f - 0.46f * cos(2.0f * (float)(M_PI) * (float)(n+1) / ((float)apertureLength - 1.0f)));
}

void findSteeredCoordinates_c(float y,float x,float r0,float steer, float & r,float & theta)
{
    float delta ;

    delta = sqrtf((2.0f*r0*cosf(steer))*(2.0f*r0*cosf(steer)) - 4.0f*(r0*r0 - x*x -y*y)) ;
    r = (- (2.0f*r0*cosf(steer)) + delta)/2.0f ;
    theta = x*(r0+r*cosf(steer)) - r*sinf(steer)*y ;
    theta /= r0*r0 + 2.0f*r*r0*cosf(steer) + r*r ;
    theta = asinf(theta);
}