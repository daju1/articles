/*
	demo-c.c
		test program for the Cuba library
		last modified 13 Mar 15 th
*/

#include "calc_RO.h"

#define pi M_PI


// p = 1

// lim_2_Rpx_Fn_integrand
// (2*n_x^2 - 2*n_x*floor(n_x) - n_x)/sqrt(n_x^2 + n_y^2 + u^2)


static inline cubareal lim_2_Rpx_Fn_integrand (cubareal u, cubareal n_x, cubareal n_y)
{
//(2*Sq(n_x) - 2*n_x*floor(n_x) - n_x)/sqrt(Sq(n_x) + Sq(n_y) + Sq(u))
	return (2*Sq(n_x) - 2*n_x*floor(n_x) - n_x)/sqrt(Sq(n_x) + Sq(n_y) +Sq(u));
}


static inline cubareal Rpy_Rpx_integrand (cubareal u, cubareal n_x, cubareal n_y, cubareal k_m)
{
return
1/4*(32*Ei(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(5.0/2))*n_x*n_y/(Cb(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Ei(k_m)) - 16*Qu(pi)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))*n_x*n_y/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Qu(k_m)) - n_x*n_y/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))))*(2*n_x - 2*floor(n_x) - 1)*(2*n_y - 2*floor(n_y) - 1)
;
}


//result_force_Fn_integrand(u, n_x, n_y, a, k_m)
// 1/4*((32*pi^8*(n_x^2 + n_y^2 + u^2)^(5/2)*n_x*n_y/((pi^4*(n_x^2 + n_y^2 + u^2)^2/k_m^4 + 1)^3*k_m^8) - 16*pi^4*sqrt(n_x^2 + n_y^2 + u^2)*n_x*n_y/((pi^4*(n_x^2 + n_y^2 + u^2)^2/k_m^4 + 1)^2*k_m^4) - n_x*n_y/((pi^4*(n_x^2 + n_y^2 + u^2)^2/k_m^4 + 1)*(n_x^2 + n_y^2 + u^2)^(3/2)))*(2*n_y - 2*floor(n_y) - 1) - 16*pi^4*(n_x^2 + n_y^2 + u^2)^(3/2)*n_x/((pi^4*(n_x^2 + n_y^2 + u^2)^2/k_m^4 + 1)^2*k_m^4) + 4*n_x/((pi^4*(n_x^2 + n_y^2 + u^2)^2/k_m^4 + 1)*sqrt(n_x^2 + n_y^2 + u^2)))*(2*n_x - 2*floor(n_x) - 1)

static inline cubareal result_force_Fn_integrand (cubareal u, cubareal n_x, cubareal n_y, cubareal k_m)
{
	return
1/4*((32*Ei(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(5.0/2))*n_x*n_y/(Cb(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Ei(k_m)) - 16*Qu(pi)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))*n_x*n_y/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Qu(k_m)) - n_x*n_y/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))))*(2*n_y - 2*floor(n_y) - 1) - 16*Qu(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))*n_x/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Qu(k_m)) + 4*n_x/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))))*(2*n_x - 2*floor(n_x) - 1)
;}

// 2*Rpx_Fn_integrand(u, n_x, n_y, 1, k_m)
// -(2*n_x - 2*floor(n_x) - 1)*(4*pi^4*(n_x^2 + n_y^2 + u^2)^(3/2)*n_x/((pi^4*(n_x^2 + n_y^2 + u^2)^2/k_m^4 + 1)^2*k_m^4) - n_x/((pi^4*(n_x^2 + n_y^2 + u^2)^2/k_m^4 + 1)*sqrt(n_x^2 + n_y^2 + u^2)))

static inline cubareal two_Rpx_Fn_integrand (cubareal u, cubareal n_x, cubareal n_y, cubareal k_m)
{
//-(2*n_x - 2*floor(n_x) - 1)*(4*Qu(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))*n_x/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Qu(k_m)) - n_x/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))))
	return
-(2*n_x - 2*floor(n_x) - 1)*(4*Qu(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))*n_x/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Qu(k_m)) - n_x/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))))
;}


// p = 2

static inline cubareal F_Rp_integrand_p2 (cubareal u, cubareal n_x, cubareal n_y, cubareal k_m)
{
	return
1/144*((6144*Sq(Ei(pi))*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(9.0/2))*Sq(n_x)*Sq(n_y)/(Fi(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Sq(Ei(k_m))) - 6144*Tw(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(5.0/2))*Sq(n_x)*Sq(n_y)/(Qu(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Tw(k_m)) - 384*Tw(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(7.0/2))*Sq(n_x)/(Qu(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Tw(k_m)) - 384*Tw(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(7.0/2))*Sq(n_y)/(Qu(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Tw(k_m)) + 960*Ei(pi)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))*Sq(n_x)*Sq(n_y)/(Cb(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Ei(k_m)) + 288*Ei(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))*Sq(n_x)/(Cb(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Ei(k_m)) + 288*Ei(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))*Sq(n_y)/(Cb(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Ei(k_m)) + 32*Ei(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(5.0/2))/(Cb(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Ei(k_m)) - 12*Qu(pi)*Sq(n_x)/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))*Qu(k_m)) - 12*Qu(pi)*Sq(n_y)/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))*Qu(k_m)) - 15*Sq(n_x)*Sq(n_y)/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(7.0/2))) - 16*Qu(pi)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Qu(k_m)) + 3*Sq(n_x)/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(5.0/2))) + 3*Sq(n_y)/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(5.0/2))) - 1/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))))*(6*Sq(n_y - floor(n_y)) - 6*n_y + 6*floor(n_y) + 1) - 768*Ei(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(5.0/2))*Sq(n_x)/(Cb(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Ei(k_m)) + 384*Qu(pi)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))*Sq(n_x)/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Qu(k_m)) + 96*Qu(pi)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))/(Sq(Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*Qu(k_m)) + 24*Sq(n_x)/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*pow((Sq(n_x) + Sq(n_y) + Sq(u)),(3.0/2))) - 24/((Qu(pi)*Sq(Sq(n_x) + Sq(n_y) + Sq(u))/Qu(k_m) + 1)*sqrt(Sq(n_x) + Sq(n_y) + Sq(u))))*(6*Sq(n_x - floor(n_x)) - 6*n_x + 6*floor(n_x) + 1)

;}

#define UP_LIMIT 1000000000000000000.0

#include <float.h>
//#define UP_LIMIT FLT_MAX



int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

    // u, n_x, n_y

    #define u     xx[0]
    #define n_x   xx[1]
    #define n_y   xx[2]

    #define f ff[0]

    // f = lim_2_Rpx_Fn_integrand (u * 1000, n_x * UP_LIMIT, n_y * UP_LIMIT );
    // f = F_Rp_integrand (u * 1000, n_x * UP_LIMIT, n_y * UP_LIMIT, 10);
    f = Rpy_Rpx_integrand (u * 1000, n_x * UP_LIMIT, n_y * UP_LIMIT, FLT_MAX);
    // f = result_force_Fn_integrand_p2 (u * 1000, n_x * UP_LIMIT, n_y * UP_LIMIT, FLT_MAX);
    // f = two_Rpx_Fn_integrand (u * 1000, n_x * UP_LIMIT, n_y * UP_LIMIT, FLT_MAX);

    return 0;
}

