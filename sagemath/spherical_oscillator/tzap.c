#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>


 
#ifdef USE_DEBUG
#define DBG_INFO(fmt, args...) \
    do \
    { \
        printf(fmt, ## args); \
        /*printf("\n");*/ \
    } \
    while(0)
#else
#define DBG_INFO(fmt, args...)
#endif

#ifdef USE_DEBUG
#define DBG_ERROR(fmt, args...) \
    do \
    { \
        printf(fmt, ## args); \
        /*printf("\n");*/ \
    } \
    while(0)
#else
#define DBG_ERROR(fmt, args...)
#endif

double t_start = 0;
 
double get_a(double t_zap, double a0)
{
	if (t_zap < t_start)
		return 0;
	return a0;
}
 
double get_v(double t_zap, double a0)
{
	if (t_zap < t_start)
		return 0;
	return a0*t_zap;
}
 
double get_s(double t_zap, double a0)
{
	if (t_zap < t_start)
	{
		/*DBG_INFO("get_s t_zap=%f a0=%f returns 0\n", t_zap, a0);*/
		return 0;
	}
	/*DBG_INFO("get_s t_zap=%f a0=%f returns %f\n", t_zap, a0, a0*t_zap*t_zap/2);*/
	return a0*t_zap*t_zap/2;
 }
 
double get_r(double t_zap, double r0, double a0)
{
	return r0+get_s(t_zap, a0);
}

double get_R(double R0, double r, double theta)
{
	double R = sqrt(R0*R0 - 2*R0*r*cos(theta) + r*r);
	return R;
}
 
double calc_tzap(double t, double R0, double r0, double a0, double theta)
{
	double epsilon = 1.0e-15;
	double t1 = t;
	double t2 = t - 2 * epsilon;
	double r, R, R_pre = DBL_MAX;
	double dR, dR_pre = DBL_MAX;
	double R_tmp;
	double c = 1.0;
	int i = 0;
	double n = 0.9;


	DBG_INFO("epsilon=%e\n", epsilon);
	DBG_INFO("t1=%f\n", t1);
	DBG_INFO("t2=%f\n", t2);

	while (fabs(t1 - t2) > epsilon)
	{
		t1 = t2;
		r = get_r(t1, r0, a0);
		R = get_R(R0, r, theta);
		t2 = t -  R / c;
		dR = c*(t-t1) - R;

		DBG_INFO("t2=%f t1=%f t=%f R=%f dR=%e dR_pre=%e ", t2, t1, t, R, dR, dR_pre);
		
		if (fabs(dR) >= fabs(dR_pre))
		{
			int j = 0;

			do
			{
				R_tmp = R_pre + n * (R - R_pre);
				DBG_INFO("R_tmp = R_pre + n * (R - R_pre);= %f ", R_tmp);

				t2 = t -  R_tmp / c;
				dR = c*(t-t1) - R_tmp;
				DBG_INFO("t2 = %f ", t2);

				n *= 0.9;
				DBG_INFO("n = %f ", n);
			}
			while (fabs(dR) > fabs(dR_pre));
	
			R = R_tmp;
		}
		dR_pre = dR;
		R_pre = R;
		
		DBG_INFO("\n");
	 }
	 
    DBG_INFO("fabs(t1 - t2) = %f \n", fabs(t1 - t2)); 
	return t2;
 }
 
float calc_tzap_float(float t, float R0, float r0, float a0, float theta)
{
	return calc_tzap(t, R0, r0, a0, theta);
}
