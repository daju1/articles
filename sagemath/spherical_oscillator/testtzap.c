#include <stdio.h>

#define pi 3.14

extern  double calc_tzap(double t, double R0, double r0, double a0, double theta);
double get_r(double t_zap, double r0, double a0);
extern double get_R(double R0, double r, double theta);

int main()
{
	double t = 5;
	double R0 = 2;
	double r0 = 1;
	double a0 = 1.0;
	double theta = pi/2;
	double t_zap = calc_tzap(t, R0, r0, a0, theta);
	printf("t_zap = %f\n", t_zap );
	
	double r_zap = get_r(t_zap, r0, a0);
	double R_zap = get_R(R0, r_zap, theta);
	double c = 1.0;
	printf("R_zap = %f c*(t-t_zap) = %f delta = %e\n", R_zap, c*(t-t_zap), R_zap - c*(t-t_zap));
	return 0;
}