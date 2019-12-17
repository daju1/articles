#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI
#include "tzap.h"
#include <assert.h>
//#define USE_DEBUG
#include "dbg_info.h"
#include "integrate.h"
#include "calc.h"

extern velocity g_c;
extern timevalue * v_t;
extern int v_n_t; // итератор полноты заполения двумерных массивов по координате времени

int do_v1_calc(charge q, mass m_pos, mass m_neg, coordinate r0_pos, coordinate r0_neg, velocity v0_pos, velocity v0_neg, power pw_pos, power pw_neg, timespan t_a0)
{
	printf("do_v1_calc\n");
	printf("q = %0.20Le\n", q);
	printf("t_a0 = %0.20Le\n", t_a0);
	printf("m_pos = %0.20Le m_neg = %0.20Le\n", m_pos, m_neg);
	printf("r0_pos = %0.20Le r0_neg = %0.20Le\n", r0_pos, r0_neg);
	printf("v0_pos = %0.20Le v0_neg = %0.20Le\n", v0_pos, v0_neg);
	printf("pw_pos = %0.20Le pw_neg = %0.20Le\n", pw_pos, pw_neg);
#ifdef USE_NORM
	// q = q_calc * k_q
	// m = m_calc * k_m
	// r = r_calc * k_r
	// t = t_calc * k_t
	// v = v_valc * k_v
	// c = c_calc * k_v -> LIGHT_VELONCITY = g_c * k_v
	// a = a_calc * k_a
	// E = E_calc * k_E

    // SI :
	// k_E = (LIGHT_VELONCITY * LIGHT_VELONCITY / 10000000.0) * k_q / (k_r * k_r)
	// k_E = SI_multiplier_E * k_q / (k_r * k_r)
	// sgs :
	// k_E = k_q / (k_r * k_r)

	// k_a = k_E * k_q  / k_m
	// k_t = k_v / k_a
	// k_v = k_t * k_a
	// k_a * k_r = k_v * k_v


	long double k_q = 1.0;
	long double k_m = 0.0001;
	long double k_r = 0.001;
	long double k_E = multiplier_E * k_q / (k_r * k_r);
	long double k_a = k_E * k_q  / k_m;
	long double k_v = sqrt(k_a * k_r);
	long double k_t = k_v  / k_a;
	g_c = LIGHT_VELONCITY / k_v;

	printf("k_q = %0.20Lf\n", k_q);
	printf("k_m = %0.20Lf\n", k_m);
	printf("k_r = %0.20Lf\n", k_r);
	printf("m_E = %0.20Lf\n", multiplier_E);
	printf("k_E = %0.20Lf\n", k_E);
	printf("k_a = %0.20Lf\n", k_a);
	printf("k_v = %0.20Lf\n", k_v);
	printf("g_c = %0.20Lf\n", g_c);
	printf("defc= %0.20f\n", LIGHT_VELONCITY);
	printf("k_t = %0.20Lf\n", k_t);

	charge q_calc = q / k_q;
	power pw_pos_calc = a0_pos / k_a;
	power pw_neg_calc = a0_neg / k_a;
	timespan t_a0_calc = t_a0 / k_t;
	mass m_pos_calc = m_pos / k_m, m_neg_calc = m_neg / k_m;
	coordinate r0_pos_calc = r0_pos / k_r, r0_neg_calc = r0_neg / k_r;
	velocity v0_pos_calc = v0_pos / k_v, v0_neg_calc = v0_neg / k_v;
#else
	charge q_calc = q;
	power pw_pos_calc = pw_pos;
	power pw_neg_calc = pw_neg;
	timespan t_a0_calc = t_a0;
	mass m_pos_calc = m_pos, m_neg_calc = m_neg;
	coordinate r0_pos_calc = r0_pos, r0_neg_calc = r0_neg;
	velocity v0_pos_calc = v0_pos, v0_neg_calc = v0_neg;
#endif

	coordinate r_finish_calc = r0_pos_calc * 10;
	distance dr_calc = r_finish_calc / 100;
	coordinate r_min_pos_calc = dr_calc;
	coordinate r_min_neg_calc = dr_calc;

	set_dr(dr_calc);
	set_r_finish(r_finish_calc);
	return do_v1_calc_priv(q_calc, m_pos_calc, m_neg_calc, r0_pos_calc, r0_neg_calc, v0_pos_calc, v0_neg_calc, pw_pos_calc, pw_neg_calc, t_a0_calc , r_min_pos_calc, r_min_neg_calc);
}

int calc_E(charge q, timevalue t, coordinate R0,
			coordinate r0_pos, coordinate r0_neg,
			velocity v0_pos, velocity v0_neg,
			acceleration a0_pos, acceleration a0_neg,
			coordinate r_min_pos, coordinate r_min_neg,
			field * E_minus_grad_phi_R0_pos, field * E_minus_1_c_dA_dt_R0_pos,
			field * E_minus_grad_phi_R0_neg, field * E_minus_1_c_dA_dt_R0_neg,
			field * E1, field * E2, field * E)
{
	int error = 0, err;
	potential phi_lw_pos, A_lw_pos;
	potential phi_lw_neg, A_lw_neg;

	printf("\n");
	error = integral_phi_and_E(+q, t, R0, r0_pos, v0_pos, a0_pos, E_minus_grad_phi_R0_pos, E_minus_1_c_dA_dt_R0_pos, r_min_pos, &phi_lw_pos, &A_lw_pos);
	printf("pos err = %d phi_lw = %Lf E1=%Lf E2 = %0.20Lf\n", error, phi_lw_pos, *E_minus_grad_phi_R0_pos, *E_minus_1_c_dA_dt_R0_pos);

	printf("\n");
	error = integral_phi_and_E(-q, t, R0, r0_neg, v0_neg, a0_neg, E_minus_grad_phi_R0_neg, E_minus_1_c_dA_dt_R0_neg, r_min_neg, &phi_lw_neg, &A_lw_neg);
	printf("neg err = %d phi_lw = %Lf E1=%Lf E2 = %0.20Lf\n", error, phi_lw_neg, *E_minus_grad_phi_R0_neg, *E_minus_1_c_dA_dt_R0_neg);

	*E1 = *E_minus_grad_phi_R0_pos + *E_minus_grad_phi_R0_neg;
	*E2 = *E_minus_1_c_dA_dt_R0_pos + *E_minus_1_c_dA_dt_R0_neg;

	*E = *E1 + *E2;
}


int do_v1_calc_priv(charge q, mass m_pos, mass m_neg, coordinate r0_pos, coordinate r0_neg, velocity v0_pos, velocity v0_neg, acceleration a0_pos, acceleration a0_neg, timespan t_a0, coordinate r_min_pos, coordinate r_min_neg)
{
	printf("sizeof(double) %ld\n", sizeof(double));
	printf("sizeof(long double) %ld\n", sizeof(long double));

	printf ("nr = %d ", get_nr());
	printf ("nt = %d\n", get_nt());
	printf ("dr = %0.20Lf ", get_dr());

	printf("do_v1_calc_priv\n");
	printf("q = %0.20Le\n", q);
	printf("t_a0 = %0.20Le\n", t_a0);
	printf("m_pos = %0.20Le m_neg = %0.20Le\n", m_pos, m_neg);
	printf("r0_pos = %0.20Le r0_neg = %0.20Le\n", r0_pos, r0_neg);
	printf("v0_pos = %0.20Le v0_neg = %0.20Le\n", v0_pos, v0_neg);
	printf("a0_pos = %0.20Le a0_neg = %0.20Le\n", a0_pos, a0_neg);

#ifdef ALGORITHM_VERSION_0
	init_array_0();
	timespan dt = 1e-3;
#endif

#ifdef ALGORITHM_VERSION_1
	init_array_1(0.0, v0_pos, r0_pos, 0.0, v0_neg, r0_neg);
	timespan dt = 1e-18;
#endif
	v_n_t = 0;

	printf("dt = %Le\n", dt);
	while (v_n_t < get_nt())
	{
		int error = 0, err;
		//potential phi_lw_pos;
		//potential phi_lw_neg;
		field E_minus_grad_phi_R0_pos, E_minus_1_c_dA_dt_R0_pos;
		field E_minus_grad_phi_R0_neg, E_minus_1_c_dA_dt_R0_neg;
		coordinate r_pos;
		coordinate r_neg;
		field E1, E2, E;

		timevalue t = v_t[v_n_t];
#ifdef ALGORITHM_VERSION_0
		printf("t = %Le v_n_t = %d\n", t, v_n_t);
		for (int v_n_r = 0; v_n_r < get_nr(); ++v_n_r)
		{
			coordinate R0 = v_n_r * get_dr();
			printf("R0 = %Lf v_n_r = %d dr = %Lf\n", R0, v_n_r, get_dr());

			calc_E(q, t, R0,
				r0_pos, r0_neg,
				v0_pos, v0_neg,
				a0_pos, a0_neg,
				r_min_pos, r_min_neg,
				&E_minus_grad_phi_R0_pos, &E_minus_1_c_dA_dt_R0_pos,
				&E_minus_grad_phi_R0_neg, &E_minus_1_c_dA_dt_R0_neg,
				&E1, &E2, &E);
			if (fabs(E) > 1e-20){
#if 1
				printf(
					"R0 = %Lf t = %Lf "
					"E = % 0.20Lf E1 = % 0.20Lf E2 = % 0.20Lf "
					//"phi_lw_pos = %Lf phi_lw_neg = %Lf"
                    "\n"
					, R0, t
					, E
					, E1, E2
					//, phi_lw_pos, phi_lw_neg
                    );
#endif
#if 1
				printf(
					"R0 = %0.10Lf t = %0.10Lf "
					"E1_pos % 0.20Lf "
					"E1_neg % 0.20Lf "
					//"phi_lw_pos = %Lf phi_lw_neg = %Lf "
					//"E2_pos % 0.20f "
					//"E2_neg % 0.20f "
					"\n"
					, R0, t
					, E_minus_grad_phi_R0_pos
					, E_minus_grad_phi_R0_neg
					//, E_minus_1_c_dA_dt_R0_pos
					//, E_minus_1_c_dA_dt_R0_neg
					//, phi_lw_pos, phi_lw_neg
					);
#endif
			}
#ifdef ALGORITHM_VERSION_1
			set_E_ex_1(v_n_t, v_n_r, E);
#endif
		}
		printf("\n");
		#endif

#ifdef ALGORITHM_VERSION_0
		acceleration a_pos = get_a(t, a0_pos);
		acceleration a_neg = get_a(t, a0_neg);

		velocity v_pos = get_v(t, v0_pos, a0_pos);
		velocity v_neg = get_v(t, v0_neg, a0_neg);

		distance s_pos = get_s(t, v0_pos, a0_pos);
		distance s_neg = get_s(t, v0_neg, a0_neg);
		distance ds = s_neg - s_pos;

		err = get_r(+q, t, r0_pos, v0_pos, a0_pos, r_min_pos, &r_pos);
		err = get_r(-q, t, r0_neg, v0_neg, a0_neg, r_min_neg, &r_neg);
#endif

#ifdef ALGORITHM_VERSION_1
		acceleration a_pos = get_a_ex1(t, +q);
		acceleration a_neg = get_a_ex1(t, -q);

		velocity v_pos = get_v_ex1(t, v0_pos, +q);
		velocity v_neg = get_v_ex1(t, v0_neg, -q);

		distance s_pos = get_s_ex1(t, v0_pos, +q);
		distance s_neg = get_s_ex1(t, v0_neg, -q);
		distance ds = s_neg - s_pos;

		err = get_r_ex1(+q, t, r0_pos, v0_pos, r_min_pos, &r_pos, 1);
		err = get_r_ex1(-q, t, r0_neg, v0_neg, r_min_neg, &r_neg, 1);
#endif
		printf("t = %Le\n"
			"a_pos = % 0.20Le, a_neg = % 0.20Le\n"
			"v_pos = % 0.20Le, v_neg = % 0.20Le\n"
			"vcpos = % 0.20Le, vcneg = % 0.20Le\n"
			"s_pos = % 0.20Le, s_neg = % 0.20Le s_neg - s_pos = % 0.20Le  \n"
			"r_pos = % 0.20Le, r_neg = % 0.20Le r_neg - r_pos = % 0.20Le \n"
			, t
			, a_pos, a_neg
			, v_pos, v_neg
			, v_pos / g_c, v_neg / g_c
			, s_pos, s_neg, (s_neg - s_pos)
			, r_pos, r_neg, (r_neg - r_pos)
			);


		//
		field E_pos, E_neg; // электрическое поле в облаасти нахождения положительной и отрицательной обкладки
		field E1_pos, E1_neg; // электрическое поле в облаасти нахождения положительной и отрицательной обкладки
		field E2_pos, E2_neg; // электрическое поле в облаасти нахождения положительной и отрицательной обкладки

		printf("\nCalc E on positive side");

		calc_E(q, t,
			r_pos, // координата наблюдения совпадает с координатой обкладки
			r0_pos, r0_neg,
			v0_pos, v0_neg,
			a0_pos, a0_neg,
			r_min_pos, r_min_neg,
			&E_minus_grad_phi_R0_pos, &E_minus_1_c_dA_dt_R0_pos,
			&E_minus_grad_phi_R0_neg, &E_minus_1_c_dA_dt_R0_neg,
			&E1_pos, &E2_pos, &E_pos);

		printf(
			"R0_pos = %0.10Lf t = %0.10Lf "
			"E1_pos %Le "
			"E1_neg %Le "
			"E1 %Le "
			"E2_pos %Le "
			"E2_neg %Le "
			"E2 %0.20Le "
			"\n"
			, r_pos, t
			, E_minus_grad_phi_R0_pos
			, E_minus_grad_phi_R0_neg
			, E_minus_grad_phi_R0_pos + E_minus_grad_phi_R0_neg

			, E_minus_1_c_dA_dt_R0_pos
			, E_minus_1_c_dA_dt_R0_neg
			, E_minus_1_c_dA_dt_R0_pos + E_minus_1_c_dA_dt_R0_neg
			);

		printf("\nCalc E on negative side");

		calc_E(q, t,
			r_neg, // координата наблюдения совпадает с координатой обкладки
			r0_pos, r0_neg,
			v0_pos, v0_neg,
			a0_pos, a0_neg,
			r_min_pos, r_min_neg,
			&E_minus_grad_phi_R0_pos, &E_minus_1_c_dA_dt_R0_pos, 
			&E_minus_grad_phi_R0_neg, &E_minus_1_c_dA_dt_R0_neg,
			&E1_neg, &E2_neg, &E_neg);

		printf(
			"R0_neg = %0.10Lf t = %0.10Lf "
			"E1_pos %Le "
			"E1_neg %Le "
			"E1 %Le "
			"E2_pos %Le "
			"E2_neg %Le "
			"E2 %0.20Le "
			"\n"
			, r_pos, t
			, E_minus_grad_phi_R0_pos
			, E_minus_grad_phi_R0_neg
			, E_minus_grad_phi_R0_pos + E_minus_grad_phi_R0_neg

			, E_minus_1_c_dA_dt_R0_pos
			, E_minus_1_c_dA_dt_R0_neg
			, E_minus_1_c_dA_dt_R0_pos + E_minus_1_c_dA_dt_R0_neg
			);

		printf(
			"E_pos = % 0.20Le, E_neg = % 0.20Le (E_pos - E_neg) = % 0.20Le\n"
			, E_pos, E_neg, (E_pos - E_neg)
			);

		if (v_n_t < get_nt() - 1)
		{
			v_t[v_n_t + 1] = v_t[v_n_t] + dt;
			timevalue t1 = v_t[v_n_t + 1];

			printf("t = %Le t1 = %Le v_n_t = %d dt = %Le\n", t, t1, v_n_t, dt);
#ifdef ALGORITHM_VERSION_1
			if (0 != set_a_ex1(t1, dt, r_neg, pw_neg, t_a0, -q, m_neg, E_neg, E1_neg, E2_neg))
			{
				error += 1;
				int * p = 0;
				*p += 1;
			}
			if (0 != set_a_ex1(t1, dt, r_pos, pw_pos, t_a0, +q, m_pos, E_pos, E1_pos, E2_pos))
			{
				error += 1;
			}

			velocity v_pos, v_neg;
			if (0 != set_v_ex1(t1, v0_neg, t_a0, -q, m_neg, &v_neg))
			{
				error += 1;
			}
			if (0 != set_v_ex1(t1, v0_pos, t_a0, +q, m_pos, &v_pos))
			{
				error += 1;
			}

			printf(
				"v_pos = % 0.20Le, v_neg = % 0.20Le\n"
				, v_pos, v_neg
				);

			if (0 == error) {
				set_s_ex1(t1, r0_pos, v0_pos, r_min_pos, +q);
				set_s_ex1(t1, r0_neg, v0_neg, r_min_neg, -q);
			}
#endif
		}

		if (0 == error)
		{
			++v_n_t;
		}
		else
		{
			v_n_t -= 2;
			if (v_n_t < 0) {
				v_n_t = 0;
			}
			else if (0 == v_n_t) {
				printf("0 == v_n_t\n");
			}
			else {
				printf("v_n_t > 0\n");
			}
			dt /= 2.0;
			printf("dt = %Le\n", dt);
		}
	}
}



