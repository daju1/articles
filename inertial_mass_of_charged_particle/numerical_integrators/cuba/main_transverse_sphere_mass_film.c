/*
    main_transverse_sphere_mass_film.c
*/

//#define INCLUDE_CUBA_H
#include "calc_RO.h"
#include "integrator_enum.h"

#include <getopt.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata);

/*********************************************************************/

#define NCOMP 9

int integrate(
    double R0,       /* Радиус сферы */
    double rho0,     /* Радиус орбиты */
    double omega,    /* Угловая скорость */
    double c,        /* Скорость света */
    int use_delay,
    int use_lorentz_factor,
    int use_lorentz_general_factor,
    int use_fermi_factor_O,
    int use_fermi_factor,
    int use_fermi_general_factor,
    int use_fast_integrand,

    int use_fermi_general_factor_O,
    int use_fermi_film_factor,
    int use_fermi_film_instant,
    int fermi_film_K,
    int use_fermi_pde_factor,
    int pde_grid_n,
    double *pde_T_grid,
    int use_rot_to_lab,

    int seed,
    int maxeval,

    Integrator integrator,

    cubareal* integral, cubareal* error, cubareal* prob);

/*********************************************************************/
/*                 Разбор параметров командной строки                */
/*********************************************************************/

typedef struct {
    double R0, rho0, c;
    double omega;  int omega_set;
    double v_c;    int v_c_set;
    double mass;   /* m: omega^2 = 5 R0 c^2 / (24 m rho0^3); m=1 или 4/3 */
    int use_delay;
    int use_lorentz_factor;
    int use_lorentz_general_factor;
    int use_fermi_factor_O;
    int use_fermi_factor;
    int use_fermi_general_factor;
    int use_fast_integrand;
    int use_fermi_general_factor_O;
    int use_fermi_film_factor;
    int use_fermi_film_instant;
    int fermi_film_K;
    int use_fermi_pde_factor;
    int pde_grid_n;
    int use_rot_to_lab;

    int seed;
    int maxeval;

    Integrator integrator;
} Options;

enum {
    OPT_R0 = 1000, OPT_RHO0, OPT_C, OPT_OMEGA, OPT_VC, OPT_MASS, OPT_INTEGRATOR,
    OPT_DELAY, OPT_LORENTZ, OPT_LORENTZ_GENERAL,
    OPT_FERMI_O, OPT_FERMI, OPT_FERMI_GENERAL, OPT_FERMI_GENERAL_O,
    OPT_FAST, OPT_FILM, OPT_FILM_INSTANT, OPT_FILM_K, OPT_PDE, OPT_PDE_N,
    OPT_ROT, OPT_SEED, OPT_MAXEVAL, OPT_HELP
};

static const struct option long_opts[] = {
    {"R0",              required_argument, NULL, OPT_R0},
    {"rho0",            required_argument, NULL, OPT_RHO0},
    {"c",               required_argument, NULL, OPT_C},
    {"omega",           required_argument, NULL, OPT_OMEGA},
    {"v-c",             required_argument, NULL, OPT_VC},
    {"mass",            required_argument, NULL, OPT_MASS},
    {"integrator",      required_argument, NULL, OPT_INTEGRATOR},
    {"delay",           optional_argument, NULL, OPT_DELAY},
    {"lorentz",         optional_argument, NULL, OPT_LORENTZ},
    {"lorentz-general", optional_argument, NULL, OPT_LORENTZ_GENERAL},
    {"fermi-O",         optional_argument, NULL, OPT_FERMI_O},
    {"fermi",           optional_argument, NULL, OPT_FERMI},
    {"fermi-general",   optional_argument, NULL, OPT_FERMI_GENERAL},
    {"fermi-general-O", optional_argument, NULL, OPT_FERMI_GENERAL_O},
    {"fast",            optional_argument, NULL, OPT_FAST},
    {"film",            optional_argument, NULL, OPT_FILM},
    {"film-instant",    optional_argument, NULL, OPT_FILM_INSTANT},
    {"film-K",          required_argument, NULL, OPT_FILM_K},
    {"pde",             optional_argument, NULL, OPT_PDE},
    {"pde-n",           required_argument, NULL, OPT_PDE_N},
    {"rot",             optional_argument, NULL, OPT_ROT},
    {"seed",            required_argument, NULL, OPT_SEED},
    {"maxeval",         required_argument, NULL, OPT_MAXEVAL},
    {"help",            no_argument,       NULL, OPT_HELP},
    {NULL, 0, NULL, 0}
};

static void die(const char *fmt, ...)
{
    va_list ap;
    fprintf(stderr, "ERROR: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n(см. --help)\n");
    exit(2);
}

static void print_help(const char *prog)
{
    printf(
"Использование: %s [параметры]\n"
"\n"
"Геометрия и движение:\n"
"  --R0=X              радиус сферы                          (по умолч. 1)\n"
"  --rho0=X            радиус орбиты                         (по умолч. 100)\n"
"  --c=X               скорость света                        (по умолч. 1)\n"
"  --omega=X           угловая скорость (высший приоритет)\n"
"  --v-c=X             v/c центра сферы: omega = X*c/rho0\n"
"  --mass=M            m (число или 4/3): omega^2 = 5*R0*c^2/(24*M*rho0^3)\n"
"                      (по умолч. 1; M=4/3 даёт вашу закомментированную формулу)\n"
"  Приоритет omega: --omega > --v-c > --mass.\n"
"\n"
"Булевы параметры: --name включает, --name=0 выключает (через пробел нельзя).\n"
"  --delay             учёт запаздывания                     (по умолч. 1)\n"
"  --fast              быстрый интегранд (f_one = 0)         (по умолч. 1)\n"
"  --lorentz           лоренцево сжатие\n"
"  --lorentz-general   лоренцево сжатие, общий вариант\n"
"\n"
"Режимы поправки Ферми (не более одного):\n"
"  --fermi-O           нерелятивистский, точка O, без полусуммы\n"
"  --fermi             нерелятивистский, точка P', полусумма 1/2\n"
"  --fermi-general     релятивистский 4-векторный, точка P', 1/2\n"
"  --fermi-general-O   релятивистский 4-векторный, точка O, без 1/2\n"
"  --film              скользящая плёнка (только при --delay)\n"
"  --film-instant      плёнка берёт ускорения в один момент\n"
"  --film-K=N          число точек Монте-Карло для плёнки (0 = 20)\n"
"  --pde               PDE-версия, сетка T (только при --delay)\n"
"  --pde-n=N           узлов сетки по оси (0 = 24)\n"
"\n"
"  --integrator=NAME   vegas | suave | divonne | cuhre       (по умолч. vegas)\n"
"  --rot               повернуть поля E в лаб. оси x,y            (по умолч. 1)\n"
"  --seed=N            seed генератора Cuba                       (по умолч. 38304)\n"
"  --maxeval=N         максимум вычислений интегранда             (по умолч. 500000)\n"
"  -h, --help          эта справка\n", prog);
}

static int parse_flag(const char *arg, const char *name)
{
    if (!arg) return 1;
    if (!strcmp(arg, "1") || !strcasecmp(arg, "on") || !strcasecmp(arg, "true"))  return 1;
    if (!strcmp(arg, "0") || !strcasecmp(arg, "off") || !strcasecmp(arg, "false")) return 0;
    die("--%s: ожидалось 0 или 1, получено '%s'", name, arg);
    return 0;
}

static double parse_double(const char *arg, const char *name)
{
    char *end = NULL;
    double v = strtod(arg, &end);
    if (end == arg || *end != '\0') die("--%s: '%s' не число", name, arg);
    return v;
}

static int parse_int(const char *arg, const char *name)
{
    char *end = NULL;
    long v = strtol(arg, &end, 10);
    if (end == arg || *end != '\0') die("--%s: '%s' не целое число", name, arg);
    return (int)v;
}

/* число или дробь вида 4/3 */
static double parse_ratio(const char *arg, const char *name)
{
    if (strchr(arg, '/')) {
        double a, b;
        if (sscanf(arg, "%lf/%lf", &a, &b) != 2 || b == 0.0)
            die("--%s: '%s' не дробь a/b", name, arg);
        return a / b;
    }
    return parse_double(arg, name);
}

static Integrator parse_integrator(const char *s)
{
    if (!strcasecmp(s, "vegas"))   return VegasIntegrator;
    if (!strcasecmp(s, "suave"))   return SuaveIntegrator;
    if (!strcasecmp(s, "divonne")) return DivonneIntegrator;
    if (!strcasecmp(s, "cuhre"))   return CuhreIntegrator;
    die("--integrator: неизвестное имя '%s'", s);
    return VegasIntegrator;
}

static const char *integrator_name(Integrator i)
{
    if (i == VegasIntegrator)   return "vegas";
    if (i == SuaveIntegrator)   return "suave";
    if (i == DivonneIntegrator) return "divonne";
    if (i == CuhreIntegrator)   return "cuhre";
    return "?";
}

static void set_defaults(Options *o)
{
    memset(o, 0, sizeof *o);
    o->R0   = 1.0;
    o->rho0 = 100.0;
    o->c    = 1.0;
    o->mass = 1.0;
    o->use_delay          = 1;
    o->use_fast_integrand = 1;
    o->use_rot_to_lab     = 1;
    o->seed               = 38304;
    o->maxeval            = 500000;
    o->integrator         = VegasIntegrator;
}

static void parse_args(int argc, char **argv, Options *o)
{
    int ch;
    while ((ch = getopt_long(argc, argv, "h", long_opts, NULL)) != -1) {
        switch (ch) {
        case OPT_R0:      o->R0   = parse_double(optarg, "R0");   break;
        case OPT_RHO0:    o->rho0 = parse_double(optarg, "rho0"); break;
        case OPT_C:       o->c    = parse_double(optarg, "c");    break;
        case OPT_OMEGA:   o->omega = parse_double(optarg, "omega"); o->omega_set = 1; break;
        case OPT_VC:      o->v_c   = parse_double(optarg, "v-c");   o->v_c_set   = 1; break;
        case OPT_MASS:    o->mass  = parse_ratio(optarg, "mass");   break;
        case OPT_INTEGRATOR: o->integrator = parse_integrator(optarg); break;

        case OPT_DELAY:           o->use_delay                  = parse_flag(optarg, "delay");           break;
        case OPT_LORENTZ:         o->use_lorentz_factor         = parse_flag(optarg, "lorentz");         break;
        case OPT_LORENTZ_GENERAL: o->use_lorentz_general_factor = parse_flag(optarg, "lorentz-general"); break;
        case OPT_FERMI_O:         o->use_fermi_factor_O         = parse_flag(optarg, "fermi-O");         break;
        case OPT_FERMI:           o->use_fermi_factor           = parse_flag(optarg, "fermi");           break;
        case OPT_FERMI_GENERAL:   o->use_fermi_general_factor   = parse_flag(optarg, "fermi-general");   break;
        case OPT_FERMI_GENERAL_O: o->use_fermi_general_factor_O = parse_flag(optarg, "fermi-general-O"); break;
        case OPT_FAST:            o->use_fast_integrand         = parse_flag(optarg, "fast");            break;
        case OPT_FILM:            o->use_fermi_film_factor      = parse_flag(optarg, "film");            break;
        case OPT_FILM_INSTANT:    o->use_fermi_film_instant     = parse_flag(optarg, "film-instant");    break;
        case OPT_PDE:             o->use_fermi_pde_factor       = parse_flag(optarg, "pde");             break;
        case OPT_FILM_K:          o->fermi_film_K = parse_int(optarg, "film-K"); break;
        case OPT_PDE_N:           o->pde_grid_n   = parse_int(optarg, "pde-n");  break;
        case OPT_ROT:             o->use_rot_to_lab             = parse_flag(optarg, "rot");  break;
        case OPT_SEED:            o->seed    = parse_int(optarg, "seed");    break;
        case OPT_MAXEVAL:         o->maxeval = parse_int(optarg, "maxeval"); break;

        case 'h':
        case OPT_HELP:
            print_help(argv[0]);
            exit(0);
        default:
            exit(2);   /* getopt_long уже вывела сообщение */
        }
    }
    if (optind < argc) die("лишний аргумент '%s'", argv[optind]);
}

static void resolve_and_validate(Options *o)
{
    if (o->R0 <= 0 || o->rho0 <= 0 || o->c <= 0) die("R0, rho0 и c должны быть > 0");
    if (o->mass <= 0) die("--mass должно быть > 0");
    if (o->fermi_film_K < 0 || o->pde_grid_n < 0) die("--film-K и --pde-n не могут быть < 0");

    if (o->omega_set)
        ;  /* omega уже задана */
    else if (o->v_c_set)
        o->omega = o->v_c * o->c / o->rho0;
    else
        o->omega = sqrt(5.0 * o->R0 * o->c * o->c / (24.0 * o->mass * o->rho0 * o->rho0 * o->rho0));

    int n_fermi = o->use_fermi_factor_O + o->use_fermi_factor
                + o->use_fermi_general_factor + o->use_fermi_general_factor_O
                + o->use_fermi_film_factor + o->use_fermi_pde_factor;
    if (n_fermi > 1)
        die("включено %d режимов Ферми; в computing_electric_field сработает только один "
            "(приоритет: film/pde > general-O > general > fermi > fermi-O). Выберите один.", n_fermi);

    if ((o->use_fermi_film_factor || o->use_fermi_pde_factor) && !o->use_delay)
        die("--film и --pde работают только в ветке с запаздыванием: добавьте --delay");
    if (o->use_fermi_film_instant && !o->use_fermi_film_factor)
        die("--film-instant имеет смысл только вместе с --film");

    /* скорость самой быстрой точки сферы (без учёта сжатия) */
    double v_max = o->omega * (o->rho0 + o->R0) / o->c;
    if (v_max >= 1.0) {
        if (o->use_lorentz_factor || o->use_lorentz_general_factor)
            die("v_max/c = %.4f >= 1 при включённом лоренцевом факторе (gamma = NaN)", v_max);
        fprintf(stderr, "WARNING: v_max/c = %.4f >= 1: часть сферы движется быстрее света\n", v_max);
    }
    if (o->seed < 0)    die("--seed должен быть >= 0");
    if (o->maxeval < 1) die("--maxeval должен быть >= 1");
}

/*********************************************************************/

int main(int argc, char **argv)
{
    int comp, fail;
    cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

    Options o;
    set_defaults(&o);
    parse_args(argc, argv, &o);
    resolve_and_validate(&o);

    double v_c     = o.omega * o.rho0 / o.c;               /* v/c центра */
    double v_c_max = o.omega * (o.rho0 + o.R0) / o.c;      /* v/c края сферы */

    printf("PARAMS integrator=%s R0=%g rho0=%g c=%g omega=%.12g v_c=%.6g v_c_max=%.6g "
           "delay=%d lorentz=%d lorentz_general=%d fermi_O=%d fermi=%d fermi_general=%d "
           "fermi_general_O=%d fast=%d film=%d film_instant=%d film_K=%d pde=%d pde_n=%d rot=%d seed=%d maxeval=%d\n",
           integrator_name(o.integrator), o.R0, o.rho0, o.c, o.omega, v_c, v_c_max,
           o.use_delay, o.use_lorentz_factor, o.use_lorentz_general_factor,
           o.use_fermi_factor_O, o.use_fermi_factor, o.use_fermi_general_factor,
           o.use_fermi_general_factor_O, o.use_fast_integrand,
           o.use_fermi_film_factor, o.use_fermi_film_instant, o.fermi_film_K,
           o.use_fermi_pde_factor, o.pde_grid_n, o.use_rot_to_lab, o.seed, o.maxeval);
    fflush(stdout);

    fail = integrate(
        o.R0,
        o.rho0,
        o.omega,
        o.c,
        o.use_delay,
        o.use_lorentz_factor,
        o.use_lorentz_general_factor,
        o.use_fermi_factor_O,
        o.use_fermi_factor,
        o.use_fermi_general_factor,
        o.use_fast_integrand,

        o.use_fermi_general_factor_O,
        o.use_fermi_film_factor,
        o.use_fermi_film_instant,
        o.fermi_film_K,
        o.use_fermi_pde_factor,
        o.pde_grid_n,
        NULL,               /* сетку T строит integrate() */
        o.use_rot_to_lab,
        o.seed,
        o.maxeval,

        o.integrator,

        integral, error, prob);

    char result_names[9][32] = {
        "observed_ratio_1_rho",
        "observed_ratio_1_phi",
        "observed_ratio_1_z  ",
        "observed_ratio_2_rho",
        "observed_ratio_2_phi",
        "observed_ratio_2_z  ",
        "observed_ratio_rho  ",
        "observed_ratio_phi  ",
        "observed_ratio_z    "};

    printf("STATUS fail=%d integrator=%s\n", fail, integrator_name(o.integrator));
    for (comp = 0; comp < NCOMP; ++comp)
        printf("RESULT: %s %.8f +- %.8f\tp = %.3f\n", result_names[comp],
          (double)integral[comp], (double)error[comp], (double)prob[comp]);
    fflush(stdout);

    return 0;   /* fail = 1 у Cuba означает «не уложился в MAXEVAL», это не ошибка запуска */
}