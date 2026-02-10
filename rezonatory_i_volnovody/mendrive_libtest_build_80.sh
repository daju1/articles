#     # Выбор флагов компиляции
#     if precision == 'float128':
#         define_flag = "-DARITH_FLOAT128"
#         extra_libs = "-lquadmath"
#     elif precision == 'mpfr_512':
#         define_flag = "-DARITH_MPFR_512"
#         extra_libs = "-lmpfr -lgmp"
#     else:  # long_double
#         define_flag = "-DARITH_LONG_DOUBLE"
#         extra_libs = "-lm"

#     if 'DEBUG' == loglevel:
#         define_flag += ' -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_DEBUG'
#     elif 'TRACE' == lglevel:
#         define_flag += ' --DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_TRACE'
#     else:
#         define_flag += '-DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_OFF'





rm mendrive_libtest.so

gcc -shared -fPIC -O3 -o mendrive_libtest.so \
        mendrive_det.c \
        mendrive_isolines.c \
        mendrive_isolines_traced.c \
        mendrive_contour_intersections.c \
        mendrive_contour_intersections_sharp.c \
        mendrive_characteristic_roots.c \
        mendrive_root.c \
        mendrive_newton.c \
        mendrive_newton_classic/mendrive_newton_classic.c \
        mendrive_linsolve.c \
        -I../vimaniks/gsl/local/include \
        -L../vimaniks/gsl/local/lib \
        -DLOGGING -DKY \
        -DARITH_LONG_DOUBLE \
        -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_TRACE \
        -Wl,-rpath='$ORIGIN/../vimaniks/gsl/local/lib' \
        -lgsl -lgslcblas -lm

gcc mendrive_libtest.c -DKY -DARITH_LONG_DOUBLE \
        -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_TRACE \
        -ldl  -lm -o mendrive_libtest
