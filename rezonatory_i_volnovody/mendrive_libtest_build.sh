# gcc -shared -fPIC -O3 -lm -lgsl -lgslcblas \
#     mendrive_det.c \
#     mendrive_isolines.c \
#     mendrive_contour_intersections.c \
#     mendrive_contour_intersections_sharp.c \
#     -o mendrive.so

gcc -shared -fPIC -O3 -o mendrive.so \
        mendrive_det.c \
        mendrive_isolines.c \
        mendrive_isolines_traced.c \
        mendrive_contour_intersections.c \
        mendrive_contour_intersections_sharp.c \
        mendrive_characteristic_roots.c \
        mendrive_root.c \
        mendrive_newton.c \
        mendrive_linsolve.c \
        -I../vimaniks/gsl/local/include \
        -L../vimaniks/gsl/local/lib \
        -DLOGGING \
        -Wl,-rpath='$ORIGIN/../vimaniks/gsl/local/lib' \
        -lgsl -lgslcblas -lm

gcc mendrive_libtest.c -ldl -o mendrive_libtest