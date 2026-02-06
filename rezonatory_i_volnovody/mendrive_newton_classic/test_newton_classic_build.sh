gcc mendrive_newton_classic.c \
    test_newton_classic.c \
    -DARITH_FLOAT128 \
    -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_TRACE \
    -ldl -lm -lquadmath -o test_newton_classic