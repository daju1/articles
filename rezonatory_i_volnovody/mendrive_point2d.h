#ifndef MENDDRIVE_POINT2D_H
#define MENDDRIVE_POINT2D_H

typedef struct {
    long double kz;
    long double sz;
} point2d_t;

typedef struct {
    long double kz;
    long double sz;
    long double cosine;
    long double sine;
    int type;
    int i;
    int refined_i;
    point2d_t refined_pt;
    #ifdef USE_REFILTER
    int is_confirmed;
    #endif
} corner2d_t;

#endif