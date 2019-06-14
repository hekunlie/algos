#include "stdio.h"
#include "stdlib.h"
#include "TridMat.h"
#include "string.h"

typedef double (*MYFUNC)( double, double );

struct CN_Data{
    double *x, t_step, *y0, *y, xmin, xmax, t_start;
    int xN;
    MYFUNC f_a;
    MYFUNC f_b;
    MYFUNC f_c;
    MYFUNC f_d;
    void (*boundary)( double *x, double *y, int N );

};

void CN_init( struct CN_Data p );
void CN_forward( double t_forward );
void CN_test();
void CN_free();
