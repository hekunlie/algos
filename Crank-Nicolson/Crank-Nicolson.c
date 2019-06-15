#include "Crank-Nicolson.h"

static double *x, *dx1, *dx2, *y, *y1, t_step,
              *y0, *A, *B, *C, *D, dx, xmax, xmin, t_cur, t_start;
static int xN, equal_space;
static MYFUNC f_a;
static MYFUNC f_b;
static MYFUNC f_c;
static MYFUNC f_d;
static void (*boundary)( double *x, double *y, int );

void CN_init( struct CN_Data p) {

    int i;
    xN = p.xN;
    t_step = p.t_step;
    f_a = p.f_a;
    f_b = p.f_b;
    f_c = p.f_c;
    f_d = p.f_d;
    boundary = p.boundary;
    y = p.y;
    t_cur = t_start = p.t_start;

    y0 = malloc( sizeof( double ) * xN );
    y1 = malloc( sizeof( double ) * xN );
    memcpy( y0, p.y0, sizeof(double) *xN );

    if ( p.x != NULL ) {
        equal_space = 0;
        x = p.x;
        dx1 = malloc( sizeof( double ) * xN );
        dx2 = malloc( sizeof( double ) * xN );

        for( i=0; i<xN-2; i++ ) {
            dx1[i] = x[i+1] - x[i];
            dx2[i] = x[i+2] - x[i];
        }
        dx1[xN-2] = x[xN-1] - x[xN-2];
    }
    else {
        equal_space = 1;
        xmax = p.xmax;
        xmin = p.xmin;
        dx = ( xmax-xmin ) / ( xN-1 );
        x = malloc( sizeof( double ) * xN );
        for( i=0; i<xN; i++ )
            x[i] = xmin + i * dx;

    }

    A = malloc( sizeof(double)*xN );
    B = malloc( sizeof(double)*xN );
    C = malloc( sizeof(double)*xN );
    D = malloc( sizeof(double)*xN );

    if ( equal_space )
        printf( "User equal space step\n" );

}

void print_ABCD() {

    int i;

    printf( "A: " );
    for( i=0; i<xN; i++ )
        printf( "%g ", A[i]);
    printf( "\n" );

    printf( "B: " );
    for( i=0; i<xN; i++ )
        printf( "%g ", B[i]);
    printf( "\n" );

    printf( "C: " );
    for( i=0; i<xN; i++ )
        printf( "%g ", C[i]);
    printf( "\n" );

    printf( "D: " );
    for( i=0; i<xN; i++ )
        printf( "%g ", D[i]);
    printf( "\n" );
}

void CN_single_forward() {

    int i;
    double a, b, c, d, tmid, h, k, h2;

    tmid = t_cur + t_step / 2.0;

    if ( equal_space ) {

        k = t_step;
        h = dx;
        h2 = h*h;

        for( i=1; i<xN-1; i++ ) {

             a = (*f_a)( x[i], tmid );
             b = (*f_b)( x[i], tmid );
             c = (*f_c)( x[i], tmid );
             d = (*f_d)( x[i], tmid );
             //printf( "%g %g %g %g\n", a, b, c, d );

             A[i] = -( 2*k*a - k*h*b );
             B[i] = 4*h2 + 4*k*a - 2*h2*k*c;
             C[i] = -(2*k*a+k*h*b);

             D[i] = ( 2*k*a-k*h*b ) * y0[i-1]
                  + ( 4*h2-4*k*a+2*h2*k*c ) * y0[i]
                  + ( 2*k*a + k*h*b ) * y0[i+1]
                  + 4*h2*k*d;

        }

    }
    else {

        for( i=1; i<xN-1; i++ ) {

             a = (*f_a)( x[i], tmid );
             b = (*f_b)( x[i], tmid );
             c = (*f_c)( x[i], tmid );
             d = (*f_d)( x[i], tmid );
             //printf( "%g %g %g %g\n", a, b, c, d );

             A[i] = -a / (dx1[i-1]*dx2[i-1]) + b / (2*dx2[i-1]);
             B[i] = 1/t_step + a / (dx1[i]*dx1[i-1]) - c/2;
             C[i] = -a / ( dx1[i]*dx2[i-1]) - b / (2*dx2[i-1]);

             D[i] = ( a / (dx1[i-1]*dx2[i-1]) - b / (2*dx2[i-1]) ) * y0[i-1]
                  + ( 1/t_step - a / (dx1[i]*dx1[i-1]) + c / 2 ) * y0[i]
                  + ( a / (dx1[i]*dx2[i-1]) + b / (2*dx2[i-1]) ) * y0[i+1]
                  + d;

        }
    }

    //print_ABCD();

    TridMat( A+1, B+1, C+1, D+1, y1+1, xN-2, 0 );
    (*boundary)( x, y1, xN );

}

void CN_forward( double t_forward ) {

    //int i;
    double dt_tmp, *y_tmp, t_end;
    t_end = t_cur + t_forward;
    printf( "[%g] --> [%g]\n", t_cur, t_end );

    while( t_cur+t_step<=t_end ){
        CN_single_forward();;
        t_cur += t_step;

        y_tmp = y0;
        y0 = y1;
        y1 = y_tmp;
    }

    if ( t_cur < t_end ){

        dt_tmp = t_step;
        t_step = t_end - t_cur;
        CN_single_forward();
        t_step = dt_tmp;

        y_tmp = y0;
        y0 = y1;
        y1 = y_tmp;
    }

    t_cur = t_end;
    memcpy( y, y0, sizeof( double ) * xN );


}

void CN_free() {

    free( y0 );
    free( y1 );

    if ( !equal_space ) {
        free( dx1 );
        free( dx2 );
    }
    else {
        free(x);
    }

    free( A );
    free( B );
    free( C );
    free( D );

}

void TridMat( double *a, double *b, double *c, double *d, double *x, int N, int flag ) {

    int i;
    double *cc, *dd, w;

    if ( flag==1 ) {

        cc = malloc( sizeof( double ) * N );
        dd = malloc( sizeof( double ) * N );

        cc[0] = c[0] / b[0];
        dd[0]  = d[0] / b[0];

        for( i=1; i<N; i++ ) {

            cc[i] = c[i] / ( b[i] - a[i] * cc[i-1] );
            dd[i] = ( d[i] - a[i] * dd[i-1]  ) / ( b[i] - a[i] * cc[i-1] );

        }

        x[N-1] = dd[N-1];
        for ( i=N-2; i>=0; i-- ) {
            x[i] = dd[i] - cc[i] * x[i+1];
        }

        free( cc );
        free( dd );

    }
    else {

        for ( i=1; i<N; i++ ) {
            w = a[i] / b[i-1];
            //printf( "[%i] b: %g, d: %g,", i, b[i], d[i] );
            b[i] = b[i] - w * c[i-1];
            d[i] = d[i] - w * d[i-1];
            //printf( "w: %g, b: %g, d: %g\n", w, b[i], d[i]  );
        }

        x[N-1] = d[N-1] / b[N-1];
        for ( i=N-2; i>=0; i-- )
            x[i] = ( d[i] - c[i] * x[i+1]) / b[i];

    }

}
