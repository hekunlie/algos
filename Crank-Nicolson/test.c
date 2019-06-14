#include "stdio.h"
#include "Crank-Nicolson.h"
#include "math.h"

void TridMat_test(){

    int i;

    double a[4] = { 0, 3, 1, 3 };
    double b[4] = { 10, 10, 7, 4 };
    double c[4] = { 2, 4, 5, 0 };
    double d[4] = { 3, 4, 5, 6 };
    double x[4];

    TridMat( a, b, c, d, x, 4, 1 );
    for ( i=0; i<4; i++ )
        printf( "%g ", x[i] );
    printf( "\n" );

    TridMat( a, b, c, d, x, 4, 0 );
    for ( i=0; i<4; i++ )
        printf( "%g ", x[i] );
    printf( "\n" );


}

double f_zeros( double x, double t ){
    return 0;
}

double f_ones( double x, double t ){
    return 1;
}

void boundary_zeros( double *x, double *y, int N ){
    y[0] = y[N-1] = 0;
}


MYFUNC f_a, f_b, f_c, f_d, f_0, f_th;
double xmax, xmin, dt, dtt, t0, t1;
int xN, log_flag;
char FileName[100];
void (*boundary)( double *x, double *y, int N );

void Crank_Nicolson_test() {

    struct CN_Data p;
    double dx, *y0, *y, t, *x;
    int i;
    char buf[100];
    FILE *fd;

    if ( f_a == NULL )
        f_a = f_zeros;

    if ( f_b == NULL )
        f_b = f_zeros;

    if ( f_c == NULL )
        f_c = f_zeros;

    if ( f_d == NULL )
        f_d = f_zeros;

    if ( f_0 == NULL )
        f_0 = f_zeros;

    if ( boundary == NULL )
        boundary = boundary_zeros;

    y0 = malloc( sizeof( double ) * xN );
    y = malloc( sizeof( double ) * xN );
    x = malloc( sizeof( double ) * xN );
    if ( log_flag ) {
        dx =  log10(xmax/xmin) / ( xN-1  );
    }
    else {
        dx =  (xmax-xmin) / ( xN-1  );
    }

    sprintf( buf, "%s.dat", FileName );
    fd = fopen( buf, "w" );
    fprintf( fd, "0 " );
    for( i=0; i<xN; i++ ) {
        if ( log_flag )
            x[i] = pow( 10, log10(xmin) + i * dx);
        else
            x[i] = xmin + i * dx;
        fprintf( fd, "%g ", x[i] );
        y0[i] = (*f_0)(x[i], 0);
    }
    fprintf( fd, "\n" );

    fprintf( fd, "0 " );
    for( i=0; i<xN; i++ ) {
        fprintf( fd, "%g ", y0[i] );
    }
    fprintf( fd, "\n" );

    p.y0 = y0;
    p.y = y;
    p.xN = xN;
    p.t_step = dt;
    p.f_a = f_a;
    p.f_b = f_b;
    p.f_c = f_c;
    p.f_d = f_d;
    p.t_start = t0;
    p.boundary = boundary;
    p.xmax = xmax;
    p.xmin = xmin;

    if ( log_flag ) {
        p.x = x;
    }
    else {
        p.x = NULL;
    }

    CN_init( p );

    for( t=t0; t<t1; t+=dtt ) {

        CN_forward( dtt );

        fprintf( fd, "%g ", t+dtt );
        for( i=0; i<xN; i++ ) {
            fprintf( fd, "%g ", y[i] );
        }
        fprintf( fd, "\n" );
    }


    CN_free();
    free( y );
    fclose( fd );

    sprintf( buf, "%s_th.dat", FileName );
    fd = fopen( buf, "w" );

    fprintf( fd, "0 " );
    for( i=0; i<xN; i++ ) {
        fprintf( fd, "%g ", x[i] );
    }
    fprintf( fd, "\n" );

    fprintf( fd, "0 " );
    for( i=0; i<xN; i++ ) {
        fprintf( fd, "%g ", y0[i] );
    }
    free( y0 );
    fprintf( fd, "\n" );

    if ( f_th == NULL ) {
        free( x );
        fclose( fd );
        return;
    }

    for( t=t0; t<t1; t+=dtt ) {
        fprintf( fd, "%g ", t+dtt );
        for( i=0; i<xN; i++ ) {
            fprintf( fd, "%g ", (*f_th)( xmin+i*dx, t+dtt ) );
        }
        fprintf( fd, "\n" );
    }

    free( x );
    fclose( fd );

}

void boundary_linear( double *x, double *y, int N ){
    y[0] = (y[2]-y[1]) / (x[2]-x[1]) * (x[0]-x[1]) + y[1];
    y[N-1] = (y[N-2]-y[N-3]) / (x[N-2]-x[N-3]) * (x[N-1]-x[N-2]) + y[N-2];
}

void boundary_10( double *x, double *y, int N ){
    y[0] = 1;
    y[N-1] = 0;
}

void boundary_ylog( double *x, double *y, int N ){
    //y[0] = 0;
    //y[N-1] = 10000;
    //
    y[0] = exp( log(y[2]/y[1]) / (x[2]-x[1]) * (x[0]-x[1]) + log(y[1]));
    y[N-1] = exp( log(y[N-2]/y[N-3]) / (x[N-2]-x[N-3]) * (x[N-1]-x[N-2]) + log(y[N-2]));
}

/***************************************************/
double f_b1( double x, double t ){
    return 1/exp(x);
}

double f_c1( double x, double t ){
    return -2;
}

double f_a11( double x, double t ){
    return x*x;
}

double f_b11( double x, double t ){
    return x+1;
}

double f_c11( double x, double t ){
    return -2;
}

/***************************************************/

double f_b2( double x, double t ){
    return -1/exp(x);
}

double f_c2( double x, double t ){
    return -2;
}

double f_a21( double x, double t ){
    return x*x;
}

double f_b21( double x, double t ){
    return x-1;
}

double f_c21( double x, double t ){
    return -2;
}

/***************************************************/

double f_c3( double x, double t ){
    return -(1+1/exp(x));
}

double f_a31( double x, double t ){
    return x*x;
}

double f_b31( double x, double t ){
    return x;
}

double f_c31( double x, double t ){
    return -(1+1/x);
}

/***************************************************/
double f_a4( double x, double t ){
    return exp(x);
}

double f_b4( double x, double t ){
    return exp(x);
}

double f_c4( double x, double t ){
    return -(2*exp(x)+1);
}

double f_a41( double x, double t ){
    return x*x*x;
}

double f_b41( double x, double t ){
    return 2*x*x;
}

double f_c41( double x, double t ){
    return -2*x-1;
}

/***************************************************/

double f_delta( double x, double t ){
    double a, x0, r;
    x = exp(x);
    a=0.01;
    x0 = 0.1;
    //printf( "[in d] %g, %g, %g, %g\n", x, r, (x-x0)/a, pow( (x-x0)/a, 2 ) );
    r =  1 / (a*sqrt(M_PI)) * exp( -pow( (x-x0)/a, 2 ) );
    return r;
}

double f_delta1( double x, double t ){
    double a, x0, r;
    a=0.01;
    x0 = 0.1;
    //printf( "[in d] %g, %g, %g, %g\n", x, r, (x-x0)/a, pow( (x-x0)/a, 2 ) );
    r =  1 / (a*sqrt(M_PI)) * exp( -pow( (x-x0)/a, 2 ) );
    return r;
}


void main() {

    int i, j, model;

    dt = 1e-4;
    t0 = 0;
    t1 = 5;
    dtt = 1;
    f_0 = NULL;
    f_th = NULL;
    xN = 1000;
    boundary = NULL;

    for( i=1; i<5; i++)
        for( j=0; j<2; j++ ) {
            if ( i!=4 )
                continue;
            model = i*10 + j;
            log_flag = j;
            printf( "Model_%i%i ...\n", i, j );
            sprintf( FileName, "Model_%i%i", i, j );
            switch ( model ){
                case 10:
                    printf( "Model 1 [equal] for a=1\n" );
                    f_a = f_ones;
                    f_b = f_b1;
                    f_c = f_c1;
                    f_d = f_delta;
                    xmin = -6;
                    xmax = 20;
                    break;

                case 11:
                    printf( "Model 1 [unequal] for a=1\n" );
                    f_a = f_a11;
                    f_b = f_b11;
                    f_c = f_c11;
                    f_d = f_delta1;
                    xmin = 1e-5;
                    xmax = 1e5;
                    break;

                case 20:
                    printf( "Model 1 [equal] for a=-1\n" );
                    f_a = f_ones;
                    f_b = f_b2;
                    f_c = f_c2;
                    f_d = f_delta;
                    xmin = -6;
                    xmax = 20;
                    break;

                case 21:
                    printf( "Model 1 [unequal] for a=-1\n" );
                    f_a = f_a21;
                    f_b = f_b21;
                    f_c = f_c21;
                    f_d = f_delta1;
                    xmin = 1e-5;
                    xmax = 1e5;
                    break;
                case 30:

                    printf( "Model 2 [euqal] for a=1\n" );
                    f_a = f_ones;
                    f_b = NULL;
                    f_c = f_c3;
                    f_d = f_delta;
                    xmin = -6;
                    xmax = 20;
                    xN = 300;
                    break;

                case 31:

                    printf( "Model 2 [unequal] for a=1\n" );
                    f_a = f_a31;
                    f_b = f_b31;
                    f_c = f_c31;
                    f_d = f_delta1;
                    xmin = 1e-5;
                    xmax = 1e5;
                    xN = 300;
                    break;

                case 40:
                    printf( "Model 3 for a=1\n" );
                    f_a = f_a4;
                    f_b = f_b4;
                    f_c = f_c4;
                    f_d = NULL;
                    f_0 = f_delta;
                    xmin = -5;
                    xmax = 20;
                    break;
                case 41:
                    printf( "Model 3 for a=1\n" );
                    f_a = f_a41;
                    f_b = f_b41;
                    f_c = f_c41;
                    f_d = NULL;
                    f_0 = f_delta1;
                    xmin = 1e-5;
                    xmax = 1e5;
                    break;

                }

                Crank_Nicolson_test();
        }

}


