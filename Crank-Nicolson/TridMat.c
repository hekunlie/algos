#include "TridMat.h"

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
