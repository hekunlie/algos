#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.style.use( 'ggplot' )



def delta_fun( x ):
    a = 0.01
    x0 = 0.1
    return 1 / ( a * np.sqrt(np.pi) ) * \
           np.exp( -np.power((x-x0)/a, 2) )


for mj in range(0,2):
    fig, axs = plt.subplots( 1, 4, figsize=(4*4, 4) )

    for mi in range(1,5):

        #if mi != 4:
        #    continue

        ax = axs[mi-1]
        d = np.loadtxt( 'Model_%i%i.dat'%( mi, mj ) )
        x = d[0, 1:]
        t = d[1:, 0]
        d = d[1:, 1:]

        if mj == 0:
            x = np.exp( x )

        print( 'x:(%g, %g)'%(x.min(), x.max()) )

        d_th = np.loadtxt( 'Model_%i%i_th.dat'%( mi, mj ) )
        m, n = d_th.shape

        flag = 1
        if m == 2:
            flag = 0

        if flag:
            x_th = d_th[0, 1:]
            t_th = d_th[1:, 0]
            d_th = d_th[1:, 1:]

        ax.plot( x, delta_fun( x ), 'k-', label='inj' )

        for i in range( len(t) ):

            y = d[i, :]
            print( 'plot t: %g, y:(%g, %g)'%(t[i], y.min(), y.max()) )
            ax.plot( x, y, label='%g'%t[i] )

            if flag:
                y_th = d_th[i, :]
                ax.plot( x_th, y_th, '-.', label='%g'%t[i] )

        ax.legend( loc='lower left' )
        ax.set_ylim( [1e-10, 100] )
        ax.set_xscale( 'log' )
        ax.set_xlim( [1e-4, 1e4] )
        ax.set_yscale( 'log' )
        ax.set_ylabel( r'$\rm f(x)$' )
        ax.set_xlabel( r'$\rm x$' )

    plt.tight_layout()
    plt.savefig( 'Model_%i.png'%mj )
