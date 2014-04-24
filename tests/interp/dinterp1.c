#include <numcalc/interp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    double xt9[9], yt9[9];
    double xt73[73], yt73[73];
    double x[541], y[541];
    /* integer hint[541]; TODO */
    integer i;

    /* xt9 is ascending order [0:360] by 45 degrees */
    for (i = 0; i < 9; i++) {
        xt9[i] = M_PI * i / 4;
        yt9[i] = sin(xt9[i]);
    }

    /* xt73 is discending order [360:0] by 5 degrees */
    for (i = 0; i < 73; i++) {
        xt73[i] = M_PI * (72 - i) / 36;
        yt73[i] = sin(xt73[i]);
    }

    /* x is [-90:450] by 1 degree */
    for (i = 0; i < 541; i++) {
        x[i] = M_PI * (i - 90) / 180;
        /* hint[i] = 0; TODO */
    }

    /*
     * Then run and output results by csv format.
     * You can check the results using a plot application.
     *
     * Example: for gnuplot
     *   $ ./a.out >stdout.txt
     *   $ gnuplot -e "plot sin(x), 'stdout.txt' index 0 with line"
     */
    dinterp1(9, xt9, yt9, 541, x, y, NULL);
    for (i = 0; i < 541; i++) {
        printf("%e, %e\n", x[i], y[i]);
    }
    printf("\n\n");

    /*
     *   $ gnuplot -e "plot sin(x), 'stdout.txt' index 1 with line"
     */
    dinterp1(73, xt73, yt73, 541, x, y, NULL);
    for (i = 0; i < 541; i++) {
        printf("%e, %e\n", x[i], y[i]);
    }

    /* TODO: test for hint */
    /* interp1(73, xt73, yt73, 541, x, y, hint); */

    return EXIT_SUCCESS;
}
