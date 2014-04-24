#include <numcalc/blas.h>
#include <stdio.h>

int main(void) {
    unsigned char c;
    size_t expected, got;

    expected = 0;
    got = 0;
    c = 0;
    while (1) {
        expected += c == 'N' || c == 'n';
        got += lsame('N', c);
        if (expected != got || c == 255) break;
        c++;
    }
    printf("expected: %ld, got: %ld\n", expected, got);

    return 0;
}
