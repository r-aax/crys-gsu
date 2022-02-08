// file: src/hw.c
#include <stdio.h>
#include <math.h>
#include <sstream>
#include "swim.h"

using namespace std;

void hi(int n1, double *a1,
        int n2, double *a2,
        int n3, double *a3)
{
    (void) n2;
    (void) n3;

    for (int i = 0; i < n1; i++)
    {
        a1[i] += 1;
        a2[i] += 2;
        a3[i] += 3;
    }
    
}
