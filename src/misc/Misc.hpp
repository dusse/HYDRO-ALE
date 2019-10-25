#ifndef Misc_hpp
#define Misc_hpp

#include <stdio.h>
#include <cmath>
#include <vector>

#define  IDX_OLD(i, j, k, n0, n1, n2) ((k)+(n2)*((j)+(n1)*((i)+(n0)*(0))))
#define  IDX(i, j, k, n0, n1, n2) ((j)+(n1)*(i)+(0)*(k)+(0)*(n2)+(n0)*(0))

#define EPSILON           (1.0E-18)
#define EPSILON_VOLUME    (1.0E-10 )
#define EPSILON_DOUBLE    (1.0E-36)
#define EPSILON_COMPARATOR    (1.0E-18)
#define PI           3.14159265358979323846 


bool areSame(double, double);

#endif /* Misc_hpp */
