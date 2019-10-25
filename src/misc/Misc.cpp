
#include "Misc.hpp"
using namespace std;

bool areSame(double a, double b)
{
    return fabs(a - b) < EPSILON_COMPARATOR;
}

