#define real double
#define real4 float4

real gamma(int subdom, real x, real y)
{
    return y*y;
}

real lambda(int subdom, real x, real y)
{
    return (real)0.5;
}

real f(int subdom, real x, real y)
{
    return (y*y-1.)*exp(x+y) + x*y*y;
}

constant real localG1[4][4] = {
    { 2, -2,  1, -1},
    {-2,  2, -1,  1},
    { 1, -1,  2, -2},
    {-1,  1, -2,  2},
};

constant real localG2[4][4] = {
    { 2,  1, -2, -1},
    { 1,  2, -1, -2},
    {-2, -1,  2,  1},
    {-1, -2,  1,  2},
};

constant real localM[4][4] = {
    {4, 2, 2, 1},
    {2, 4, 1, 2},
    {2, 1, 4, 2},
    {1, 2, 2, 4},
};
