#define real double

#ifdef __MESA__
#include "FEMShared.cl"
#else
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

/*
void __attribute__((always_inline))
atomicAdd_g2(volatile global double *addr, const float val)
{
    union {
        ulong u64;
        double f64;
    } next, expected, current;
    current.f64 = *addr;
    do {
        expected.f64 = current.f64;
        next.f64 = expected.f64 + val;
        current.u64 = atom_cmpxchg(
            (volatile global ulong *)addr,
            expected.u64, next.u64);
    } while( current.u64 != expected.u64 );
}
*/

#if false
void __attribute__((always_inline))
atomicAdd_gf(volatile global float *addr, const float val)
{
    union {
        uint u32;
        float f32;
    } next, expected, current;
    current.f32 = *addr;
    do {
        next.f32 = (expected.f32 = current.f32) + val;
        current.u32 = atom_cmpxchg(
            (volatile global uint *)addr,
            expected.u32, next.u32);
    } while( current.u32 != expected.u32 );
}
#else 
void __attribute__((always_inline))
atomicAdd_g(volatile __global double *addr, double val)
{
    union {
        unsigned long u64;
        double f64;
    } next, expected, current;
    current.f64 = *addr;
    do {
        expected.f64 = current.f64;
        next.f64 = expected.f64 + val;
        current.u64 = atom_cmpxchg(
            (volatile global unsigned long *)addr,
            expected.u64, next.u64);
    } while( current.u64 != expected.u64 );
}
#endif

#endif


// test
kernel void
add_to_zero(volatile global real *zero)
{
    int y = get_global_id(1);
    int x = get_global_id(0);

    real tmp = ((x+y) % 2) == 0 ? 1 : -1;
    atomicAdd_g(zero, tmp);
}

void __attribute__((always_inline))
fill_mat(real mat[4][4], real l_avg, real g_avg, real hx, real hy)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mat[i][j] = l_avg / 6 * (hy / hx * localG1[i][j] + hx / hy * localG2[i][j])
                     + g_avg / 36 * hx * hy * localM[i][j];
        }
    }
}

kernel void global_matrix_compose(
    // матрица
    volatile global real *ld3,
    volatile global real *ld2,
    volatile global real *ld1,
    volatile global real *ld0,
    
    volatile global real *di,
    
    volatile global real *rd0,
    volatile global real *rd1,
    volatile global real *rd2,
    volatile global real *rd3,

    volatile global real *b,

    const int n,
    const int gap,
    // сетка
    global const real *axis_x,
    const int xn,
    global const real *axis_y,
    const int yn
)
{
    uint xi = get_global_id(0);
    uint yi = get_global_id(1);

    if (xi >= xn - 1 || yi >= yn - 1) return;

    // TODO: изменить позже
    int subDom = 0;

    real x0 = axis_x[xi];
    real x1 = axis_x[xi + 1];
    real y0 = axis_y[yi];
    real y1 = axis_y[yi + 1];

    real hy = y1 - y0;
    real hx = x1 - x0;

    real f1 = f(subDom, x0, y0);
    real f2 = f(subDom, x1, y0);
    real f3 = f(subDom, x0, y1);
    real f4 = f(subDom, x1, y1);

    private real localB[4];
    localB[0] = hx * hy / 36. * (4. * f1 + 2. * f2 + 2. * f3 +      f4);
    localB[1] = hx * hy / 36. * (2. * f1 + 4. * f2 +      f3 + 2. * f4);
    localB[2] = hx * hy / 36. * (2. * f1 +      f2 + 4. * f3 + 2. * f4);
    localB[3] = hx * hy / 36. * (     f1 + 2. * f2 + 2. * f3 + 4. * f4);

    real l_avg = lambda(subDom, x0, y0) + lambda(subDom, x1, y0)
               + lambda(subDom, x0, y1) + lambda(subDom, x0, y1);
    l_avg /= 4.;

    real g_avg = gamma(subDom, x0, y0) + gamma(subDom, x1, y0)
               + gamma(subDom, x0, y1) + gamma(subDom, x1, y1);
    g_avg /= 4.;

    int anchor = yi * xn + xi;

    real mat[4][4];
    fill_mat(mat, l_avg, g_avg, hx, hy);

    atomicAdd_g(&di[anchor],  mat[0][0]);
    atomicAdd_g(&ld0[anchor], mat[0][1]);
    atomicAdd_g(&rd0[anchor], mat[0][1]);
    atomicAdd_g(&ld2[anchor], mat[0][2]);
    atomicAdd_g(&rd2[anchor], mat[0][2]);
    atomicAdd_g(&ld3[anchor], mat[0][3]);
    atomicAdd_g(&rd3[anchor], mat[0][3]);

    atomicAdd_g(&di[anchor + 1], mat[1][1]);
    atomicAdd_g(&ld1[anchor + 1], mat[1][2]);
    atomicAdd_g(&rd1[anchor + 1], mat[1][2]);
    atomicAdd_g(&ld2[anchor + 1], mat[1][3]);
    atomicAdd_g(&rd2[anchor + 1], mat[1][3]);

    int a2 = anchor + xn;
    atomicAdd_g(&di[a2], mat[2][2]);
    atomicAdd_g(&ld0[a2], mat[2][3]);
    atomicAdd_g(&rd0[a2], mat[2][3]);

    atomicAdd_g(&di[a2 + 1], mat[3][3]);

    atomicAdd_g(&b[anchor], localB[0]);
    atomicAdd_g(&b[anchor + 1], localB[1]);
    atomicAdd_g(&b[a2], localB[2]);
    atomicAdd_g(&b[a2 + 1], localB[3]);

    for (int i = 0; i < 4; i++)
    {
        // b[m[i]] += localB[i];
    }
}
