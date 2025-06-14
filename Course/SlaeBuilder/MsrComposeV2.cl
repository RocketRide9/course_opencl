#define real double

#ifdef false
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
#endif


int qfind(global const int *where, int what , int start, int end)
{
    int beg = start;
    while (beg < end)
    {
        int mid = (beg + end) / 2;
        if (what > where[mid])
        {
            beg = mid + 1;
        }
        else
        {
            end = mid;
        }
    }

    return beg; 
}

int __attribute__((always_inline))
lfind(global const int *where, int what, int start)
{
    while (where[start] != what) start++;
    return start;
}

real __attribute__((always_inline)) 
get_l_avg (
    int dom,
    real x0, real y0,
    real x1, real y1
) {
    real res = lambda(dom, x0, y0)
             + lambda(dom, x1, y0)
             + lambda(dom, x0, y1)
             + lambda(dom, x1, y1);

    return res / 4;
}

real __attribute__((always_inline))
get_g_avg (
    int dom,
    real x0, real y0,
    real x1, real y1
) {
    real res = gamma(dom, x0, y0)
             + gamma(dom, x1, y0)
             + gamma(dom, x0, y1)
             + gamma(dom, x1, y1);

    return res / 4;
}

kernel void global_matrix_compose_v2(
    // матрица
    global real *mat,
    global real *di,
    global real *b,
    global const int *ia,
    global const int *ja,
    const int n,
    // сетка
    global const real *axis_x,
    const int xn,
    global const real *axis_y,
    const int yn
) {
    int xi = get_global_id(0);
    int yi = get_global_id(1);

    if (xi >= xn || yi >= yn) return;

    int targetNode = yi * xn + xi;
    // TODO: основано на нескольких предположениях об области решения
    int dom1 = (xi-1 >= 0  && yi-1 >= 0 ) ? 0 : -1;
    int dom2 = (xi+1 <  xn && yi-1 >= 0 ) ? 0 : -1;
    int dom3 = (xi-1 >= 0  && yi+1 <  yn) ? 0 : -1;
    int dom4 = (xi+1 <  xn && yi+1 <  yn) ? 0 : -1;

    // printf("xi, yi = %d, %d\n", xi, yi);
    // printf("xn, yn = %d, %d\n", xn, yn);
    // printf("dom = %d, %d, %d, %d\n", dom1, dom2, dom3, dom4);

    // номера текущего узла и соседних с ним сверху и снизу 
    // в порядке снизу вверх
    // могут выходить за пределы нумерации
    int r1 = targetNode;
    int r0 = r1 - xn;
    int r2 = r1 + xn;

    // номера этих узлов в строке матрицы, относящейся к текущему узлу
    int mr0 = -1;
    int mr1 = -1;
    int mr2 = -1;
    
    int beg = ia[targetNode];
    int bound = ia[targetNode + 1] - 1;
    
    #if 1
    // TODO: можно немного ускорить поиск за счёт перемещения левой границы поиска
    // но возможно быстрее будет сделать простой линейный поиск вместо "быстрого"
    // int curr = beg;
    if (r0 >= 0)
    {
        mr0 = lfind(ja, r0, beg) - beg;
    }
    // else оставить -1, так как узел вышел за пределы матрицы

    if (r1 % xn == 0)
    {
        mr1 = lfind(ja, r1+1, beg) - 1 - beg;
    }
    else
    {
        mr1 = lfind(ja, r1-1, beg) - beg;
    }

    if (r2 < n)
    {
        mr2 = lfind(ja, r2, beg) - beg;
    }
    #endif
    
    // printf("mr = %d, %d, %d\n", mr0, mr1, mr2);

    real matc[8] = {};
    
    // return;
    real dic = 0;
    real bc = 0;

    // есть в каждом случае
    real x1 = axis_x[xi];
    real y1 = axis_y[yi];
    
    if (dom1 != -1)
    {
        real x0 = axis_x[xi-1];
        real y0 = axis_y[yi-1];
        real l_avg = get_l_avg(dom1, x0, y0, x1, y1);
        real g_avg = get_g_avg(dom1, x0, y0, x1, y1);
        real hx0 = x1 - x0;
        real hy0 = y1 - y0;

        matc[mr0 - 1] += l_avg / 6 * (hy0 / hx0 * localG1[3][0] + hx0 / hy0 * localG2[3][0])
            + g_avg / 36 * hx0 * hy0 * localM[3][0];
        matc[mr0]     += l_avg / 6 * (hy0 / hx0 * localG1[3][1] + hx0 / hy0 * localG2[3][1])
            + g_avg / 36 * hx0 * hy0 * localM[3][1];
        matc[mr1]     += l_avg / 6 * (hy0 / hx0 * localG1[3][2] + hx0 / hy0 * localG2[3][2])
            + g_avg / 36 * hx0 * hy0 * localM[3][2];
        dic           += l_avg / 6 * (hy0 / hx0 * localG1[3][3] + hx0 / hy0 * localG2[3][3])
            + g_avg / 36 * hx0 * hy0 * localM[3][3];

        real f1 = f(dom1, x0, y0);
        real f2 = f(dom1, x1, y0);
        real f3 = f(dom1, x0, y1);
        real f4 = f(dom1, x1, y1);

        bc += hx0 * hy0 / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);
    }

    if (dom2 != -1)
    {
        real x2 = axis_x[xi+1];
        real y0 = axis_y[yi-1];
        real l_avg = get_l_avg(dom2, x1, y0, x2, y1);
        real g_avg = get_g_avg(dom2, x1, y0, x2, y1);
        real hx1 = x2 - x1;
        real hy0 = y1 - y0;

        matc[mr0]     += l_avg / 6 * (hy0 / hx1 * localG1[2][0] + hx1 / hy0 * localG2[2][0])
            + g_avg / 36 * hx1 * hy0 * localM[2][0];
        matc[mr0 + 1] += l_avg / 6 * (hy0 / hx1 * localG1[2][1] + hx1 / hy0 * localG2[2][1])
            + g_avg / 36 * hx1 * hy0 * localM[2][1];
        dic           += l_avg / 6 * (hy0 / hx1 * localG1[2][2] + hx1 / hy0 * localG2[2][2])
            + g_avg / 36 * hx1 * hy0 * localM[2][2];
        matc[mr1 + 1] += l_avg / 6 * (hy0 / hx1 * localG1[2][3] + hx1 / hy0 * localG2[2][3])
            + g_avg / 36 * hx1 * hy0 * localM[2][3];

        real f1 = f(dom2, x1, y0);
        real f2 = f(dom2, x2, y0);
        real f3 = f(dom2, x1, y1);
        real f4 = f(dom2, x2, y1);

        bc += hx1 * hy0 / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
    }

    if (dom3 != -1)
    {
        real x0 = axis_x[xi-1];
        real y2 = axis_y[yi+1];
        real l_avg = get_l_avg(dom3, x0, y1, x1, y2);
        real g_avg = get_g_avg(dom3, x0, y1, x1, y2);
        real hx0 = x1 - x0;
        real hy1 = y2 - y1;

        matc[mr1]     += l_avg / 6 * (hy1 / hx0 * localG1[1][0] + hx0 / hy1 * localG2[1][0])
            + g_avg / 36 * hx0 * hy1 * localM[1][0];
        dic           += l_avg / 6 * (hy1 / hx0 * localG1[1][1] + hx0 / hy1 * localG2[1][1])
            + g_avg / 36 * hx0 * hy1 * localM[1][1];
        matc[mr2 - 1] += l_avg / 6 * (hy1 / hx0 * localG1[1][2] + hx0 / hy1 * localG2[1][2])
            + g_avg / 36 * hx0 * hy1 * localM[1][2];
        matc[mr2]     += l_avg / 6 * (hy1 / hx0 * localG1[1][3] + hx0 / hy1 * localG2[1][3])
            + g_avg / 36 * hx0 * hy1 * localM[1][3];

        real f1 = f(dom3, x0, y1);
        real f2 = f(dom3, x1, y1);
        real f3 = f(dom3, x0, y2);
        real f4 = f(dom3, x1, y2);

        bc += hx0 * hy1 / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
    }

    if (dom4 != -1)
    {
        real x2 = axis_x[xi+1];
        real y2 = axis_y[yi+1];
        real l_avg = get_l_avg(dom4, x1, y1, x2, y2);
        real g_avg = get_g_avg(dom4, x1, y1, x2, y2);
        real hx1 = x2 - x1;
        real hy1 = y2 - y1;

        dic           += l_avg / 6 * (hy1 / hx1 * localG1[0][0] + hx1 / hy1 * localG2[0][0])
            + g_avg / 36 * hx1 * hy1 * localM[0][0];
        matc[mr1 + 1] += l_avg / 6 * (hy1 / hx1 * localG1[0][1] + hx1 / hy1 * localG2[0][1])
            + g_avg / 36 * hx1 * hy1 * localM[0][1];
        matc[mr2]     += l_avg / 6 * (hy1 / hx1 * localG1[0][2] + hx1 / hy1 * localG2[0][2])
            + g_avg / 36 * hx1 * hy1 * localM[0][2];
        matc[mr2 + 1] += l_avg / 6 * (hy1 / hx1 * localG1[0][3] + hx1 / hy1 * localG2[0][3])
            + g_avg / 36 * hx1 * hy1 * localM[0][3];

        real f1 = f(dom4, x1, y1);
        real f2 = f(dom4, x2, y1);
        real f3 = f(dom4, x1, y2);
        real f4 = f(dom4, x2, y2);

        bc += hx1 * hy1 / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
    }
    
#if 1
    di[targetNode] = dic;
    b[targetNode] = bc;

    // перемещение строки в матрицу
    for (int i = beg; i <= bound; i++)
    {
        mat[i] = matc[i-beg];
    }
#endif
    
#if 0
    real store = dic;
    store += bc;
    for (int i = beg; i <= bound; i++)
    {
        store += matc[i-beg];
    }
    di[targetNode] = store;
#endif
}
