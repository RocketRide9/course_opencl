#define real float
#define real4 float4

kernel void function(
    // матрица
    global const real *mat,
    global const real *di,
    global const real *b,
    global const int *ia,
    global const int *ja,
    const int n,
    // сетка
    global const real *X,
    const int xn,
    global const real *Y,
    const int yn
)
{
    uint row = get_global_id(0);
    uint col = get_global_id(1);
    uint idx = row * xn + col;
    
    
}
