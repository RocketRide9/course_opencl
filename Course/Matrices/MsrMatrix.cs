#define HOST_PARALLEL
using Real = double;

using System.Collections.Concurrent;

using Types;

namespace Matrices;

public class MsrMatrix : Matrix
{
    public Real[] Elems = [];
    public Real[] Di = [];
    public int[] Ia = [];
    public int[] Ja = [];

    int Matrix.Size => Di.Length;
    Span<Real> Matrix.Di => Di;

    public SparkAlgos.Types.Matrix GetComputeMatrix()
    {
        return new SparkAlgos.Matrices.MsrMatrix(new()
        {
            Elems = Elems,
            Di = Di,
            Ia = Ia,
            Ja = Ja,
        });
    }

    // TODO: можно переписать на Array<> так как эта функция теперь здесь
    unsafe void Matrix.Mul(ReadOnlySpan<Real> vec, Span<Real> res)
    {
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, Ia.Length);
        fixed(Real* _p_mat = Elems)
        fixed(Real* _p_di = Di)
        fixed(int* _p_ia = Ia)
        fixed(int* _p_ja = Ja)
        fixed(Real* _p_v = vec)
        fixed(Real* _p_res = res)
        {
            var p_mat = _p_mat;
            var p_di = _p_di;
            var p_ia = _p_ia;
            var p_ja = _p_ja;
            var p_v = _p_v;
            var p_res = _p_res;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    int start = p_ia[i];
                    int stop = p_ia[i + 1];
                    Real dot = p_di[i] * p_v[i];
                    for (int a = start; a < stop; a++)
                    {
                        dot += p_mat[a] * p_v[p_ja[a]];
                    }
                    p_res[i] = dot;
                }
            });
        }
#else
        fixed(Real* p_mat = Elems)
        fixed(Real* p_di = Di)
        fixed(int* p_ia = Ia)
        fixed(int* p_ja = Ja)
        fixed(Real* p_v = vec)
        fixed(Real* p_res = res)
        for (int i = 0; i < Ia.Length - 1; i++)
        {
            int start = p_ia[i];
            int stop = p_ia[i + 1];
            Real dot = p_di[i] * p_v[i];
            for (int a = start; a < stop; a++)
            {
                dot += p_mat[a] * p_v[p_ja[a]];
            }
            p_res[i] = dot;
        }
#endif
    }
}
