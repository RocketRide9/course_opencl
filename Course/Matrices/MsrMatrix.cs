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

    public int Size => Di.Length;
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

    public IEnumerable<Real> FlatNonZero()
    {
        for (int i = 0; i < Ia.Length - 1; i++)
        {
            int start = Ia[i];
            int stop = Ia[i + 1];

            int curr_a = start;

            for (_ = 0; curr_a < stop; curr_a++)
            {
                if (Ja[curr_a] > i) {
                    break;
                } else {
                    if (Elems[curr_a] != 0) yield return Elems[curr_a];
                }
            }
            yield return Di[i];
            for (_ = 0; curr_a < stop; curr_a++)
            {
                if (Elems[curr_a] != 0) yield return Elems[curr_a];
            }
        }
    }

    public unsafe void Mul(ReadOnlySpan<Real> vec, Span<Real> res)
    {
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, Ia.Length - 1);
        fixed(Real* _p_v = vec)
        fixed(Real* _p_res = res)
        {
            var p_v = _p_v;
            var p_res = _p_res;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    int start = Ia[i];
                    int stop = Ia[i + 1];
                    Real dot = Di[i] * p_v[i];
                    for (int a = start; a < stop; a++)
                    {
                        dot += Elems[a] * p_v[Ja[a]];
                    }
                    p_res[i] = dot;
                }
            });
        }
#else
        fixed(Real* p_v = vec)
        fixed(Real* p_res = res)
        for (int i = 0; i < Ia.Length - 1; i++)
        {
            int start = Ia[i];
            int stop = Ia[i + 1];
            Real dot = Di[i] * p_v[i];
            for (int a = start; a < stop; a++)
            {
                dot += Elems[a] * p_v[Ja[a]];
            }
            p_res[i] = dot;
        }
#endif
    }
}
