using Real = double;

using System.Collections.Concurrent;

using Types;

namespace Matrices;

public class SymDiagMatrix : Matrix
{
    // diagonal
    public Real[] d3 = [];
    public Real[] d2 = [];
    public Real[] d1 = [];
    public Real[] d0 = [];
    // Основная диагональ
    public Real[] Di = [];

    // d0 находятся "вплотную" к основной диагонали
    // d1, d2, d3 находятся стоят "вплотную" друг к другу
    // d1 смещена на Gap элементов от d0.
    // Например, если они находятся вплотную друг к друг,
    // то Gap == 1.
    public int Gap;
    
    public int Size => Di.Length;
    Span<Real> Matrix.Di => Di;

    public SparkAlgos.Types.Matrix GetComputeMatrix()
    {
        return new SparkAlgos.Matrices.SymDiagMatrix(new(){
            d3 = d3,
            d2 = d2,
            d1 = d1,
            d0 = d0,
            Di = Di,
            Gap = Gap,
        });
    }
    
    public IEnumerable<Real> FlatNonZero()
    {
        for (int i = 0; i < Size; i++)
        {
            int t = i - 3 - Gap;
            if (t >= 0 && d3[t] != 0) yield return d3[t];
            t = i - 2 - Gap;
            if (t >= 0 && d2[t] != 0) yield return d2[t];
            t = i - 1 - Gap;
            if (t >= 0 && d1[t] != 0) yield return d1[t];
            t = i - 1;
            if (t >= 0 && d0[t] != 0) yield return d0[t];

            yield return Di[i];

            t = i+1;
            if (t < Size && d0[i] != 0) yield return d0[i];
            t = i+1+Gap;
            if (t < Size && d1[i] != 0) yield return d1[i];
            t = i+2+Gap;
            if (t < Size && d2[i] != 0) yield return d2[i];
            t = i+3+Gap;
            if (t < Size && d3[i] != 0) yield return d3[i];   
        }
    }

    public unsafe void Mul(ReadOnlySpan<Real> vec, Span<Real> res)
    {
        var partitioner = Partitioner.Create(0, Size);
        fixed(Real* _p_v = vec)
        fixed(Real* _p_res = res)
        {
        var p_v = _p_v;
        var p_res = _p_res;
        Parallel.ForEach(partitioner, (range, state) =>
        {
            for (int i = range.Item1; i < range.Item2; i++)
            {
                Real dot = 0;
                
                int t = i - 3 - Gap;
                if (t >= 0) dot += d3[t] * p_v[t];
                t = i - 2 - Gap;
                if (t >= 0) dot += d2[t] * p_v[t];
                t = i - 1 - Gap;
                if (t >= 0) dot += d1[t] * p_v[t];
                t = i - 1;
                if (t >= 0) dot += d0[t] * p_v[t];
                
                dot += Di[i] * p_v[i];
    
                t = i+1;
                if (t < Size) dot += d0[i] * p_v[t];
                t = i+1+Gap;
                if (t < Size) dot += d1[i] * p_v[t];
                t = i+2+Gap;
                if (t < Size) dot += d2[i] * p_v[t];
                t = i+3+Gap;
                if (t < Size) dot += d3[i] * p_v[t];
                
                p_res[i] = dot;
            }
        });
        }
    }
}
