using Real = double;

using System.Collections.Concurrent;

using Types;

namespace Matrices;

public class DiagMatrix : Matrix
{
    // Left diagonal
    public Real[] Ld3 = [];
    public Real[] Ld2 = [];
    public Real[] Ld1 = [];
    public Real[] Ld0 = [];
    // Основная диагональ
    public Real[] Di = [];
    // Right diagonal
    public Real[] Rd0 = [];
    public Real[] Rd1 = [];
    public Real[] Rd2 = [];
    public Real[] Rd3 = [];

    // Ld0 и Rd0 (*d0) находятся "вплотную" к основной диагонали
    // *d1, *d2, *d3 находятся стоят "вплотную" друг к другу
    // *d1 смещена на Gap элементов от *d0.
    // Например, если они находятся вплотную друг к друг,
    // то Gap == 1.
    public int Gap;
    
    public int Size => Di.Length;
    Span<Real> Matrix.Di => Di;

    public SparkAlgos.Types.Matrix GetComputeMatrix()
    {
        return new SparkAlgos.Matrices.DiagMatrix(new(){
            Ld3 = Ld3,
            Ld2 = Ld2,
            Ld1 = Ld1,
            Ld0 = Ld0,
            Di = Di,
            Rd0 = Rd0,
            Rd1 = Rd1,
            Rd2 = Rd2,
            Rd3 = Rd3,
            Gap = Gap,
        });
    }
    
    public IEnumerable<Real> FlatNonZero()
    {
        for (int i = 0; i < Size; i++)
        {
            int t = i - 3 - Gap;
            if (t >= 0 && Ld3[t] != 0) yield return Ld3[t];
            t = i - 2 - Gap;
            if (t >= 0 && Ld2[t] != 0) yield return Ld2[t];
            t = i - 1 - Gap;
            if (t >= 0 && Ld1[t] != 0) yield return Ld1[t];
            t = i - 1;
            if (t >= 0 && Ld0[t] != 0) yield return Ld0[t];

            yield return Di[i];

            t = i+1;
            if (t < Size && Rd0[i] != 0) yield return Rd0[i];
            t = i+1+Gap;
            if (t < Size && Rd1[i] != 0) yield return Rd1[i];
            t = i+2+Gap;
            if (t < Size && Rd2[i] != 0) yield return Rd2[i];
            t = i+3+Gap;
            if (t < Size && Rd3[i] != 0) yield return Rd3[i];   
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
                if (t >= 0) dot += Ld3[t] * p_v[t];
                t = i - 2 - Gap;
                if (t >= 0) dot += Ld2[t] * p_v[t];
                t = i - 1 - Gap;
                if (t >= 0) dot += Ld1[t] * p_v[t];
                t = i - 1;
                if (t >= 0) dot += Ld0[t] * p_v[t];
                
                dot += Di[i] * p_v[i];
    
                t = i+1;
                if (t < Size) dot += Rd0[i] * p_v[t];
                t = i+1+Gap;
                if (t < Size) dot += Rd1[i] * p_v[t];
                t = i+2+Gap;
                if (t < Size) dot += Rd2[i] * p_v[t];
                t = i+3+Gap;
                if (t < Size) dot += Rd3[i] * p_v[t];
                
                p_res[i] = dot;
            }
        });
        }
    }
}
