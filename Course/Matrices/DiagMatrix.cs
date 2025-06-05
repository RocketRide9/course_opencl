using Real = double;

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
    
    int Matrix.Size => Di.Length;
    Span<double> Matrix.Di => Di;

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

    void Matrix.Mul(ReadOnlySpan<double> vec, Span<double> res)
    {
        for (int i = 0; i < vec.Length; i++)
        {
            Real dot = 0;
            
            int t = i - 3 - Gap;
            if (t >= 0) dot += Ld3[t] * vec[t];
            t = i - 2 - Gap;
            if (t >= 0) dot += Ld2[t] * vec[t];
            t = i - 1 - Gap;
            if (t >= 0) dot += Ld1[t] * vec[t];
            t = i - 1;
            if (t >= 0) dot += Ld0[t] * vec[t];
            
            dot += Di[i] * vec[i];

            dot += Rd0[i] * vec[i+1];
            dot += Rd1[i] * vec[i+1+Gap];
            dot += Rd2[i] * vec[i+2+Gap];
            dot += Rd3[i] * vec[i+3+Gap];
            
            res[i] = dot;
        }
    }
}
