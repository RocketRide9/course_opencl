/*
using Real = double;

using System.Globalization;

using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Running;

namespace Benchmarks;
#if false
[SimpleJob(RuntimeMoniker.Net90)]
public class BenchBicgStabPure
{
    ProblemLine prob;
    BiCGStabPure solver;
    Real[] x0 = [];
    Real[] res = [];
    Slae2 slae;
    
    [Benchmark]
    public void Primary()
    {
        Array.Fill(x0, 0);
        solver.Solve(slae.AsRef(), x0);
    }
    
    [GlobalSetup]
    public void Setup()
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../../../../../InputRect4x5");

        x0 = [.. Enumerable.Repeat((Real)0.7, prob.femSlae.B.Length)];
        res = [.. Enumerable.Repeat((Real)0.7, prob.femSlae.B.Length)];

        slae = prob.femSlae;

        solver = new BiCGStabPure(prob.ProblemParams.maxIter, prob.ProblemParams.eps);
        solver.AllocateTemps(x0.Length);
        
        Console.WriteLine($"Размерность матрицы: {x0.Length}, mat.Length: {slae.Mat.Length}");
    }
}
#endif
[SimpleJob(RuntimeMoniker.Net90)]
public class BenchMsrMul
{
    Real[] x0;
    Real[] res;
    Types.Matrix slae;

    static void MSRMulSpans(
        ReadOnlySpan<Real> mat,
        ReadOnlySpan<Real> di,
        ReadOnlySpan<int> ia,
        ReadOnlySpan<int> ja,
        int n,
        ReadOnlySpan<Real> v,
        Span<Real> res)
    {
        for (int i = 0; i < ia.Length - 1; i++)
        {
            int start = ia[i];
            int stop = ia[i + 1];
            Real dot = di[i] * v[i];
            for (int a = start; a < stop; a++)
            {
                dot += mat[a] * v[ja[a]];
            }
            res[i] = dot;
        }
    }

    static unsafe void MSRMulSpansUnsafe(
        ReadOnlySpan<Real> mat,
        ReadOnlySpan<Real> di,
        ReadOnlySpan<int> ia,
        ReadOnlySpan<int> ja,
        int n,
        ReadOnlySpan<Real> v,
        Span<Real> res)
    {
        fixed(Real* p_mat = mat)
        fixed(Real* p_di = di)
        fixed(int* p_ia = ia)
        fixed(int* p_ja = ja)
        fixed(Real* p_v = v)
        fixed(Real* p_res = res)
        for (int i = 0; i < ia.Length - 1; i++)
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
    }
    
    static void MSRMulHybrid(
        Real[] mat,
        ReadOnlySpan<Real> di,
        ReadOnlySpan<int> ia,
        int[] ja,
        int n,
        Real[] v,
        Span<Real> res)
    {
        for (int i = 0; i < ia.Length - 1; i++)
        {
            int start = ia[i];
            int stop = ia[i + 1];
            Real dot = di[i] * v[i];
            for (int a = start; a < stop; a++)
            {
                dot += mat[a] * v[ja[a]];
            }
            res[i] = dot;
        }
    }
    
    static void MSRMulArrays(
        Real[] mat,
        Real[] di,
        int[] ia,
        int[] ja,
        int n,
        Real[] v,
        Real[] res)
    {
        for (int i = 0; i < ia.Length - 1; i++)
        {
            int start = ia[i];
            int stop = ia[i + 1];
            Real dot = di[i] * v[i];
            for (int a = start; a < stop; a++)
            {
                dot += mat[a] * v[ja[a]];
            }
            res[i] = dot;
        }
    }
    
    static void MSRMulSpansSlice(
        ReadOnlySpan<Real> mat,
        ReadOnlySpan<Real> di,
        ReadOnlySpan<int> ia,
        ReadOnlySpan<int> ja,
        int n,
        ReadOnlySpan<Real> v,
        Span<Real> res)
    {
        for (int i = 0; i < ia.Length - 1; i++)
        {
            int start = ia[i];
            int stop = ia[i + 1];
            var len = stop - start;
            Real dot = di[i] * v[i];
            var row = mat.Slice(start, len);
            var jarow = ja.Slice(start, len);
            for (int a = 0; a < len; a++)
            {
                dot += row[a] * v[jarow[a]];
            }
            res[i] = dot;
        }
    }
    
    
    [Benchmark(Baseline = true)]
    public void Array() => MSRMulArrays(
        slae.Elems,
        slae.Di,
        slae.Ia,
        slae.Ja,
        slae.Di.Length,
        x0,
        res
    );
    
    [Benchmark]
    public void Hybrid() => MSRMulHybrid(
        slae.Elems,
        slae.Di,
        slae.Ia,
        slae.Ja,
        slae.Di.Length,
        x0,
        res
    );
    
    [Benchmark]
    public void Span() => MSRMulSpans(
        slae.Elems,
        slae.Di,
        slae.Ia,
        slae.Ja,
        slae.Di.Length,
        x0,
        res
    );
    
    [Benchmark]
    public void SpanSlice() => MSRMulSpansSlice(
        slae.Elems,
        slae.Di,
        slae.Ia,
        slae.Ja,
        slae.Di.Length,
        x0,
        res
    );

    [Benchmark]
    public void SpanUnsafe() => MSRMulSpansUnsafe(
        slae.Elems,
        slae.Di,
        slae.Ia,
        slae.Ja,
        slae.Di.Length,
        x0,
        res
    );
    

    [GlobalSetup]
    public void Setup()
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../../../../../InputRect4x5");

        x0 = [.. Enumerable.Repeat((Real)0.7, prob.b.Length)];
        res = [.. Enumerable.Repeat((Real)0.7, prob.b.Length)];

        slae = prob.matrix;
        
        Console.WriteLine($"Размерность матрицы: {x0.Length}, mat.Length: {slae.Elems.Length}");
        
        // var (rr, _, iter) = solver.Solve(prob.femSlae, x0);
    }
}
 */
