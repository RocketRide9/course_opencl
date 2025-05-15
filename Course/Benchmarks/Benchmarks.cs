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
    Slae2 slae;

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
        slae.Mat,
        slae.Di,
        slae.Ia,
        slae.Ja,
        slae.Di.Length,
        x0,
        res
    );
    
    [Benchmark]
    public void Hybrid() => MSRMulHybrid(
        slae.Mat,
        slae.Di,
        slae.Ia,
        slae.Ja,
        slae.Di.Length,
        x0,
        res
    );
    
    [Benchmark]
    public void Span() => MSRMulSpans(
        slae.Mat,
        slae.Di,
        slae.Ia,
        slae.Ja,
        slae.Di.Length,
        x0,
        res
    );
    
    [Benchmark]
    public void SpanSlice() => MSRMulSpansSlice(
        slae.Mat,
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

        x0 = [.. Enumerable.Repeat((Real)0.7, prob.femSlae.B.Length)];
        res = [.. Enumerable.Repeat((Real)0.7, prob.femSlae.B.Length)];

        slae = prob.femSlae;
        
        Console.WriteLine($"Размерность матрицы: {x0.Length}, mat.Length: {slae.Mat.Length}");
        
        // var (rr, _, iter) = solver.Solve(prob.femSlae, x0);
    }
}
