using Real = double;

using System.Globalization;

using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Running;

namespace Benchmarks;

[SimpleJob(RuntimeMoniker.Net90)]
public class BicgStabPure
{
    ProblemLine prob;
    
    [Benchmark]
    public void Primary()
    {
        prob.SolveBiCGStabPure();
    }
    
    [GlobalSetup]
    public void Setup()
    {
        // TODO: preallocate tmp arrays
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var task = new TaskRect4x5();
        prob = new ProblemLine(task, "../../../../../../../InputRect4x5");
    }
}
