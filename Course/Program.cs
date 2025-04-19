using SparkCL;
using System.Globalization;
using System.Diagnostics;
using Real = float;

class Course
{
    static void Main(string[] args)
    {
        Core.Init();
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

        // SolveAndExportSomeSlae();
        TestConvergence();
        // TestHostVSOpenCLOnce();

        // TestAtomicAdd();
    }

    static void TestAtomicAdd()
    {
        var zero = new ComputeBuffer<Real>([0]);
        var prog = new Program("Kernels.clcpp");
        var kernAddZeros = prog.GetKernel(
            "add_to_zero",
            new(1024, 1024),
            new(2, 2)
        );
        kernAddZeros.SetArg(0, zero);
        kernAddZeros.Execute();
        var zeroHost = new Real[1];
        zero.ReadTo(zeroHost);

        Console.WriteLine("This is zero: {0}!", zeroHost[0]);
    }
    
    static void TestHostVSOpenCLOnce()
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");
        
        var sw = new Stopwatch();
        { // CSharp pure
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStabPure();
            sw.Stop();
            var err = prob.Lebeg2Err(ans.AsSpan());
            Console.WriteLine($"{err} {iters} (discrep: {rr})  {sw.ElapsedMilliseconds}мс");
            sw.Reset();
        }
        return;
        { // OpenCL
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStab();
            sw.Stop();
            var err = prob.Lebeg2Err(ans.AsSpan());
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"{err} {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();
            prob.femSlae.MeshDouble();
        }
    }
    
    static void SolveAndExportSomeSlae()
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");
        
        var sw = new Stopwatch();
        sw.Start();
        var (ans, iters, rr) = prob.SolveBiCGStabPure();
        sw.Stop();
        var err = prob.Lebeg2Err(ans.AsSpan());
        Console.WriteLine($"{err} {iters} (discrep: {rr})  {sw.ElapsedMilliseconds}мс");
        sw.Reset();

        prob.Serialize();
    }
    
    [deprecate]
    static void TestConvergence()
    {
        var sw = new Stopwatch();
        var task = new TaskRect4x5();
        
        for (int g = 0; g < 2; g++)
        {
            ProblemLine prob;
            prob = new ProblemLine(task, "../../../InputRect4x5");

#if false
            for (int i = 0; i < 5; i++)
            {
                sw.Start();
                var (ans, iters, rr) = prob.SolveBiCGStab();
                sw.Stop();
                var err = prob.Lebeg2Err(ans.AsSpan());
                var (ioTime, kernTime) = Core.MeasureTime();
                ioTime /= (ulong)1e+6;
                kernTime /= (ulong)1e+6;
                Console.WriteLine($"{err} {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
                sw.Reset();
                prob.femSlae.MeshDouble();
            }
            Console.WriteLine();
#endif
            prob = new ProblemLine(task, "../../../InputRect4x5");
            for (int i = 0; i < 5; i++)
            {
                sw.Start();
                var (ans, iters, rr) = prob.SolveBiCGStabPure();
                sw.Stop();
                var err = prob.Lebeg2Err(ans.AsSpan());
                Console.WriteLine($"{err} {iters} (discrep: {rr})  {sw.ElapsedMilliseconds}мс");
                sw.Reset();
                prob.femSlae.MeshDouble();
            }
            Console.WriteLine();
            
            
            Console.WriteLine();
        }
    }
}
