using SparkCL;
using System.Globalization;
using System.Diagnostics;
using Real = double;

using BenchmarkDotNet.Running;

class Course
{
    static void Main(string[] args)
    {
        // Core.Init();
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

        // var summaries = BenchmarkSwitcher.FromAssembly(typeof(Benchmarks.BicgStabPure).Assembly).RunAll();
        SolveAndExportSomeSlae();
        // TestConvergence();
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
            prob.MeshDouble();
        }
    }

    static void SolveAndExportSomeSlae()
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");
        prob.slaeBuilder.GlobalMatrixImpl = GlobalMatrixImplType.Host;

        var sw = new Stopwatch();
        sw.Start();
        var (ans, iters, rr) = prob.SolveBiCGStabPure();
        sw.Stop();
        var err = prob.Lebeg2Err(ans.AsSpan());
        Console.WriteLine($"(iters {iters}) (discrep {rr}) (time {sw.ElapsedMilliseconds}ms)");
        sw.Reset();

        prob.Serialize();
    }

    static void TestConvergence()
    {
        var sw_bicg = new Stopwatch();
        var sw_glob = new Stopwatch();
        var task = new TaskRect4x5();

        void _TestBicgOclBatch(ProblemLine prob)
        {
            Console.WriteLine("n sw_glob err iters discrep total_time(ms) kerns(ms) pci_io(ms)");
            for (int i = 0; i < 4; i++)
            {
                // TODO: подразумевается подсчёт времени, потраченного на сборку матрицы
                // но считается дробление, которое включает в себя ещё другие операции
                sw_glob.Start();
                prob.MeshDouble();
                sw_glob.Stop();

                sw_bicg.Start();
                Core.ResetTime();
                var (ans, iters, rr) = prob.SolveBiCGStab();
                sw_bicg.Stop();
                var err = prob.Lebeg2Err(ans.AsSpan());

                var (ioTime, kernTime) = Core.MeasureTime();
                ioTime /= (ulong)1e+6;
                kernTime /= (ulong)1e+6;
                Console.Write($"{prob.femSlae.B.Length} {sw_glob.ElapsedMilliseconds} ");
                Console.WriteLine($"{err} {iters} {rr} {sw_bicg.ElapsedMilliseconds} {kernTime} {ioTime}");

                sw_bicg.Reset();
                sw_glob.Reset();
            }
        }
        void _TestBicgHostBatch(ProblemLine prob)
        {
            Console.WriteLine("n sw_glob err iters discrep total_time(ms)");
            for (int i = 0; i < 4; i++)
            {
                sw_glob.Start();
                prob.MeshDouble();
                sw_glob.Stop();

                sw_bicg.Start();
                var (ans, iters, rr) = prob.SolveBiCGStabPure();
                sw_bicg.Stop();
                var err = prob.Lebeg2Err(ans.AsSpan());

                Console.Write($"{prob.femSlae.B.Length} {sw_glob.ElapsedMilliseconds} ");
                Console.WriteLine($"{err} {iters} {rr} {sw_bicg.ElapsedMilliseconds}");

                sw_bicg.Reset();
                sw_glob.Reset();
            }
        }

        Console.WriteLine("Сборка на хосте: ");
        for (int g = 0; g < 2; g++)
        {
            ProblemLine prob;

            prob = new ProblemLine(task, "../../../InputRect4x5");
            _TestBicgOclBatch(prob);
            Console.WriteLine();

#if false
            prob = new ProblemLine(task, "../../../InputRect4x5");
            _TestBicgHostBatch(prob);
            Console.WriteLine();
#endif

            Console.WriteLine();
        }
        Console.WriteLine("Сборка на GPU: ");
        for (int g = 0; g < 2; g++)
        {
            ProblemLine prob;

            prob = new ProblemLine(task, "../../../InputRect4x5");
            prob.slaeBuilder.GlobalMatrixImpl = GlobalMatrixImplType.OpenCL;
            _TestBicgOclBatch(prob);
            Console.WriteLine();

#if false
            prob = new ProblemLine(task, "../../../InputRect4x5");
            prob.femSlae.GlobalMatrixImpl = GlobalMatrixImplType.OpenCL;
            _TestBicgHostBatch(prob);
            Console.WriteLine();
#endif

            Console.WriteLine();
        }
    }
}
