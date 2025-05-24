using SparkCL;
using System.Globalization;
using System.Diagnostics;
using Real = double;

using BenchmarkDotNet.Running;

class Course
{
    static void Main(string[] args)
    {
        Core.Init();
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

        // var summaries = BenchmarkSwitcher.FromAssembly(typeof(Benchmarks.BenchMsrMul).Assembly).RunAll();
        // var summaries = BenchmarkSwitcher.FromAssembly(typeof(Benchmarks.BenchBicgStabPure).Assembly).RunAll();
        // SolveAndExportSomeSlae();
        // TestConvergence();
        // TestParallelBuildV2();
        // TestHostVSOpenCLOnce();
        // TestAtomicAdd();
        // TestBuildImpls();
        
        /*
        Console.WriteLine("Waiting for debugger to attach");
        while (!Debugger.IsAttached)
        {
            Thread.Sleep(100);
        }
        Console.WriteLine("Debugger attached");
        Debugger.Break();
        */
        
        TestBuildImplsIncreasing();
    }

    static void TestParallelBuildV2()
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");

        var sw = new Stopwatch();
        { // Classic build
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStab();
            sw.Stop();
            var err = prob.Lebeg2Err(ans);
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"(err {err}) {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();
        }
        
        { // ParallelV2 build
            prob.buildType = GlobalMatrixImplType.HostV2;
            prob.Rebuild();
            
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStab();
            sw.Stop();
            var err = prob.Lebeg2Err(ans);
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"(err {err}) {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();
        }

        // return;
        { // ParallelOclV2 build
            prob.buildType = GlobalMatrixImplType.OpenCLV2;
            prob.Rebuild();
            
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStab();
            sw.Stop();
            var err = prob.Lebeg2Err(ans);
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"(err {err}) {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();
        }
    }
    
    static void TestAtomicAdd()
    {
        var zero = new ComputeBuffer<Real>([0], BufferFlags.OnDevice);
        var prog = new Program("Kernels.clcpp");
        var kernAddZeros = prog.GetKernel(
            "add_to_zero",
            new(1024, 1024),
            new(2, 2)
        );
        kernAddZeros.SetArg(0, zero);
        kernAddZeros.Execute();
        Span<Real> zeroHost = stackalloc Real[1];
        zero.DeviceReadTo(zeroHost);

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
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс");
            sw.Reset();
        }

        { // OpenCL
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStab();
            sw.Stop();
            var err = prob.Lebeg2Err(ans);
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"(err {err}) {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();
        }
    }

    static void SolveAndExportSomeSlae()
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");
        // TODO: 
        // prob.slaeBuilder.GlobalMatrixImpl = GlobalMatrixImplType.Host;

        for (int i = 0; i < 4; i++)
        {
            var sw = new Stopwatch();
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStabPure();
            sw.Stop();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) (iters {iters}) (discrep {rr}) (time {sw.ElapsedMilliseconds}ms)");
            sw.Reset();
        }
        prob.Serialize();
    }

    static void TestBuildImplsIncreasing()
    {
        const int REPEAT_COUNT = 3;
        const int REFINE_COUNT = 6;
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");
        
        
        Console.WriteLine("Classic Parallel:");
        for (int re = 0; re < REFINE_COUNT; re++)
        {
            prob.MeshDouble();
            Console.WriteLine($"n: {prob.femSlae.B.Length}");
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // Classic Parallel build
                var sw = Stopwatch.StartNew();
                prob.buildType = GlobalMatrixImplType.HostParallel;
                prob.Rebuild();
                sw.Stop();
                Console.WriteLine($"(build time {sw.ElapsedMilliseconds}ms)");
    
                // var (ans, iters, rr) = prob.SolveBiCGStab();
                // var err = prob.Lebeg2Err(ans);
                // Console.WriteLine($"(err {err}) {iters} (discrep: {rr})");
            }
        }
        
        
        prob = new ProblemLine(task, "../../../InputRect4x5");
        Console.WriteLine("OpenCL V2:");
        for (int re = 0; re < REFINE_COUNT; re++)
        {
            prob.MeshDouble();
            Console.WriteLine($"n: {prob.femSlae.B.Length}");
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // ParallelOclV2 build
                var sw = Stopwatch.StartNew();
                prob.buildType = GlobalMatrixImplType.OpenCLV2;
                prob.Rebuild();
                sw.Stop();
                Console.WriteLine($"(build time {sw.ElapsedMilliseconds}ms)");

                // var (ans, iters, rr) = prob.SolveBiCGStab();
                // var err = prob.Lebeg2Err(ans);
                // Console.WriteLine($"(err {err}) {iters} (discrep: {rr})");
            }
        }

        prob = new ProblemLine(task, "../../../InputRect4x5");
        Console.WriteLine("OpenCL:");
        for (int re = 0; re < REFINE_COUNT; re++)
        {
            prob.MeshDouble();
            Console.WriteLine($"n: {prob.femSlae.B.Length}");
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // ParallelOclV2 build
                var sw = Stopwatch.StartNew();
                prob.buildType = GlobalMatrixImplType.OpenCL;
                prob.Rebuild();
                sw.Stop();
                Console.WriteLine($"(build time: {sw.ElapsedMilliseconds}ms)");

                // var (ans, iters, rr) = prob.SolveBiCGStab();
                // var err = prob.Lebeg2Err(ans);
                // Console.WriteLine($"(err {err}) {iters} (discrep: {rr})");
            }
        }
    }
    
    static void TestBuildImpls()
    {
        const int REPEAT_COUNT = 5;
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");

        // норма 1706 итераций на решение слау
        /*
        Console.WriteLine("Classic Parallel:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // Classic Parallel build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.HostParallel;
            prob.Rebuild();
            sw.Stop();
            Console.WriteLine($"(build time {sw.ElapsedMilliseconds}ms)");

            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) {iters} (discrep: {rr})");
        }
        */
        
        Console.WriteLine("Classic:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // Classic Parallel build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.Host;
            prob.Rebuild();
            sw.Stop();
            Console.WriteLine($"(build time {sw.ElapsedMilliseconds}ms)");

            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) (iters {iters}) (discrep {rr})");
        }
        /*
        */
        Console.WriteLine("Host V2:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // ParallelV2 build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.HostV2;
            prob.Rebuild();
            sw.Stop();
            Console.WriteLine($"(build time {sw.ElapsedMilliseconds}ms)");
            
            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) (iters {iters}) (discrep {rr})");
        }
        Console.WriteLine("OpenCL V2:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // ParallelOclV2 build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.OpenCLV2;
            prob.Rebuild();
            sw.Stop();
            Console.WriteLine($"(build time {sw.ElapsedMilliseconds}ms)");
            
            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) (iters {iters}) (discrep: {rr})");
        }

        /*
        Console.WriteLine("OpenCL:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // ParallelOclV2 build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.OpenCL;
            prob.Rebuild();
            sw.Stop();
            Console.WriteLine($"(build time: {sw.ElapsedMilliseconds}ms)");
            
            // var (ans, iters, rr) = prob.SolveBiCGStab();
            // var err = prob.Lebeg2Err(ans);
            // Console.WriteLine($"(err {err}) {iters} (discrep: {rr})");
        }
        */
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
                var err = prob.Lebeg2Err(ans);

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
                var err = prob.Lebeg2Err(ans);

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

#if false
            prob = new ProblemLine(task, "../../../InputRect4x5");
            _TestBicgOclBatch(prob);
            Console.WriteLine();
#endif

#if true
            prob = new ProblemLine(task, "../../../InputRect4x5");
            _TestBicgHostBatch(prob);
            Console.WriteLine();
#endif

            Console.WriteLine();
        }
        return;
        Console.WriteLine("Сборка на GPU: ");
        for (int g = 0; g < 2; g++)
        {
            ProblemLine prob;

            prob = new ProblemLine(task, "../../../InputRect4x5");
            // TODO: 
            // prob.slaeBuilder.GlobalMatrixImpl = GlobalMatrixImplType.OpenCL;
            _TestBicgOclBatch(prob);
            Console.WriteLine();

#if false
            prob = new ProblemLine(task, "../../../InputRect4x5");
            prob.slaeBuilder.GlobalMatrixImpl = GlobalMatrixImplType.OpenCL;
            _TestBicgHostBatch(prob);
            Console.WriteLine();
#endif

            Console.WriteLine();
        }
    }
}
