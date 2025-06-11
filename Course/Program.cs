#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Globalization;
using System.Diagnostics;

using SparkCL;

using SlaeBuilder;
using Matrices;
using System.Drawing;

class Course
{
    static void Main(string[] args)
    {
        Core.Init();
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

        // var summaries = BenchmarkSwitcher.FromAssembly(typeof(Benchmarks.BenchMsrMul).Assembly).RunAll();
        // var summaries = BenchmarkSwitcher.FromAssembly(typeof(Benchmarks.BenchBicgStabPure).Assembly).RunAll();
        // SolveAndExportSomeSlae();
        // VerifyBuilds<SymDiagSlaeBuilder, DiagSlaeBuilder>();
        // return;
        for (int i = 0; i < 5; i++)
        {
            TestHostVSOpenCLOnce<SymDiagSlaeBuilder>();
        }
        // TestConvergence<DiagSlaeBuilder>();
        // TestConvergence<SymDiagSlaeBuilder>();
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
        
        // TestBuildImplsIncreasing();
    }
    
    static void VerifyBuilds<T1, T2>()
    where T1: ISlaeBuilder
    where T2: ISlaeBuilder
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");
        prob.buildType = GlobalMatrixImplType.HostParallel;

        prob.Build<T1>();
        var matrix1 = prob.matrix;
        var b1 = prob.b;

        prob.Build<T2>();
        var matrix2 = prob.matrix;
        var b2 = prob.b;

        if (matrix1.FlatNonZero().SequenceEqual(matrix2.FlatNonZero()))
        {
            if (b1.SequenceEqual(b2))
            {
                Console.WriteLine("Слау одинаковые");   
            } else {
                Console.WriteLine("Правый части разные");
                Console.WriteLine(string.Join(',',  b1));
                Console.WriteLine(string.Join(',',  b2));    
            }
        } else {
            Console.WriteLine("Матрицы разные");
            Console.WriteLine(string.Join(',', matrix1.FlatNonZero()));
            Console.WriteLine(string.Join(',', matrix2.FlatNonZero()));
        }
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
            prob.Build<MsrSlaeBuilder>();
            
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
            prob.Build<MsrSlaeBuilder>();
            
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
        var prog = new ComputeProgram("Kernels.clcpp");
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

    static void TestHostVSOpenCLOnce<T>()
    where T: ISlaeBuilder
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");

        prob.Build<T>();

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

    /*
    static void SolveAndExportSomeSlae()
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");

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
    */

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
            Console.WriteLine($"n: {prob.b.Length}");
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // Classic Parallel build
                var sw = Stopwatch.StartNew();
                prob.buildType = GlobalMatrixImplType.HostParallel;
                prob.Build<MsrSlaeBuilder>();
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
            Console.WriteLine($"n: {prob.b.Length}");
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // ParallelOclV2 build
                var sw = Stopwatch.StartNew();
                prob.buildType = GlobalMatrixImplType.OpenCLV2;
                prob.Build<MsrSlaeBuilder>();
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
            Console.WriteLine($"n: {prob.b.Length}");
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // ParallelOclV2 build
                var sw = Stopwatch.StartNew();
                prob.buildType = GlobalMatrixImplType.OpenCL;
                prob.Build<MsrSlaeBuilder>();
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

        // норма 1706 итераций на решение слау n=?
        // 385 итерация на решение слау n = 2^16
        
        Console.WriteLine("Classic Parallel:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // Classic Parallel build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.HostParallel;
            prob.Build<MsrSlaeBuilder>();
            sw.Stop();
            Console.WriteLine($"(combined build time {sw.ElapsedMilliseconds}ms)");

            /*
            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) {iters} (discrep: {rr})");
            */
        }
        
        Console.WriteLine("Classic:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // Classic Parallel build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.Host;
            prob.Build<MsrSlaeBuilder>();
            sw.Stop();
            Console.WriteLine($"(combined build time {sw.ElapsedMilliseconds}ms)");
            /* 
            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) (iters {iters}) (discrep {rr})");
            */
        }
        /*
        Console.WriteLine("Host V2:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // ParallelV2 build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.HostV2;
            prob.Build<MsrSlaeBuilder>();
            sw.Stop();
            Console.WriteLine($"(combined build time {sw.ElapsedMilliseconds}ms)");
            
            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) (iters {iters}) (discrep {rr})");
        }
        */
        /*
        Console.WriteLine("OpenCL V2:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // ParallelOclV2 build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.OpenCLV2;
            prob.Build<MsrSlaeBuilder>();
            sw.Stop();
            Console.WriteLine($"(combined build time {sw.ElapsedMilliseconds}ms)");
            
            
            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) (iters {iters}) (discrep: {rr})");
            
        }

        
        Console.WriteLine("OpenCL:");
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // ParallelOclV2 build
            var sw = Stopwatch.StartNew();
            prob.buildType = GlobalMatrixImplType.OpenCL;
            prob.Build<MsrSlaeBuilder>();
            sw.Stop();
            Console.WriteLine($"(build time: {sw.ElapsedMilliseconds}ms)");
            
            var (ans, iters, rr) = prob.SolveBiCGStab();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"(err {err}) {iters} (discrep: {rr})");
        }
        */
    }

    static void TestConvergence<T>()
    where T: ISlaeBuilder
    {
        var sw_bicg = new Stopwatch();
        var sw_glob = new Stopwatch();
        var task = new TaskRect4x5();

        int REFINE_COUNT = 4;
        int REPEAT_COUNT = 2;

        void _TestBicgOclBatch(ProblemLine prob)
        {
            Console.WriteLine("n sw_glob err iters discrep total_time(ms) kerns(ms) pci_io(ms)");
            int i = 0;
            while (true)
            {
                sw_glob.Start();
                prob.Build<T>();
                sw_glob.Stop();

                sw_bicg.Start();
                Core.ResetTime();
                var (ans, iters, rr) = prob.SolveBiCGStab();
                sw_bicg.Stop();
                var err = prob.Lebeg2Err(ans);

                var (ioTime, kernTime) = Core.MeasureTime();
                ioTime /= (ulong)1e+6;
                kernTime /= (ulong)1e+6;
                Console.Write($"{prob.matrix.Size} {sw_glob.ElapsedMilliseconds} ");
                Console.WriteLine($"{err} {iters} {rr} {sw_bicg.ElapsedMilliseconds} {kernTime} {ioTime}");

                sw_bicg.Reset();
                sw_glob.Reset();

                if (i >= REFINE_COUNT) break;

                prob.MeshDouble();

                i++;
            }
        }
        void _TestBicgHostBatch(ProblemLine prob)
        {
            Console.WriteLine("n sw_glob err iters discrep total_time(ms)");
            int i = 0;
            while (true)
            {
                sw_glob.Start();
                prob.Build<T>();
                sw_glob.Stop();

                sw_bicg.Start();
                var (ans, iters, rr) = prob.SolveBiCGStabPure();
                sw_bicg.Stop();
                var err = prob.Lebeg2Err(ans);

                Console.Write($"{prob.matrix.Size} {sw_glob.ElapsedMilliseconds} ");
                Console.WriteLine($"{err} {iters} {rr} {sw_bicg.ElapsedMilliseconds}");

                sw_bicg.Reset();
                sw_glob.Reset();

                if (i >= REFINE_COUNT) break;

                prob.MeshDouble();

                i++;
            }
        }

        Console.WriteLine("Сборка на хосте: ");
        for (int g = 0; g < REPEAT_COUNT; g++)
        {
            ProblemLine prob;

#if true
            prob = new ProblemLine(task, "../../../InputRect4x5");
            prob.buildType = GlobalMatrixImplType.HostParallel;
            _TestBicgOclBatch(prob);
            Console.WriteLine();
#endif

#if true
            prob = new ProblemLine(task, "../../../InputRect4x5");
            prob.buildType = GlobalMatrixImplType.HostParallel;
            _TestBicgHostBatch(prob);
            Console.WriteLine();
#endif

            Console.WriteLine();
        }
        return;
        Console.WriteLine("Сборка на GPU: ");
        for (int g = 0; g < REPEAT_COUNT; g++)
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
