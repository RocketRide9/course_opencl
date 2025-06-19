#define SPARKCL_COLLECT_TIME

#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif


using System.Globalization;
using System.Diagnostics;

using SparkCL;
using SparkAlgos.SlaeSolver;
using SlaeSolver;

using SlaeBuilder;
using Matrices;


class Course
{
    const string SRC_DIR = "../../../";
    static void Main(string[] args)
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var unixMs = new DateTimeOffset(DateTime.Now).ToUnixTimeMilliseconds();
        Directory.CreateDirectory(SRC_DIR + "measurements");
        Trace.Listeners.Add(new TextWriterTraceListener(new StreamWriter(SRC_DIR + "measurements/" + unixMs + ".txt")));
        Trace.AutoFlush = true;
        
        var sw = Stopwatch.StartNew();
        Core.Init();
        Trace.WriteLine($"SparkCL Init: {sw.ElapsedMilliseconds}ms");

        // var summaries = BenchmarkSwitcher.FromAssembly(typeof(Benchmarks.BenchMsrMul).Assembly).RunAll();
        // var summaries = BenchmarkSwitcher.FromAssembly(typeof(Benchmarks.BenchBicgStabPure).Assembly).RunAll();

        // BuildImplsXMatrices();
        // MatricesXSolvers();
        Table1();
        Table2();

        /*
        Console.WriteLine("Waiting for debugger to attach");
        while (!Debugger.IsAttached)
        {
            Thread.Sleep(100);
        }
        Console.WriteLine("Debugger attached");
        Debugger.Break();
        */
    }
#if SPARKCL_COLLECT_TIME
    static void OCLTraceTime()
    {
        var (ioTime, kernTime) = Core.MeasureTime();
        ioTime /= (ulong)1e+6;
        kernTime /= (ulong)1e+6;
        Trace.WriteLine($"OpenCL Tracing: {kernTime}ms + {ioTime}ms");
        Core.ResetTime();
    }
#endif
    
    static void Table1()
    {
        const int REPEAT_COUNT = 3;
        const int REFINE_COUNT = 6;
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
        prob.MeshRefine(new()
        {
            XSplitCount = [64],
            YSplitCount = [64],
            XStretchRatio = [1],
            YStretchRatio = [1],
        });

        Console.WriteLine("Table 1:");

        for (int r = 0; r < REFINE_COUNT; r++)
        {
            Console.WriteLine($"n = {prob.Mesh.nodesCount}");
            Trace.WriteLine($"n = {prob.Mesh.nodesCount}");
    
            prob.buildType = GlobalMatrixImplType.HostParallel;
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // Element Parallel build
                prob.Build<MsrSlaeBuilder>();
                prob.Build<SymDiagSlaeBuilder>();
                prob.Build<DiagSlaeBuilder>();
                // SolveOCL<BicgStabOCL>(prob);
            }
    
            Trace.WriteLine("");
            
#if true
            prob.buildType = GlobalMatrixImplType.OpenCL;
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // Element OpenCL build
#if SPARKCL_COLLECT_TIME
                Core.ResetTime();
#endif
                prob.Build<MsrSlaeBuilder>();
#if SPARKCL_COLLECT_TIME
                OCLTraceTime();
#endif
                prob.Build<SymDiagSlaeBuilder>();
#if SPARKCL_COLLECT_TIME
                OCLTraceTime();
#endif
                prob.Build<DiagSlaeBuilder>();
#if SPARKCL_COLLECT_TIME
                OCLTraceTime();
#endif
                // SolveOCL<CgmOCL>(prob);
            }
            
            Trace.WriteLine("");
#endif
            prob.buildType = GlobalMatrixImplType.OpenCLV2;
            for (int i = 0; i < REPEAT_COUNT; i++)
            { // Node OpenCL build
#if SPARKCL_COLLECT_TIME
                Core.ResetTime();
#endif
                prob.Build<MsrSlaeBuilder>();
#if SPARKCL_COLLECT_TIME
                OCLTraceTime();
#endif
                // SolveOCL<CgmOCL>(prob);
            }
            
            Trace.WriteLine("");

            prob.MeshDouble();
        }
    }
    
    static void Table2()
    {
        const int REPEAT_COUNT = 3;
        const int REFINE_COUNT = 5;
        
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5", GlobalMatrixImplType.HostParallel);
        prob.MeshRefine(new()
        {
            XSplitCount = [64],
            YSplitCount = [64],
            XStretchRatio = [1],
            YStretchRatio = [1],
        });

        Console.WriteLine("Table 2:");

        for (int r = 0; r < REFINE_COUNT; r++)
        {
            Trace.WriteLine($"n = {prob.Mesh.nodesCount}");
            Console.WriteLine($"n = {prob.Mesh.nodesCount}");
            
            for (int i = 0; i < REPEAT_COUNT; i++)
            {
                prob.Build<SymDiagSlaeBuilder>();
                SolveHost<BicgStabHost>(prob);
                SolveOCL<BicgStabOCL>(prob);
                SolveHost<CgmHost>(prob);
                SolveOCL<CgmOCL>(prob);

                Trace.WriteLine("");

                prob.Build<MsrSlaeBuilder>();
                SolveHost<BicgStabHost>(prob);
                SolveOCL<BicgStabOCL>(prob);

                Trace.WriteLine("");
            }
            
            prob.MeshDouble();
        }
    }
    
    static void MatricesXSolvers()
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
        prob.buildType = GlobalMatrixImplType.OpenCLV2;

        Trace.WriteLine($"n = {prob.Mesh.nodesCount}");

        for (int i = 0; i < 5; i++)
        {
            // prob.Build<DiagSlaeBuilder>();
            // // SolveHost<BicgStabHost>(prob);
            // SolveOCL<BicgStabOCL>(prob);
            // // SolveHost<CgmHost>(prob);
            // SolveOCL<CgmOCL>(prob);
            
            // Trace.WriteLine("");
            
            // prob.Build<SymDiagSlaeBuilder>();
            // // SolveHost<BicgStabHost>(prob);
            // SolveOCL<BicgStabOCL>(prob);
            // // SolveHost<CgmHost>(prob);
            // SolveOCL<CgmOCL>(prob);
            
            // Trace.WriteLine("");
            
            prob.Build<MsrSlaeBuilder>();
            // SolveHost<BicgStabHost>(prob);
            SolveOCL<BicgStabOCL>(prob);
            
            Trace.WriteLine("");
        }
    }
    
    static void VerifyBuilds<T1, T2>()
    where T1: ISlaeBuilder
    where T2: ISlaeBuilder
    {
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
        prob.buildType = GlobalMatrixImplType.HostParallel;

        prob.Build<T1>();
        var matrix1 = prob.matrix;
        var b1 = prob.b;

        prob.Build<T2>();
        var matrix2 = prob.matrix;
        var b2 = prob.b;

        VerifySlaes(matrix1, b1, matrix2, b2);
    }
    
    static void VerifySlaes(Types.Matrix m1, Real[] b1, Types.Matrix m2, Real[] b2)
    {
        if (m1.FlatNonZero().SequenceEqual(m2.FlatNonZero()))
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
            Console.WriteLine(string.Join(',', m1.FlatNonZero()));
            Console.WriteLine(string.Join(',', m2.FlatNonZero()));
        }
    }
    
    static void SolveOCL<T>(ProblemLine prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        Trace.WriteLine("SolveOCL");
        Trace.Indent();
#if SPARKCL_COLLECT_TIME
        Core.ResetTime();
#endif
        var sw = Stopwatch.StartNew();
        var (ans, iters, rr) = prob.SolveOCL<T>();
        sw.Stop();
        var err = prob.Lebeg2Err(ans);
#if SPARKCL_COLLECT_TIME
        var (ioTime, kernTime) = Core.MeasureTime();
        ioTime /= (ulong)1e+6;
        kernTime /= (ulong)1e+6;
#endif
        Console.WriteLine($"(err {err}) (iters {iters}) (discrep: {rr})");
        Trace.Unindent();
        Trace.Write($"Solver total: {sw.ElapsedMilliseconds}мс");
#if SPARKCL_COLLECT_TIME
        Trace.Write($": {kernTime}мс + {ioTime}мс");
#endif
        Trace.WriteLine("");
    }
    
    static void SolveHost<T>(ProblemLine prob)
    where T: SlaeSolver.ISlaeSolver
    {
        Trace.WriteLine("Setting up solver");
        Trace.Indent();
        var sw = Stopwatch.StartNew();
        var (ans, iters, rr) = prob.SolveHost<T>();
        sw.Stop();
        var err = prob.Lebeg2Err(ans);
        Console.WriteLine($"(err {err}) (iters {iters}) (discrep: {rr})");
        Trace.WriteLine($"Solver total: {sw.ElapsedMilliseconds}мс");
        Trace.Unindent();
    }

    static void TestBuildImplsIncreasing()
    {
        const int REPEAT_COUNT = 3;
        const int REFINE_COUNT = 6;
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
        
        
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
        
        
        prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
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

        prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
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
    
    static void BuildImplsXMatrices()
    {
        const int REPEAT_COUNT = 5;
        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");

        Trace.WriteLine($"n = {prob.Mesh.nodesCount}");
        
        prob.buildType = GlobalMatrixImplType.HostParallel;
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // Classic Parallel build
            // prob.Build<MsrSlaeBuilder>();
            prob.Build<SymDiagSlaeBuilder>();
            prob.Build<DiagSlaeBuilder>();
            // SolveOCL<BicgStabOCL>(prob);
        }

        Trace.WriteLine("");

// эта реализация сборки нужна только для проверки правильности
#if false
        prob.buildType = GlobalMatrixImplType.Host;
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // Classic build
            prob.Build<MsrSlaeBuilder>();
            prob.Build<SymDiagSlaeBuilder>();
            prob.Build<DiagSlaeBuilder>();
            // SolveOCL<BicgStabOCL>(prob);
        }
        
        Trace.WriteLine("");
#endif
        
        prob.buildType = GlobalMatrixImplType.OpenCL;
        for (int i = 0; i < REPEAT_COUNT; i++)
        { // Classic OpenCL build
            // prob.Build<MsrSlaeBuilder>();
            prob.Build<DiagSlaeBuilder>();
            prob.Build<SymDiagSlaeBuilder>();
            // SolveOCL<CgmOCL>(prob);
        }
        
        Trace.WriteLine("");
    }

    static void TestConvergence<T>()
    where T: ISlaeBuilder
    {
        var task = new TaskRect4x5();

        int REFINE_COUNT = 4;
        int REPEAT_COUNT = 2;

        void _TestBicgOclBatch(ProblemLine prob)
        {
            Console.WriteLine("n err iters discrep kerns(ms) pci_io(ms)");
            int i = 0;
            while (true)
            {
                prob.Build<T>();

#if SPARKCL_COLLECT_TIME
                Core.ResetTime();
#endif
                var (ans, iters, rr) = prob.SolveOCL<BicgStabOCL>();
                var err = prob.Lebeg2Err(ans);

#if SPARKCL_COLLECT_TIME
                var (ioTime, kernTime) = Core.MeasureTime();
                ioTime /= (ulong)1e+6;
                kernTime /= (ulong)1e+6;
#endif
                Console.Write($"{prob.matrix.Size} {err} {iters} {rr}");
#if SPARKCL_COLLECT_TIME
                Console.WriteLine($" {kernTime} {ioTime}");
#endif
                if (i >= REFINE_COUNT) break;

                prob.MeshDouble();

                i++;
            }
        }
        void _TestBicgHostBatch(ProblemLine prob)
        {
            Console.WriteLine("n err iters discrep");
            int i = 0;
            while (true)
            {
                prob.Build<T>();

                var (ans, iters, rr) = prob.SolveHost<BicgStabHost>();
                var err = prob.Lebeg2Err(ans);

                Console.WriteLine($"{prob.matrix.Size} {err} {iters} {rr}");

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
            prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
            prob.buildType = GlobalMatrixImplType.HostParallel;
            _TestBicgOclBatch(prob);
            Console.WriteLine();
#endif

#if true
            prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
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

            prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
            prob.buildType = GlobalMatrixImplType.OpenCL;
            _TestBicgOclBatch(prob);
            Console.WriteLine();

#if false
            prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
            prob.buildType = GlobalMatrixImplType.OpenCL;
            _TestBicgHostBatch(prob);
            Console.WriteLine();
#endif

            Console.WriteLine();
        }
    }
}
