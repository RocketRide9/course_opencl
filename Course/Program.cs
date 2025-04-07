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

        var sw = new Stopwatch();

        var task = new TaskRect4x5();

#if false
        for (int g = 0; g < 3; g++)
        {
            var prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStabMkl();
            sw.Stop();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"{err} {iters} discrep: {rr}  {sw.ElapsedMilliseconds}мс");
            sw.Reset();

            prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            (ans, iters, rr) = prob.SolveBiCGStab();
            sw.Stop();
            err = prob.Lebeg2Err(ans);
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"{err} {iters} discrep: {rr}  {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();

            Console.WriteLine();
        }
        return;
#endif
        for (int g = 0; g < 2; g++)
        {
            ProblemLine prob;
            Real err;
            prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            var (ans, iters, rr) = prob.SolveBiCGStab();
            sw.Stop();
            err = prob.Lebeg2Err(ans.AsSpan());
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"{err} {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();

            for (int i = 0; i < 4; i++)
            {
                prob.femSlae.MeshDouble();
                sw.Start();
                (ans, iters, rr) = prob.SolveBiCGStab();
                sw.Stop();
                err = prob.Lebeg2Err(ans.AsSpan());
                (ioTime, kernTime) = Core.MeasureTime();
                ioTime /= (ulong)1e+6;
                kernTime /= (ulong)1e+6;
                Console.WriteLine($"{err} {iters} (discrep: {rr}) {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
                sw.Reset();
            }
            
            Console.WriteLine();
#if false
            prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            var (ans1, iters1, rr1) = prob.SolveBiCGStabMkl();
            sw.Stop();
            err = prob.Lebeg2Err(ans1.AsSpan());
            Console.WriteLine($"{err} {iters1} (discrep: {rr1}) {sw.ElapsedMilliseconds}мс");
            sw.Reset();
    
            for (int i = 0; i < 4; i++)
            {
                prob.femSlae.MeshDouble();
                sw.Start();
                (ans1, iters1, rr1) = prob.SolveBiCGStabMkl();
                sw.Stop();
                err = prob.Lebeg2Err(ans1.AsSpan());
                Console.WriteLine($"{err} {iters1} (discrep: {rr1})  {sw.ElapsedMilliseconds}мс");
                sw.Reset();
            }
            
            Console.WriteLine();
            Console.WriteLine();
#endif
        }
    }
}
