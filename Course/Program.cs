using SparkCL;
using System.Globalization;
using System.Diagnostics;

class Course
{
    static void Main(string[] args)
    {
        Core.Init();
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

        var sw = new Stopwatch();

        var task = new TaskRect4x5();

        for (int g = 0; g < 3; g++)
        {
            var prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            var (ans, iters) = prob.SolveBiCGStabMkl();
            sw.Stop();
            var err = prob.Lebeg2Err(ans);
            Console.WriteLine($"{err} {iters} {sw.ElapsedMilliseconds}мс");
            sw.Reset();

            prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            (ans, iters) = prob.SolveBiCGStab();
            sw.Stop();
            err = prob.Lebeg2Err(ans);
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"{err} {iters} {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();

            Console.WriteLine();
        }

        return;
        for (int g = 0; g < 7; g++)
        {
            
            var prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            var (ans, iters) = prob.SolveBiCGStab();
            sw.Stop();
            var err = prob.Lebeg2Err(ans);
            var (ioTime, kernTime) = Core.MeasureTime();
            ioTime /= (ulong)1e+6;
            kernTime /= (ulong)1e+6;
            Console.WriteLine($"{err} {iters} {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
            sw.Reset();

            for (int i = 0; i < 4; i++)
            {
                prob.femSlae.MeshDouble();
                sw.Start();
                (ans, iters) = prob.SolveBiCGStab();
                sw.Stop();
                (ioTime, kernTime) = Core.MeasureTime();
                ioTime /= (ulong)1e+6;
                kernTime /= (ulong)1e+6;
                Console.WriteLine($"{err} {iters} {sw.ElapsedMilliseconds}мс: {kernTime}мс + {ioTime}мс");
                sw.Reset();
            }
            
            Console.WriteLine();
    
            prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            (ans, iters) = prob.SolveBiCGStabMkl();
            sw.Stop();
            Console.WriteLine($"{err} {iters} {sw.ElapsedMilliseconds}мс");
            sw.Reset();
    
            for (int i = 0; i < 4; i++)
            {
                prob.femSlae.MeshDouble();
                sw.Start();
                (ans, iters) = prob.SolveBiCGStabMkl();
                sw.Stop();
                Console.WriteLine($"{err} {iters} {sw.ElapsedMilliseconds}мс");
                sw.Reset();
            }
            
            Console.WriteLine();
            Console.WriteLine();
        }
    }
}
