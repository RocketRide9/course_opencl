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

        for (int g = 0; g < 7; g++)
        {
            
            var prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            var info = prob.SolveMCG();
            sw.Stop();
            Console.WriteLine($"{info.err} {info.iters} {sw.ElapsedMilliseconds}мс: {info.kern}мс + {info.io}мс");
            sw.Reset();

            for (int i = 0; i < 4; i++)
            {
                prob.femSlae.MeshDouble();
                sw.Start();
                info = prob.SolveMCG();
                sw.Stop();
                Console.WriteLine($"{info.err} {info.iters} {sw.ElapsedMilliseconds}мс: {info.kern}мс + {info.io}мс");
                sw.Reset();
            }
            
            Console.WriteLine();
    
            prob = new ProblemLine(task, "../../../InputRect4x5");
            sw.Start();
            var info_mkl = prob.SolveMcgMkl();
            sw.Stop();
            Console.WriteLine($"{info_mkl.err} {info_mkl.iters} {sw.ElapsedMilliseconds}мс");
            sw.Reset();
    
            for (int i = 0; i < 4; i++)
            {
                prob.femSlae.MeshDouble();
                sw.Start();
                info_mkl = prob.SolveMcgMkl();
                sw.Stop();
                Console.WriteLine($"{info_mkl.err} {info_mkl.iters} {sw.ElapsedMilliseconds}мс");
                sw.Reset();
            }
            
            Console.WriteLine();
            Console.WriteLine();
        }
    }
}
