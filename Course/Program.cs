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

        var prob = new ProblemLine(task, "../../../InputRect4x5");
        sw.Start();
        prob.SolveMCG();

        for (int i = 0; i < 4; i++)
        {
            prob.femSlae.MeshDouble();
            prob.SolveMCG();
        }
        sw.Stop();
        var timeOcl = sw.ElapsedMilliseconds;
        sw.Reset();
        Console.WriteLine();

        prob = new ProblemLine(task, "../../../InputRect4x5");
        sw.Start();
        prob.SolveMcgMkl();

        for (int i = 0; i < 4; i++)
        {
            prob.femSlae.MeshDouble();
            prob.SolveMcgMkl();
        }
        sw.Stop();
        var timeMkl = sw.ElapsedMilliseconds;
        sw.Reset();
        Console.WriteLine($"OCL time: {timeOcl}");
        Console.WriteLine($"Mkl time: {timeMkl}");
    }
}
