using SparkCL;
using System.Globalization;

class Course
{
    static void Main(string[] args)
    {
        Core.Init();
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

        var task = new TaskBook();
        var prob = new ProblemLine(task, "../../../InputBook");
        prob.SolveMCG();
    }
}
