using SparkCL;
using System.Globalization;

class Course
{
    static void Main(string[] args)
    {
        Core.Init();
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

        var task = new TaskRect4x5();
        var prob = new ProblemLine(task, "../../../InputRect4x5");
        prob.SolveMCG();

        for (int i = 0; i < 4; i++)
        {
            prob.femSlae.MeshDouble();
            prob.SolveMCG();
        }
    }
}
