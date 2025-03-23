using System.Text.Json;
using Real = float;
using SparkAlgos;

class ProblemLine {
    ProblemParams problemParams;
    RefineParams refineParams;
    ComputationalDomain computationalDomain;
    BoundaryCondition[] boundaryConditions;
    public FEMSlae femSlae;
    
    int[] XMonitor = [];
    int[] YMonitor = [];
    
    TaskFuncs _funcs;

    void Repurpose (TaskFuncs taskFunctions, string taskFolder)
    {
        computationalDomain = ReadDomains(taskFolder);
        boundaryConditions = ReadConditions(taskFolder);

        var mesh = new RectMesh(
            computationalDomain.xAxis,
            computationalDomain.yAxis,
            computationalDomain.subDomains,
            boundaryConditions
        );

        femSlae = new FEMSlae(mesh, taskFunctions, refineParams);
    }
    
    // folder - директория с условиями задачи
    public ProblemLine(TaskFuncs taskFunctions, string taskFolder)
    {
        var json = File.ReadAllText("ProblemParams.json");
        problemParams = JsonSerializer.Deserialize<ProblemParams>(json)!;
        
        json = File.ReadAllText(Path.Combine(taskFolder, "RefineParams.json"));
        refineParams = JsonSerializer.Deserialize<RefineParams>(json)!;

        _funcs = taskFunctions;
        Repurpose(taskFunctions, taskFolder);
    }
    
    // сохранить узлы текущего разбиения как узлы наблюдения
    // void MonitorNodesFix()
    // {
    //     // в результате выполнения этой функции у мониторных и текущих
    //     // узлов будет соотношение 1к1
    //     XMonitor = Enumerable.Range(0, slae.mesh.X.Count).ToArray();
    //     YMonitor = Enumerable.Range(0, slae.mesh.Y.Count).ToArray();
    // }
    
    /* Взятие нормы погрешности в пространстве Лебега 2.
        Интеграл считается методом прямоугольников  */
    Real Lebeg2Err (SparkOCL.Array<Real> q)
    {
        var mesh = femSlae.Mesh;
        Real sum = 0;
        for (int yi = 0; yi < mesh.Y.Count - 1; yi++)
        {
            for (int xi = 0; xi < mesh.X.Count - 1; xi++)
            {
                Real x0 = mesh.X[xi];
                Real y0 = mesh.Y[yi];
                Real hx = mesh.X[xi + 1] - x0;
                Real hy = mesh.Y[yi + 1] - y0;
                var subdom = femSlae.GetSubdomNumAtElCoords(xi, yi);

                if (subdom.HasValue)
                {
                    Real u = femSlae.ResultAt(q, (Real)(x0 + hx / 2d), (Real)(y0 + hy / 2d));
                    Real u_true = femSlae.AnswerAt((Real)(x0 + hx / 2d), (Real)(y0 + hy / 2d));
                    Real func = u_true - u;
                    sum += hx * hy * func * func;
                }
            }
        }
        return (Real)Math.Sqrt(sum);
    }
    
    public void SolveMCG ()
    {
        var x0 = new SparkOCL.Array<Real>(Enumerable.Repeat((Real)0, femSlae.Slae.B.Count).ToArray());
        var slae = femSlae.Slae;
        var solver = new BicgStab(
            slae.Mat,
            slae.Di,
            slae.B,
            slae.Ia,
            slae.Ja,
            x0,
            problemParams.maxIter,
            problemParams.eps
        );
        var info = solver.Solve();
        x0 = info.ans;

        Real max_err = 0;

        var mesh = femSlae.Mesh;
        int i = 0;
        for (int row = 0; row < mesh.Y.Count; row++)
        {
            for (int col = 0; col < mesh.X.Count; col++)
            {
                var ans = femSlae.AnswerAt(mesh.X[col], mesh.Y[row]);
                var err = Math.Abs(ans - x0[i]);
                if (err > max_err)
                {
                    max_err = err;
                }
                i++;
            }
        }
        //Console.WriteLine($"Max error: {max_err}");
        //Console.WriteLine($"Info: {info}");
        Console.Write($"Error: {Lebeg2Err(x0)}, Iters: {info.iter}, ");
        var (ioTime, kernTime) = SparkCL.Core.MeasureTime();
        Console.WriteLine($"IO: {ioTime / 1e6}мс, Kernels: {kernTime / 1e6}мс");
    }

    public void SolveMcgMkl ()
    {
        var x0 = new SparkOCL.Array<Real>(Enumerable.Repeat((Real)0, femSlae.Slae.B.Count).ToArray());
        var slae = femSlae.Slae;
        var solver = new BiCGStabMkl(
            slae.Mat,
            slae.Di,
            slae.B,
            slae.Ia,
            slae.Ja,
            x0,
            problemParams.maxIter,
            problemParams.eps
        );
        var info = solver.Solve();
        x0 = info.ans;

        Real max_err = 0;

        var mesh = femSlae.Mesh;
        int i = 0;
        for (int row = 0; row < mesh.Y.Count; row++)
        {
            for (int col = 0; col < mesh.X.Count; col++)
            {
                var ans = femSlae.AnswerAt(mesh.X[col], mesh.Y[row]);
                var err = Math.Abs(ans - x0[i]);
                if (err > max_err)
                {
                    max_err = err;
                }
                i++;
            }
        }
        //Console.WriteLine($"Max error: {max_err}");
        //Console.WriteLine($"Info: {info}");
        Console.WriteLine($"Error: {Lebeg2Err(x0)}, Iters: {info.iter}");
    }

    static ComputationalDomain ReadDomains(string taskFolder)
    {
        var file = new StreamReader(Path.Combine(taskFolder, "ComputationalDomain.txt"));
        ComputationalDomain res;
        
        res.xAxis = file.ReadLine()!.Split().Select(Real.Parse).ToArray();
        
        res.yAxis = file.ReadLine()!.Split().Select(Real.Parse).ToArray();
        
        var domains_num = uint.Parse(file.ReadLine()!.Trim());
        res.subDomains = new Subdomain[domains_num];
        for (int i = 0; i < domains_num; i++)
        {
            var parts = file.ReadLine()!.Trim().Split().Select(int.Parse).ToArray();
            res.subDomains[i] = new Subdomain
            {
                Num = parts[0] - 1,
                X1 = parts[1] - 1,
                X2 = parts[2] - 1,
                Y1 = parts[3] - 1,
                Y2 = parts[4] - 1
            };
        }

        return res;
    }

    static BoundaryCondition[] ReadConditions(string taskFolder)
    {
        var file = new StreamReader(Path.Combine(taskFolder, "BoundaryConditions.txt"));
        
        var condsNum = uint.Parse(file.ReadLine()!.Trim());
        var res = new BoundaryCondition[condsNum];
        
        for (int i = 0; i < condsNum; i++)
        {
            var numbers = file.ReadLine()!.Trim().Split().Select(int.Parse).ToArray();
            res[i] = new BoundaryCondition
            {
                Num = numbers[0] - 1,
                Type = numbers[1],
                X1 = numbers[2] - 1,
                X2 = numbers[3] - 1,
                Y1 = numbers[4] - 1,
                Y2 = numbers[5] - 1,
            };
        }
        
        return res;
    }
}
