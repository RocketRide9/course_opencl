using System.Text.Json;
using Real = float;
using SparkAlgos;

class ProblemLine {
    ProblemParams problemParams;
    RefineParams refineParams;
    ComputationalDomain computationalDomain;
    BoundaryCondition[] boundaryConditions;
    FEMSlae femSlae;
    
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
        var json = File.ReadAllText(Path.Combine(taskFolder, "ProblemParams.json"));
        problemParams = JsonSerializer.Deserialize<ProblemParams>(json)!;
        
        json = File.ReadAllText(Path.Combine(taskFolder, "RefineParams.json"));
        refineParams = JsonSerializer.Deserialize<RefineParams>(json)!;

        _funcs = taskFunctions;
        Repurpose(taskFunctions, taskFolder);
    }
    
    
    
    // int? GetSubdomAtNode (int x, int y)
    // {
    //     foreach (var a in computationalDomain.Subdomains)
    //     {
    //         if (x >= IXw[a.X1] && x <= IXw[a.X2] &&
    //             y >= IYw[a.Y1] && y <= IYw[a.Y2]
    //         ) {
    //             return a.Num;
    //         }
    //     }
        
    //     return null;
    // }

    // (int xi, int yi) GetElCoordsAtPoint (Real x, Real y)
    // {
    //     int? xi = null;
    //     int? yi = null;
    //     for (int i = 0; i < X.Count - 1; i++)
    //     {
    //         if (X[i] <= x && x <= X[i+1])
    //         {
    //             xi = i;
    //             break;
    //         }
    //     }
    //     for (int i = 0; i < Y.Count - 1; i++)
    //     {
    //         if (Y[i] <= y && y <= Y[i+1])
    //         {
    //             yi = i;
    //             break;
    //         }
    //     }
        
    //     if (xi == null || yi == null)
    //     {
    //         throw new Exception("Point is outside the domain");
    //     }

    //     return ((int)xi, (int)yi);
    // }
    
    // /* Нахождение значения полученного рещения в вещественной точке */
    // Real ResultAtPoint(Real x, Real y, Real[] q)
    // {
    //     Real result = 0;
    //     if (!(
    //           X[0] <= x && x <= X.Last() &&
    //           Y[0] <= y && y <= Y.Last()
    //     )) {
    //         throw new Exception("Функция вне области вычислений");
    //     }

    //     var (xi, yi) = GetElCoordsAtPoint(x, y);

    //     Real hx = X[xi + 1] - X[xi];
    //     Real hy = Y[yi + 1] - Y[yi];

    //     var sub_dom = GetSubdomNumAtElCoords(xi, yi);
    //     if (sub_dom != null)
    //     {
    //         var m = new int[4];
    //         m[0] = yi * X.Count + xi;
    //         m[1] = m[0] + 1;
    //         m[2] = (yi + 1) * X.Count + xi;
    //         m[3] = m[2] + 1;
            
    //         result = (
    //               q[m[0]] * (X[xi + 1] - x) * (Y[yi + 1] - y)
    //             + q[m[1]] * (x - X[xi])     * (Y[yi + 1] - y)
    //             + q[m[2]] * (X[xi + 1] - x) * (y - Y[yi])
    //             + q[m[3]] * (x - X[xi])     * (y - Y[yi])
    //             ) /hx/hy;
    //     }
    //     return result;
    // }

    // /* Взятие нормы погрешности в пространстве Лебега 2.
    //     Интеграл считается методом прямоугольников  */
    // Real IntegrateLebeg2 (Real[] q)
    // {
    //     Real sum = 0;
    //     for (int yi = 0; yi < Y.Count - 1; yi++)
    //     {
    //         for (int xi = 0; xi < X.Count - 1; xi++)
    //         {
    //             Real x0 = X[xi];
    //             Real y0 = Y[yi];
    //             Real hx = X[xi + 1] - x0;
    //             Real hy = Y[yi + 1] - y0;
    //             var subdom = GetSubdomNumAtElCoords(xi, yi);

    //             if (subdom != null)
    //             {
    //                 Real u = ResultAtPoint(x0 + hx / 2d, y0 + hy / 2d, q);
    //                 Real u_true = funcs.Answer((int)subdom, x0 + hx / 2d, y0 + hy / 2d);
    //                 Real func = u_true - u;
    //                 sum += hx * hy * func * func;
    //             }
    //         }
    //     }
    //     return Math.Sqrt(sum);
    // }
    
    // сохранить узлы текущего разбиения как узлы наблюдения
    // void MonitorNodesFix()
    // {
    //     // в результате выполнения этой функции у мониторных и текущих
    //     // узлов будет соотношение 1к1
    //     XMonitor = Enumerable.Range(0, slae.mesh.X.Count).ToArray();
    //     YMonitor = Enumerable.Range(0, slae.mesh.Y.Count).ToArray();
    // }
    
    Real FirstStepSize(Real stretch, int seg_count, Real gap)
    {
        Real sum;
        if (stretch != 1d)
        {
            sum = (Real)((1 - Math.Pow(stretch, seg_count)) / ((Real)1 - stretch));
        } else {
            sum = seg_count;
        }

        return gap / sum;
    }
    
    public void SolveMCG ()
    {
        var x0 = new SparkCL.Memory<float>(Enumerable.Repeat(0f, femSlae.slae.B.Count).ToArray());
        var slae = femSlae.slae;
        var solver = new BicgStab(
            slae.Mat,
            slae.Di,
            slae.B,
            slae.Ia,
            slae.Ja,
            x0,
            problemParams.maxIter
        );
        var info = solver.Solve();
        solver.X.Read();

        var max_err = 0d;

        var mesh = femSlae.mesh;
        int i = 0;
        for (int row = 0; row < mesh.Y.Count; row++)
        {
            for (int col = 0; col < mesh.X.Count; col++)
            {
                var ans = femSlae.AnswerAt(mesh.X[col], mesh.Y[row]);
                var err = Math.Abs(ans - solver.X[i]);
                if (err > max_err)
                {
                    max_err = err;
                }
                i++;
            }
        }
        Console.WriteLine($"Max error: {max_err}");
    }
    
    ComputationalDomain ReadDomains(string taskFolder)
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
    
    BoundaryCondition[] ReadConditions(string taskFolder)
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
