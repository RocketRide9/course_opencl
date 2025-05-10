using System.Text.Json;
using Real = double;
using SparkAlgos;

class ProblemLine {
    ProblemParams _problemParams;
    RefineParams _refineParams;
    RectMesh _mesh;
    ComputationalDomain _computationalDomain;
    BoundaryCondition[] _boundaryConditions;

    public Slae2 femSlae;
    public FEMSlaeBuilder slaeBuilder;

    int[] XMonitor = [];
    int[] YMonitor = [];

    TaskFuncs _funcs;

    // void Repurpose (TaskFuncs taskFunctions, string taskFolder)
    // {
    //     computationalDomain = ReadDomains(taskFolder);
    //     boundaryConditions = ReadConditions(taskFolder);

    //     var mesh = new RectMesh(
    //         computationalDomain.xAxis,
    //         computationalDomain.yAxis,
    //         computationalDomain.subDomains,
    //         boundaryConditions
    //     );

    //     femSlae = new FEMSlae(mesh, taskFunctions, refineParams);
    // }

    // folder - директория с условиями задачи
    public ProblemLine(TaskFuncs taskFunctions, string taskFolder)
    {
        _funcs = taskFunctions;

        var json = File.ReadAllText("ProblemParams.json");
        _problemParams = JsonSerializer.Deserialize<ProblemParams>(json)!;

        json = File.ReadAllText(Path.Combine(taskFolder, "RefineParams.json"));
        _refineParams = JsonSerializer.Deserialize<RefineParams>(json)!;

        _computationalDomain = ReadDomains(taskFolder);
        _boundaryConditions = ReadConditions(taskFolder);

        _mesh = new RectMesh(
            _computationalDomain.xAxis,
            _computationalDomain.yAxis,
            _computationalDomain.subDomains,
            _boundaryConditions
        );

        _mesh.Refine(_refineParams);

        slaeBuilder = new FEMSlaeBuilder(_mesh, taskFunctions);
        femSlae = slaeBuilder.Build();
    }

    void MeshRefine(RefineParams refineParams)
    {
        _mesh.Refine(refineParams);
    }

    public void MeshDouble()
    {
        _mesh.RefineDiv2();
    }


    public Real AnswerAt (Real x, Real y)
    {
        var num = _mesh.GetSubdomNumAtPoint(x, y);
        if (num.HasValue)
        {
            return _funcs.Answer(num.Value, x, y);
        }
        else
        {
            return 0;
        }
    }

    public Real ResultAt(Span<Real> q, Real x, Real y)
    {
        var X = _mesh.X;
        var Y = _mesh.Y;
        Real result = 0;

        var (xi, yi) = _mesh.GetElCoordsAtPoint(x, y);

        Real hx = X[xi + 1] - X[xi];
        Real hy = Y[yi + 1] - Y[yi];

        var subdom = _mesh.GetSubdomNumAtElCoords(xi, yi);
        if (subdom.HasValue)
        {
            Span<int> m = stackalloc int[4];
            m[0] = yi * X.Length + xi;
            m[1] = m[0] + 1;
            m[2] = (yi + 1) * X.Length + xi;
            m[3] = m[2] + 1;

            result =
            ( q[m[0]] * (X[xi + 1] - x) * (Y[yi + 1] - y)
            + q[m[1]] * (x - X[xi])     * (Y[yi + 1] - y)
            + q[m[2]] * (X[xi + 1] - x) * (y - Y[yi])
            + q[m[3]] * (x - X[xi])     * (y - Y[yi])) /hx/hy;
        }
        return result;
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
    public Real Lebeg2Err (Span<Real> q)
    {
        var mesh = _mesh;
        Real sum = 0;
        for (int yi = 0; yi < mesh.Y.Length - 1; yi++)
        {
            for (int xi = 0; xi < mesh.X.Length - 1; xi++)
            {
                Real x0 = mesh.X[xi];
                Real y0 = mesh.Y[yi];
                Real hx = mesh.X[xi + 1] - x0;
                Real hy = mesh.Y[yi + 1] - y0;
                var subdom = _mesh.GetSubdomNumAtElCoords(xi, yi);

                if (subdom.HasValue)
                {
                    Real u = ResultAt(q, (Real)(x0 + hx / 2d), (Real)(y0 + hy / 2d));
                    Real u_true = AnswerAt((Real)(x0 + hx / 2d), (Real)(y0 + hy / 2d));
                    Real func = u_true - u;
                    sum += hx * hy * func * func;
                }
            }
        }
        return (Real)Math.Sqrt(sum);
    }

    public void Serialize() => femSlae.Serialize();

    public (SparkCL.Accessor<Real> ans, int iters, Real rr) SolveBiCGStab ()
    {
        var x0 = Enumerable.Repeat((Real)0, femSlae.B.Length).ToArray();
        var slae = femSlae;
        var solver = new BicgStab(
            slae.Mat,
            slae.Di,
            slae.B,
            slae.Ia,
            slae.Ja,
            x0,
            _problemParams.maxIter,
            _problemParams.eps
        );
        var (ans, rr, _, iter) = solver.Solve();

        return (ans, iter, rr);
    }

    public (Real[] ans, int iters, Real rr) SolveBiCGStabMkl ()
    {
        Real[] x0 = [.. Enumerable.Repeat((Real)0, femSlae.B.Length)];
        var slae = femSlae;
        var solver = new BiCGStabMkl(
            slae.Mat,
            slae.Di,
            slae.B,
            slae.Ia,
            slae.Ja,
            x0,
            _problemParams.maxIter,
            _problemParams.eps
        );
        var (ans, rr, _, iter) = solver.Solve();

        return (ans, iter, rr);
    }

    public (Real[] ans, int iters, Real rr) SolveBiCGStabPure ()
    {
        Real[] x0 = [.. Enumerable.Repeat((Real)0, femSlae.B.Length)];
        var solver = new BiCGStabPure(
            _problemParams.maxIter,
            _problemParams.eps
        );
        
        var (rr, _, iter) = solver.Solve(femSlae, x0);

        return (x0, iter, rr);
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
