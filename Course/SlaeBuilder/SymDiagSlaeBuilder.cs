using Real = double;

using System.Collections.Concurrent;
using System.Diagnostics;

using SparkCL;
using OCLHelper;

using Matrices;
using Types;

namespace SlaeBuilder;
using static Shared;

class SymDiagSlaeBuilder : ISlaeBuilder
{
    SymDiagMatrix _matrix;
    Real[] _b = [];

    readonly RectMesh _mesh;
    public RectMesh Mesh { get => _mesh; }
    public GlobalMatrixImplType GlobalMatrixImpl { get; set; } = GlobalMatrixImplType.Host;

    readonly TaskFuncs _funcs;

    public SymDiagSlaeBuilder(RectMesh mesh, TaskFuncs funcs)
    {
        _mesh = mesh;
        _matrix = new SymDiagMatrix();
        _funcs = funcs;
    }

    public static ISlaeBuilder Construct(RectMesh mesh, TaskFuncs funcs)
        => new SymDiagSlaeBuilder(mesh, funcs);

    public (Matrix, Real[]) Build()
    {
        Trace.WriteLine($"SymDiag Builder: {GlobalMatrixImpl}");
        
        Trace.Indent();
        var sw = Stopwatch.StartNew();
        GlobalMatrixInit();
        Trace.WriteLine($"Init: {sw.ElapsedMilliseconds}ms");

        sw.Restart();
        GlobalMatrixBuild();
        Trace.WriteLine($"Build: {sw.ElapsedMilliseconds}ms");

        sw.Restart();
        BoundaryConditionsApply();
        Trace.WriteLine($"Conds: {sw.ElapsedMilliseconds}");
        Trace.Unindent();

        return (_matrix, _b);
    }

    void GlobalMatrixInit()
    {
        GlobalMatrixPortraitCompose();

        _matrix.Di = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _matrix.d0 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.d1 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.d2 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.d3 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _b = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
    }

    void BoundaryConditionType1Apply(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        void SharedBody(int m, Real b)
        {
            _b[m] = b;
            _matrix.Di[m] = 1;

            int t;

#if true
            /* Гауссово исключение столбца */
            t = m-3-_matrix.Gap;
            if (t >= 0) _b[t] -= b * _matrix.d3[t];
            t = m-2-_matrix.Gap;
            if (t >= 0) _b[t] -= b * _matrix.d2[t];
            t = m-1-_matrix.Gap;
            if (t >= 0) _b[t] -= b * _matrix.d1[t];
            t = m-1;
            if (t >= 0) _b[t] -= b * _matrix.d0[t];
            
            t = m+1;
            if (t < _matrix.Size) _b[t] -= b * _matrix.d0[m];
            t = m+1+_matrix.Gap;
            if (t < _matrix.Size) _b[t] -= b * _matrix.d1[m];
            t = m+2+_matrix.Gap;
            if (t < _matrix.Size) _b[t] -= b * _matrix.d2[m];
            t = m+3+_matrix.Gap;
            if (t < _matrix.Size) _b[t] -= b * _matrix.d3[m];
#endif

            /* Обнуление строки и столбца */
            _matrix.d3[m] = 0;
            _matrix.d2[m] = 0;
            _matrix.d1[m] = 0;
            _matrix.d0[m] = 0;

            t = m - 3 - _matrix.Gap;
            if (t >= 0) _matrix.d3[t] = 0;
            t = m - 2 - _matrix.Gap;
            if (t >= 0) _matrix.d2[t] = 0;
            t = m - 1 - _matrix.Gap;
            if (t >= 0) _matrix.d1[t] = 0;
            t = m - 1;
            if (t >= 0) _matrix.d0[t] = 0;
        }
        
        if (x1 == x2)
        {
            for (int yi = y1; yi <= y2; yi++)
            {
                var m = yi * _mesh.X.Length + x1;
                var b = _funcs.Ug(num, _mesh.X[x1], _mesh.Y[yi]);
                SharedBody(m, b);
            }
        }
        else if (y1 == y2)
        {
            for (int xi = x1; xi <= x2; xi++)
            {
                var m = y1 * _mesh.X.Length + xi;
                var b = _funcs.Ug(num, _mesh.X[xi], _mesh.Y[y1]);
                SharedBody(m, b);
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }

    void BoundaryConditionType2Apply(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        if (x1 == x2)
        {
            for (int yi = y1; yi < y2; yi++)
            {
                var h = _mesh.Y[yi + 1] - _mesh.Y[yi];

                Real k1 = _funcs.Theta(num, _mesh.X[x1], _mesh.Y[yi]); // aka theta1
                Real k2 = _funcs.Theta(num, _mesh.X[x2], _mesh.Y[yi + 1]);
                int n1 = yi * _mesh.X.Length + x1;
                int n2 = (yi + 1) * _mesh.X.Length + x2;

                _b[n1] += h * (2 * k1 + k2) / 6;
                _b[n2] += h * (k1 + 2 * k2) / 6;
            }
        }
        else if (y1 == y2)
        {
            for (int xi = x1; xi < x2; xi++)
            {
                var h = _mesh.X[xi + 1] - _mesh.X[xi];

                Real k1 = _funcs.Theta(num, _mesh.X[xi], _mesh.Y[y1]); // aka theta1
                Real k2 = _funcs.Theta(num, _mesh.X[xi + 1], _mesh.Y[y2]);
                int n1 = y1 * _mesh.X.Length + xi;
                int n2 = y2 * _mesh.X.Length + xi + 1;

                _b[n1] += h * (2 * k1 + k2) / 6;
                _b[n2] += h * (k1 + 2 * k2) / 6;
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }

    void BoundaryConditionType3Apply(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        var localB = new Real[2]; // 'hat B'
        var localA = new Real[2, 2]; // 'hat A'
        Real h;
        if (x1 == x2)
        {
            for (int yi = y1; yi < y2; yi++)
            {
                h = _mesh.Y[yi + 1] - _mesh.Y[yi];
                localA[0, 0] = localA[1, 1] = _funcs.Beta(num) * h / 3;
                localA[0, 1] = localA[1, 0] = _funcs.Beta(num) * h / 6;

                Real k1 = _funcs.uBeta(num, _mesh.X[x1], _mesh.Y[yi]);
                Real k2 = _funcs.uBeta(num, _mesh.X[x2], _mesh.Y[yi + 1]);
                localB[0] = h * _funcs.Beta(num) * (2 * k1 + k2) / 6;
                localB[1] = h * _funcs.Beta(num) * (k1 + 2 * k2) / 6;

                var m = new int[2];
                m[0] = yi * _mesh.X.Length + x1;
                m[1] = (yi + 1) * _mesh.X.Length + x2;

                _b[m[0]] += localB[0];
                _b[m[1]] += localB[1];

                _matrix.Di[m[0]] += localA[0, 0];
                _matrix.Di[m[1]] += localA[1, 1];

                _matrix.d2[m[0]] += localA[0, 1];
            }
        }
        else if (y1 == y2)
        {
            for (int xi = x1; xi < x2; xi++)
            {
                h = _mesh.X[xi + 1] - _mesh.X[xi];
                localA[0, 0] = localA[1, 1] = _funcs.Beta(num) * h / 3;
                localA[0, 1] = localA[1, 0] = _funcs.Beta(num) * h / 6;

                Real k1 = _funcs.uBeta(num, _mesh.X[xi], _mesh.Y[y1]);
                Real k2 = _funcs.uBeta(num, _mesh.X[xi + 1], _mesh.Y[y2]);
                localB[0] = h * _funcs.Beta(num) * (2 * k1 + k2) / 6;
                localB[1] = h * _funcs.Beta(num) * (k1 + 2 * k2) / 6;

                var m = new int[2];
                m[0] = y1 * _mesh.X.Length + xi;
                m[1] = y2 * _mesh.X.Length + xi + 1;

                _b[m[0]] += localB[0];
                _b[m[1]] += localB[1];

                _matrix.Di[m[0]] += localA[0, 0];
                _matrix.Di[m[1]] += localA[1, 1];

                _matrix.d0[m[0]] += localA[0, 1];
            }
        }
        else
        {
            throw new ArgumentException("Странное краевое условие");
        }
    }

    void BoundaryConditionsApply()
    {
        var bc_type1 = new List<BoundaryCondition>();

        foreach (var bc in _mesh.BoundaryConditions)
        {
            var type = bc.Type;

            switch (type)
            {
                case 1:
                    /* К.у. первого рода будут применены последними */
                    bc_type1.Add(bc);
                    break;
                case 2:
                    BoundaryConditionType2Apply(bc);
                    break;
                case 3:
                    BoundaryConditionType3Apply(bc);
                    break;

                default:
                    throw new Exception("Странный тип краевого условия");
            }
        }

        foreach (var b1 in bc_type1)
        {
            BoundaryConditionType1Apply(b1);
        }
    }

    void GlobalMatrixPortraitCompose()
    {
        // в случае диагональной матрицы это просто Gap
        _matrix.Gap = _mesh.X.Length - 2;
    }

    void GlobalMatrixBuild()
    {
        switch (GlobalMatrixImpl)
        {
            case GlobalMatrixImplType.OpenCL:
                GlobalMatrixBuildImplOcl();
                break;
            case GlobalMatrixImplType.Host:
                GlobalMatrixBuildImplHost();
                break;
            case GlobalMatrixImplType.HostParallel:
                GlobalMatrixBuildImplHostParallel();
                break;
            case GlobalMatrixImplType.HostV2:
                throw new NotImplementedException();
            case GlobalMatrixImplType.OpenCLV2:
                throw new NotImplementedException();
            default:
                throw new InvalidOperationException();
        }
    }

    void GlobalMatrixBuildImplHost()
    {
        // csharp не нравится stackalloc в циклах
        Span<Real> localB = stackalloc Real[4];

        for (int yi = 0; yi < _mesh.Y.Length - 1; yi++)
        {
            for (int xi = 0; xi < _mesh.X.Length - 1; xi++)
            {
                var subDom = _mesh.GetSubdomNumAtElCoords(xi, yi);

                if (!subDom.HasValue) continue;

                Real x0 = _mesh.X[xi];
                Real x1 = _mesh.X[xi + 1];
                Real y0 = _mesh.Y[yi];
                Real y1 = _mesh.Y[yi + 1];

                Real GetGammaAverage()
                {
                    Real res = _funcs.Gamma(subDom.Value, x0, y0)
                            + _funcs.Gamma(subDom.Value, x1, y0)
                            + _funcs.Gamma(subDom.Value, x0, y1)
                            + _funcs.Gamma(subDom.Value, x1, y1);

                    return res / 4;
                }

                Real GetLamdaAverage()
                {
                    Real res = _funcs.Lambda(subDom.Value, x0, y0)
                            + _funcs.Lambda(subDom.Value, x1, y0)
                            + _funcs.Lambda(subDom.Value, x0, y1)
                            + _funcs.Lambda(subDom.Value, x1, y1);

                    return res / 4;
                }

                Real hy = y1 - y0;
                Real hx = x1 - x0;
                // Заменить на интеграл от биквадратичного разложения
                Real l_avg = GetLamdaAverage();
                Real g_avg = GetGammaAverage();

                // мне понравилось это слово для обозначения левого нижнего узла элемента
                int anchor = yi * _mesh.X.Length + xi;

                Real Mat(int i, int j)
                {
                    return l_avg / 6 * (hy / hx * LocalG1[i, j] + hx / hy * LocalG2[i, j])
                         + g_avg / 36 * hx * hy * LocalM[i, j];
                }

                _matrix.Di[anchor] += Mat(0, 0);
                _matrix.d0[anchor] += Mat(0, 1);
                _matrix.d2[anchor] += Mat(0, 2);
                _matrix.d3[anchor] += Mat(0, 3);

                _matrix.Di[anchor + 1] += Mat(1, 1);
                _matrix.d1[anchor + 1] += Mat(1, 2);
                _matrix.d2[anchor + 1] += Mat(1, 3);

                var a2 = anchor + _mesh.X.Length;
                _matrix.Di[a2] += Mat(2, 2);
                _matrix.d0[a2] += Mat(2, 3);

                _matrix.Di[a2 + 1] += Mat(3, 3);

                /* правая часть */
                Real f1 = _funcs.F(subDom.Value, x0, y0);
                Real f2 = _funcs.F(subDom.Value, x1, y0);
                Real f3 = _funcs.F(subDom.Value, x0, y1);
                Real f4 = _funcs.F(subDom.Value, x1, y1);

                localB[0] = hx * hy / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
                localB[1] = hx * hy / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
                localB[2] = hx * hy / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
                localB[3] = hx * hy / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);

                _b[anchor] += localB[0];
                _b[anchor + 1] += localB[1];
                _b[a2] += localB[2];
                _b[a2 + 1] += localB[3];
            }
        }

        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < _matrix.Di.Length; i++)
        {
            if (_matrix.Di[i] == 0)
            {
                _matrix.Di[i] = 1;
            }
        }
    }

    void GlobalMatrixBuildImplHostParallel()
    {
        // csharp не нравится stackalloc в циклах
        var kernTime = Stopwatch.StartNew();
        var part_y = Partitioner.Create(0, _mesh.Y.Length - 1);

        Parallel.ForEach(part_y, (range, state) =>
        {
            Span<Real> localB = stackalloc Real[4];
            for (int yi = range.Item1; yi < range.Item2; yi++)
            {
                for (int xi = 0; xi < _mesh.X.Length - 1; xi++)
                {
                    var subDom = _mesh.GetSubdomNumAtElCoords(xi, yi);

                    if (!subDom.HasValue) continue;

                    Real x0 = _mesh.X[xi];
                    Real x1 = _mesh.X[xi + 1];
                    Real y0 = _mesh.Y[yi];
                    Real y1 = _mesh.Y[yi + 1];

                    Real GetGammaAverage()
                    {
                        Real res = _funcs.Gamma(subDom.Value, x0, y0)
                                + _funcs.Gamma(subDom.Value, x1, y0)
                                + _funcs.Gamma(subDom.Value, x0, y1)
                                + _funcs.Gamma(subDom.Value, x1, y1);

                        return res / 4;
                    }

                    Real GetLamdaAverage()
                    {
                        Real res = _funcs.Lambda(subDom.Value, x0, y0)
                                + _funcs.Lambda(subDom.Value, x1, y0)
                                + _funcs.Lambda(subDom.Value, x0, y1)
                                + _funcs.Lambda(subDom.Value, x1, y1);

                        return res / 4;
                    }

                    Real hy = y1 - y0;
                    Real hx = x1 - x0;
                    // Заменить на интеграл от биквадратичного разложения
                    Real l_avg = GetLamdaAverage();
                    Real g_avg = GetGammaAverage();

                    // мне понравилось это слово для обозначения левого нижнего узла элемента
                    int anchor = yi * _mesh.X.Length + xi;

                    Real Mat(int i, int j)
                        => l_avg / 6 * (hy / hx * LocalG1[i, j] + hx / hy * LocalG2[i, j])
                         + g_avg / 36 * hx * hy * LocalM[i, j];


                    Add(ref _matrix.Di[anchor], Mat(0, 0));
                    Add(ref _matrix.d0[anchor], Mat(0, 1));
                    Add(ref _matrix.d2[anchor], Mat(0, 2));
                    Add(ref _matrix.d3[anchor], Mat(0, 3));

                    Add(ref _matrix.Di[anchor + 1], Mat(1, 1));
                    Add(ref _matrix.d1[anchor + 1], Mat(1, 2));
                    Add(ref _matrix.d2[anchor + 1], Mat(1, 3));

                    var a2 = anchor + _mesh.X.Length;
                    Add(ref _matrix.Di[a2], Mat(2, 2));
                    Add(ref _matrix.d0[a2], Mat(2, 3));

                    Add(ref _matrix.Di[a2 + 1], Mat(3, 3));

                    /* правая часть */
                    Real f1 = _funcs.F(subDom.Value, x0, y0);
                    Real f2 = _funcs.F(subDom.Value, x1, y0);
                    Real f3 = _funcs.F(subDom.Value, x0, y1);
                    Real f4 = _funcs.F(subDom.Value, x1, y1);

                    localB[0] = hx * hy / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
                    localB[1] = hx * hy / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
                    localB[2] = hx * hy / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
                    localB[3] = hx * hy / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);

                    Add(ref _b[anchor], localB[0]);
                    Add(ref _b[anchor + 1], localB[1]);
                    Add(ref _b[a2], localB[2]);
                    Add(ref _b[a2 + 1], localB[3]);
                }
            }
        }
        );
        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < _matrix.Di.Length; i++)
        {
            if (_matrix.Di[i] == 0)
            {
                _matrix.Di[i] = 1;
            }
        }
    }

    static ComputeProgram? progCompose = null;
    void GlobalMatrixBuildImplOcl ()
    {
        Trace.Indent();
        var sw = Stopwatch.StartNew();
        
        if (progCompose == null)
        {
            progCompose = new ComputeProgram("SlaeBuilder/SymDiagCompose.cl");
        }
        var kernCompose = progCompose.GetKernel(
            "global_matrix_compose",
            globalWork: new(NDRange.PaddedTo(_mesh.X.Length, 4), NDRange.PaddedTo(_mesh.Y.Length, 4)),
            localWork: new(4, 4)
        );
        
        Trace.WriteLine($"Kernel prepare: {sw.ElapsedMilliseconds}");
        sw.Restart();
        
        var d3     = new ComputeBuffer<Real>(_matrix.d3, BufferFlags.OnDevice);
        var d2     = new ComputeBuffer<Real>(_matrix.d2, BufferFlags.OnDevice);
        var d1     = new ComputeBuffer<Real>(_matrix.d1, BufferFlags.OnDevice);
        var d0     = new ComputeBuffer<Real>(_matrix.d0, BufferFlags.OnDevice);
        var di     = new ComputeBuffer<Real>(_matrix.Di, BufferFlags.OnDevice);
        var b      = new ComputeBuffer<Real>(_b        , BufferFlags.OnDevice);
        var x_axis = new ComputeBuffer<Real>(_mesh.X   , BufferFlags.OnDevice);
        var y_axis = new ComputeBuffer<Real>(_mesh.Y   , BufferFlags.OnDevice);
        
        Trace.WriteLine($"Transfer Host->Device: {sw.ElapsedMilliseconds}");
        sw.Restart();
        
        kernCompose.SetArg(0, d3);
        kernCompose.SetArg(1, d2);
        kernCompose.SetArg(2, d1);
        kernCompose.SetArg(3, d0);
        kernCompose.SetArg(4, di);
        kernCompose.SetArg(5, b);
        kernCompose.SetArg(6, _matrix.Size);
        kernCompose.SetArg(7, _matrix.Gap);
        kernCompose.SetArg(8, x_axis);
        kernCompose.SetArg(9, x_axis.Length);
        kernCompose.SetArg(10, y_axis);
        kernCompose.SetArg(11, y_axis.Length);
        
        Trace.WriteLine($"Setargs: {sw.ElapsedMilliseconds}");
        sw.Restart();
        
        kernCompose.Execute();
        
        Trace.WriteLine($"Build time: {sw.ElapsedMilliseconds}ms");
        sw.Restart();
        
        d3.DeviceReadTo(_matrix.d3);
        d2.DeviceReadTo(_matrix.d2);
        d1.DeviceReadTo(_matrix.d1);
        d0.DeviceReadTo(_matrix.d0);
        di.DeviceReadTo(_matrix.Di);
        b.DeviceReadTo(_b);
        
        Trace.WriteLine($"Transfer Device->Host: {sw.ElapsedMilliseconds}");
        sw.Restart();
        
        // После сборки матрицы надо нулевые диагональные элементы заменить на 1
        for (int i = 0; i < _matrix.Size; i++)
        {
            if (_matrix.Di[i] == 0)
            {
                _matrix.Di[i] = 1;
            }
        }
        
        Trace.WriteLine($"0->1 on diag: {sw.ElapsedMilliseconds}ms");
        Trace.Unindent();
    }
}
