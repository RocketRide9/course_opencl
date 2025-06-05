using Real = double;

using System.Collections.Concurrent;
using System.Diagnostics;

namespace SlaeBuilder;
using static Shared;
using Matrices;

class DiagSlaeBuilder
{
    DiagMatrix _slae;
    Real[] _b = [];

    readonly RectMesh _mesh;
    public RectMesh Mesh { get => _mesh; }
    public GlobalMatrixImplType GlobalMatrixImpl { get; set; } = GlobalMatrixImplType.Host;

    readonly TaskFuncs _funcs;

    public DiagSlaeBuilder(RectMesh mesh, TaskFuncs funcs)
    {
        _mesh = mesh;
        _slae = new DiagMatrix();
        _funcs = funcs;
    }

    public (DiagMatrix, Real[]) Build()
    {
        GlobalMatrixInit();
        GlobalMatrixBuild();
        BoundaryConditionsApply();

        return (_slae, _b);
    }

    void GlobalMatrixInit()
    {
        GlobalMatrixPortraitCompose();

        _slae.Di = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _slae.Ld0 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _slae.Ld1 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _slae.Ld2 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _slae.Ld3 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _slae.Rd0 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _slae.Rd1 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _slae.Rd2 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _slae.Rd3 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

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

        if (x1 == x2)
        {
            for (int yi = y1; yi <= y2; yi++)
            {
                var m = yi * _mesh.X.Length + x1;
                _b[m] = _funcs.Ug(num, _mesh.X[x1], _mesh.Y[yi]);
                _slae.Di[m] = 1;

                /* Обнуление строки */
                _slae.Rd3[m] = 0;
                _slae.Rd2[m] = 0;
                _slae.Rd1[m] = 0;
                _slae.Rd0[m] = 0;

                if (m - 1 >= 0) _slae.Ld0[m - 1] = 0;
                if (m - 1 - _slae.Gap >= 0) _slae.Ld1[m - 1 - _slae.Gap] = 0;
                if (m - 2 - _slae.Gap >= 0) _slae.Ld2[m - 2 - _slae.Gap] = 0;
                if (m - 3 - _slae.Gap >= 0) _slae.Ld3[m - 3 - _slae.Gap] = 0;

            }
        }
        else if (y1 == y2)
        {
            for (int xi = x1; xi <= x2; xi++)
            {
                var m = y1 * _mesh.X.Length + xi;
                _b[m] = _funcs.Ug(num, _mesh.X[xi], _mesh.Y[y1]);
                _slae.Di[m] = 1;

                /* Обнуление строки */
                _slae.Rd3[m] = 0;
                _slae.Rd2[m] = 0;
                _slae.Rd1[m] = 0;
                _slae.Rd0[m] = 0;

                if (m - 1 >= 0) _slae.Ld0[m - 1] = 0;
                if (m - 1 - _slae.Gap >= 0) _slae.Ld1[m - 1 - _slae.Gap] = 0;
                if (m - 2 - _slae.Gap >= 0) _slae.Ld2[m - 2 - _slae.Gap] = 0;
                if (m - 3 - _slae.Gap >= 0) _slae.Ld3[m - 3 - _slae.Gap] = 0;
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

                _slae.Di[m[0]] += localA[0, 0];
                _slae.Di[m[1]] += localA[1, 1];

                _slae.Rd2[m[0]] += localA[0, 1];
                _slae.Ld2[m[1] - 2 - _slae.Gap] += localA[1, 0];
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

                _slae.Di[m[0]] += localA[0, 0];
                _slae.Di[m[1]] += localA[1, 1];

                _slae.Ld0[m[0]] += localA[0, 1];
                _slae.Ld0[m[1] - 1] += localA[1, 0];
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
        _slae.Gap = _mesh.X.Length - 2;
    }

#if false
    void GlobalMatrixBuildImplHostV2 ()
    {
        Real GetGammaAverage (int dom, Real x0, Real y0, Real x1, Real y1)
        {
            Real res = _funcs.Gamma(dom, x0, y0)
                       + _funcs.Gamma(dom, x1, y0)
                       + _funcs.Gamma(dom, x0, y1)
                       + _funcs.Gamma(dom, x1, y1);

            return res / 4;
        }

        Real GetLamdaAverage (int dom, Real x0, Real y0, Real x1, Real y1)
        {
            Real res = _funcs.Lambda(dom, x0, y0)
                       + _funcs.Lambda(dom, x1, y0)
                       + _funcs.Lambda(dom, x0, y1)
                       + _funcs.Lambda(dom, x1, y1);

            return res / 4;
        }
        
        // не функция из csharp потому что мне её ещё нужно самому реализовать в OpenCL
        static int QFind<T> (T[] @where, int start, int end, T what)
        where T: unmanaged, INumber<T>
        {
            int beg = start;
            while (beg < end)
            {
                int mid = (beg + end) / 2;
                if (what > where[mid])
                {
                    beg = mid + 1;
                }
                else
                {
                    end = mid;
                }
            }

            if (where[beg] != what)
            {
                throw new Exception("Quick search failed");
            }

            return beg;
        }
        
        static int LFind<T> (T[] where, T what, int start)
        where T: unmanaged, INumber<T>
        {
            while (@where[start] != what) start++;
            return start;
        }

        Span<Real> matc = stackalloc Real[8];
        for (int yi = 0; yi < _mesh.Y.Length; yi++)
        {
            for (int xi = 0; xi < _mesh.X.Length; xi++)
            {
                int targetNode = yi * _mesh.X.Length + xi;
                var dom1 = _mesh.GetSubdomNumAtElCoords(xi - 1, yi - 1);
                var dom2 = _mesh.GetSubdomNumAtElCoords(xi, yi - 1);
                var dom3 = _mesh.GetSubdomNumAtElCoords(xi - 1, yi);
                var dom4 = _mesh.GetSubdomNumAtElCoords(xi, yi);

                // номера текущего узла и соседних с ним сверху и снизу 
                // в порядке снизу вверх
                // могут выходить за пределы нумерации
                int r1 = targetNode;
                int r0 = r1 - _mesh.X.Length;
                int r2 = r1 + _mesh.X.Length;

                // номера этих узлов в строке матрицы, относящейся к текущему узлу
                int mr0 = -1;
                int mr1 = -1;
                int mr2 = -1;
                
                int beg = _slae.Ia[targetNode];
                int bound = _slae.Ia[targetNode + 1] - 1;

                // TODO: можно немного ускорить поиск за счёт перемещения левой границы поиска
                // но возможно быстрее будет сделать простой линейный поиск вместо "быстрого"
                // int curr = beg;
                if (r0 >= 0)
                {
                    mr0 = LFind(_slae.Ja, r0, beg) - beg;
                }
                // else оставить -1, так как узел вышел за пределы матрицы

                if (r1 % _mesh.X.Length == 0)
                {
                    mr1 = LFind(_slae.Ja, r1+1, beg) - 1 - beg;
                }
                else
                {
                    mr1 = LFind(_slae.Ja, r1-1, beg) - beg;
                }

                if (r2 < _slae.Di.Length)
                {
                    mr2 = LFind(_slae.Ja, r2, beg) - beg;
                }

                // Console.WriteLine("mr = {0}, {1}, {2}", mr0, mr1, mr2);
                // Console.WriteLine("xi, yi = {0}, {1}", xi, yi);
                // Console.WriteLine("dom = {0}, {1}, {2}, {3}", dom1, dom2, dom3, dom4);

                for (int i = 0; i < 8; i++)
                {
                    matc[i] = 0;
                }
                
                Real dic = 0;
                Real bc = 0;

                Real x1 = _mesh.X[xi];
                Real y1 = _mesh.Y[yi];

                if (dom1.HasValue)
                {
                    var x0 = _mesh.X[xi-1];
                    var y0 = _mesh.Y[yi-1];
                    var l_avg = GetLamdaAverage(dom1.Value, x0, y0, x1, y1);
                    var g_avg = GetGammaAverage(dom1.Value, x0, y0, x1, y1);
                    var hx0 = x1 - x0;
                    var hy0 = y1 - y0;

                    matc[mr0 - 1] += l_avg / 6 * (hy0 / hx0 * _localG1[3, 0] + hx0 / hy0 * _localG2[3, 0])
                        + g_avg / 36 * hx0 * hy0 * _localM[3, 0];
                    matc[mr0] += l_avg / 6 * (hy0 / hx0 * _localG1[3, 1] + hx0 / hy0 * _localG2[3, 1])
                        + g_avg / 36 * hx0 * hy0 * _localM[3, 1];
                    matc[mr1]     += l_avg / 6 * (hy0 / hx0 * _localG1[3, 2] + hx0 / hy0 * _localG2[3, 2])
                        + g_avg / 36 * hx0 * hy0 * _localM[3, 2];
                    dic += l_avg / 6 * (hy0 / hx0 * _localG1[3, 3] + hx0 / hy0 * _localG2[3, 3])
                        + g_avg / 36 * hx0 * hy0 * _localM[3, 3];

                    Real f1 = _funcs.F(dom1.Value, x0, y0);
                    Real f2 = _funcs.F(dom1.Value, x1, y0);
                    Real f3 = _funcs.F(dom1.Value, x0, y1);
                    Real f4 = _funcs.F(dom1.Value, x1, y1);

                    bc += hx0 * hy0 / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);
                }

                // continue;

                if (dom2.HasValue)
                {
                    var x2 = _mesh.X[xi+1];
                    var y0 = _mesh.Y[yi-1];
                    var l_avg = GetLamdaAverage(dom2.Value, x1, y0, x2, y1);
                    var g_avg = GetGammaAverage(dom2.Value, x1, y0, x2, y1);
                    var hx1 = x2 - x1;
                    var hy0 = y1 - y0;

                    matc[mr0] += l_avg / 6 * (hy0 / hx1 * _localG1[2, 0] + hx1 / hy0 * _localG2[2, 0])
                        + g_avg / 36 * hx1 * hy0 * _localM[2, 0];
                    matc[mr0 + 1] += l_avg / 6 * (hy0 / hx1 * _localG1[2, 1] + hx1 / hy0 * _localG2[2, 1])
                        + g_avg / 36 * hx1 * hy0 * _localM[2, 1];
                    dic += l_avg / 6 * (hy0 / hx1 * _localG1[2, 2] + hx1 / hy0 * _localG2[2, 2])
                        + g_avg / 36 * hx1 * hy0 * _localM[2, 2];
                    matc[mr1 + 1] += l_avg / 6 * (hy0 / hx1 * _localG1[2, 3] + hx1 / hy0 * _localG2[2, 3])
                        + g_avg / 36 * hx1 * hy0 * _localM[2, 3];

                    Real f1 = _funcs.F(dom2.Value, x1, y0);
                    Real f2 = _funcs.F(dom2.Value, x2, y0);
                    Real f3 = _funcs.F(dom2.Value, x1, y1);
                    Real f4 = _funcs.F(dom2.Value, x2, y1);

                    bc += hx1 * hy0 / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
                }

                if (dom3.HasValue)
                {
                    var x0 = _mesh.X[xi-1];
                    var y2 = _mesh.Y[yi+1];
                    var l_avg = GetLamdaAverage(dom3.Value, x0, y1, x1, y2);
                    var g_avg = GetGammaAverage(dom3.Value, x0, y1, x1, y2);
                    var hx0 = x1 - x0;
                    var hy1 = y2 - y1;

                    matc[mr1] += l_avg / 6 * (hy1 / hx0 * _localG1[1, 0] + hx0 / hy1 * _localG2[1, 0])
                        + g_avg / 36 * hx0 * hy1 * _localM[1, 0];
                    dic += l_avg / 6 * (hy1 / hx0 * _localG1[1, 1] + hx0 / hy1 * _localG2[1, 1])
                        + g_avg / 36 * hx0 * hy1 * _localM[1, 1];
                    matc[mr2 - 1] += l_avg / 6 * (hy1 / hx0 * _localG1[1, 2] + hx0 / hy1 * _localG2[1, 2])
                        + g_avg / 36 * hx0 * hy1 * _localM[1, 2];
                    matc[mr2] += l_avg / 6 * (hy1 / hx0 * _localG1[1, 3] + hx0 / hy1 * _localG2[1, 3])
                        + g_avg / 36 * hx0 * hy1 * _localM[1, 3];

                    Real f1 = _funcs.F(dom3.Value, x0, y1);
                    Real f2 = _funcs.F(dom3.Value, x1, y1);
                    Real f3 = _funcs.F(dom3.Value, x0, y2);
                    Real f4 = _funcs.F(dom3.Value, x1, y2);

                    bc += hx0 * hy1 / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
                }

                if (dom4.HasValue)
                {
                    var x2 = _mesh.X[xi+1];
                    var y2 = _mesh.Y[yi+1];
                    var l_avg = GetLamdaAverage(dom4.Value, x1, y1, x2, y2);
                    var g_avg = GetGammaAverage(dom4.Value, x1, y1, x2, y2);
                    var hx1 = x2 - x1;
                    var hy1 = y2 - y1;

                    dic += l_avg / 6 * (hy1 / hx1 * _localG1[0, 0] + hx1 / hy1 * _localG2[0, 0])
                        + g_avg / 36 * hx1 * hy1 * _localM[0, 0];
                    matc[mr1 + 1] += l_avg / 6 * (hy1 / hx1 * _localG1[0, 1] + hx1 / hy1 * _localG2[0, 1])
                        + g_avg / 36 * hx1 * hy1 * _localM[0, 1];
                    matc[mr2] += l_avg / 6 * (hy1 / hx1 * _localG1[0, 2] + hx1 / hy1 * _localG2[0, 2])
                        + g_avg / 36 * hx1 * hy1 * _localM[0, 2];
                    matc[mr2 + 1] += l_avg / 6 * (hy1 / hx1 * _localG1[0, 3] + hx1 / hy1 * _localG2[0, 3])
                        + g_avg / 36 * hx1 * hy1 * _localM[0, 3];

                    Real f1 = _funcs.F(dom4.Value, x1, y1);
                    Real f2 = _funcs.F(dom4.Value, x2, y1);
                    Real f3 = _funcs.F(dom4.Value, x1, y2);
                    Real f4 = _funcs.F(dom4.Value, x2, y2);

                    bc += hx1 * hy1 / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
                }

                _slae.Di[targetNode] = dic;
                B[targetNode]  = bc;

                // перемещение строки в матрицу
                for (int i = beg; i <= bound; i++)
                {
                    _slae.Mat[i] = matc[i-beg];
                }
            }
        }
    }
#endif

    void GlobalMatrixBuild()
    {
        switch (GlobalMatrixImpl)
        {
            case GlobalMatrixImplType.OpenCL:
                // GlobalMatrixBuildImplOcl();
                break;
            case GlobalMatrixImplType.Host:
                GlobalMatrixBuildImplHost();
                break;
            case GlobalMatrixImplType.HostParallel:
                GlobalMatrixBuildImplHostParallel();
                break;
            case GlobalMatrixImplType.HostV2:
                // GlobalMatrixBuildImplHostV2();
                break;
            case GlobalMatrixImplType.OpenCLV2:
                // GlobalMatrixBuildImplOclV2();
                break;
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

                _slae.Di[anchor] += Mat(0, 0);
                _slae.Ld0[anchor] += Mat(0, 1);
                _slae.Rd0[anchor] += Mat(0, 1);
                _slae.Ld2[anchor] += Mat(0, 2);
                _slae.Rd2[anchor] += Mat(0, 2);
                _slae.Ld3[anchor] += Mat(0, 3);
                _slae.Rd3[anchor] += Mat(0, 3);

                _slae.Di[anchor + 1] += Mat(1, 1);
                _slae.Ld1[anchor + 1] += Mat(1, 2);
                _slae.Rd1[anchor + 1] += Mat(1, 2);
                _slae.Ld2[anchor + 1] += Mat(1, 3);
                _slae.Rd2[anchor + 1] += Mat(1, 3);

                var a2 = anchor + _slae.Gap + 1;
                _slae.Di[a2] += Mat(2, 2);
                _slae.Ld0[a2] += Mat(2, 3);
                _slae.Rd0[a2] += Mat(2, 3);

                _slae.Di[a2 + 1] += Mat(3, 3);

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
                _b[anchor + _slae.Gap + 1] += localB[2];
                _b[anchor + _slae.Gap + 2] += localB[3];
            }
        }

        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < _slae.Di.Length; i++)
        {
            if (_slae.Di[i] == 0)
            {
                _slae.Di[i] = 1;
            }
        }
    }

    // https://stackoverflow.com/a/16893641
    public static double Add(ref double location1, double value)
    {
        double newCurrentValue = location1; // non-volatile read, so may be stale
        while (true)
        {
            double currentValue = newCurrentValue;
            double newValue = currentValue + value;
            newCurrentValue = Interlocked.CompareExchange(ref location1, newValue, currentValue);
            if (newCurrentValue.Equals(currentValue))
            {
                return newValue;
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


                    Add(ref _slae.Di[anchor], Mat(0, 0));
                    Add(ref _slae.Ld0[anchor], Mat(0, 1));
                    Add(ref _slae.Rd0[anchor], Mat(0, 1));
                    Add(ref _slae.Ld2[anchor], Mat(0, 2));
                    Add(ref _slae.Rd2[anchor], Mat(0, 2));
                    Add(ref _slae.Ld3[anchor], Mat(0, 3));
                    Add(ref _slae.Rd3[anchor], Mat(0, 3));

                    Add(ref _slae.Di[anchor + 1], Mat(1, 1));
                    Add(ref _slae.Ld1[anchor + 1], Mat(1, 2));
                    Add(ref _slae.Rd1[anchor + 1], Mat(1, 2));
                    Add(ref _slae.Ld2[anchor + 1], Mat(1, 3));
                    Add(ref _slae.Rd2[anchor + 1], Mat(1, 3));

                    var a2 = anchor + _slae.Gap + 1;
                    Add(ref _slae.Di[a2], Mat(2, 2));
                    Add(ref _slae.Ld0[a2], Mat(2, 3));
                    Add(ref _slae.Rd0[a2], Mat(2, 3));

                    Add(ref _slae.Di[a2 + 1], Mat(3, 3));

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

        Console.WriteLine($"Build time: {kernTime.ElapsedMilliseconds}ms");
        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < _slae.Di.Length; i++)
        {
            if (_slae.Di[i] == 0)
            {
                _slae.Di[i] = 1;
            }
        }
    }

#if false
    void GlobalMatrixBuildImplOcl ()
    {
        var prog = new Program("Kernels.cl");
        var kernCompose = prog.GetKernel(
            "global_matrix_compose",
            globalWork: new(PaddedTo(Mesh.X.Length, 4), PaddedTo(Mesh.Y.Length, 4)),
            localWork: new(4, 4)
        );
        var mat    = new ComputeBuffer<Real>(_slae.Mat, BufferFlags.OnDevice);
        var di     = new ComputeBuffer<Real>(_slae.Di , BufferFlags.OnDevice);
        var b      = new ComputeBuffer<Real>(B  , BufferFlags.OnDevice);
        var ia     = new ComputeBuffer<int> (_slae.Ia , BufferFlags.OnDevice);
        var ja     = new ComputeBuffer<int> (_slae.Ja , BufferFlags.OnDevice);
        var x_axis = new ComputeBuffer<Real>(Mesh.X   , BufferFlags.OnDevice);
        var y_axis = new ComputeBuffer<Real>(Mesh.Y   , BufferFlags.OnDevice);

        kernCompose.SetArg(0, mat);
        kernCompose.SetArg(1, di);
        kernCompose.SetArg(2, b);
        kernCompose.SetArg(3, ia);
        kernCompose.SetArg(4, ja);
        kernCompose.SetArg(5, di.Length);
        kernCompose.SetArg(6, x_axis);
        kernCompose.SetArg(7, x_axis.Length);
        kernCompose.SetArg(8, y_axis);
        kernCompose.SetArg(9, y_axis.Length);

        var kernTime = Stopwatch.StartNew();
        kernCompose.Execute();
        Console.WriteLine($"Build time: {kernTime.ElapsedMilliseconds}ms");

        mat.DeviceReadTo(_slae.Mat);
        di.DeviceReadTo(_slae.Di);
        b.DeviceReadTo(B);

        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < _slae.Di.Length; i++)
        {
            if (_slae.Di[i] == 0)
            {
                _slae.Di[i] = 1;
            }
        }
    }
#endif

#if false
    void GlobalMatrixBuildImplOclV2 ()
    {
        var prog = new Program("ComposeV2.cl");
        var kernCompose = prog.GetKernel(
            "global_matrix_compose_v2",
            globalWork: new(PaddedTo(Mesh.X.Length, 4), PaddedTo(Mesh.Y.Length, 4)),
            localWork: new(4, 4)
        );
        var mat    = new ComputeBuffer<Real>(_slae.Mat, BufferFlags.OnDevice);
        var di     = new ComputeBuffer<Real>(_slae.Di , BufferFlags.OnDevice);
        var b      = new ComputeBuffer<Real>(B  , BufferFlags.OnDevice);
        var ia     = new ComputeBuffer<int> (_slae.Ia , BufferFlags.OnDevice);
        var ja     = new ComputeBuffer<int> (_slae.Ja , BufferFlags.OnDevice);
        var x_axis = new ComputeBuffer<Real>(Mesh.X   , BufferFlags.OnDevice);
        var y_axis = new ComputeBuffer<Real>(Mesh.Y   , BufferFlags.OnDevice);

        kernCompose.SetArg(0, mat);
        kernCompose.SetArg(1, di);
        kernCompose.SetArg(2, b);
        kernCompose.SetArg(3, ia);
        kernCompose.SetArg(4, ja);
        kernCompose.SetArg(5, di.Length);
        kernCompose.SetArg(6, x_axis);
        kernCompose.SetArg(7, x_axis.Length);
        kernCompose.SetArg(8, y_axis);
        kernCompose.SetArg(9, y_axis.Length);

        var kernTime = Stopwatch.StartNew();
        kernCompose.Execute();
        Console.WriteLine($"Build time: {kernTime.ElapsedMilliseconds}ms");

        mat.DeviceReadTo(_slae.Mat);
        di.DeviceReadTo(_slae.Di);
        b.DeviceReadTo(B);

        // TODO: это легко сделать в OpenCL
        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < _slae.Di.Length; i++)
        {
            if (_slae.Di[i] == 0)
            {
                _slae.Di[i] = 1;
            }
        }
    }
#endif
}
