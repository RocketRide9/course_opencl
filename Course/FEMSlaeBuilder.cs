using System.Collections.Concurrent;
using System.Numerics;
using SparkCL;
using Real = double;

public enum GlobalMatrixImplType
{
    OpenCL,
    OpenCLV2,
    Host,
    HostParallel,
    HostV2
}

class FEMSlaeBuilder
{
    Slae2 _slae;
    readonly RectMesh _mesh;
    public RectMesh Mesh { get => _mesh; }
    public GlobalMatrixImplType GlobalMatrixImpl { get; set; } = GlobalMatrixImplType.Host;

    readonly Real[,] _localG1 = {
        { 2, -2,  1, -1},
        {-2,  2, -1,  1},
        { 1, -1,  2, -2},
        {-1,  1, -2,  2},
    };
    readonly Real[,] _localG2 = {
        { 2,  1, -2, -1},
        { 1,  2, -1, -2},
        {-2, -1,  2,  1},
        {-1, -2,  1,  2},
    };
    readonly Real[,] _localM = {
        {4, 2, 2, 1},
        {2, 4, 1, 2},
        {2, 1, 4, 2},
        {1, 2, 2, 4},
    };

    readonly TaskFuncs _funcs;

    public FEMSlaeBuilder(RectMesh mesh, TaskFuncs funcs)
    {
        _mesh = mesh;
        _slae = new Slae2();
        _funcs = funcs;
    }

    public Slae2 Build()
    {
        GlobalMatrixInit();
        GlobalMatrixBuild();
        BoundaryConditionsApply();

        return _slae;
    }

    void GlobalMatrixInit()
    {
        GlobalMatrixPortraitCompose();

        _slae.Mat = Enumerable.Repeat((Real)0, _slae.Ja.Length).ToArray();
        _slae.Di  = Enumerable.Repeat((Real)0, _slae.Ia.Length - 1).ToArray();
        _slae.B   = Enumerable.Repeat((Real)0, _slae.Ia.Length - 1).ToArray();
    }

    void BoundaryConditionType1Apply(BoundaryCondition bc)
    {
        /* учитывание разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        Span<Real> localB = stackalloc Real[2]; // 'hat B'

        int a1 = x1;
        int a2 = x2;
        int b1 = y1;
        int b2 = y2;

        ref int e1 = ref a1;
        ref int e2 = ref a2;
        int upperBound;

        if (y1 == y2)
        {
            e1 = ref a1;
            e2 = ref a2;
            upperBound = x2;
        } else if (x1 == x2) {
            e1 = ref b1;
            e2 = ref b2;
            upperBound = y2;
        } else {
            throw new ArgumentException("Странное краевое условие");
        }

        for (e2 = e1 + 1; e2 <= upperBound; e2++)
        {
            localB[0] = _funcs.Ug(num, _mesh.X[a1], _mesh.Y[b1]);
            localB[1] = _funcs.Ug(num, _mesh.X[a2], _mesh.Y[b2]);

            // номера узлов, через которые проъодит первое краевое условие
            Span<int> m = [
                b1 * _mesh.X.Length + a1,
                b2 * _mesh.X.Length + a2
            ];
            _slae.B[m[0]] = localB[0];
            _slae.B[m[1]] = localB[1];

            _slae.Di[m[0]] = 1;
            _slae.Di[m[1]] = 1;

            /* Обнуление строки */
            for (int idx = 0; idx < 2; idx++)
            {
                int ig0 = _slae.Ia[m[idx]];
                int ig1 = _slae.Ia[m[idx]+1];
                for (int i = ig0; i < ig1; i++)
                {
                    _slae.Mat[i] = 0;
                }
            }
            e1 = e2;
        }
    }

    void BoundaryConditionType2Apply(BoundaryCondition bc)
    {
        /* учитывание разбиения сетки */
        int x1 = _mesh.XAfterGridInit(bc.X1);
        int x2 = _mesh.XAfterGridInit(bc.X2);
        int y1 = _mesh.YAfterGridInit(bc.Y1);
        int y2 = _mesh.YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        Span<Real> localB = stackalloc Real[2];

        int a1 = x1;
        int a2 = x2;
        int b1 = y1;
        int b2 = y2;

        ref int e1 = ref a1;
        ref int e2 = ref a2;
        int upperBound;

        if (y1 == y2)
        {
            e1 = ref a1;
            e2 = ref a2;
            upperBound = x2;
        } else if (x1 == x2) {
            e1 = ref b1;
            e2 = ref b2;
            upperBound = y2;
        } else {
            throw new ArgumentException("Странное краевое условие");
        }

        for (e2 = e1 + 1; e2 <= upperBound; e2++)
        {
            /* Формула опирается на предположение что одна из разностей
            равна нулю */
            Real h = (_mesh.X[a2] - _mesh.X[a1]) + (_mesh.Y[b2] - _mesh.Y[b1]);

            Real k1 = _funcs.Theta(num, _mesh.X[a1], _mesh.Y[b1]); // aka theta1
            Real k2 = _funcs.Theta(num, _mesh.X[a2], _mesh.Y[b2]);
            localB[0] = h * (2 * k1 + k2) / 6;
            localB[1] = h * (k1 + 2 * k2) / 6;

            int node1_num = b1 * _mesh.X.Length + a1;
            int node2_num = b2 * _mesh.X.Length + a2;
            _slae.B[node1_num] += localB[0];
            _slae.B[node2_num] += localB[1];

            e1 = e2;
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

        int a1 = x1;
        int a2 = x2;
        int b1 = y1;
        int b2 = y2;

        ref int e1 = ref a1;
        ref int e2 = ref a2;
        int upperBound;

        if (y1 == y2)
        {
            e1 = ref a1;
            e2 = ref a2;
            upperBound = x2;
        } else if (x1 == x2) {
            e1 = ref b1;
            e2 = ref b2;
            upperBound = y2;
        } else {
            throw new ArgumentException("Странное краевое условие");
        }

        for (e2 = e1 + 1; e2 <= upperBound; e2++)
        {
            Real h = _mesh.X[a2] - _mesh.X[a1] + _mesh.Y[b2] - _mesh.Y[b1];

            localA[0, 0] = localA[1, 1] = _funcs.Beta(num) * h / 3;
            localA[0, 1] = localA[1, 0] = _funcs.Beta(num) * h / 6;

            Real k1 = _funcs.uBeta(num, _mesh.X[a1], _mesh.Y[b1]);
            Real k2 = _funcs.uBeta(num, _mesh.X[a2], _mesh.Y[b2]);
            localB[0] = h * _funcs.Beta(num) * (2  * k1 + k2) / 6;
            localB[1] = h * _funcs.Beta(num) * (k1 + 2  * k2) / 6;

            var m = new int[2];
            m[0] = b1 * _mesh.X.Length + a1;
            m[1] = b2 * _mesh.X.Length + a2;

            _slae.B[m[0]] += localB[0];
            _slae.B[m[1]] += localB[1];

            _slae.Di[m[0]] += localA[0, 0];
            _slae.Di[m[1]] += localA[1, 1];

            for (int i = 0; i < 2; i++)
            {
                int beg = _slae.Ia[m[i]];
                for (int j = 0; j < 2; j++)
                {
                    // TODO: возможно есть лучше способ пропускать диагональные
                    // элементы
                    if (i == j)
                    {
                        continue;
                    }
                    int end = _slae.Ia[m[i] + 1] - 1;
                    while (beg < end)
                    {
                        int mid = (beg + end) / 2;
                        if (m[j] > _slae.Ja[mid])
                        {
                            beg = mid + 1;
                        }
                        else
                        {
                            end = mid;
                        }
                    }

                    if (_slae.Ja[beg] != m[j])
                    {
                        throw new Exception("Quick search failed");
                    }

                    _slae.Mat[beg] += localA[i, j];
                    beg++;
                }
            }
            e1 = e2;
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
        int numberOfUnknowns(int i) => 4;
        int idxOfUnknown (int i, int j)
        {
            int x0 = i % (_mesh.X.Length - 1);
            int y0 = i / (_mesh.X.Length - 1);

            if (j == 0 || j == 1)
            {
                return y0 * _mesh.X.Length + x0 + j;
            } else if (j == 2 || j == 3)
            {
                return (y0 + 1) * _mesh.X.Length + x0 + (j - 2);
            } else {
                throw new ArgumentException("Странная координата конечного элемента");
            }
        }


        HashSet<int>[] list = new HashSet<int>[Mesh.nodesCount];
        for (int i = 0; i < list.Length; i++)
        {
            list[i] = [];
        }

        /* цикл по всем конечным элементам */
        for (int ielem = 0; ielem < Mesh.feCount; ielem++)
        {
            /* цикл по всем узлам данного к.э. */
            for (int idx0 = 0; idx0 < numberOfUnknowns(ielem); idx0++)
            {
                /* цикл по узлам, соседним с idx0 */
                for (int idx1 = 0; idx1 < numberOfUnknowns(ielem); idx1++)
                {
                    if (idx0 == idx1) continue;
                    /* нахождение глобальных номеров локальных узлов */
                    int k1 = idxOfUnknown(ielem, idx0);
                    int k2 = idxOfUnknown(ielem, idx1);
                    /* */
                    /* заносим в list[k2] номера всех узлов, "соседних" с ним
                        по сути здесь k2 - номер строки, а k1 - номер ненулевого
                        элемента в глобальной матрице*/
                    list[k2].Add(k1);
                }
            }
        }

        _slae.Ia = new int[list.Length + 1];
        _slae.Ia[0] = 0;
        /* формирование массивов ig jg по списку list */
        for (int i = 1; i < _slae.Ia.Length; i++)
        {
            _slae.Ia[i] = _slae.Ia[i-1] + list[i-1].Count;
        }
        _slae.Ja = new int[_slae.Ia[list.Length]];
        for (var i = 0; i < list.Length; i++)
        {
            var row = list[i].Order().ToArray();
            for (int j = _slae.Ia[i]; j < _slae.Ia[i+1]; j++)
            {
                _slae.Ja[j] = row[j-_slae.Ia[i]];
            }
        }
    }

    void GlobalMatrixBuildImplHostV2 ()
    {
        Real GetGammaAverage (int dom, int x0, int y0)
        {
            Real res = _funcs.Gamma(dom, _mesh.X[x0], _mesh.Y[y0])
                       + _funcs.Gamma(dom, _mesh.X[x0 + 1], _mesh.Y[y0])
                       + _funcs.Gamma(dom, _mesh.X[x0], _mesh.Y[y0 + 1])
                       + _funcs.Gamma(dom, _mesh.X[x0 + 1], _mesh.Y[y0 + 1]);

            return res / 4;
        }

        Real GetLamdaAverage (int dom, int x0, int y0)
        {
            Real res = _funcs.Lambda(dom, _mesh.X[x0], _mesh.Y[y0])
                       + _funcs.Lambda(dom, _mesh.X[x0 + 1], _mesh.Y[y0])
                       + _funcs.Lambda(dom, _mesh.X[x0], _mesh.Y[y0 + 1])
                       + _funcs.Lambda(dom, _mesh.X[x0 + 1], _mesh.Y[y0 + 1]);

            return res / 4;
        }
        
        // не функция из csharp потому что мне её ещё нужно самому реализовать в OpenCL
        static int QFind<T> (T[] where, int start, int end, T what)
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
                var r = new int[3];
                r[1] = targetNode;
                r[0] = r[1] - _mesh.X.Length;
                r[2] = r[1] + _mesh.X.Length;

                // номера этих узлов в строке матрицы, относящейся к текущему узлу
                int[] mr = [-1, -1, -1];
                int beg = _slae.Ia[targetNode];
                int bound = _slae.Ia[targetNode + 1] - 1;
                if (r[0] >= 0)
                {
                    mr[0] = QFind(_slae.Ja, beg, bound, r[0]);
                    beg = mr[0];
                }
                // else оставить -1, так как узел вышел за пределы матрицы

                if (r[1] % _mesh.X.Length == 0)
                {
                    mr[1] = QFind(_slae.Ja, beg, bound, r[1]+1) - 1;
                }
                else
                {
                    mr[1] = QFind(_slae.Ja, beg, bound, r[1]-1);
                }

                if (r[2] < _slae.Di.Length)
                {
                    mr[2] = QFind(_slae.Ja, beg, bound, r[2]);
                }

                // Console.WriteLine("mr = {0}, {1}, {2}", mr[0], mr[1], mr[2]);
                // Console.WriteLine("xi, yi = {0}, {1}", xi, yi);
                // Console.WriteLine("dom = {0}, {1}, {2}, {3}", dom1, dom2, dom3, dom4);

                if (dom1.HasValue)
                {
                    var l_avg = GetLamdaAverage(dom1.Value, xi - 1, yi - 1);
                    var g_avg = GetGammaAverage(dom1.Value, xi - 1, yi - 1);
                    var hx0 = _mesh.X[xi] - _mesh.X[xi - 1];
                    var hy0 = _mesh.Y[yi] - _mesh.Y[yi - 1];

                    _slae.Mat[mr[0] - 1] += l_avg / 6 * (hy0 / hx0 * _localG1[3, 0] + hx0 / hy0 * _localG2[3, 0])
                        + g_avg / 36 * hx0 * hy0 * _localM[3, 0];
                    _slae.Mat[mr[0]]     += l_avg / 6 * (hy0 / hx0 * _localG1[3, 1] + hx0 / hy0 * _localG2[3, 1])
                        + g_avg / 36 * hx0 * hy0 * _localM[3, 1];
                    _slae.Mat[mr[1]]     += l_avg / 6 * (hy0 / hx0 * _localG1[3, 2] + hx0 / hy0 * _localG2[3, 2])
                        + g_avg / 36 * hx0 * hy0 * _localM[3, 2];
                    _slae.Di[targetNode] += l_avg / 6 * (hy0 / hx0 * _localG1[3, 3] + hx0 / hy0 * _localG2[3, 3])
                        + g_avg / 36 * hx0 * hy0 * _localM[3, 3];

                    var x0 = _mesh.X[xi-1];
                    var x1 = _mesh.X[xi];
                    var y0 = _mesh.Y[yi-1];
                    var y1 = _mesh.Y[yi];
                    Real f1 = _funcs.F(dom1.Value, x0, y0);
                    Real f2 = _funcs.F(dom1.Value, x1, y0);
                    Real f3 = _funcs.F(dom1.Value, x0, y1);
                    Real f4 = _funcs.F(dom1.Value, x1, y1);

                    _slae.B[targetNode] += hx0 * hy0 / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);
                }

                // continue;

                if (dom2.HasValue)
                {
                    var l_avg = GetLamdaAverage(dom2.Value, xi, yi - 1);
                    var g_avg = GetGammaAverage(dom2.Value, xi, yi - 1);
                    var hx1 = _mesh.X[xi + 1] - _mesh.X[xi];
                    var hy0 = _mesh.Y[yi] - _mesh.Y[yi - 1];

                    _slae.Mat[mr[0]] += l_avg / 6 * (hy0 / hx1 * _localG1[2, 0] + hx1 / hy0 * _localG2[2, 0])
                        + g_avg / 36 * hx1 * hy0 * _localM[2, 0];
                    _slae.Mat[mr[0] + 1] += l_avg / 6 * (hy0 / hx1 * _localG1[2, 1] + hx1 / hy0 * _localG2[2, 1])
                        + g_avg / 36 * hx1 * hy0 * _localM[2, 1];
                    _slae.Di[targetNode] += l_avg / 6 * (hy0 / hx1 * _localG1[2, 2] + hx1 / hy0 * _localG2[2, 2])
                        + g_avg / 36 * hx1 * hy0 * _localM[2, 2];
                    _slae.Mat[mr[1] + 1] += l_avg / 6 * (hy0 / hx1 * _localG1[2, 3] + hx1 / hy0 * _localG2[2, 3])
                        + g_avg / 36 * hx1 * hy0 * _localM[2, 3];

                    var x0 = _mesh.X[xi];
                    var x1 = _mesh.X[xi+1];
                    var y0 = _mesh.Y[yi-1];
                    var y1 = _mesh.Y[yi];
                    Real f1 = _funcs.F(dom2.Value, x0, y0);
                    Real f2 = _funcs.F(dom2.Value, x1, y0);
                    Real f3 = _funcs.F(dom2.Value, x0, y1);
                    Real f4 = _funcs.F(dom2.Value, x1, y1);

                    _slae.B[targetNode] += hx1 * hy0 / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
                }

                if (dom3.HasValue)
                {
                    var l_avg = GetLamdaAverage(dom3.Value, xi - 1, yi);
                    var g_avg = GetGammaAverage(dom3.Value, xi - 1, yi);
                    var hx0 = _mesh.X[xi] - _mesh.X[xi - 1];
                    var hy1 = _mesh.Y[yi + 1] - _mesh.Y[yi];

                    _slae.Mat[mr[1]] += l_avg / 6 * (hy1 / hx0 * _localG1[1, 0] + hx0 / hy1 * _localG2[1, 0])
                        + g_avg / 36 * hx0 * hy1 * _localM[1, 0];
                    _slae.Di[targetNode] += l_avg / 6 * (hy1 / hx0 * _localG1[1, 1] + hx0 / hy1 * _localG2[1, 1])
                        + g_avg / 36 * hx0 * hy1 * _localM[1, 1];
                    _slae.Mat[mr[2] - 1] += l_avg / 6 * (hy1 / hx0 * _localG1[1, 2] + hx0 / hy1 * _localG2[1, 2])
                        + g_avg / 36 * hx0 * hy1 * _localM[1, 2];
                    _slae.Mat[mr[2]] += l_avg / 6 * (hy1 / hx0 * _localG1[1, 3] + hx0 / hy1 * _localG2[1, 3])
                        + g_avg / 36 * hx0 * hy1 * _localM[1, 3];

                    var x0 = _mesh.X[xi-1];
                    var x1 = _mesh.X[xi];
                    var y0 = _mesh.Y[yi];
                    var y1 = _mesh.Y[yi+1];
                    Real f1 = _funcs.F(dom3.Value, x0, y0);
                    Real f2 = _funcs.F(dom3.Value, x1, y0);
                    Real f3 = _funcs.F(dom3.Value, x0, y1);
                    Real f4 = _funcs.F(dom3.Value, x1, y1);

                    _slae.B[targetNode] += hx0 * hy1 / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
                }

                if (dom4.HasValue)
                {
                    var l_avg = GetLamdaAverage(dom4.Value, xi, yi);
                    var g_avg = GetGammaAverage(dom4.Value, xi, yi);
                    var hy1 = _mesh.Y[yi + 1] - _mesh.Y[yi];
                    var hx1 = _mesh.X[xi + 1] - _mesh.X[xi];

                    _slae.Di[targetNode] += l_avg / 6 * (hy1 / hx1 * _localG1[0, 0] + hx1 / hy1 * _localG2[0, 0])
                        + g_avg / 36 * hx1 * hy1 * _localM[0, 0];
                    _slae.Mat[mr[1] + 1] += l_avg / 6 * (hy1 / hx1 * _localG1[0, 1] + hx1 / hy1 * _localG2[0, 1])
                        + g_avg / 36 * hx1 * hy1 * _localM[0, 1];
                    _slae.Mat[mr[2]] += l_avg / 6 * (hy1 / hx1 * _localG1[0, 2] + hx1 / hy1 * _localG2[0, 2])
                        + g_avg / 36 * hx1 * hy1 * _localM[0, 2];
                    _slae.Mat[mr[2] + 1] += l_avg / 6 * (hy1 / hx1 * _localG1[0, 3] + hx1 / hy1 * _localG2[0, 3])
                        + g_avg / 36 * hx1 * hy1 * _localM[0, 3];

                    var x0 = _mesh.X[xi];
                    var x1 = _mesh.X[xi+1];
                    var y0 = _mesh.Y[yi];
                    var y1 = _mesh.Y[yi+1];
                    Real f1 = _funcs.F(dom4.Value, x0, y0);
                    Real f2 = _funcs.F(dom4.Value, x1, y0);
                    Real f3 = _funcs.F(dom4.Value, x0, y1);
                    Real f4 = _funcs.F(dom4.Value, x1, y1);

                    _slae.B[targetNode] += hx1 * hy1 / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
                }
            }
        }
    }

    void GlobalMatrixBuild()
    {
        switch (GlobalMatrixImpl) {
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
                GlobalMatrixBuildImplHostV2();
                break;
            case GlobalMatrixImplType.OpenCLV2:
                GlobalMatrixBuildImplOclV2();
                break;
            default:
                throw new InvalidOperationException();
        }
    }

    void GlobalMatrixBuildImplHost()
    {
        // csharp не нравится stackalloc в циклах
        Span<Real> localB = stackalloc Real[4];
        Span<int> m = stackalloc int[4];

        for (int yi = 0; yi < _mesh.Y.Length - 1; yi++)
        {
            for (int xi = 0; xi < _mesh.X.Length - 1; xi++)
            {
                var subDom = _mesh.GetSubdomNumAtElCoords(xi, yi);

                m[0] = yi * _mesh.X.Length + xi;
                m[1] = m[0] + 1;
                m[2] = (yi + 1) * _mesh.X.Length + xi;
                m[3] = m[2] + 1;

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

                Real f1 = _funcs.F(subDom.Value, x0, y0);
                Real f2 = _funcs.F(subDom.Value, x1, y0);
                Real f3 = _funcs.F(subDom.Value, x0, y1);
                Real f4 = _funcs.F(subDom.Value, x1, y1);

                localB[0] = hx * hy / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
                localB[1] = hx * hy / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
                localB[2] = hx * hy / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
                localB[3] = hx * hy / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);

                /* нахождение в ja индексов элементов в al/au, куда
                    нужно добавить элементы локальных матриц */
                for (int i = 0; i < 4; i++)
                {
                    var v2 = l_avg / 6 * (hy / hx * _localG1[i, i] + hx / hy * _localG2[i, i])
                        + g_avg / 36 * hx * hy * _localM[i, i];
                    _slae.Di[m[i]] += v2;

                    int beg = _slae.Ia[m[i]];
                    for (int j = 0; j < 4; j++)
                    {
                        // TODO: пропуск
                        if (i == j)
                        {
                            continue;
                        }
                        int end = _slae.Ia[m[i] + 1] - 1;
                        while (beg < end)
                        {
                            int mid = (beg + end) / 2;
                            if (m[j] > _slae.Ja[mid])
                            {
                                beg = mid + 1;
                            }
                            else
                            {
                                end = mid;
                            }
                        }

                        if (_slae.Ja[beg] != m[j])
                        {
                            throw new Exception("Quick search failed");
                        }

                        v2 = l_avg / 6 * (hy / hx * _localG1[i, j] + hx / hy * _localG2[i, j])
                            + g_avg / 36 * hx * hy * _localM[i, j];
                        _slae.Mat[beg] += v2;
                        beg++;
                    }
                }

                /* добавление локальной правой части в слау */
                for (int i = 0; i < 4; i++)
                {
                    _slae.B[m[i]] += localB[i];
                }
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

        var part_y = Partitioner.Create(0, _mesh.Y.Length - 1);

        Parallel.ForEach(part_y, (range, state) =>
        {
            Span<Real> localB = stackalloc Real[4];
            Span<int> m = stackalloc int[4];
            for (int yi = range.Item1; yi < range.Item2; yi++)
            {
                for (int xi = 0; xi < _mesh.X.Length - 1; xi++)
                {
                    var subDom = _mesh.GetSubdomNumAtElCoords(xi, yi);

                    m[0] = yi * _mesh.X.Length + xi;
                    m[1] = m[0] + 1;
                    m[2] = (yi + 1) * _mesh.X.Length + xi;
                    m[3] = m[2] + 1;

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

                    Real f1 = _funcs.F(subDom.Value, x0, y0);
                    Real f2 = _funcs.F(subDom.Value, x1, y0);
                    Real f3 = _funcs.F(subDom.Value, x0, y1);
                    Real f4 = _funcs.F(subDom.Value, x1, y1);

                    localB[0] = hx * hy / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
                    localB[1] = hx * hy / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
                    localB[2] = hx * hy / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
                    localB[3] = hx * hy / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);

                    /* нахождение в ja индексов элементов в al/au, куда
                        нужно добавить элементы локальных матриц */
                    for (int i = 0; i < 4; i++)
                    {
                        var v2 = l_avg / 6 * (hy / hx * _localG1[i, i] + hx / hy * _localG2[i, i])
                            + g_avg / 36 * hx * hy * _localM[i, i];
                        Add(ref this._slae.Di[m[i]], v2);
                        // Slae.Di[m[i]] += v2;

                        int beg = this._slae.Ia[m[i]];
                        for (int j = 0; j < 4; j++)
                        {
                            // TODO: пропуск
                            if (i == j)
                            {
                                continue;
                            }
                            int end = this._slae.Ia[m[i] + 1] - 1;
                            while (beg < end)
                            {
                                int mid = (beg + end) / 2;
                                if (m[j] > this._slae.Ja[mid])
                                {
                                    beg = mid + 1;
                                }
                                else
                                {
                                    end = mid;
                                }
                            }

                            if (this._slae.Ja[beg] != m[j])
                            {
                                throw new Exception("Quick search failed");
                            }

                            v2 = l_avg / 6 * (hy / hx * _localG1[i, j] + hx / hy * _localG2[i, j])
                                + g_avg / 36 * hx * hy * _localM[i, j];
                            Add(ref this._slae.Mat[beg], v2);
                            // Slae.Mat[beg] += v2;
                            beg++;
                        }
                    }

                    /* добавление локальной правой части в слау */
                    for (int i = 0; i < 4; i++)
                    {
                        Add(ref this._slae.B[m[i]], localB[i]);
                        // Slae.B[m[i]] += localB[i];
                    }
                }
            }
        }
        );

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

    static nuint PaddedTo(int initial, int multiplier)
    {
        if (initial % multiplier == 0)
        {
            return (nuint) initial;
        } else {
            return (nuint) ( (initial / multiplier + 1 ) * multiplier );
        }
    }

    void GlobalMatrixBuildImplOcl ()
    {
        var prog = new Program("Kernels.clcpp");
        var kernCompose = prog.GetKernel(
            "global_matrix_compose",
            globalWork: new(PaddedTo(Mesh.X.Length, 4), PaddedTo(Mesh.Y.Length, 4)),
            localWork: new(4, 4)
        );
        var mat    = new ComputeBuffer<Real>(_slae.Mat, BufferFlags.OnDevice);
        var di     = new ComputeBuffer<Real>(_slae.Di , BufferFlags.OnDevice);
        var b      = new ComputeBuffer<Real>(_slae.B  , BufferFlags.OnDevice);
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

        kernCompose.Execute();

        mat.DeviceReadTo(_slae.Mat);
        di.DeviceReadTo(_slae.Di);
        b.DeviceReadTo(_slae.B);

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

    void GlobalMatrixBuildImplOclV2 ()
    {
        var prog = new Program("Kernels.clcpp");
        var kernCompose = prog.GetKernel(
            "global_matrix_compose_v2",
            globalWork: new(PaddedTo(Mesh.X.Length, 4), PaddedTo(Mesh.Y.Length, 4)),
            localWork: new(4, 4)
        );
        var mat    = new ComputeBuffer<Real>(_slae.Mat, BufferFlags.OnDevice);
        var di     = new ComputeBuffer<Real>(_slae.Di , BufferFlags.OnDevice);
        var b      = new ComputeBuffer<Real>(_slae.B  , BufferFlags.OnDevice);
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

        kernCompose.Execute();

        mat.DeviceReadTo(_slae.Mat);
        di.DeviceReadTo(_slae.Di);
        b.DeviceReadTo(_slae.B);

        // TODO: это легко сделать с OpenCL
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
}
