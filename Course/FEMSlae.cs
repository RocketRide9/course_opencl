using Real = float;

class FEMSlae
{
    Slae2 _slae;
    public Slae2 Slae { get => _slae; }
    RectMesh _mesh;
    public RectMesh Mesh { get => _mesh; }

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
    
    // Перестроить СЛАУ с разбитой сеткой
    void MeshRefine(RefineParams refineParams)
    {
        _mesh.Refine(refineParams);
        GlobalMatrixInit();
        GlobalMatrixBuild();
        BoundaryConditionsApply();
    }
    
    public FEMSlae(RectMesh mesh, TaskFuncs funcs, RefineParams? refineParams)
    {
        _mesh = mesh;
        _slae = new Slae2();
        _funcs = funcs;

        if (refineParams != null)
        {
            _mesh.Refine(refineParams.Value);
        }
        
        GlobalMatrixInit();
        GlobalMatrixBuild();
        BoundaryConditionsApply();
    }

    public void MeshDouble()
    {
        Mesh.RefineDiv2();

        GlobalMatrixInit();
        GlobalMatrixBuild();
        BoundaryConditionsApply();
    }

    void GlobalMatrixInit()
    {
        GlobalMatrixPortraitCompose();  

        _slae.Mat = Enumerable.Repeat((Real)0, Slae.Ja.Length).ToArray();
        _slae.Di = Enumerable.Repeat((Real)0, Slae.Ia.Length - 1).ToArray();
        _slae.B =  Enumerable.Repeat((Real)0, Slae.Ia.Length - 1).ToArray();
    }
    
    /* Перевод координаты x до разбития в координату после разбития расчётной
        области */
    int XAfterGridInit (int x)
    {
        return _mesh.IXw[x];
    }

    /* См. выше */
    int YAfterGridInit (int y)
    {
        return _mesh.IYw[y];
    }
    
    void BoundaryConditionType1Apply(BoundaryCondition bc)
    {
        /* учитывание разбиения сетки */
        int x1 = XAfterGridInit(bc.X1);
        int x2 = XAfterGridInit(bc.X2);
        int y1 = YAfterGridInit(bc.Y1);
        int y2 = YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        var localB = new Real[2]; // 'hat B'
        
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
            var m = new int[2];
            m[0] = b1 * _mesh.X.Count + a1;
            m[1] = b2 * _mesh.X.Count + a2;

            Slae.B[m[0]] = localB[0];
            Slae.B[m[1]] = localB[1];

            Slae.Di[m[0]] = 1;
            Slae.Di[m[1]] = 1;

            /* Обнуление строки */
            for (int idx = 0; idx < 2; idx++)
            {
                int ig0 = Slae.Ia[m[idx]];
                int ig1 = Slae.Ia[m[idx]+1];
                for (int i = ig0; i < ig1; i++)
                {
                    Slae.Mat[i] = 0;
                }
            }
            e1 = e2;
        }
    }
    
    void BoundaryConditionType2Apply(BoundaryCondition bc)
    {
        /* учитывание разбиения сетки */
        int x1 = XAfterGridInit(bc.X1);
        int x2 = XAfterGridInit(bc.X2);
        int y1 = YAfterGridInit(bc.Y1);
        int y2 = YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        var localB = new Real[2];
        
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

            int node1_num = b1 * _mesh.X.Count + a1;
            int node2_num = b2 * _mesh.X.Count + a2;
            Slae.B[node1_num] += localB[0];
            Slae.B[node2_num] += localB[1];

            e1 = e2;
        }
    }
    
    void BoundaryConditionType3Apply(BoundaryCondition bc)
    {
        /* учёт разбиения сетки */
        int x1 = XAfterGridInit(bc.X1);
        int x2 = XAfterGridInit(bc.X2);
        int y1 = YAfterGridInit(bc.Y1);
        int y2 = YAfterGridInit(bc.Y2);
        /*  */

        var num = bc.Num;

        var localB = new Real[2]; // 'hat B'
        var localA = new Real[2 ,2]; // 'hat A'
        
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
            m[0] = b1 * _mesh.X.Count + a1;
            m[1] = b2 * _mesh.X.Count + a2;

            Slae.B[m[0]] += localB[0];
            Slae.B[m[1]] += localB[1];

            Slae.Di[m[0]] += localA[0, 0];
            Slae.Di[m[1]] += localA[1, 1];

            for (int i = 0; i < 2; i++)
            {
                int beg = Slae.Ia[m[i]];
                for (int j = 0; j < 2; j++)
                {
                    // TODO: возможно есть лучше способ пропускать диагональные
                    // элементы
                    if (i == j)
                    {
                        continue;
                    }
                    int end = Slae.Ia[m[i] + 1] - 1;
                    while (beg < end)
                    {
                        int mid = (beg + end) / 2;
                        if (m[j] > Slae.Ja[mid])
                        {
                            beg = mid + 1;
                        }
                        else
                        {
                            end = mid;
                        }
                    }

                    if (Slae.Ja[beg] != m[j])
                    {
                        throw new Exception("Quick search failed");
                    }

                    Slae.Mat[beg] += localA[i, j];
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
            int x0 = i % (_mesh.X.Count - 1);
            int y0 = i / (_mesh.X.Count - 1);
            
            if (j == 0 || j == 1)
            {
                return y0 * _mesh.X.Count + x0 + j;
            } else if (j == 2 || j == 3)
            {
                return (y0 + 1) * _mesh.X.Count + x0 + (j - 2);
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
        Slae.Ia[0] = 0;
        /* формирование массивов ig jg по списку list */
        for (int i = 1; i < Slae.Ia.Length; i++)
        {
            Slae.Ia[i] = Slae.Ia[i-1] + list[i-1].Count;
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
    
    public int? GetSubdomNumAtElCoords (int x1, int y1)
    {
        foreach (var a in _mesh.SubDomains)
        {
            if (x1 >= _mesh.IXw[a.X1] && x1 < _mesh.IXw[a.X2] &&
                y1 >= _mesh.IYw[a.Y1] && y1 < _mesh.IYw[a.Y2]
            ) {
                return a.Num;
            }
        }
        
        return null;
    }

    int? GetSubdomNumAtPoint (Real x1, Real y1)
    {
        foreach (var a in _mesh.SubDomains)
        {
            if (x1 >= _mesh.Xw[a.X1] && x1 <= _mesh.Xw[a.X2] &&
                y1 >= _mesh.Yw[a.Y1] && y1 <= _mesh.Yw[a.Y2]
            ) {
                return a.Num;
            }
        }
        
        return null;
    }
    
    (int xi, int yi) GetElCoordsAtPoint(Real x, Real y)
    {
        int xi = -1;
        int yi = -1;
        for (int i = 0; i < _mesh.X.Count; i++)
        {
            if (_mesh.X[i] <= x && x <= _mesh.X[i+1])
            {
                xi = i;
                break;
            }
        }
        
        for (int i = 0; i < _mesh.Y.Count; i++)
        {
            if (_mesh.Y[i] <= y && y <= _mesh.Y[i+1])
            {
                yi = i;
                break;
            }
        }

        if (xi < 0 || yi < 0)
        {
            throw new Exception("Bad");
        }
        return (xi, yi);
    }
    
    public Real AnswerAt (Real x, Real y)
    {
        var num = GetSubdomNumAtPoint(x, y);
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

        var (xi, yi) = GetElCoordsAtPoint(x, y);

        Real hx = X[xi + 1] - X[xi];
        Real hy = Y[yi + 1] - Y[yi];

        var subdom = GetSubdomNumAtElCoords(xi, yi);
        if (subdom.HasValue)
        {
            var m = new int[4];
            m[0] = yi * X.Count + xi;
            m[1] = m[0] + 1;
            m[2] = (yi + 1) * X.Count + xi;
            m[3] = m[2] + 1;
            
            result =
            ( q[m[0]] * (X[xi + 1] - x) * (Y[yi + 1] - y)
            + q[m[1]] * (x - X[xi])     * (Y[yi + 1] - y)
            + q[m[2]] * (X[xi + 1] - x) * (y - Y[yi])
            + q[m[3]] * (x - X[xi])     * (y - Y[yi])) /hx/hy;
        }
        return result;
    }

    void GlobalMatrixBuildParallel ()
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

        for (int yi = 0; yi < _mesh.Y.Count - 1; yi++)
        {
            for (int xi = 0; xi < _mesh.X.Count - 1; xi++)
            {
                int targetNode = yi * _mesh.X.Count + xi;
                var dom1 = GetSubdomNumAtElCoords(xi-1, yi-1);
                var dom2 = GetSubdomNumAtElCoords(xi, yi-1);
                var dom3 = GetSubdomNumAtElCoords(xi-1, yi);
                var dom4 = GetSubdomNumAtElCoords(xi, yi);

                var r = new int[3];
                r[1] = targetNode - 1;
                r[0] = r[1] - _mesh.X.Count;
                r[2] = r[1] + _mesh.X.Count;

                var mr = new int[3];
                int beg = Slae.Ia[targetNode];
                int bound = Slae.Ia[targetNode + 1] - 1;
                for (int i = 0; i < 3; i++)
                {
                    int end = bound;
                    while (beg < end)
                    {
                        int mid = (beg + end) / 2;
                        if (r[i] > Slae.Ja[mid])
                        {
                            beg = mid + 1;
                        }
                        else
                        {
                            end = mid;
                        }
                    }

                    if (Slae.Ja[beg] != r[i])
                    {
                        throw new Exception("Quick search failed");
                    }

                    mr[i] = beg;
                    beg++;
                }

                var hx0 = _mesh.X[targetNode] - _mesh.X[targetNode - 1];
                var hx1 = _mesh.X[targetNode + 1] - _mesh.X[targetNode];
                var hy0 = _mesh.Y[targetNode] - _mesh.Y[r[0] + 1];
                var hy1 = _mesh.Y[r[2] + 1] - _mesh.Y[targetNode];

                if (dom1.HasValue)
                {
                    var lambda = GetLamdaAverage(dom1.Value, xi, yi);
                    var gamma = GetGammaAverage(dom1.Value, xi, yi);
                    // -1 1 -2 2
                    Slae.Mat[mr[0]] += lambda / 6 * (-hy0 / hx0 - hx1 / hy0);
                    Slae.Mat[mr[0] + 1] += lambda / 6 * (hy0 / hx0 - 2*hx1 / hy0);
                    Slae.Mat[mr[1]] += lambda / 6 * (-2*hy0 / hx0 + hx1 / hy0);
                    Slae.Mat[mr[1] + 1] += lambda / 6 * (2*hy0 / hx0 - 2*hx1 / hy0);
                }

            }
        }
    }

    void GlobalMatrixBuild ()
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

        for (int yi = 0; yi < _mesh.Y.Count - 1; yi++)
        {
            for (int xi = 0; xi < _mesh.X.Count - 1; xi++)
            {
                var subDom = GetSubdomNumAtElCoords(xi, yi);

                var m = new int[4];
                m[0] = yi * _mesh.X.Count + xi;
                m[1] = m[0] + 1;
                m[2] = (yi + 1) * _mesh.X.Count + xi;
                m[3] = m[2] + 1;

                if (subDom == null) continue;

                Real x0 = _mesh.X[xi];
                Real x1 = _mesh.X[xi + 1];
                Real y0 = _mesh.Y[yi];
                Real y1 = _mesh.Y[yi + 1];

                Real hy = y1 - y0;
                Real hx = x1 - x0;
                // Заменить на интеграл от биквадратичного разложения
                Real l_avg = GetLamdaAverage(subDom.Value, xi, yi);
                Real g_avg = GetGammaAverage(subDom.Value, xi, yi);
                var localB = LocalBVector(subDom.Value, x0, x1, y0, y1);

                /* нахождение в ja индексов элементов в al/au, куда
                    нужно добавить элементы локальных матриц */
                for (int i = 0; i < 4; i++)
                {
                    var v2 = l_avg/6 * (hy/hx * _localG1[i, i] + hx/hy * _localG2[i, i])
                        + g_avg/36 * hx*hy * _localM[i, i];
                    Slae.Di[m[i]] += v2;

                    int beg = Slae.Ia[m[i]];
                    for (int j = 0; j < 4; j++)
                    {
                        // TODO: пропуск
                        if (i == j)
                        {
                            continue;
                        }
                        int end = Slae.Ia[m[i] + 1] - 1;
                        while (beg < end)
                        {
                            int mid = (beg + end) / 2;
                            if (m[j] > Slae.Ja[mid])
                            {
                                beg = mid + 1;
                            }
                            else
                            {
                                end = mid;
                            }
                        }

                        if (Slae.Ja[beg] != m[j])
                        {
                            throw new Exception("Quick search failed");
                        }

                        v2 = l_avg/6 * (hy/hx * _localG1[i, j] + hx/hy * _localG2[i, j])
                            + g_avg/36 * hx*hy * _localM[i, j];
                        Slae.Mat[beg] += v2;
                        beg++;
                    }
                }

                /* добавление локальной правой части в слау */
                for (int i = 0; i < 4; i++)
                {
                    Slae.B[m[i]] += localB[i];
                }
            }
        }

        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < Slae.Di.Length; i++)
        {
            if (Slae.Di[i] == 0)
            {
                Slae.Di[i] = 1;
            }
        }
    }
    Real[] LocalBVector (
        int subDom,
        Real x0, Real x1,
        Real y0, Real y1
    ) {
        var localB = new Real[4];

        Real hy = y1 - y0;
        Real hx = x1 - x0;

        Real f1 = _funcs.F(subDom, x0, y0);
        Real f2 = _funcs.F(subDom, x1, y0);
        Real f3 = _funcs.F(subDom, x0, y1);
        Real f4 = _funcs.F(subDom, x1, y1);

        localB[0] = hx * hy / 36 * (4 * f1 + 2 * f2 + 2 * f3 +     f4);
        localB[1] = hx * hy / 36 * (2 * f1 + 4 * f2 +     f3 + 2 * f4);
        localB[2] = hx * hy / 36 * (2 * f1 +     f2 + 4 * f3 + 2 * f4);
        localB[3] = hx * hy / 36 * (    f1 + 2 * f2 + 2 * f3 + 4 * f4);

        return localB;
    }
}
