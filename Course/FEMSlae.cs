using Real = float;

class FEMSlae
{
    Slae _slae;
    public Slae slae { get => _slae; }
    RectMesh _mesh;
    public RectMesh mesh { get => _mesh; }

    TaskFuncs _funcs;
    
    void MeshRefine(RefineParams refineParams)
    {
        _mesh.Refine(refineParams);
    }
    
    public FEMSlae(RectMesh mesh, TaskFuncs funcs, RefineParams? refineParams)
    {
        _mesh = mesh;
        _slae = new Slae();
        _funcs = funcs;

        if (refineParams != null)
        {
            MeshRefine(refineParams.Value);
        }
        
        GlobalMatrixInit();
        GlobalMatrixBuild();
        BoundaryConditionsApply();
    }
    
    void GlobalMatrixInit()
    {
        GlobalMatrixPortraitCompose();

        _slae.Mat = new SparkCL.Memory<Real>(Enumerable.Repeat((Real)0, slae.Ja.Count).ToArray());
        _slae.Di = new SparkCL.Memory<Real>(Enumerable.Repeat((Real)0, slae.Ia.Count - 1).ToArray());
        _slae.B =  new SparkCL.Memory<Real>(Enumerable.Repeat((Real)0, slae.Ia.Count - 1).ToArray());
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

            slae.B[m[0]] = localB[0];
            slae.B[m[1]] = localB[1];

            slae.Di[m[0]] = 1;
            slae.Di[m[1]] = 1;

            /* Обнуление строки */
            for (int idx = 0; idx < 2; idx++)
            {
                int ig0 = slae.Ia[m[idx]];
                int ig1 = slae.Ia[m[idx]+1];
                for (int i = ig0; i < ig1; i++)
                {
                    slae.Mat[i] = 0;
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
            slae.B[node1_num] += localB[0];
            slae.B[node2_num] += localB[1];

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

            slae.B[m[0]] += localB[0];
            slae.B[m[1]] += localB[1];

            slae.Di[m[0]] += localA[0, 0];
            slae.Di[m[1]] += localA[1, 1];

            for (int i = 0; i < 2; i++)
            {
                int beg = slae.Ia[m[i]];
                for (int j = 0; j < 2; j++)
                {
                    // TODO: возможно есть лучше способ пропускать диагональные
                    // элементы
                    if (i == j)
                    {
                        continue;
                    }
                    int end = slae.Ia[m[i] + 1] - 1;
                    while (beg < end)
                    {
                        int mid = (beg + end) / 2;
                        if (m[j] > slae.Ja[mid])
                        {
                            beg = mid + 1;
                        }
                        else
                        {
                            end = mid;
                        }
                    }

                    if (slae.Ja[beg] != m[j])
                    {
                        throw new Exception("Quick search failed");
                    }

                    slae.Mat[beg] += localA[i, j];
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
        
        
        HashSet<int>[] list = new HashSet<int>[mesh.nodesCount];
        for (int i = 0; i < list.Length; i++)
        {
            list[i] = [];
        }

        /* цикл по всем конечным элементам */
        for (int ielem = 0; ielem < mesh.feCount; ielem++)
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

        _slae.Ia = new SparkCL.Memory<int>(list.Length + 1);
        slae.Ia[0] = 0;
        /* формирование массивов ig jg по списку list */
        for (int i = 1; i < slae.Ia.Count; i++)
        {
            slae.Ia[i] = slae.Ia[i-1] + list[i-1].Count;
        }
        _slae.Ja = new SparkCL.Memory<int>(_slae.Ia[list.Length]);
        for (var i = 0; i < list.Length; i++)
        {
            var row = list[i].Order().ToArray();
            for (int j = _slae.Ia[i]; j < _slae.Ia[i+1]; j++)
            {
                _slae.Ja[j] = row[j-_slae.Ia[i]];
            }
        }
    }
    
    int? GetSubdomNumAtElCoords (int x1, int y1)
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
                
                if (subDom != null)
                {
                    Real x0 = _mesh.X[xi];
                    Real x1 = _mesh.X[xi + 1];
                    Real y0 = _mesh.Y[yi];
                    Real y1 = _mesh.Y[yi + 1];
                    
                    Real h_y = y1 - y0;
                    Real h_x = x1 - x0;
                    // Заменить на интеграл от биквадратичного разложения
                    Real gamma_average = GetGammaAverage((int)subDom, xi, yi);
                    Real lamda_average = GetLamdaAverage((int)subDom, xi, yi);
                    var localG = LocalGMatrix(h_x, h_y, lamda_average);
                    var localM = LocalMMatrix(h_x, h_y, gamma_average);
                    var localB = LocalBVector((int)subDom, x0, x1, y0, y1);

                    /* нахождение в ja индексов элементов в al/au, куда
                        нужно добавить элементы локальных матриц */
                    for (int i = 0; i < 4; i++)
                    {
                        slae.Di[m[i]] += localG[i, i] + localM[i, i];
                        int beg = slae.Ia[m[i]];
                        for (int j = 0; j < 4; j++)
                        {
                            // TODO: пропуск
                            if (i == j)
                            {
                                continue;
                            }
                            int end = slae.Ia[m[i] + 1] - 1;
                            while (beg < end)
                            {
                                int mid = (beg + end) / 2;
                                if (m[j] > slae.Ja[mid])
                                {
                                    beg = mid + 1;
                                }
                                else
                                {
                                    end = mid;
                                }
                            }

                            if (slae.Ja[beg] != m[j])
                            {
                                throw new Exception("Quick search failed");
                            }

                            slae.Mat[beg] += localG[i, j] + localM[i, j];
                            beg++;
                        }
                    }
                    
                    /* добавление локальной правой части в слау */
                    for (int i = 0; i < 4; i++)
                    {
                        slae.B[m[i]] += localB[i];
                    }
                }
            }
        }

        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < slae.Di.Count; i++)
        {
            if (slae.Di[i] == 0)
            {
                slae.Di[i] = 1;
            }
        }
    }
    
    Real[,] LocalGMatrix(Real hx, Real hy, Real gamma)
    {
        Real[,] G = new Real[4, 4];
        
        Real GElem (Real k1, Real k2)
        {
            return gamma * (hy / hx / k1 + hx / hy / k2);
        }
                  G[0, 0] = GElem ( 3,  3);
        G[1, 0] = G[0, 1] = GElem (-3,  6);
        G[2, 0] = G[0, 2] = GElem ( 6, -3);
        G[3, 0] = G[0, 3] = GElem (-6, -6);
        
                  G[1, 1] = GElem ( 3,  3);
        G[2, 1] = G[1, 2] = GElem (-6, -6);
        G[3, 1] = G[1, 3] = GElem ( 6, -3);
        
                  G[2, 2] = GElem ( 3,  3);
        G[3, 2] = G[2, 3] = GElem (-3,  6);
    
                  G[3, 3] = GElem ( 3,  3);
        
        return G;
    }

    Real[,] LocalMMatrix(Real hx, Real hy, Real gamma)
    {
        Real[,] M = new Real[4, 4];
        
        Real MElem (Real k)
        {
            return gamma * hx * hy * k / 36;
        }
        
                  M[0, 0] = MElem (4);
        M[1, 0] = M[0, 1] = MElem (2);
        M[2, 0] = M[0, 2] = MElem (2);
        M[3, 0] = M[0, 3] = MElem (1);
            
                  M[1, 1] = MElem (4);
        M[2, 1] = M[1, 2] = MElem (1);
        M[3, 1] = M[1, 3] = MElem (2);
            
                  M[2, 2] = MElem (4);
        M[3, 2] = M[2, 3] = MElem (2);
            
                  M[3, 3] = MElem (4);

        return M;
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
