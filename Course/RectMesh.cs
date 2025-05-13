using Real = double;

public class RectMesh
{
    Subdomain[] _subDomains;
    public Subdomain[] SubDomains { get => _subDomains; }

    BoundaryCondition[] _boundaryConditions;
    public BoundaryCondition[] BoundaryConditions { get => _boundaryConditions; }

    RefineParams _refineParams;

    public Real[] Xw;
    public Real[] Yw;
    public Real[] X { get; private set; }
    public Real[] Y { get; private set; }
    public int[] IXw { get; private set; }
    public int[] IYw { get; private set; }

    public int nodesCount;
    public int feCount;

    public RectMesh(
        Real[] xAxis, Real[] yAxis,
        Subdomain[] subDomains,
        BoundaryCondition[] boundaryConditions
    ) {
        Xw = xAxis;
        Yw = yAxis;
        _subDomains = subDomains;
        _boundaryConditions = boundaryConditions;

        X = (Real[])Xw.Clone();
        Y = (Real[])Yw.Clone();
        IXw = Enumerable.Range(0, X.Length).ToArray();
        IYw = Enumerable.Range(0, Y.Length).ToArray();

        nodesCount = X.Length * Y.Length;
        feCount = (X.Length - 1) * (Y.Length - 1);
        
        _refineParams = new()
        {
            XSplitCount = Enumerable.Repeat(1, X.Length - 1).ToArray(),
            XStretchRatio = Enumerable.Repeat((Real)1, X.Length - 1).ToArray(),

            YSplitCount = Enumerable.Repeat(1, Y.Length - 1).ToArray(),
            YStretchRatio = Enumerable.Repeat((Real)1, Y.Length - 1).ToArray(),
        };
    }

    public int? GetSubdomNumAtElCoords (int x1, int y1)
    {
        foreach (var a in SubDomains)
        {
            if (x1 >= IXw[a.X1] && x1 < IXw[a.X2] &&
                y1 >= IYw[a.Y1] && y1 < IYw[a.Y2]
            ) {
                return a.Num;
            }
        }

        return null;
    }

    public int? GetSubdomNumAtPoint (Real x1, Real y1)
    {
        foreach (var a in SubDomains)
        {
            if (x1 >= Xw[a.X1] && x1 <= Xw[a.X2] &&
                y1 >= Yw[a.Y1] && y1 <= Yw[a.Y2]
            ) {
                return a.Num;
            }
        }

        return null;
    }

    public (int xi, int yi) GetElCoordsAtPoint(Real x, Real y)
    {
        int xi = -1;
        int yi = -1;
        for (int i = 0; i < X.Length; i++)
        {
            if (X[i] <= x && x <= X[i+1])
            {
                xi = i;
                break;
            }
        }

        for (int i = 0; i < Y.Length; i++)
        {
            if (Y[i] <= y && y <= Y[i+1])
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

    /* Перевод координаты x до разбития в координату после разбития расчётной
        области */
    public int XAfterGridInit (int x)
    {
        return IXw[x];
    }

    /* См. выше */
    public int YAfterGridInit (int y)
    {
        return IYw[y];
    }

    public void RefineDiv2()
    {
        for (int i = 0; i < _refineParams.XSplitCount.Length; i++)
        {
            _refineParams.XSplitCount[i] *= 2;
            _refineParams.XStretchRatio[i] = (Real)Math.Sqrt(_refineParams.XStretchRatio[i]);
        }
        for (int i = 0; i < _refineParams.YSplitCount.Length; i++)
        {
            _refineParams.YSplitCount[i] *= 2;
            _refineParams.YStretchRatio[i] = (Real)Math.Sqrt(_refineParams.YStretchRatio[i]);
        }

        Refine(_refineParams);
    }

    static Real FirstStepSize(Real stretch, int seg_count, Real gap)
    {
        Real sum;
        if (stretch != 1d)
        {
            sum = (Real)(1 - Math.Pow(stretch, seg_count)) / (1 - stretch);
        } else {
            sum = seg_count;
        }

        return gap / sum;
    }

    public void Refine(RefineParams refineParams)
    {
        _refineParams = refineParams;

        { // ось X
            var xLength = _refineParams.XSplitCount.Sum() + 1;
            IXw = new int[Xw.Length];
            X = new Real[xLength];

            IXw[0] = 0;
            X[0] = Xw[0];
            int xCount = 1;
            for (int i = 1; i < Xw.Length; i++)
            {
                Real gap = Xw[i] - Xw[i - 1];

                int seg_count = _refineParams.XSplitCount[i - 1];
                Real stretch = _refineParams.XStretchRatio[i - 1];

                var step = FirstStepSize(stretch, seg_count, gap);
                var step_n = step;
                var stretch_n = stretch;
                int idx = xCount - 1;
                for (int j = 0; j < seg_count - 1; j++)
                {
                    X[xCount] = X[idx] + step_n;
                    xCount++;
                    stretch_n *= stretch;
                    if (stretch != 1d)
                    {
                        step_n = step * (stretch_n - 1) / (stretch - 1);
                    } else {
                        step_n = step * (j + 2);
                    }
                }
                IXw[i] = xCount;
                X[xCount] = Xw[i];
                xCount++;
            }
        }

        { // ось Y
            var yLength = _refineParams.YSplitCount.Sum() + 1;
            IYw = new int[Yw.Length];
            Y = new Real[yLength];
            IYw[0] = 0;
            Y[0] = Yw[0];
            int yCount = 1;
            for (int i = 1; i < Yw.Length; i++)
            {
                Real gap = Yw[i] - Yw[i - 1];

                int seg_count = _refineParams.YSplitCount[i - 1];
                Real stretch = _refineParams.YStretchRatio[i - 1];

                var step = FirstStepSize(stretch, seg_count, gap);
                var step_n = step;
                var stretch_n = stretch;
                int idx = yCount - 1;
                for (int j = 0; j < seg_count - 1; j++)
                {
                    Y[yCount] = Y[idx] + step_n;
                    yCount++;
                    stretch_n *= stretch;
                    if (stretch != 1d)
                    {
                        step_n = step * (stretch_n - 1) / (stretch - 1);
                    } else {
                        step_n = step * (j + 2);
                    }
                }
                IYw[i] = yCount;
                Y[yCount] = Yw[i];
                yCount++;
            }
        }

        nodesCount = X.Length * Y.Length;
        feCount = (X.Length - 1) * (Y.Length - 1);
    }
}
