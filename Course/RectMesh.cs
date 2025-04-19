using Real = float;

public class RectMesh
{
    Subdomain[] _subDomains;
    public Subdomain[] SubDomains { get => _subDomains; }

    BoundaryCondition[] _boundaryConditions;
    public BoundaryCondition[] BoundaryConditions { get => _boundaryConditions; }
    
    RefineParams? _refineParams;
    public RefineParams? RefineParams { get => _refineParams; }
    
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
    }

    public void RefineDiv2()
    {
        var rparams = RefineParams.Value;
        for (int i = 0; i < rparams.XSplitCount.Length; i++)
        {
            rparams.XSplitCount[i] *= 2;
            rparams.XStretchRatio[i] = (Real)Math.Sqrt(rparams.XStretchRatio[i]);
        }
        for (int i = 0; i < rparams.YSplitCount.Length; i++)
        {
            rparams.YSplitCount[i] *= 2;
            rparams.YStretchRatio[i] = (Real)Math.Sqrt(rparams.YStretchRatio[i]);
        }

        Refine(rparams);
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
            var xLength = _refineParams.Value.XSplitCount.Sum() + 1;
            IXw = new int[Xw.Length];
            X = new Real[xLength];
    
            IXw[0] = 0;
            X[0] = Xw[0];
            int xCount = 1;
            for (int i = 1; i < Xw.Length; i++)
            {
                Real gap = Xw[i] - Xw[i - 1];
                
                int seg_count = _refineParams.Value.XSplitCount[i - 1];
                Real stretch = _refineParams.Value.XStretchRatio[i - 1];
                
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
            var yLength = _refineParams.Value.YSplitCount.Sum() + 1;
            IYw = new int[Yw.Length];
            Y = new Real[yLength];
            IYw[0] = 0;
            Y[0] = Yw[0];
            int yCount = 1;        
            for (int i = 1; i < Yw.Length; i++)
            {
                Real gap = Yw[i] - Yw[i - 1];
                
                int seg_count = _refineParams.Value.YSplitCount[i - 1];
                Real stretch = _refineParams.Value.YStretchRatio[i - 1];
                
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
