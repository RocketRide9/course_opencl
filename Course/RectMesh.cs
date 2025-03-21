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
    public List<Real> X { get; }
    public List<Real> Y { get; }
    public List<int> IXw { get; }
    public List<int> IYw { get; }
    
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

        X = new(Xw);
        Y = new(Yw);
        IXw = Enumerable.Range(0, X.Count).ToList();
        IYw = Enumerable.Range(0, Y.Count).ToList();
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

        IXw.Clear();
        X.Clear();
        IXw.Add(X.Count);
        X.Add(Xw[0]);
        for (int i = 1; i < Xw.Length; i++)
        {
            Real gap = Xw[i] - Xw[i - 1];
            
            int seg_count = _refineParams.Value.XSplitCount[i - 1];
            Real stretch = _refineParams.Value.XStretchRatio[i - 1];
            
            Real step = FirstStepSize(stretch, seg_count, gap);
            Real step_n = step;
            Real stretch_n = stretch;
            int idx = X.Count - 1;
            for (int j = 0; j < seg_count - 1; j++)
            {
                X.Add(X[idx] + step_n);
                stretch_n *= stretch;
                if (stretch != 1d)
                {
                    step_n = step * (stretch_n - 1) / (stretch - 1);
                } else {
                    step_n = step * (j + 2);
                }
            }
            IXw.Add(X.Count);
            X.Add(Xw[i]);
        }
        
        IYw.Clear();
        Y.Clear();
        IYw.Add(Y.Count);
        Y.Add(Yw[0]);
        for (int i = 1; i < Yw.Length; i++)
        {
            Real gap = Yw[i] - Yw[i - 1];
            
            int seg_count = _refineParams.Value.YSplitCount[i - 1];
            Real stretch = _refineParams.Value.YStretchRatio[i - 1];
            
            Real step = FirstStepSize(stretch, seg_count, gap);
            Real step_n = step;
            Real stretch_n = stretch;
            int idx = Y.Count - 1;
            for (int j = 0; j < seg_count - 1; j++)
            {
                Y.Add(Y[idx] + step_n);
                stretch_n *= stretch;
                if (stretch != 1d)
                {
                    step_n = step * (stretch_n - 1) / (stretch - 1);
                } else {
                    step_n = step * (j + 2);
                }
            }
            IYw.Add(Y.Count);
            Y.Add(Yw[i]);
        }
        
        nodesCount = X.Count * Y.Count;
        feCount = (X.Count - 1) * (Y.Count - 1);
    }
}
