using Real = float;

public struct ProblemParams 
{
    public Real eps { get; set; }
    public int maxIter { get; set; }
}

public struct Subdomain 
{
    // номера подобласти начинаются с 0
    public int Num;
    public int X1;
    public int X2;
    public int Y1;
    public int Y2;
}

public struct RefineParams
{
    public int[] XSplitCount { get; set; }
    public Real[] XStretchRatio { get; set; }
    
    public int[] YSplitCount { get; set; }
    public Real[] YStretchRatio { get; set; }
}

public struct BoundaryCondition
{
    // нач. с 0
    public int Num;
    // тип краевого условия (первый, второй ...) нумеруется с 1
    public int Type;
    public int X1;
    public int X2;
    public int Y1;
    public int Y2;
}

public struct Slae2
{
    public Real[] Mat;
    public Real[] Di;
    public Real[] B;
    public int[] Ia;
    public int[] Ja;
}

public struct Slae1
{
    public SparkOCL.Array<Real> Mat;
    public SparkOCL.Array<Real> Di;
    public SparkOCL.Array<Real> B;
    public SparkOCL.Array<int> Ia;
    public SparkOCL.Array<int> Ja;
}

public struct ComputationalDomain
{
    public Real[] xAxis;
    public Real[] yAxis;
    public Subdomain[] subDomains;
}

interface IElement
{
    int NumberOfDofs { get; }
    Real[,] LocalGMatrix (Real hx, Real hy, Real gamma);
    Real[,] LocalMMatrix (Real hx, Real hy, Real gamma);
}
