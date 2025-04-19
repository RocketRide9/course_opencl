using System.Numerics;
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

    void ArraySerialize(Real[] arr, string fileName)
    {
        var stream = new StreamWriter(fileName);
        stream.Write(
            string.Join(
                "\n",
                arr.Select(
                    e =>
                    {
                        var bytes = BitConverter.GetBytes(e);
                        var i = BitConverter.ToInt64(bytes, 0);
                        return "0x" + i.ToString("X");
                    }
                )
            )
        );
        stream.Close();
    }
    
    void ArraySerialize(int[] arr, string fileName)
    {
        var stream = new StreamWriter(fileName);
        stream.Write(
            string.Join(
                "\n",
                arr
            )
        );
        stream.Close();
    }
    
    public void Serialize()
    {
        ArraySerialize(Mat, "mat.txt");
        ArraySerialize(Di, "di.txt");
        Console.WriteLine($"{Di[4500]}");
        ArraySerialize(B, "b.txt");
        ArraySerialize(Ia, "ia.txt");
        ArraySerialize(Ja, "ja.txt");
    }
}

public struct Slae1
{
    public SparkOCL.DeprecatedArray<Real> Mat;
    public SparkOCL.DeprecatedArray<Real> Di;
    public SparkOCL.DeprecatedArray<Real> B;
    public SparkOCL.DeprecatedArray<int> Ia;
    public SparkOCL.DeprecatedArray<int> Ja;
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
