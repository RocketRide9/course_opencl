using Real = double;

using System.Numerics;

using Types;

namespace SlaeBuilder;

public interface ISlaeBuilder
{
    RectMesh Mesh { get; }
    GlobalMatrixImplType GlobalMatrixImpl { get; set; }
    static abstract ISlaeBuilder Construct(RectMesh mesh, TaskFuncs funcs);
    (Matrix, Real[]) Build();
}

public enum GlobalMatrixImplType
{
    OpenCL,
    OpenCLV2,
    Host,
    HostParallel,
    HostV2
}

static class Shared {
    public static readonly Real[,] LocalG1 = {
        { 2, -2,  1, -1},
        {-2,  2, -1,  1},
        { 1, -1,  2, -2},
        {-1,  1, -2,  2},
    };
    public static readonly Real[,] LocalG2 = {
        { 2,  1, -2, -1},
        { 1,  2, -1, -2},
        {-2, -1,  2,  1},
        {-1, -2,  1,  2},
    };
    public static readonly Real[,] LocalM = {
        {4, 2, 2, 1},
        {2, 4, 1, 2},
        {2, 1, 4, 2},
        {1, 2, 2, 4},
    };
    
    // https://stackoverflow.com/a/16893641
    public static Real Add(ref Real location1, Real value)
    {
        Real newCurrentValue = location1; // non-volatile read, so may be stale
        while (true)
        {
            Real currentValue = newCurrentValue;
            Real newValue = currentValue + value;
            newCurrentValue = Interlocked.CompareExchange(ref location1, newValue, currentValue);
            if (newCurrentValue.Equals(currentValue))
            {
                return newValue;
            }
        }
    }
    
    // не функция из csharp потому что мне её ещё нужно самому реализовать в OpenCL
    public static int QFind<T> (T[] @where, int start, int end, T what)
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
    
    public static int LFind<T> (T[] where, T what, int start)
    where T: unmanaged, INumber<T>
    {
        while (@where[start] != what) start++;
        return start;
    }
}
