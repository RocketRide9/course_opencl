using Real = double;

namespace Types;
public interface Matrix
{
    int Size { get; }
    // TODO: нужно для предобуславливания.
    // Надо придумать что-то более разумное
    Span<Real> Di { get; }

    SparkAlgos.Types.Matrix GetComputeMatrix();
    void Mul(ReadOnlySpan<Real> vec, Span<Real> res);
    // не нулевый, потому что так проще
    IEnumerable<Real> FlatNonZero();
}
