using Real = double;

interface TaskFuncs
{
    string Description { get; }

    Real Ug(int bcNum, Real x, Real y);
    Real Theta(int bcNum, Real x, Real y);
    Real Beta(int bcNum);
    Real uBeta(int bcNum, Real x, Real y);

    Real Lambda(int subdom, Real x, Real y);
    Real Gamma(int subdom, Real x, Real y);
    Real F(int subdom, Real x, Real y);
    Real Answer(int subdom, Real x, Real y);
}
