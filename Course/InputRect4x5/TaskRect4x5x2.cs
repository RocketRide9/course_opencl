using Real = double;

class TaskRect4x5x2 : TaskFuncs
{
    public string Description => "Прямоугольник 4на5 x2";

    public Real Answer(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real) x*x + y*y,
            _ => throw new ArgumentException("Неверный номер подобласти"),
        };
    }

    public Real Beta(int bcNum)
    {
        return bcNum switch
        {
            0 => (Real)0.5,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real F(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real)(x*x + y*y - 2),
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Gamma(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 1,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Lambda(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real)0.5,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Theta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 4,
            1 => -1,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real uBeta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 2*y + x*x + 25,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Ug(int bcNum, Real x, Real y)
    {
        return Answer(bcNum, x, y);
    }
}
