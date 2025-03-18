using Real = float;

class TaskBook : TaskFuncs
{
    public string Description => "Область из книги с примером";

    public Real Answer(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => x,
            1 or 2 => (Real) (1.8 + 0.1*x),
            _ => throw new ArgumentException("Неверный номер подобласти"),
        };
    }

    public Real Beta(int bcNum)
    {
        return bcNum switch
        {
            0 => 1,
            1 => 2,
            2 => (Real)0.5,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real F(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 2*x,
            1 => (Real) (1.8 + 0.1*x),
            2 => 0,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Gamma(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 2,
            1 => 1,
            2 => 0,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Lambda(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 1,
            1 => 10,
            2 => 10,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Theta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 1,
            1 => 0,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real uBeta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => x,
            1 => (Real) (1.8 + 0.1*x),
            2 => -1,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Ug(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 2,
            1 => (Real) (0.1*x + 1.8),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }    

}
