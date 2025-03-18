using Real = float;

class TaskSquare : TaskFuncs
{
    public string Description => "Квадрат";

    public Real Answer(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => x + y,
            _ => throw new ArgumentException("Неверный номер подобласти"),
        };
    }

    public Real Beta(int bcNum)
    {
        return bcNum switch
        {
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real F(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 6*(x + y),
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Gamma(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 36,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Lambda(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 6,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Theta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => -6,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real uBeta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Ug(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => x + 1,
            1 => 16 + y,
            2 => x + 25,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }    

}
