using Real = double;

class TaskRect4x5 : TaskFuncs
{
    public string Description => "Прямоугольник 4на5 с экспонентой";

    public Real Answer(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real)Math.Exp(x+y) + x,
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
            0 => (Real)((y*y-1)*Math.Exp(x+y) + x*y*y),
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Gamma(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => y*y,
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
            0 => (Real)(Math.Exp(4+y) + 1)/2,
            1 => -(Real)(Math.Exp(1+y) + 1)/2,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real uBeta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => (Real)(2*Math.Exp(x+5) + x),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Ug(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => (Real)(Math.Exp(x+y) + x),
            1 => (Real)(Math.Exp(x+y) + x),
            2 => (Real)(Math.Exp(x+y) + x),
            3 => (Real)(Math.Exp(x+y) + x),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

}
