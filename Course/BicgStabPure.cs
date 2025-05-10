// #define HOST_PARALLEL
#define USE_BLAS

using Quasar.Native;
using SparkAlgos;
using System.Collections.Concurrent;
using System.Runtime.CompilerServices;
using Real = double;

public class BiCGStabPure
{
    int _maxIter;
    Real _eps;

    int _n = 0; // размерность СЛАУ
    Real[] r;
    Real[] di_inv;
    Real[] y;
    Real[] z;
    Real[] ks;
    Real[] kt;
    Real[] r_hat;
    Real[] p;
    Real[] nu;
    Real[] h;
    Real[] s;
    Real[] t;

    public BiCGStabPure(
        int maxIter,
        Real eps)
    {
        _maxIter = maxIter;
        _eps = eps;

        r =      [];
        di_inv = [];
        y =      [];
        z =      [];
        ks =     [];
        kt =     [];
        r_hat =  [];
        p =      [];
        nu =     [];
        h =      [];
        s =      [];
        t =      [];
    }

//     [MethodImpl(MethodImplOptions.AggressiveInlining)]
//     public static void MyFor(int i0, int i1, Action<int> iteration)
//     {
// #if HOST_PARALLEL
//         var partitioner = System.Collections.Concurrent.Partitioner.Create(i0, i1);
//         Parallel.ForEach(partitioner, (range, state) =>
//         {
//             for (int i = range.Item1; i < range.Item2; i++)
//             {
//                 iteration(i);
//             }
//         });
// #else
//         for (int i = i0; i < i1; i++)
//         {
//             iteration(i);
//         }        
// #endif
//     }
    
    // y *= x
    static void Vmul(Real[] y, Real[] x)
    {
        if (x.Length != y.Length)
        {
            throw new ArgumentException("Vectors must have the same length");
        }

        for (int i = 0; i < y.Length; i++)
        {
            y[i] *= x[i];
        }
        // MyFor(0, y.Length, (i) =>
        // {
        //     y[i] *= x[i];
        // });
    }
    // y = y*(-1/2)
    static void Rsqrt(Real[] y)
    {
        for (int i = 0; i < y.Length; i++)
        {
            y[i] = (Real)(1 / Math.Sqrt(y[i]));
        }
        // MyFor(0, y.Length, (i) =>
        // {
        //     y[i] = (Real)(1 / Math.Sqrt(y[i]));
        // });
    }
    // y += alpha*x
    static void Axpy(Real alpha, Real[] x, Real[] y)
    {
        if (x.Length != y.Length)
        {
            throw new ArgumentException("Vectors must have the same length");
        }
        #if USE_BLAS
            BLAS.axpy(x.Length, alpha, x, y);
        #else
            for (int i = 0; i < y.Length; i++)
            {
                y[i] += (Real)(alpha * x[i]);
            }
        #endif
    }
    // x·y
    static Real Dot(Real[] x, Real[] y)
    {
        if (x.Length != y.Length)
        {
            throw new ArgumentException("Vectors must have the same length");
        }
        #if USE_BLAS
            return (Real)BLAS.dot(x.Length, x, y);
        #else
            Real sum = 0;
            for (int i = 0; i < y.Length; i++)
            {
                sum += x[i] * y[i];
            }
            return sum;
        #endif
    }
    // y_i = alpha * y[i]
    static void Scale(Real alpha, Real[] y)
    {
        #if USE_BLAS
            BLAS.scal(y.Length, alpha, y);
        #else
            for (int i = 0; i < y.Length; i++)
            {
                y[i] *= alpha;
            }
        #endif
    }

    public static void MSRMul(
        Real[] mat,
        Real[] di,
        int[] ia,
        int[] ja,
        int n,
        Real[] v,
        Real[] res)
    {
        for (int i = 0; i < ia.Length - 1; i++)
        {
            int start = ia[i];
            int stop = ia[i + 1];
            Real dot = di[i] * v[i];
            for (int a = start; a < stop; a++)
            {
                dot += mat[a] * v[ja[a]];
            }
            res[i] = dot;
        }
    }

    // x используется как начальное приближение, туда же попадёт ответ
    public (Real rr, Real pp, int iter) Solve(Slae2 slae, Real[] x)
    {
        if (x.Length != _n)
        {
            _n = x.Length;

            r       = new Real[_n];
            r_hat   = new Real[_n];
            p       = new Real[_n];
            nu      = new Real[_n];
            h       = new Real[_n];
            s       = new Real[_n];
            t       = new Real[_n];
            di_inv  = new Real[_n];
            y       = new Real[_n];
            z       = new Real[_n];
            ks      = new Real[_n];
            kt      = new Real[_n];
        }

        var _mat = slae.Mat;
        var _di  = slae.Di;
        var _b   = slae.B;
        var _ia  = slae.Ia;
        var _ja  = slae.Ja;

        // precond
        _di.CopyTo(di_inv, 0);
        Rsqrt(di_inv);
        // 1.
        MSRMul(_mat, _di, _ia, _ja, _n, x, t);
        _b.CopyTo(r, 0);
        Axpy(-1, t, r);
        // BLAS.axpy(_n, -1, t, r);
        // 2.
        r.CopyTo(r_hat, 0);
        // 3.
        Real pp = Dot(r, r); // r_hat * r
        // 4.
        r.CopyTo(p, 0);

        int iter = 0;
        Real rr;
        for (; iter < _maxIter; iter++)
        {
            // 1.
            p.CopyTo(y, 0);
            Vmul(y, di_inv);
            Vmul(y, di_inv);

            // 2.
            MSRMul(_mat, _di, _ia, _ja, _n, y, nu);

            // 3.
            Real rnu = Dot(r_hat, nu);
            Real alpha = pp / rnu;

            // 4.
            x.CopyTo(h, 0);
            Axpy(alpha, y, h);
            // BLAS.axpy(_n, alpha, y, h);

            // 5.
            r.CopyTo(s, 0);
            Axpy(-alpha, nu, s);
            // BLAS.axpy(_n, -alpha, nu, s);

            // 6.
            Real ss = Dot(s, s);
            if (ss < _eps)
            {
                h.CopyTo(x, 0);
                // _x.Dispose();
                // _x = h;
                break;
            }

            // 7.
            s.CopyTo(ks, 0);
            Vmul(ks, di_inv);
            ks.CopyTo(z, 0);
            Vmul(z, di_inv);

            // 8.
            MSRMul(_mat, _di, _ia, _ja, _n, z, t);

            // 9.
            t.CopyTo(kt, 0);
            Vmul(kt, di_inv);

            Real ts = Dot(ks, kt);
            Real tt = Dot(kt, kt);
            Real w = ts / tt;

            // 10.
            h.CopyTo(x, 0);
            Axpy(w, z, x);
            // BLAS.axpy(_n, w, z, _x);

            // 11.
            s.CopyTo(r, 0);
            Axpy(-w, t, r);
            // BLAS.axpy(_n, -w, t, r);

            // 12.
            rr = Dot(r, r);
            if (rr < _eps)
            {
                break;
            }

            // 13-14
            Real pp1 = Dot(r, r_hat);
            Real beta = (pp1 / pp) * (alpha / w);

            // 15.
            Axpy(-w, nu, p);
            // BLAS.axpy(_n, -w, nu, p);
            Scale(beta, p);
            // BLAS.scal(_n, beta, p);
            // BLAS.axpy(_n, 1, r, p);
            Axpy(1, r, p);

            pp = pp1;
        }

        MSRMul(_mat, _di, _ia, _ja, _n, x, t);
        _b.CopyTo(r, 0);
        Axpy(-1, t, r);
        // BLAS.axpy(_x.Length, -1, t, r);
        rr = Dot(r, r);

        return (rr, pp, iter);
    }

    #if false
    public void SolveAndBreakdown()
    {
        var sw_host = new Stopwatch();
        sw_host.Start();
        var (rr, pp, iter) = Solve();
        sw_host.Stop();

        var x = slae.x;
        Real max_err = Math.Abs(x[0] - slae.ans[0]);
        for (int i = 0; i < (int)x.Length; i++)
        {
            var err = Math.Abs(x[i] - slae.ans[i]);
            if (err > max_err)
            {
                max_err = err;
            }
        }

        Console.WriteLine("Решение с MKL");
        Console.WriteLine($"rr = {rr}");
        Console.WriteLine($"pp = {pp}");
        Console.WriteLine($"max err. = {max_err}");
        Console.WriteLine($"Итераций: {iter}");
        Console.WriteLine($"Вычисления на хосте: {sw_host.ElapsedMilliseconds}мс");
    }
    #endif
}
