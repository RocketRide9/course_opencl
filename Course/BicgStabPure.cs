#define HOST_PARALLEL
#define USE_BLAS

using Quasar.Native;
using SparkAlgos;
using System.Collections.Concurrent;
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

    // y *= x
    static unsafe void Vmul(Span<Real> y, ReadOnlySpan<Real> x)
    {
        if (x.Length != y.Length)
        {
            throw new ArgumentException("Vectors must have the same length");
        }
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, y.Length);
        fixed(Real* _p_y = y)
        fixed(Real* _p_x = x)
        {
            var p_y = _p_y;
            var p_x = _p_x;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    p_y[i] *= p_x[i];
                }
            });

        }
#else
        for (int i = 0; i < y.Length; i++)
        {
            y[i] *= x[i];
        }
#endif
    }

    // y = y*(-1/2)
    static unsafe void Rsqrt(Span<Real> y)
    {
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, y.Length);
        fixed(Real* _p_y = y)
        {
            var p_y = _p_y;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    p_y[i] = (Real)(1 / Math.Sqrt(p_y[i]));
                }
            });

        }
#else
        for (int i = 0; i < y.Length; i++)
        {
            y[i] = (Real)(1 / Math.Sqrt(y[i]));
        }
#endif
    }

    // y += alpha*x
    static void Axpy(Real alpha, ReadOnlySpan<Real> x, Span<Real> y)
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
    static Real Dot(ReadOnlySpan<Real> x, ReadOnlySpan<Real> y)
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
    static void Scale(Real alpha, Span<Real> y)
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

    public static unsafe void MSRMul(
        ReadOnlySpan<Real> mat,
        ReadOnlySpan<Real> di,
        ReadOnlySpan<int> ia,
        ReadOnlySpan<int> ja,
        int n,
        ReadOnlySpan<Real> v,
        Span<Real> res)
    {
#if HOST_PARALLEL
        var partitioner = Partitioner.Create(0, ia.Length);
        fixed(Real* _p_mat = mat)
        fixed(Real* _p_di = di)
        fixed(int* _p_ia = ia)
        fixed(int* _p_ja = ja)
        fixed(Real* _p_v = v)
        fixed(Real* _p_res = res)
        {
            var p_mat = _p_mat;
            var p_di = _p_di;
            var p_ia = _p_ia;
            var p_ja = _p_ja;
            var p_v = _p_v;
            var p_res = _p_res;
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    int start = p_ia[i];
                    int stop = p_ia[i + 1];
                    Real dot = p_di[i] * p_v[i];
                    for (int a = start; a < stop; a++)
                    {
                        dot += p_mat[a] * p_v[p_ja[a]];
                    }
                    p_res[i] = dot;
                }
            });
        }
#else
        fixed(Real* p_mat = mat)
        fixed(Real* p_di = di)
        fixed(int* p_ia = ia)
        fixed(int* p_ja = ja)
        fixed(Real* p_v = v)
        fixed(Real* p_res = res)
        for (int i = 0; i < ia.Length - 1; i++)
        {
            int start = p_ia[i];
            int stop = p_ia[i + 1];
            Real dot = p_di[i] * p_v[i];
            for (int a = start; a < stop; a++)
            {
                dot += p_mat[a] * p_v[p_ja[a]];
            }
            p_res[i] = dot;
        }
#endif
    }

    // Выделить память для временных массивов
    // n - длина каждого массива
    public void AllocateTemps(int n)
    {
        if (n != _n)
        {
            _n = n;

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
    }

    // x используется как начальное приближение, туда же попадёт ответ
    public (Real rr, Real pp, int iter) Solve(SlaeRef slae, Span<Real> x)
    {
        AllocateTemps(x.Length);

        var _mat = slae.Mat;
        var _di  = slae.Di;
        var _b   = slae.B;
        var _ia  = slae.Ia;
        var _ja  = slae.Ja;

        var r       = this.r     .AsSpan();
        var r_hat   = this.r_hat .AsSpan();
        var p       = this.p     .AsSpan();
        var nu      = this.nu    .AsSpan();
        var h       = this.h     .AsSpan();
        var s       = this.s     .AsSpan();
        var t       = this.t     .AsSpan();
        var di_inv  = this.di_inv.AsSpan();
        var y       = this.y     .AsSpan();
        var z       = this.z     .AsSpan();
        var ks      = this.ks    .AsSpan();
        var kt      = this.kt    .AsSpan();
        
        
        // precond
        _di.CopyTo(di_inv);
        Rsqrt(di_inv);
        // 1.
        MSRMul(_mat, _di, _ia, _ja, _n, x, t);
        _b.CopyTo(r);
        Axpy(-1, t, r);
        // BLAS.axpy(_n, -1, t, r);
        // 2.
        r.CopyTo(r_hat);
        // 3.
        Real pp = Dot(r, r); // r_hat * r
        // 4.
        r.CopyTo(p);

        int iter = 0;
        Real rr;
        for (; iter < _maxIter; iter++)
        {
            // 1.
            p.CopyTo(y);
            Vmul(y, di_inv);
            Vmul(y, di_inv);

            // 2.
            MSRMul(_mat, _di, _ia, _ja, _n, y, nu);

            // 3.
            Real rnu = Dot(r_hat, nu);
            Real alpha = pp / rnu;

            // 4.
            x.CopyTo(h);
            Axpy(alpha, y, h);
            // BLAS.axpy(_n, alpha, y, h);

            // 5.
            r.CopyTo(s);
            Axpy(-alpha, nu, s);
            // BLAS.axpy(_n, -alpha, nu, s);

            // 6.
            Real ss = Dot(s, s);
            if (ss < _eps)
            {
                h.CopyTo(x);
                // _x.Dispose();
                // _x = h;
                break;
            }

            // 7.
            s.CopyTo(ks);
            Vmul(ks, di_inv);
            ks.CopyTo(z);
            Vmul(z, di_inv);

            // 8.
            MSRMul(_mat, _di, _ia, _ja, _n, z, t);

            // 9.
            t.CopyTo(kt);
            Vmul(kt, di_inv);

            Real ts = Dot(ks, kt);
            Real tt = Dot(kt, kt);
            Real w = ts / tt;

            // 10.
            h.CopyTo(x);
            Axpy(w, z, x);
            // BLAS.axpy(_n, w, z, _x);

            // 11.
            s.CopyTo(r);
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
        _b.CopyTo(r);
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
