#define HOST_PARALLEL

using Real = double;
using Quasar.Native;

public class BiCGStabMkl
{
    Real[] _mat;
    Real[] _di;
    Real[] _b;
    int[] _ia;
    int[] _ja;

    int _maxIter;
    Real _eps;
    Real[] _x;

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

    public BiCGStabMkl(
        Real[] Mat,
        Real[] Di,
        Real[] B,
        int[] Ia,
        int[] Ja,

        Real[] x0,
        int maxIter,
        Real eps)
    {
        _maxIter = maxIter;
        _eps = eps;

        _mat = Mat;
        _di = Di;
        _b = B;
        _ia = Ia;
        _ja = Ja;

        _x = x0;

        var zeros = Enumerable.Repeat((Real)0, _b.Length).ToArray();
        r       = [.. zeros];
        r_hat   = [.. zeros];
        p       = [.. zeros];
        nu      = [.. zeros];
        h       = [.. zeros];
        s       = [.. zeros];
        t       = [.. zeros];
        di_inv  = [.. zeros];
        y       = [.. zeros];
        z       = [.. zeros];
        ks      = [.. zeros];
        kt      = [.. zeros];
    }

    // y *= x
    static void Vmul(Real[] y, Real[] x)
    {
        var partitioner = System.Collections.Concurrent.Partitioner.Create(0, y.Length);
        Parallel.ForEach(partitioner, (range, state) =>
        {
            for (int i = range.Item1; i < range.Item2; i++)
            {
                y[i] *= x[i];
            }
        });
    }
    // y = y*(-1/2)
    static void Rsqrt(Real[] y)
    {
        var partitioner = System.Collections.Concurrent.Partitioner.Create(0, y.Length);
        Parallel.ForEach(partitioner, (range, state) =>
        {
            for (int i = range.Item1; i < range.Item2; i++)
            {
                y[i] = (Real)(1 / Math.Sqrt(y[i]));
            }
        });
    }

    public (Real[] ans, Real rr, Real pp, int iter) Solve()
    {
        // precond
        _di.CopyTo(di_inv, 0);
        Rsqrt(di_inv);
        // 1.
        MSRMul(_mat, _di, _ia, _ja, _x.Length, _x, t);
        _b.CopyTo(r, 0);
        BLAS.axpy(_x.Length, -1, t.AsSpan(), r.AsSpan());
        // 2.
        r.CopyTo(r_hat, 0);
        // 3.
        Real pp = (Real)BLAS.dot(_x.Length, r.AsSpan(), r.AsSpan()); // r_hat * r
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
            MSRMul(_mat, _di, _ia, _ja, _x.Length, y, nu);

            // 3.
            Real rnu = (Real)BLAS.dot(_x.Length, r_hat.AsSpan(), nu.AsSpan());
            Real alpha = pp / rnu;

            // 4.
            _x.CopyTo(h, 0);
            BLAS.axpy(_x.Length, alpha, y.AsSpan(), h.AsSpan());

            // 5.
            r.CopyTo(s, 0);
            BLAS.axpy(_x.Length, -alpha, nu.AsSpan(), s.AsSpan());

            // 6.
            Real ss = (Real)BLAS.dot(_x.Length, s.AsSpan(), s.AsSpan());
            if (ss < _eps)
            {
                h.CopyTo(_x, 0);
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
            MSRMul(_mat, _di, _ia, _ja, _x.Length, z, t);

            // 9.
            t.CopyTo(kt, 0);
            Vmul(kt, di_inv);

            Real ts = (Real)BLAS.dot(_x.Length, ks.AsSpan(), kt.AsSpan());
            Real tt = (Real)BLAS.dot(_x.Length, kt.AsSpan(), kt.AsSpan());
            Real w = ts / tt;

            // 10.
            h.CopyTo(_x, 0);
            BLAS.axpy(_x.Length, w, z.AsSpan(), _x.AsSpan());

            // 11.
            s.CopyTo(r, 0);
            BLAS.axpy(_x.Length, -w, t.AsSpan(), r.AsSpan());

            // 12.
            rr = (Real)BLAS.dot(_x.Length, r.AsSpan(), r.AsSpan());
            if (rr < _eps)
            {
                break;
            }

            // 13-14
            Real pp1 = (Real)BLAS.dot(_x.Length, r.AsSpan(), r_hat.AsSpan());
            Real beta = (pp1 / pp) * (alpha / w);

            // 15.
            BLAS.axpy(_x.Length, -w, nu.AsSpan(), p.AsSpan());
            BLAS.scal(_x.Length, beta, p.AsSpan());
            BLAS.axpy(_x.Length, 1, r.AsSpan(), p.AsSpan());

            pp = pp1;
        }

        MSRMul(_mat, _di, _ia, _ja, _x.Length, _x, t);
        _b.CopyTo(r, 0);
        BLAS.axpy(_x.Length, -1, t.AsSpan(), r.AsSpan());
        rr = (Real)BLAS.dot(r.Length, r.AsSpan(), r.AsSpan());

        return (_x, rr, pp, iter);
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

    public static void MSRMul(
        Real[] mat,
        Real[] di,
        int[] ia,
        int[] ja,
        int n,
        Real[] v,
        Real[] res)
    {
        var partitioner = System.Collections.Concurrent.Partitioner.Create(0, n);
        Parallel.ForEach(partitioner, (range, state) =>
        {
            for (int i = range.Item1; i < range.Item2; i++)
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
        });
    }
}
