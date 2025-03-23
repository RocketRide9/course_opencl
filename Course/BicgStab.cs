#define HOST_PARALLEL

using Real = float;
using Quasar.Native;


public class BiCGStabMkl
{
    SparkOCL.Array<Real> _mat;
    SparkOCL.Array<Real> _di;
    SparkOCL.Array<Real> _b;
    SparkOCL.Array<int> _ia;
    SparkOCL.Array<int> _ja;

    int _maxIter;
    Real _eps;
    SparkOCL.Array<Real> _x;
    
    SparkOCL.Array<Real> r;
    SparkOCL.Array<Real> di_inv;
    SparkOCL.Array<Real> y;
    SparkOCL.Array<Real> z;
    SparkOCL.Array<Real> ks;
    SparkOCL.Array<Real> kt;
    SparkOCL.Array<Real> r_hat;
    SparkOCL.Array<Real> p;
    SparkOCL.Array<Real> nu;
    SparkOCL.Array<Real> h;
    SparkOCL.Array<Real> s;
    SparkOCL.Array<Real> t;
    
    public BiCGStabMkl(
        SparkOCL.Array<Real> Mat,
        SparkOCL.Array<Real> Di,
        SparkOCL.Array<Real> B,
        SparkOCL.Array<int> Ia,
        SparkOCL.Array<int> Ja,

        SparkOCL.Array<Real> x0,
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
        
        var zeros = Enumerable.Repeat((Real)0, _b.Count).ToArray();
        r       = new SparkOCL.Array<Real>(zeros);
        r_hat   = new SparkOCL.Array<Real>(zeros);
        p       = new SparkOCL.Array<Real>(zeros);
        nu      = new SparkOCL.Array<Real>(zeros);
        h       = new SparkOCL.Array<Real>(zeros);
        s       = new SparkOCL.Array<Real>(zeros);
        t       = new SparkOCL.Array<Real>(zeros);
        di_inv  = new SparkOCL.Array<Real>(zeros);
        y       = new SparkOCL.Array<Real>(zeros);
        z       = new SparkOCL.Array<Real>(zeros);
        ks      = new SparkOCL.Array<Real>(zeros);
        kt      = new SparkOCL.Array<Real>(zeros);
    }

    // y *= x
    static void Vmul(SparkOCL.Array<Real> y, SparkOCL.Array<Real> x)
    {
        for (int i = 0; i < y.Count; i++)
        {
            y[i] *= x[i];
        }
    }
    // y = y*(-1/2)
    static void Rsqrt(SparkOCL.Array<Real> y)
    {
        for (int i = 0; i < y.Count; i++)
        {
            y[i] = (Real)( 1 / Math.Sqrt(y[i]) );
        }
    }
    
    public (SparkOCL.Array<Real> ans, Real rr, Real pp, int iter) Solve()
    {
        // precond
        _di.CopyTo(di_inv);
        Rsqrt(di_inv);
        // 1.
        MSRMul(_mat, _di, _ia, _ja, _x.Count, _x, t);
        _b.CopyTo(r);
        BLAS.axpy(_x.Count, -1, t.AsSpan(), r.AsSpan());
        // 2.
        r.CopyTo(r_hat);
        // 3.
        Real pp = (Real)BLAS.dot(_x.Count, r.AsSpan(), r.AsSpan()); // r_hat * r
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
            MSRMul(_mat, _di, _ia, _ja, _x.Count, y, nu);
            
            // 3.
            Real rnu = (Real)BLAS.dot(_x.Count, r_hat.AsSpan(), nu.AsSpan());
            Real alpha = pp / rnu;

            // 4.
            _x.CopyTo(h);
            BLAS.axpy(_x.Count, alpha, y.AsSpan(), h.AsSpan());
            
            // 5.
            r.CopyTo(s);
            BLAS.axpy(_x.Count, -alpha, nu.AsSpan(), s.AsSpan());

            // 6.
            Real ss = (Real)BLAS.dot(_x.Count, s.AsSpan(), s.AsSpan());
            if (ss < _eps)
            {
                _x.Dispose();
                _x = h;
                break;
            }
            
            // 7.
            s.CopyTo(ks);
            Vmul(ks, di_inv);
            ks.CopyTo(z);
            Vmul(z, di_inv);

            // 8.
            MSRMul(_mat, _di, _ia, _ja, _x.Count, z, t);

            // 9.
            t.CopyTo(kt);
            Vmul(kt, di_inv);
            
            Real ts = (Real)BLAS.dot(_x.Count, ks.AsSpan(), kt.AsSpan());
            Real tt = (Real)BLAS.dot(_x.Count, kt.AsSpan(), kt.AsSpan());
            Real w = ts / tt;

            // 10.
            h.CopyTo(_x);
            BLAS.axpy(_x.Count, w, z.AsSpan(), _x.AsSpan());

            // 11.
            s.CopyTo(r);
            BLAS.axpy(_x.Count, -w, t.AsSpan(), r.AsSpan());

            // 12.
            rr = (Real)BLAS.dot(_x.Count, r.AsSpan(), r.AsSpan());
            if (rr < _eps)
            {
                break;
            }
            
            // 13-14
            Real pp1 = (Real)BLAS.dot(_x.Count, r.AsSpan(), r_hat.AsSpan());
            Real beta = (pp1 / pp) * (alpha / w);
            
            // 15.
            BLAS.axpy(_x.Count, -w, nu.AsSpan(), p.AsSpan());
            BLAS.scal(_x.Count, beta, p.AsSpan());
            BLAS.axpy(_x.Count, 1, r.AsSpan(), p.AsSpan());

            pp = pp1;
        }

        MSRMul(_mat, _di, _ja, _ia, _x.Count, _x, t);
        _b.CopyTo(r);
        BLAS.axpy(_x.Count, -1, t.AsSpan(), r.AsSpan());
        rr = (Real)BLAS.dot(r.Count, r.AsSpan(), r.AsSpan());

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
            SparkOCL.Array<Real> mat,
            SparkOCL.Array<Real> di,
            SparkOCL.Array<int> aptr,
            SparkOCL.Array<int> jptr,
            int n,
            SparkOCL.Array<Real> v,
            SparkOCL.Array<Real> res)
        {
            MyFor(0, n, i => {
                int start = aptr[i];
                int stop = aptr[i + 1];
                Real dot = di[i] * v[i];
                for (int a = start; a < stop; a++)
                {
                    dot += mat[a] * v[jptr[a]];
                }
                res[i] = dot;
            });
        }
        
        public static void MyFor(int i0, int i1, Action<int> iteration)
        {
#if HOST_PARALLEL
            var partitioner = System.Collections.Concurrent.Partitioner.Create(i0, i1);
            Parallel.ForEach(partitioner, (range, state) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    iteration(i);
                }
            });
#else
            for (int i = i0; i < i1; i++)
            {
                iteration(i);
            }        
#endif
        }
}
