using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AdvUtils.Regression
{
    class SingularValueDecomposition
    {
        private Double[,] u; // left singular vectors
        private Double[,] v; // right singular vectors
        private Double[] s;  // singular values
        private int m;
        private int n;

        private int[] si; // sorting order

        private const Double eps = 2 * (1.11022302462515654042E-16);
        private const Double tiny = 1.493221789605150e-300;



        /// <summary>
        ///   Returns the singularity threshold.
        /// </summary>
        ///
        public Double Threshold
        {
            get { return Double.Epsilon * System.Math.Max(m, n) * s[0]; }
        }

        /// <summary>Constructs a new singular value decomposition.</summary>
        /// <param name="value">
        ///   The matrix to be decomposed.</param>
        public SingularValueDecomposition(Double[,] value)
        {
            if (value == null)
            {
                throw new ArgumentNullException("value", "Matrix cannot be null.");
            }

            //m should not less than n
            Double[,] a;
            m = value.GetLength(0); // rows
            n = value.GetLength(1); // cols

            if (m < n)
            {
                throw new System.Exception("rows should not less than cols in value matrix");
            }

            // Proceed anyway
            a = (Double[,])value.Clone();


            int nu = System.Math.Min(m, n);
            int ni = System.Math.Min(m + 1, n);
            s = new Double[ni];
            u = new Double[m, nu];
            v = new Double[n, n];
            Double[] e = new Double[n];
            Double[] work = new Double[m];

            // Will store ordered sequence of indices after sorting.
            si = new int[ni]; for (int i = 0; i < ni; i++) si[i] = i;


            // Reduce A to bidiagonal form, storing the diagonal elements in s and the super-diagonal elements in e.
            int nct = System.Math.Min(m - 1, n);
            int nrt = System.Math.Max(0, System.Math.Min(n - 2, m));
            for (int k = 0; k < System.Math.Max(nct, nrt); k++)
            {
                if (k < nct)
                {
                    // Compute the transformation for the k-th column and place the k-th diagonal in s[k].
                    // Compute 2-norm of k-th column without under/overflow.
                    s[k] = 0;
                    for (int i = k; i < m; i++)
                    {
                        s[k] = Hypotenuse(s[k], a[i, k]);
                    }

                    if (s[k] != 0)
                    {
                        if (a[k, k] < 0)
                            s[k] = -s[k];

                        for (int i = k; i < m; i++)
                            a[i, k] /= s[k];

                        a[k, k] += 1;
                    }

                    s[k] = -s[k];
                }

                for (int j = k + 1; j < n; j++)
                {
                    if ((k < nct) & (s[k] != 0))
                    {
                        // Apply the transformation.
                        Double t = 0;
                        for (int i = k; i < m; i++)
                        {
                            t += a[i, k] * a[i, j];
                        }

                        t = -t / a[k, k];

                        for (int i = k; i < m; i++)
                        {
                            a[i, j] += t * a[i, k];
                        }
                    }

                    // Place the k-th row of A into e for the subsequent calculation of the row transformation.
                    e[j] = a[k, j];
                }


                if (k < nct)
                {
                    // Place the transformation in U for subsequent back
                    // multiplication.
                    for (int i = k; i < m; i++)
                        u[i, k] = a[i, k];
                }

                if (k < nrt)
                {
                    // Compute the k-th row transformation and place the k-th super-diagonal in e[k].
                    // Compute 2-norm without under/overflow.
                    e[k] = 0;
                    for (int i = k + 1; i < n; i++)
                        e[k] = Hypotenuse(e[k], e[i]);

                    if (e[k] != 0)
                    {
                        if (e[k + 1] < 0)
                            e[k] = -e[k];

                        for (int i = k + 1; i < n; i++)
                            e[i] /= e[k];

                        e[k + 1] += 1;
                    }

                    e[k] = -e[k];
                    if ((k + 1 < m) & (e[k] != 0))
                    {
                        // Apply the transformation.
                        for (int i = k + 1; i < m; i++)
                            work[i] = 0;

                        int k1 = k + 1;
                        for (int i = k1; i < m; i++)
                        {
                            for (int j = k1; j < n; j++)
                            {
                                work[i] += e[j] * a[i, j];
                            }
                        }

                        for (int j = k1; j < n; j++)
                        {
                            Double t = -e[j] / e[k1];
                            for (int i = k1; i < m; i++)
                            {
                                a[i, j] += t * work[i];
                            }
                        }
                    }

                    // Place the transformation in V for subsequent back multiplication.
                    for (int i = k + 1; i < n; i++)
                        v[i, k] = e[i];
                }
            }

            // Set up the final bidiagonal matrix or order p.
            int p = System.Math.Min(n, m + 1);
            if (nct < n) s[nct] = a[nct, nct];
            if (m < p) s[p - 1] = 0;
            if (nrt + 1 < p) e[nrt] = a[nrt, p - 1];
            e[p - 1] = 0;

            //generate U.
            for (int j = nct; j < nu; j++)
            {
                for (int i = 0; i < m; i++)
                    u[i, j] = 0;
                u[j, j] = 1;
            }

            for (int k = nct - 1; k >= 0; k--)
            {
                if (s[k] != 0)
                {
                    for (int j = k + 1; j < nu; j++)
                    {
                        Double t = 0;
                        for (int i = k; i < m; i++)
                        {
                            t += u[i, k] * u[i, j];
                        }

                        t = -t / u[k, k];

                        for (int i = k; i < m; i++)
                        {
                            u[i, j] += t * u[i, k];
                        }
                    }

                    for (int i = k; i < m; i++)
                    {
                        u[i, k] = -1.0 * u[i, k];
                    }

                    u[k, k] = 1 + u[k, k];
                    for (int i = 0; i < k - 1; i++)
                        u[i, k] = 0;
                }
                else
                {
                    for (int i = 0; i < m; i++)
                        u[i, k] = 0;
                    u[k, k] = 1;
                }
            }


            //generate V.
            for (int k = n - 1; k >= 0; k--)
            {
                if ((k < nrt) & (e[k] != 0))
                {
                    // TODO: The following is a pseudo correction to make SVD
                    //  work on matrices with n > m (less rows than columns).

                    // For the proper correction, compute the decomposition of the
                    //  transpose of A and swap the left and right eigenvectors

                    // Original line:
                    //   for (int j = k + 1; j < nu; j++)
                    // Pseudo correction:
                    //   for (int j = k + 1; j < n; j++)

                    for (int j = k + 1; j < n; j++) // pseudo-correction
                    {
                        Double t = 0;
                        for (int i = k + 1; i < n; i++)
                        {
                            t += v[i, k] * v[i, j];
                        }

                        t = -t / v[k + 1, k];

                        for (int i = k + 1; i < n; i++)
                        {
                            v[i, j] += t * v[i, k];
                        }
                    }
                }

                for (int i = 0; i < n; i++)
                {
                    v[i, k] = 0;
                }
                v[k, k] = 1;
            }

            // Main iteration loop for the singular values.
            int pp = p - 1;
            int iter = 0;
            while (p > 0)
            {
                int k, kase;

                // Here is where a test for too many iterations would go.

                // This section of the program inspects for
                // negligible elements in the s and e arrays.  On
                // completion the variables kase and k are set as follows.

                // kase = 1     if s(p) and e[k-1] are negligible and k<p
                // kase = 2     if s(k) is negligible and k<p
                // kase = 3     if e[k-1] is negligible, k<p, and
                //              s(k), ..., s(p) are not negligible (qr step).
                // kase = 4     if e(p-1) is negligible (convergence).

                for (k = p - 2; k >= -1; k--)
                {
                    if (k == -1)
                        break;

                    if (System.Math.Abs(e[k]) <=
                       tiny + eps * (System.Math.Abs(s[k]) + System.Math.Abs(s[k + 1])))
                    {
                        e[k] = 0;
                        break;
                    }
                }

                if (k == p - 2)
                {
                    kase = 4;
                }
                else
                {
                    int ks;
                    for (ks = p - 1; ks >= k; ks--)
                    {
                        if (ks == k)
                            break;

                        Double t = (ks != p ? System.Math.Abs(e[ks]) : 0) +
                                   (ks != k + 1 ? System.Math.Abs(e[ks - 1]) : 0);
                        if (System.Math.Abs(s[ks]) <= tiny + eps * t)
                        {
                            s[ks] = 0;
                            break;
                        }
                    }

                    if (ks == k)
                        kase = 3;
                    else if (ks == p - 1)
                        kase = 1;
                    else
                    {
                        kase = 2;
                        k = ks;
                    }
                }

                k++;

                // Perform the task indicated by kase.
                switch (kase)
                {
                    // Deflate negligible s(p).
                    case 1:
                        {
                            Double f = e[p - 2];
                            e[p - 2] = 0;
                            for (int j = p - 2; j >= k; j--)
                            {
                                Double t = Hypotenuse(s[j], f);
                                Double cs = s[j] / t;
                                Double sn = f / t;
                                s[j] = t;
                                if (j != k)
                                {
                                    f = -sn * e[j - 1];
                                    e[j - 1] = cs * e[j - 1];
                                }

                                for (int i = 0; i < n; i++)
                                {
                                    t = cs * v[i, j] + sn * v[i, p - 1];
                                    v[i, p - 1] = -sn * v[i, j] + cs * v[i, p - 1];
                                    v[i, j] = t;
                                }
                            }
                        }
                        break;

                    // Split at negligible s(k).
                    case 2:
                        {
                            Double f = e[k - 1];
                            e[k - 1] = 0;
                            for (int j = k; j < p; j++)
                            {
                                Double t = Hypotenuse(s[j], f);
                                Double cs = s[j] / t;
                                Double sn = f / t;
                                s[j] = t;
                                f = -sn * e[j];
                                e[j] = cs * e[j];

                                for (int i = 0; i < m; i++)
                                {
                                    t = cs * u[i, j] + sn * u[i, k - 1];
                                    u[i, k - 1] = -sn * u[i, j] + cs * u[i, k - 1];
                                    u[i, j] = t;
                                }
                            }
                        }
                        break;

                    // Perform one qr step.
                    case 3:
                        {
                            // Calculate the shift.
                            Double scale = System.Math.Max(System.Math.Max(System.Math.Max(System.Math.Max(System.Math.Abs(s[p - 1]), System.Math.Abs(s[p - 2])), System.Math.Abs(e[p - 2])), System.Math.Abs(s[k])), System.Math.Abs(e[k]));
                            Double sp = s[p - 1] / scale;
                            Double spm1 = s[p - 2] / scale;
                            Double epm1 = e[p - 2] / scale;
                            Double sk = s[k] / scale;
                            Double ek = e[k] / scale;
                            Double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2;
                            Double c = (sp * epm1) * (sp * epm1);
                            double shift = 0;

                            if ((b != 0) | (c != 0))
                            {
                                if (b < 0)
                                    shift = -System.Math.Sqrt(b * b + c);
                                else
                                    shift = System.Math.Sqrt(b * b + c);

                                shift = c / (b + shift);
                            }

                            Double f = (sk + sp) * (sk - sp) + (Double)shift;
                            Double g = sk * ek;

                            // Chase zeros.
                            for (int j = k; j < p - 1; j++)
                            {
                                Double t = Hypotenuse(f, g);
                                Double cs = f / t;
                                Double sn = g / t;
                                if (j != k) e[j - 1] = t;
                                f = cs * s[j] + sn * e[j];
                                e[j] = cs * e[j] - sn * s[j];
                                g = sn * s[j + 1];
                                s[j + 1] = cs * s[j + 1];

                                for (int i = 0; i < n; i++)
                                {
                                    /*t = cs * v[i, j] + sn * v[i, j + 1];
                                    v[i, j + 1] = -sn * v[i, j] + cs * v[i, j + 1];
                                    v[i, j] = t;*/

                                    Double vij = v[i, j]; // *vj;
                                    Double vij1 = v[i, j + 1]; // *vj1;

                                    t = cs * vij + sn * vij1;
                                    v[i, j + 1] = -sn * vij + cs * vij1;
                                    v[i, j] = t;
                                }

                                t = Hypotenuse(f, g);
                                cs = f / t;
                                sn = g / t;
                                s[j] = t;
                                f = cs * e[j] + sn * s[j + 1];
                                s[j + 1] = -sn * e[j] + cs * s[j + 1];
                                g = sn * e[j + 1];
                                e[j + 1] = cs * e[j + 1];

                                if (j < m - 1)
                                {
                                    for (int i = 0; i < m; i++)
                                    {
                                        /* t = cs * u[i, j] + sn * u[i, j + 1];
                                         u[i, j + 1] = -sn * u[i, j] + cs * u[i, j + 1];
                                         u[i, j] = t;*/

                                        Double uij = u[i, j]; // *uj;
                                        Double uij1 = u[i, j + 1]; // *uj1;

                                        t = cs * uij + sn * uij1;

                                        u[i, j + 1] = -sn * uij + cs * uij1;
                                        u[i, j] = t;
                                    }
                                }

                            }

                            e[p - 2] = f;
                            iter = iter + 1;
                        }
                        break;

                    // Convergence.
                    case 4:
                        {
                            // Make the singular values positive.
                            if (s[k] <= 0)
                            {
                                s[k] = (s[k] < 0 ? -s[k] : 0);

                                for (int i = 0; i <= pp; i++)
                                    v[i, k] = -v[i, k];

                            }

                            // Order the singular values.
                            while (k < pp)
                            {
                                if (s[k] >= s[k + 1])
                                    break;

                                Double t = s[k];
                                s[k] = s[k + 1];
                                s[k + 1] = t;

                                int ti = si[k];
                                si[k] = si[k + 1];
                                si[k + 1] = ti;

                                if (k < n - 1)
                                {
                                    for (int i = 0; i < n; i++)
                                    {
                                        t = v[i, k + 1];
                                        v[i, k + 1] = v[i, k];
                                        v[i, k] = t;
                                    }
                                }

                                if (k < m - 1)
                                {
                                    for (int i = 0; i < m; i++)
                                    {
                                        t = u[i, k + 1];
                                        u[i, k + 1] = u[i, k];
                                        u[i, k] = t;
                                    }
                                }

                                k++;
                            }

                            iter = 0;
                            p--;
                        }
                        break;
                }
            }
        }


     

        /// <summary>
        ///   Solves a linear equation system of the form Ax = b.
        /// </summary>
        /// <param name="value">The b from the equation Ax = b.</param>
        /// <returns>The x from equation Ax = b.</returns>
        public Double[] Solve(Double[] value)
        {
            // Additionally an important property is that if there does not exists a solution
            // when the matrix A is singular but replacing 1/Li with 0 will provide a solution
            // that minimizes the residue |AX -Y|. SVD finds the least squares best compromise
            // solution of the linear equation system. Interestingly SVD can be also used in an
            // over-determined system where the number of equations exceeds that of the parameters.

            // L is a diagonal matrix with non-negative matrix elements having the same
            // dimension as A, Wi ? 0. The diagonal elements of L are the singular values of matrix A.

            //singularity threshold
            Double e = this.Threshold;

            var Y = value;

            // Create L*, which is a diagonal matrix with elements
            //    L*i = 1/Li  if Li = e, else 0, 
            // where e is the so-called singularity threshold.

            // In other words, if Li is zero or close to zero (smaller than e),
            // one must replace 1/Li with 0. The value of e depends on the precision
            // of the hardware. This method can be used to solve linear equations
            // systems even if the matrices are singular or close to singular.


            int scols = s.Length;

            var Ls = new Double[scols, scols];
            for (int i = 0; i < s.Length; i++)
            {
                if (System.Math.Abs(s[i]) <= e)
                    Ls[i, i] = 0;
                else Ls[i, i] = 1 / s[i];
            }

            //(V x L*) x Ut x Y
            var VL = v.Multiply(Ls);

            //(V x L* x Ut) x Y
            int urows = u.GetLength(0);
            int vrows = v.GetLength(0);
            var VLU = new Double[vrows, urows];
            for (int i = 0; i < vrows; i++)
            {
                for (int j = 0; j < urows; j++)
                {
                    Double sum = 0;
                    for (int k = 0; k < scols; k++)
                        sum += VL[i, k] * u[j, k];
                    VLU[i, j] = sum;
                }
            }

            //(V x L* x Ut x Y)
            return VLU.Multiply(Y);
        }

        
        public double Hypotenuse(double a, double b)
        {
            double r = 0.0;
            double absA = System.Math.Abs(a);
            double absB = System.Math.Abs(b);

            if (absA > absB)
            {
                r = b / a;
                r = absA * System.Math.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = absB * System.Math.Sqrt(1 + r * r);
            }

            return r;
        }
    }
}
