using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AdvUtils.Regression
{
    public class SingularPolynomialRegression
    {
        double[] coefficients;

        public double[] Coefficients
        {
            get { return coefficients; }
        }

        public void Regress(double[] arrX, double[] arrY, int dimension)
        {
            int n = dimension + 1;
            double[,] Guass = new double[n, n + 1];
            for (int i = 0; i < n; i++)
            {
                int j;
                for (j = 0; j < n; j++)
                {
                    //Compute sum of x^(i+j)
                    Guass[i, j] = SumArr(arrX, i + j);
                }

                //Compute sum of (x^i)*y
                Guass[i, j] = SumArr(arrX, i, arrY, 1);
            }

            //y=a0+a1*x+a2*x*x+...+aN * x^N
            //a0 = coefficients[0]
            //a1 = coefficients[1]
            //...
            //aN = coefficients[N]

            coefficients = ComputGauss(Guass, n);

        }

        public static double SumArr(double[] arr, int n)
        {
            double s = 0;
            for (int i = 0; i < arr.Length; i++)
            {
                if (arr[i] != 0 || n != 0)
                    s = s + Math.Pow(arr[i], n);
                else
                    s = s + 1;
            }
            return s;
        }
        public static double SumArr(double[] arr1, int n1, double[] arr2, int n2)
        {
            double s = 0;
            for (int i = 0; i < arr1.Length; i++)
            {
                if ((arr1[i] != 0 || n1 != 0) && (arr2[i] != 0 || n2 != 0))
                    s = s + Math.Pow(arr1[i], n1) * Math.Pow(arr2[i], n2);
                else
                    s = s + 1;
            }
            return s;

        }

        public static double[] ComputGauss(double[,] Guass, int n)
        {
            int i, j;
            int k, m;
            double temp;
            double max;
            double s;
            double[] x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            for (j = 0; j < n; j++)
            {
                max = 0;

                k = j;
                for (i = j; i < n; i++)
                {
                    if (Math.Abs(Guass[i, j]) > max)
                    {
                        max = Guass[i, j];
                        k = i;
                    }
                }

                if (k != j)
                {
                    for (m = j; m < n + 1; m++)
                    {
                        temp = Guass[j, m];
                        Guass[j, m] = Guass[k, m];
                        Guass[k, m] = temp;

                    }
                }

                if (0 == max)
                {
                    return x;
                }


                for (i = j + 1; i < n; i++)
                {
                    s = Guass[i, j];
                    for (m = j; m < n + 1; m++)
                    {
                        Guass[i, m] = Guass[i, m] - Guass[j, m] * s / (Guass[j, j]);
                    }
                }
            }

            for (i = n - 1; i >= 0; i--)
            {
                s = 0;
                for (j = i + 1; j < n; j++)
                {
                    s = s + Guass[i, j] * x[j];
                }

                x[i] = (Guass[i, n] - s) / Guass[i, i];

            }

            return x;
        }
    }
}
