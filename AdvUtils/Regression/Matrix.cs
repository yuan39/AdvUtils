using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AdvUtils.Regression
{
    public static class Matrix
    {
        public static double[,] Multiply(this double[,] a, double[,] b)
        {
            double[,] result = new double[a.GetLength(0), b.GetLength(1)];

            int n = a.GetLength(1);
            int m = result.GetLength(0); //a.GetLength(0);
            int p = result.GetLength(1); //b.GetLength(1);

            double[] Bcolj = new double[n];
            for (int j = 0; j < p; j++)
            {
                for (int k = 0; k < Bcolj.Length; k++)
                    Bcolj[k] = b[k, j];

                for (int i = 0; i < m; i++)
                {
                    double s = 0;
                    for (int k = 0; k < Bcolj.Length; k++)
                        s += a[i, k] * Bcolj[k];
                    result[i, j] = s;
                }
            }

            return result;
        }

        /// <summary>
        ///   Gets the maximum element in a vector.
        /// </summary>
        /// 
        public static double Max(double[] values)
        {
            double imax = 0;
            double max = values[0];
            for (int i = 1; i < values.Length; i++)
            {
                if (values[i].CompareTo(max) > 0)
                {
                    max = values[i];
                    imax = i;
                }
            }
            return max;
        }

        public static double[] Multiply(this double[,] matrix, double[] columnVector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            if (cols != columnVector.Length)
                throw new Exception(
                    "Vector must have the same length as columns in the matrix.");

            double[] r = new double[rows];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columnVector.Length; j++)
                    r[i] += matrix[i, j] * columnVector[j];

            return r;
        }
    }

}
