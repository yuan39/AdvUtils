using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AdvUtils.Regression
{
    public class IterativeReweightedLeastSquares
    {
        private LogisticRegression regression;

        private int parameterCount;

        private double[,] hessian;
        private double[] gradient;
        private double[] previous;

        private bool computeStandardErrors = true;
        private SingularValueDecomposition decomposition;


        /// <summary>
        ///   Gets the previous values for the coefficients which were
        ///   in place before the last learning iteration was performed.
        /// </summary>
        ///
        public double[] Previous { get { return previous; } }

        /// <summary>
        ///   Gets the current values for the coefficients.
        /// </summary>
        ///
        public double[] Solution { get { return regression.Coefficients; } }

        /// <summary>
        ///   Gets the Hessian matrix computed in
        ///   the last Newton-Raphson iteration.
        /// </summary>
        ///
        public double[,] Hessian { get { return hessian; } }

        /// <summary>
        ///   Gets the Gradient vector computed in
        ///   the last Newton-Raphson iteration.
        /// </summary>
        ///
        public double[] Gradient { get { return gradient; } }

        /// <summary>
        ///   Gets the total number of parameters in the model.
        /// </summary>
        ///
        public int Parameters { get { return parameterCount; } }


        /// <summary>
        ///   Gets or sets a value indicating whether standard
        ///   errors should be computed in the next iteration.
        /// </summary>
        /// <value>
        /// 	<c>true</c> to compute standard errors; otherwise, <c>false</c>.
        /// </value>
        ///
        public bool ComputeStandardErrors
        {
            get { return computeStandardErrors; }
            set { computeStandardErrors = value; }
        }


        /// <summary>
        ///   Constructs a new Iterative Reweighted Least Squares.
        /// </summary>
        ///
        /// <param name="regression">The regression to estimate.</param>
        ///
        public IterativeReweightedLeastSquares(LogisticRegression regression)
        {
            this.regression = regression;

            this.parameterCount = regression.Coefficients.Length;

            this.hessian = new double[parameterCount, parameterCount];
            this.gradient = new double[parameterCount];
        }

        /// <summary>
        ///   Runs one iteration of the Reweighted Least Squares algorithm.
        /// </summary>
        /// <param name="input">The input data.</param>
        /// <param name="outputs">The outputs associated with each input vector.</param>
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        ///
        public double Run(double[][] input, double[][] outputs)
        {
            if (outputs[0].Length != 1)
                throw new ArgumentException("Function must have a single output.", "outputs");

            double[] output = new double[outputs.Length];
            for (int i = 0; i < outputs.Length; i++)
                output[i] = outputs[i][0];

            return Run(input, output);
        }

        /// <summary>
        ///   Runs one iteration of the Reweighted Least Squares algorithm.
        /// </summary>
        /// <param name="inputs">The input data.</param>
        /// <param name="outputs">The outputs associated with each input vector.</param>
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        ///
        public double Run(double[][] inputs, double[] outputs)
        {
            // Regress using Iteratively Reweighted Least Squares estimation.

            // References:
            //  - Bishop, Christopher M.; Pattern Recognition
            //    and Machine Learning. Springer; 1st ed. 2006.


            // Initial definitions and memory allocations
            int N = inputs.Length;

            double[][] design = new double[N][];
            double[] errors = new double[N];
            double[] weights = new double[N];
            double[] coefficients = this.regression.Coefficients;
            double[] deltas;

            // Compute the regression matrix
            for (int i = 0; i < inputs.Length; i++)
            {
                double[] row = design[i] = new double[parameterCount];

                row[0] = 1; // for intercept
                for (int j = 0; j < inputs[i].Length; j++)
                    row[j + 1] = inputs[i][j];
            }


            // Compute errors and weighing matrix
            for (int i = 0; i < inputs.Length; i++)
            {
                double y = regression.Compute(inputs[i]);

                // Calculate error vector
                errors[i] = y - outputs[i];

                // Calculate weighting matrix
                weights[i] = y * (1.0 - y);
            }


            // Reset Hessian matrix and gradient
            for (int i = 0; i < gradient.Length; i++)
            {
                gradient[i] = 0;
                for (int j = 0; j < gradient.Length; j++)
                    hessian[i, j] = 0;
            }


            // (Re-) Compute error gradient
            for (int j = 0; j < design.Length; j++)
                for (int i = 0; i < gradient.Length; i++)
                    gradient[i] += design[j][i] * errors[j];

            // (Re-) Compute weighted "Hessian" matrix
            for (int k = 0; k < weights.Length; k++)
            {
                double[] rk = design[k];

                for (int j = 0; j < rk.Length; j++)
                    for (int i = 0; i < rk.Length; i++)
                        hessian[j, i] += rk[i] * rk[j] * weights[k];
            }


            // The SVD is very stable, but is quite expensive, being on average
            // about 10-15 times more expensive than LU decomposition. There are
            // other ways to avoid a singular Hessian. For a very interesting
            // reading on the subject, please see:
            //
            //  - Jeff Gill & Gary King, "What to Do When Your Hessian Is Not Invertible",
            //    Sociological Methods & Research, Vol 33, No. 1, August 2004, 54-87.
            //    Available in: http://gking.harvard.edu/files/help.pdf
            //

            // Moreover, the computation of the inverse is optional, as it will
            // be used only to compute the standard errors of the regression.


            // Hessian Matrix is singular, try pseudo-inverse solution
            decomposition = new SingularValueDecomposition(hessian);
            deltas = decomposition.Solve(gradient);


            previous = (double[])coefficients.Clone();

            // Update coefficients using the calculated deltas
            for (int i = 0; i < coefficients.Length; i++)
                coefficients[i] -= deltas[i];

            // Return the relative maximum parameter change
            for (int i = 0; i < deltas.Length; i++)
                deltas[i] = Math.Abs(deltas[i]) / Math.Abs(previous[i]);

            return Matrix.Max(deltas);
        }


    }
}
