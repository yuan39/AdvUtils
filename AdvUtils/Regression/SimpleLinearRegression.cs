using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AdvUtils.Regression
{
    public class SimpleLinearRegression
    {
        private MultipleLinearRegression regression;

        /// <summary>
        ///   Creates a new Simple Linear Regression of the form y = Ax + B.
        /// </summary>
        /// 
        public SimpleLinearRegression()
        {
            this.regression = new MultipleLinearRegression(2);
        }

        /// <summary>
        ///   Angular coefficient (Slope).
        /// </summary>
        /// 
        public double Slope
        {
            get { return regression.Coefficients[1]; }
        }

        /// <summary>
        ///   Linear coefficient (Intercept).
        /// </summary>
        /// 
        public double Intercept
        {
            get { return regression.Coefficients[0]; }
        }


        /// <summary>
        ///   Performs the regression using the input and output
        ///   data, returning the sum of squared errors of the fit.
        /// </summary>
        /// 
        /// <param name="inputs">The input data.</param>
        /// <param name="outputs">The output data.</param>
        /// <returns>The regression Sum-of-Squares error.</returns>
        /// 
        public double Regress(double[] inputs, double[] outputs)
        {
            if (inputs.Length != outputs.Length)
                throw new ArgumentException("Number of input and output samples does not match", "outputs");

            double[][] X = new double[inputs.Length][];

            for (int i = 0; i < inputs.Length; i++)
            {
                // b[0]*1 + b[1]*inputs[i]
                X[i] = new double[] { 1.0, inputs[i] };
            }

            return regression.Regress(X, outputs);
        }

        /// <summary>
        ///   Computes the regression for a single input.
        /// </summary>
        /// 
        /// <param name="input">The input value.</param>
        /// <returns>The calculated output.</returns>
        /// 
        public double Compute(double input)
        {
            return Slope * input + Intercept;
        }
    }
}
