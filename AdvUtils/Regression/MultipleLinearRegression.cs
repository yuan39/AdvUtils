using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AdvUtils.Regression
{
    public class MultipleLinearRegression
    {
        private double[] coefficients;
        private bool addIntercept;

        /// <summary>
        ///   Creates a new Multiple Linear Regression.
        /// </summary>
        /// 
        /// <param name="inputs">The number of inputs for the regression.</param>
        /// <param name="intercept">Whether add intercept into coefficients</param>
        public MultipleLinearRegression(int inputs, bool intercept = false)
        {
            if (intercept) inputs++;
            this.coefficients = new double[inputs];
            this.addIntercept = intercept;
        }


        /// <summary>
        ///   Gets the coefficients used by the regression model. If the model
        ///   contains an intercept term, it will be in the end of the vector.
        /// </summary>
        /// 
        public double[] Coefficients
        {
            get { return coefficients; }
        }


        /// <summary>
        ///   Performs the regression using the input vectors and output
        ///   data, returning the sum of squared errors of the fit.
        /// </summary>
        /// 
        /// <param name="inputs">The input vectors to be used in the regression.</param>
        /// <param name="outputs">The output values for each input vector.</param>
        /// <returns>The Sum-Of-Squares error of the regression.</returns>
        /// 
        public virtual double Regress(double[][] inputs, double[] outputs)
        {
            if (inputs.Length != outputs.Length)
                throw new ArgumentException("Number of input and output samples does not match", "outputs");

            double[,] X;

            return regress(inputs, outputs, out X);
        }

    
        private double regress(double[][] inputs, double[] outputs, out double[,] X)
        {
            if (inputs.Length != outputs.Length)
                throw new ArgumentException("Number of input and output samples does not match", "outputs");

            int rows = inputs.Length;     // inputs
            int cols = inputs[0].Length;        // points

            if (addIntercept)
            {
                X = new double[rows, cols + 1];
                for (int i = 0; i < rows; i++)
                    X[i, cols] = 1;
            }
            else
            {
                X = new double[rows, cols];
            }

            for (int i = 0; i < inputs.Length; i++)
                for (int j = 0; j < inputs[i].Length; j++)
                    X[i, j] = inputs[i][j];


            // Solve V*C = B to find C (the coefficients)
            coefficients = new SingularValueDecomposition(X).Solve(outputs);

            // Calculate Sum-Of-Squares error
            double error = 0.0;
            double e;
            for (int i = 0; i < outputs.Length; i++)
            {
                e = outputs[i] - Compute(inputs[i]);
                error += e * e;
            }

            return error;
        }
        /// <summary>
        ///   Computes the Multiple Linear Regression for an input vector.
        /// </summary>
        /// 
        /// <param name="input">The input vector.</param>
        /// 
        /// <returns>The calculated output.</returns>
        /// 
        public double Compute(double[] input)
        {
            double output = 0.0;

            for (int i = 0; i < input.Length; i++)
                output += coefficients[i] * input[i];

            if (addIntercept) output += coefficients[input.Length];

            return output;
        }
    }
}
