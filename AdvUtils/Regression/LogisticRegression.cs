using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AdvUtils.Regression
{
    public class LogisticRegression
    {
        private double[] coefficients;
        private double[] standardErrors;


        //---------------------------------------------


        #region Constructor
        /// <summary>
        ///   Creates a new Logistic Regression Model.
        /// </summary>
        /// 
        /// <param name="inputs">The number of input variables for the model.</param>
        /// 
        public LogisticRegression(int inputs)
        {
            this.coefficients = new double[inputs + 1];
            this.standardErrors = new double[inputs + 1];
        }

        /// <summary>
        ///   Creates a new Logistic Regression Model.
        /// </summary>
        /// 
        /// <param name="inputs">The number of input variables for the model.</param>
        /// <param name="intercept">The starting intercept value. Default is 0.</param>
        /// 
        public LogisticRegression(int inputs, double intercept)
            : this(inputs)
        {
            this.coefficients[0] = intercept;
        }
        #endregion


        //---------------------------------------------


        #region Properties
        /// <summary>
        ///   Gets the coefficient vector, in which the
        ///   first value is always the intercept value.
        /// </summary>
        /// 
        public double[] Coefficients
        {
            get { return coefficients; }
        }

        /// <summary>
        ///   Gets the standard errors associated with each
        ///   cofficient during the model estimation phase.
        /// </summary>
        /// 
        public double[] StandardErrors
        {
            get { return standardErrors; }
        }

        /// <summary>
        ///   Gets the number of inputs handled by this model.
        /// </summary>
        /// 
        public int Inputs
        {
            get { return coefficients.Length - 1; }
        }
        #endregion


        //---------------------------------------------


        #region Public Methods
        /// <summary>
        ///   Computes the model output for the given input vector.
        /// </summary>
        /// 
        /// <param name="input">The input vector.</param>
        /// <returns>The output value.</returns>
        /// 
        public double Compute(double[] input)
        {
            double logit = coefficients[0];

            for (int i = 1; i < coefficients.Length; i++)
                logit += input[i - 1] * coefficients[i];

            return Logistic(logit);
        }

        /// <summary>
        ///   Computes the model output for each of the given input vectors.
        /// </summary>
        /// 
        /// <param name="input">The array of input vectors.</param>
        /// <returns>The array of output values.</returns>
        /// 
        public double[] Compute(double[][] input)
        {
            double[] output = new double[input.Length];

            for (int i = 0; i < input.Length; i++)
                output[i] = Compute(input[i]);

            return output;
        }


        /// <summary>
        ///   Gets the Odds Ratio for a given coefficient.
        /// </summary>
        /// <remarks>
        ///   The odds ratio can be computed raising euler's number
        ///   (e ~~ 2.71) to the power of the associated coefficient.
        /// </remarks>
        /// <param name="index">
        ///   The coefficient's index. The first value
        ///   (at zero index) is the intercept value.
        /// </param>
        /// <returns>
        ///   The Odds Ratio for the given coefficient.
        /// </returns>
        /// 
        public double GetOddsRatio(int index)
        {
            return Math.Exp(coefficients[index]);
        }

        /// <summary>
        ///   Gets the 95% confidence interval for the
        ///   Odds Ratio for a given coefficient.
        /// </summary>
        /// 
        /// <param name="index">
        ///   The coefficient's index. The first value
        ///   (at zero index) is the intercept value.
        /// </param>
        /// 
        //public DoubleRange GetConfidenceInterval(int index)
        //{
        //    double coeff = coefficients[index];
        //    double error = standardErrors[index];

        //    double upper = coeff + 1.9599 * error;
        //    double lower = coeff - 1.9599 * error;

        //    DoubleRange ci = new DoubleRange(Math.Exp(lower), Math.Exp(upper));

        //    return ci;
        //}

        /// <summary>
        ///   Gets the Wald Test for a given coefficient.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Wald statistical test is a test for a model parameter in which
        ///   the estimated parameter θ is compared with another proposed parameter
        ///   under the assumption that the difference between them will be approximately
        ///   normal. There are several problems with the use of the Wald test. Please
        ///   take a look on substitute tests based on the log-likelihood if possible.
        /// </remarks>
        /// 
        /// <param name="index">
        ///   The coefficient's index. The first value
        ///   (at zero index) is the intercept value.
        /// </param>
        /// 
        //public WaldTest GetWaldTest(int index)
        //{
        //    return new WaldTest(coefficients[index], 0.0, standardErrors[index]);
        //}


        /// <summary>
        ///   Gets the Log-Likelihood for the model.
        /// </summary>
        /// 
        /// <param name="input">A set of input data.</param>
        /// <param name="output">A set of output data.</param>
        /// <returns>
        ///   The Log-Likelihood (a measure of performance) of
        ///   the model calculated over the given data sets.
        /// </returns>
        /// 
        public double GetLogLikelihood(double[][] input, double[] output)
        {
            double sum = 0;

            for (int i = 0; i < input.Length; i++)
            {
                double actualOutput = Compute(input[i]);
                double expectedOutput = output[i];

                if (actualOutput != 0)
                    sum += expectedOutput * Math.Log(actualOutput);

                if (actualOutput != 1)
                    sum += (1 - expectedOutput) * Math.Log(1 - actualOutput);

#if DEBUG
                if (Double.IsNaN(sum))
                    throw new Exception();
#endif
            }

            return sum;
        }

        /// <summary>
        ///   Gets the Deviance for the model.
        /// </summary>
        /// 
        /// <remarks>
        ///   The deviance is defined as -2*Log-Likelihood.
        /// </remarks>
        /// 
        /// <param name="input">A set of input data.</param>
        /// <param name="output">A set of output data.</param>
        /// <returns>
        ///   The deviance (a measure of performance) of the model
        ///   calculated over the given data sets.
        /// </returns>
        /// 
        public double GetDeviance(double[][] input, double[] output)
        {
            return -2.0 * GetLogLikelihood(input, output);
        }

        /// <summary>
        ///   Gets the Log-Likelihood Ratio between two models.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Log-Likelihood ratio is defined as 2*(LL - LL0).
        /// </remarks>
        /// 
        /// <param name="input">A set of input data.</param>
        /// <param name="output">A set of output data.</param>
        /// <param name="regression">Another Logistic Regression model.</param>
        /// <returns>The Log-Likelihood ratio (a measure of performance
        /// between two models) calculated over the given data sets.</returns>
        /// 
        public double GetLogLikelihoodRatio(double[][] input, double[] output, LogisticRegression regression)
        {
            return 2.0 * (this.GetLogLikelihood(input, output) - regression.GetLogLikelihood(input, output));
        }


        /// <summary>
        ///   The likelihood ratio test of the overall model, also called the model chi-square test.
        /// </summary>
        /// 
        /// <remarks>
        ///   <para>
        ///   The Chi-square test, also called the likelihood ratio test or the log-likelihood test
        ///   is based on the deviance of the model (-2*log-likelihood). The log-likelihood ratio test 
        ///   indicates whether there is evidence of the need to move from a simpler model to a more
        ///   complicated one (where the simpler model is nested within the complicated one).</para>
        ///   <para>
        ///   The difference between the log-likelihood ratios for the researcher's model and a
        ///   simpler model is often called the "model chi-square".</para>
        /// </remarks>
        /// 
        //public ChiSquareTest ChiSquare(double[][] input, double[] output)
        //{
        //    double y0 = output.Count(y => y == 0.0);
        //    double y1 = output.Length - y0;

        //    LogisticRegression regression = new LogisticRegression(Inputs, Math.Log(y1 / y0));

        //    double ratio = GetLogLikelihoodRatio(input, output, regression);
        //    return new ChiSquareTest(ratio, coefficients.Length - 1);
        //}



        /// <summary>
        ///   Creates a new LogisticRegression that is a copy of the current instance.
        /// </summary>
        /// 
        public object Clone()
        {
            var regression = new LogisticRegression(coefficients.Length);
            regression.coefficients = (double[])this.coefficients.Clone();
            regression.standardErrors = (double[])this.standardErrors.Clone();
            return regression;
        }
        #endregion


        //---------------------------------------------


        #region Static Methods
        /// <summary>
        ///   The Logistic function.
        /// </summary>
        /// 
        /// <param name="value">The logit parameter.</param>
        /// 
        public static double Logistic(double value)
        {
            return 1.0 / (1.0 + System.Math.Exp(-value));
        }
        #endregion


        public double Regress(double[][] inputs, double[] outputs, int maxIterations = 100000, double limit = 10e-5)
        {
            IterativeReweightedLeastSquares learning = new IterativeReweightedLeastSquares(this);

            int iteration = 0;
            double delta = 0.0;

            do // learning iterations until convergence
            {
                delta = learning.Run(inputs, outputs);
                iteration++;

                Console.WriteLine("Iteration {0}, Delta {1}", iteration, delta);
                for (int i = 0; i < Coefficients.Length; i++)
                {
                    Console.WriteLine("Coefficient {0}: {1}", i, Coefficients[i]);
                }
                Console.WriteLine();

            } while (iteration < maxIterations && delta >= limit);

            return delta;
        }
    }
}
