using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using AdvUtils.Regression;

namespace TestLogisticRegression
{
    class Program
    {
        public static bool LoadCorpusFromFile(string strFileName, out double[][] inputs, out double[] outputs)
        {
            StreamReader sr = new StreamReader(strFileName);
            List<double[]> inputList = new List<double[]>();
            while (sr.EndOfStream == false)
            {
                string strLine = sr.ReadLine();
                string[] items = strLine.Split('\t');
                double[] featureArray = new double[items.Length];
                for (int i = 0; i < items.Length; i++)
                {
                    featureArray[i] = double.Parse(items[i]);
                }
                inputList.Add(featureArray);
            }
            sr.Close();

            inputs = new double[inputList.Count][];
            outputs = new double[inputList.Count];

            for (int i = 0; i < inputList.Count; i++)
            {
                inputs[i] = new double[inputList[i].Length - 1];
                outputs[i] = inputList[i][inputList[i].Length - 1];
                for (int j = 0; j < inputList[i].Length - 1; j++)
                {
                    inputs[i][j] = inputList[i][j];
                }
            }

            return false;
        }

        static void Normalize(double[] outputs)
        {
            double maxResult = 0.0;
            for (int i = 0; i < outputs.Length; i++)
            {
                if (maxResult < outputs[i])
                {
                    maxResult = outputs[i];
                }
            }
            if (maxResult > 1.0)
            {
                Console.WriteLine("Normalizing output value...");
                for (int i = 0; i < outputs.Length; i++)
                {
                    outputs[i] = outputs[i] / maxResult;
                }
            }
        }

        static void Main(string[] args)
        {
            if (args.Length != 2)
            {
                Console.WriteLine("TestLogisticRegression [train corpus file name] [test corpus file name]");
                return;
            }

            double[][] inputs;
            double[] outputs;

            Console.WriteLine("Loading train corpus...");
            //Load training corpus
            LoadCorpusFromFile(args[0], out inputs, out outputs);

            //Try to normalize output value
            Normalize(outputs);

            Console.WriteLine("Logistic regression...");
            LogisticRegression lr = new LogisticRegression(inputs[0].Length);
            double error = lr.Regress(inputs, outputs);

            Console.WriteLine("Parameter list:");
            for (int i = 0; i < lr.Coefficients.Length; i++)
            {
                Console.WriteLine("Coefficient {0}: {1}", i, lr.Coefficients[i]);
            }
            Console.WriteLine("Delta: {0}", error);
            Console.WriteLine();

            //Load test corpus
            Console.WriteLine("Testing regress result:");
            LoadCorpusFromFile(args[1], out inputs, out outputs);
            //Try to normalize output value
            Normalize(outputs);

            for (int i = 0; i < outputs.Length; i++)
            {
                StringBuilder sb = new StringBuilder();
                double output = lr.Compute(inputs[i]);
                for (int j = 0; j < inputs[i].Length; j++)
                {
                    sb.Append(inputs[i][j].ToString() + " ");
                }
                sb.Append(outputs[i] + " RV:" + output);
                sb.Append(" Err:" + (output - outputs[i]).ToString());

                Console.WriteLine(sb.ToString());
            }
        }
    }
}
