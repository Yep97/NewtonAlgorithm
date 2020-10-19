using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using org.mariuszgromada.math.mxparser;
using AwokeKnowing.GnuplotCSharp;
using System.Numerics;
using System.IO;
using MathNet.Numerics;
using Accord.Math;
using System.Collections;
using OpenTK.Graphics.OpenGL;
using System.Xml;
using HelixToolkit.Wpf;
using System.Windows.Media;
using System.Windows.Controls;

namespace NewtonAlg
{
    class NonlinearFunction
    {
        public string filename = @"E:\Users\User\source\repos\NewtonAlg\Text.txt";
        Function Fx;
        Argument x1 = new Argument("x1");
        Argument x2 = new Argument("x2");
        Argument x3 = new Argument("x3");
        Argument x4 = new Argument("x4");
        Argument x5 = new Argument("x5");

        Expression e1;
        protected int roundAccuracy = 3;
        protected double E = 0.0001;
        public int numberOfXs = 2;
        protected double epsilon = 0.001;
        public int numberOfIterations = 50;

        public string argument1;
        public string argument2;
        public string argument3;
        public string argument4;
        public string argument5;

        public bool wasFunctionChanged = true;

        public void LoadArguments(double[] functionArguments)
        {
            if (numberOfXs >= 1)
                x1.setArgumentValue(functionArguments[0]);
            if (numberOfXs >= 2)
                x2.setArgumentValue(functionArguments[1]);
            if (numberOfXs >= 3)
                x3.setArgumentValue(functionArguments[2]);
            if (numberOfXs >= 4)
                x4.setArgumentValue(functionArguments[3]);
            if (numberOfXs >= 5)
                x5.setArgumentValue(functionArguments[4]);
        }

        public void ParseFunction(string functionString)
        {
            Fx = new Function(functionString);
        }

        public string MakeFunction()
        {
            string func = "Fx(";

            if (numberOfXs >= 1)
                func += "x1";
            if (numberOfXs >= 2)
                func += ", x2";
            if (numberOfXs >= 3)
                func += ", x3";
            if (numberOfXs >= 4)
                func += ", x4";
            if (numberOfXs >= 5)
                func += ", x5";

            func += ")";
            return func;

        }

        public double CalculateFunction()
        {
            if (numberOfXs == 1)
                e1 = new Expression(MakeFunction(), Fx, x1);
            else if (numberOfXs == 2)
                e1 = new Expression(MakeFunction(), Fx, x1, x2);
            else if (numberOfXs == 3)
                e1 = new Expression(MakeFunction(), Fx, x1, x2, x3);
            else if (numberOfXs == 4)
                e1 = new Expression(MakeFunction(), Fx, x1, x2, x3, x4);
            else if (numberOfXs == 5)
                e1 = new Expression(MakeFunction(), Fx, x1, x2, x3, x4, x5);
            return (e1.calculate());
        }

        public void PrepareSaveFile(double[] xBounds, double[] yBounds, double step)
        {
            int sizeX = (int)((Math.Abs(xBounds[0]) + Math.Abs(xBounds[1])) / step) + 1;
            int sizeY = (int)((Math.Abs(yBounds[0]) + Math.Abs(yBounds[1])) / step) + 1;
            double[] X = new double[sizeX];
            double[] Y = new double[sizeY];
            double[] tmpValue = new double[numberOfXs];
            for (int i = 0; i < numberOfXs; i++)
                tmpValue[i] = 0;
            X[0] = xBounds[0];
            Y[0] = yBounds[0];
            for (int i = 1; i < sizeX; i++)
                X[i] = Math.Round(X[i - 1] + step, roundAccuracy, MidpointRounding.ToEven);
            for (int i = 1; i < sizeY; i++)
                Y[i] = Math.Round(Y[i - 1] + step, roundAccuracy, MidpointRounding.ToEven);
            double[,] Fx = new double[sizeX, sizeY];
            for (int i = 0; i < sizeX; i++)
            {
                for (int j = 0; j < sizeY; j++)
                {
                    tmpValue[0] = X[i];
                    tmpValue[1] = Y[j];
                    LoadArguments(tmpValue);
                    Fx[i, j] = CalculateFunction();
                }
            }
            SaveToFile(X, Y, Fx);
        }

        public void MakePlot(double[] xBound, double[] yBound)
        {
            string xr = "xr[" + xBound[0] + ":" + xBound[1] + "]";
            string yr = "yr[" + yBound[0] + ":" + yBound[1] + "]";
            GnuPlot.Set("cntrparam levels 30", "isosamples 40", xr, yr, "decimalsign locale");
            GnuPlot.SPlot(filename, "with pm3d");
        }

        public void SaveToFile(double[] X, double[] Y, double[,] Z)
        {
            string tempfile = Path.GetTempFileName();
            string[] lines = new string[X.Length * Y.Length];
            int k = 0;
            for (int i = 0; i < X.Length; i++)
            {
                for (int j = 0; j < Y.Length; j++)
                {
                    lines[k] = X[i] + " " + Y[j] + " " + Z[i,j];
                    k++;
                }
            }

            File.WriteAllLines(filename, lines);

            k = 0;
            using (var writer = new StreamWriter(tempfile))
            using (var reader = new StreamReader(filename))
            {
                writer.WriteLine("# X Y Z");
                while (!reader.EndOfStream)
                {
                    k++;
                    writer.WriteLine(reader.ReadLine());
                    if (k == X.Length)
                    {
                        writer.WriteLine("");
                        k = 0;
                    }
                }
            }
            File.Copy(tempfile, filename, true);
        }

        public double[] GradientInPoint(double[] X)
        {
            double tmpFunction;
            double[] Xtmp = new double[numberOfXs];
            Array.Copy(X, Xtmp, numberOfXs);
            LoadArguments(X);
            double Fxy = CalculateFunction();
            double[] gradient = new double[numberOfXs];
            for (int i = 0; i < numberOfXs; i++)
            {
                Xtmp[i] += E;
                LoadArguments(Xtmp);
                tmpFunction = CalculateFunction();
                gradient[i] = Math.Round((tmpFunction - Fxy) / E, roundAccuracy, MidpointRounding.ToEven);
                Array.Copy(X, Xtmp, numberOfXs);
            }
            return gradient;
        }

        public double[,] HMatrixInPoint(double[] X)
        {
            double[,] HMatrixValues = new double[numberOfXs, numberOfXs];
            double[] tmpFunction = new double[4];
            double[] Xtmp = new double[numberOfXs];
            double[,] DerivativeMatrix = { { 1, 1 }, { 1, -1 }, { -1, 1 }, { -1, -1 } };
            Array.Copy(X, Xtmp, numberOfXs);
            LoadArguments(X);
            double Fxy = CalculateFunction();
            for (int i = 0; i < numberOfXs; i++)
            {
                for (int j = 0; j < numberOfXs; j++)
                {
                    if (i == j)
                    {
                        Xtmp[i] += E;
                        LoadArguments(Xtmp);
                        tmpFunction[0] = CalculateFunction();
                        Xtmp[i] -= 2 * E;
                        LoadArguments(Xtmp);
                        tmpFunction[1] = CalculateFunction();
                        HMatrixValues[i, j] = Math.Round((tmpFunction[0] - 2 * Fxy + tmpFunction[1]) / (E * E), roundAccuracy, MidpointRounding.ToEven);
                        Array.Copy(X, Xtmp, numberOfXs);
                    }
                    else if (i < j)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            Xtmp[i] += E * DerivativeMatrix[k, 0];
                            Xtmp[j] += E * DerivativeMatrix[k, 1];

                            LoadArguments(Xtmp);
                            tmpFunction[k] = CalculateFunction();
                            Array.Copy(X, Xtmp, numberOfXs);
                        }
                        HMatrixValues[j, i] = HMatrixValues[i, j] = Math.Round((tmpFunction[0] - tmpFunction[1] - tmpFunction[2] + tmpFunction[3]) / (4 * E * E), roundAccuracy, MidpointRounding.ToEven);
                    }
                }
            }
            return HMatrixValues;
        }

        public double DirectionalDerivative(double[] ksi, double[] gradient)
        {
            double lengthOfKsi = 0;
            double[] directionalVector = new double[numberOfXs];
            double direcrionalDerivative = 0;
            lengthOfKsi = Math.Sqrt(ScalarProduct(ksi));
            for (int i = 0; i < numberOfXs; i++)
            {
                directionalVector[i] = (1 / lengthOfKsi) * ksi[i];
            }
            for (int i = 0; i < numberOfXs; i++)
            {
                direcrionalDerivative += gradient[i] * directionalVector[i];
            }
            return direcrionalDerivative;
        }
        public double ScalarProduct(double[] ksi)
        {
            double vectorLength = 0;
            for (int i = 0; i < numberOfXs; i++)
            {
                vectorLength += ksi[i] * ksi[i];
            }
            return vectorLength;
        }

        public double[] FunctionArguments()
        {
            double[] functionArguments = new double[numberOfXs];
            if (numberOfXs >= 1)
            {
                functionArguments[0] = 1;
                argument1 = string.Format("{0}", functionArguments[0]);
            }
            if (numberOfXs >= 2)
            {
                functionArguments[1] = 1;
                argument2 = string.Format("{0}", functionArguments[1]);
            }
            if (numberOfXs >= 3)
            {
                functionArguments[2] = 1;
                argument3 = string.Format("{0}", functionArguments[2]);
            }
            if (numberOfXs >= 4)
            {
                functionArguments[3] = 1;
                argument4 = string.Format("{0}", functionArguments[3]);
            }
            if (numberOfXs >= 5)
            {
                functionArguments[4] = 1;
                argument5 = string.Format("{0}", functionArguments[4]);
            }
            LoadArguments(functionArguments);
            return functionArguments;

        }

    }

}