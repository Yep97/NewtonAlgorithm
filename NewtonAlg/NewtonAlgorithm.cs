using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using org.mariuszgromada.math.mxparser;
using AwokeKnowing.GnuplotCSharp;
using System.Numerics;
using Accord.Math;
using System.IO;
using System.Windows.Interop;
using System.Windows.Controls;
using Accord;
using System.Windows.Documents;
using OpenTK.Graphics.OpenGL;

namespace NewtonAlg
{
    class NewtonAlgorithm:NonlinearFunction
    {
        public List<double[]> x = new List<double[]>();
        List<double> fxy = new List<double>();
        List<double[]> gradient = new List<double[]>();
        List<double[,]> hessian = new List<double[,]>();
        List<double[]> ksi = new List<double[]>();
        List<double> criterionGradient = new List<double>();
        List<double> criterionX = new List<double>();
        List<double> criterionFx = new List<double>();

        public string criterionText;
        public string optimalX;
        public string optimalFunction;


        public void CalculateNewtonAlgorithm(double[] xStart)
        {
            double[] xTmp = new double[numberOfXs];
            double directionalDerivative;
            double t;
            double fTmp;
            double beta;
            double alfa;

            int i = 0;                                                                     // i = 0
            bool lesserLoop;
            bool greaterLoop = true;
            bool tmpLoop = true;
            double[,] hessianTmp = new double[numberOfXs, numberOfXs];

            // Punkt pierwszy
            // Wybierz punkt startowy
            x.Add(xStart);                                                                  // x_i = x0
            LoadArguments(x[i]);                                                            // Wczytaj x_i
            fxy.Add(CalculateFunction());                                                   // Oblicz funkcje w punkcie x_i
            gradient.Add(GradientInPoint(x[i]));                                            // Oblicz gradient w punkcie x_i
            criterionGradient.Add(1);
            criterionX.Add(1);
            criterionFx.Add(1);
            //hessian.Add(HMatrixInPoint(x[i]).Inverse());

            hessian.Add(new double[numberOfXs, numberOfXs]);                                // Dodaj hesjan h_i

            for (int l = 0; l < numberOfXs; l++)                                            // Hesjan h_0 = macierz jednostkowa
            {
                for (int p = 0; p < numberOfXs; p++)
                {
                    if (l == p)
                        hessian[0][l, p] = 1;
                    else
                        hessian[0][l, p] = 0;
                }
            }

            // Punkt drugi
            while (greaterLoop)
            {
                ksi.Add(new double[numberOfXs]);                                            // Dodaj nowy kierunek poszukiwań                

                for (int l = 0; l < numberOfXs; l++)                                        // Wyznacz kierunek poszukiwań ksi_i = -hessian(x_i)*gradient(x_i)
                {
                    for (int p = 0; p < numberOfXs; p++)
                    {
                        ksi[i][l] += -hessian[i][l, p] * gradient[i][p];
                    }
                }

                // Sprawdź kryteria zbieżności
                if (i != 0)
                {
                    criterionGradient.Add(ScalarProduct(gradient[i]));
                    criterionX.Add(Math.Sqrt(ScalarProduct(x[i].Subtract(x[i - 1]))));
                    criterionFx.Add(Math.Abs(fxy[i] - fxy[i - 1]));

                    if (criterionGradient[i] <= epsilon)                              // Kryterium <gradient(x_i),gradient(x_i)> <= epsilon
                    {
                        criterionText = "Iloczyn skalarny gradientu";
                        return;
                    }
                    else if (criterionX[i] <= epsilon)  // Kryterium ||x_i - x_i-1|| <= epsilon
                    {
                        criterionText = "Długość pomiędzy x_n i x_n-1";
                        return;
                    }
                    else if (criterionFx[i] <= epsilon)                      // Kryterium |f(x_i) - f(x_i-1)| <= epsilon
                    {
                        criterionText = "Wartość bezwzględna z różnicy funkcji f(x_n) i f(x_n-1)";
                        return;
                    }
                    else if (i == numberOfIterations)                                       // Kryterium liczba iteracji
                    {
                        criterionText = "Liczba iteracji";
                        return;
                    }
                }
                alfa = 1;

                while (tmpLoop)
                {
                    LoadArguments(x[i].Add(ksi[i].Multiply(alfa)));
                    if (fxy[i] < CalculateFunction())
                    {
                        alfa *= 0.8;
                    }
                    else
                        tmpLoop = false;
                }
                tmpLoop = true;

                beta = 0.0001;                                                                 // Współczynnik beta
                directionalDerivative = DirectionalDerivative(ksi[i], gradient[i]);         // Oblicz pochodną kierunkową p funkcji w punkcie x_i w kierunku ksi_i
                lesserLoop = true;
                t = 0.5;

                // Algorytm bisekcji z testem dwuskośnym Goldsteina
                while (lesserLoop)
                {
                    Array.Copy(x[i].Add(ksi[i].Multiply(alfa)), xTmp, numberOfXs);      // Oblicz punkt xTmp = x_i + ksi_i * tau
                    LoadArguments(xTmp);                                                // Wczytaj punkt xTmp = x_i + ksi_i * tau
                    fTmp = CalculateFunction();                                         // Oblicz wartość funkcji w tym punkcie xTmp
                    if (fTmp > (fxy[i] + (beta * directionalDerivative * alfa)))        // Czy fTmp < f(x_i) + beta * p * tau
                        alfa = t * alfa;                                                // tauL = tau
                    else                                                                // Wyznaczono minimum w kierunku
                    {
                        x.Add(new double[numberOfXs]);
                        Array.Copy(xTmp, x[i + 1], numberOfXs);
                        fxy.Add(fTmp);
                        i++;                                                            // i + 1
                        gradient.Add(GradientInPoint(x[i]));                            // Dodaj gradient w punkcie x_i do listy
                        hessianTmp = HMatrixInPoint(x[i]);
                        if (hessianTmp.Determinant() != 0)
                            hessian.Add(HMatrixInPoint(x[i]).Inverse());
                        else
                        {
                            hessian.Add(new double[numberOfXs, numberOfXs]);
                            Array.Copy(hessian[0], hessian[i], numberOfXs * numberOfXs);
                        }
                        lesserLoop = false;                                             // Zakończ pętle lesserLoop
                    }
                }
            }
            //return x;
        }

        public void ListIterationsValues(double[] xStart, Table myTable)
        {
            CalculateNewtonAlgorithm(xStart);
            int id = x.Count() - 1;
   
            int cols = 5 + numberOfXs;
            int rows = id + 1;

            for (int c = 0; c < cols; c++)
                myTable.Columns.Add(new TableColumn());

            for (int r = 0; r <= rows; r++)
            {
                TableRow tr = new TableRow();

                if (r == 0)
                {
                    tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("iter.")))));
                    for (int i = 0; i < numberOfXs; i++)
                        tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("x{0}", (i + 1).ToString())))));
                    tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("F(x)")))));
                    tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("kryt. grad.")))));
                    tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("kryt. x")))));
                    tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("kryt. fx")))));
                }
                else
                {
                    tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("0{0}", (r - 1).ToString())))));
                    for (int i = 0; i < numberOfXs; i++)
                        tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("{0}", Math.Round(x[r - 1][i], roundAccuracy, MidpointRounding.ToEven).ToString("0.000"))))));

                    tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("{0}", Math.Round(fxy[r - 1], roundAccuracy, MidpointRounding.ToEven).ToString("0.000"))))));
                    if (r == 1)
                    {
                        for (int i = 0; i < 3; i++)
                            tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("Brak")))));
                    }
                    else
                    {
                        tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("{0}", Math.Round(criterionGradient[r - 1], roundAccuracy, MidpointRounding.ToEven).ToString("0.000"))))));
                        tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("{0}", Math.Round(criterionX[r - 1], roundAccuracy, MidpointRounding.ToEven).ToString("0.000"))))));
                        tr.Cells.Add(new TableCell(new Paragraph(new Run(string.Format("{0}", Math.Round(criterionFx[r - 1], roundAccuracy, MidpointRounding.ToEven).ToString("0.000"))))));
                    }
                }
                TableRowGroup trg = new TableRowGroup();
                trg.Rows.Add(tr);
                myTable.RowGroups.Add(trg);
            }

            for (int i = 0; i < numberOfXs; i++)
                optimalX += string.Format("{0} ", Math.Round(x[id][i], roundAccuracy, MidpointRounding.ToEven).ToString("0.000"));
            optimalFunction = string.Format("{0}", Math.Round(fxy[id], roundAccuracy, MidpointRounding.ToEven).ToString("0.000"));
        }

        public void MakeContour(double[] xBound, double[] yBound)
        {
            int id = x.Count() - 1;
            double[] x1plot = new double[id + 1];
            double[] x2plot = new double[id + 1];
            for (int i = 0; i <= id; i++)
            {
                x1plot[i] = x[i][0];
                x2plot[i] = x[i][1];
            }
            string xr = "xr[" + xBound[0] + ":" + xBound[1] + "]";
            string yr = "yr[" + yBound[0] + ":" + yBound[1] + "]";
            GnuPlot.Unset("key");                                                      //hide the key or legend
            GnuPlot.HoldOn();

            GnuPlot.Plot(x1plot, x2plot, "with linespoints linestyle 1");
            GnuPlot.Set("cntrparam levels 30", "isosamples 40", xr, yr, "decimalsign locale"); //notice cntrparam levels (# height levels)
            GnuPlot.Contour(filename, "lc rgb 'blue'");
        }

        public void ResetAlgotithm()
        {
            x.Clear();
            fxy.Clear();
            gradient.Clear();
            hessian.Clear();
            ksi.Clear();
            criterionGradient.Clear();
            criterionX.Clear();
            criterionFx.Clear();
            optimalX = "";
            optimalFunction = "";
            criterionText = "";
        }

        public void ResetTable (Table myTable)
        {
            myTable.RowGroups.Clear();
        }
    }

}
