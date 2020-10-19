using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using org.mariuszgromada.math.mxparser;
using AwokeKnowing.GnuplotCSharp;

namespace NewtonAlg
{
    /// <summary>
      /// Logika interakcji dla klasy MainWindow.xaml
      /// </summary>
    public partial class MainWindow : Window
    {
        NewtonAlgorithm newtonAlgorithm = new NewtonAlgorithm();
        public MainWindow()
        {
            InitializeComponent();
            string equals = " = ";

            text1.Text = "Funkcja:";
            text2.Text = "Punkt startowy:";
            text4.Text = "Przedział wykresu:";

            //string function = "x1^4+x2^4-0.62*x1^2-0.62*x2^2";
            string function = "100*(x2-x1^2)^2+(1-x1)^2";
            //string function = "4*x1^2 - 2.1*x1^4 + 1/3 * x1^6 + x1*x2 - 4*x2^2 + 4*x2^4";
            functionText.Text = newtonAlgorithm.MakeFunction() + equals;
            textfun.Text = function;
            string realTextFun = newtonAlgorithm.MakeFunction() + equals + function;

            double[] functionArguments = newtonAlgorithm.FunctionArguments();

            newtonAlgorithm.LoadArguments(functionArguments);
            newtonAlgorithm.ParseFunction(realTextFun);

            double[] X = new double[2];
            X[0] = -2; X[1] = 2;
            double[] Y = new double[2];
            Y[0] = -2; Y[1] = 2;
            double step = 0.1;
            
            startAxisX.Text = X[0].ToString();
            endAxisX.Text = X[1].ToString();
            startAxisY.Text = Y[0].ToString();
            endAxisY.Text = Y[1].ToString();
            stepText.Text = step.ToString();

            newtonAlgorithm.ListIterationsValues(functionArguments,myTable);
            criterionText.Text = newtonAlgorithm.criterionText;
            optimalXText.Text = newtonAlgorithm.optimalX;
            optimalFunctionText.Text = newtonAlgorithm.optimalFunction;
            argumentText1.Text = newtonAlgorithm.argument1;
            argumentText2.Text = newtonAlgorithm.argument2;
            argumentText3.Text = newtonAlgorithm.argument3;
            argumentText4.Text = newtonAlgorithm.argument4;
            argumentText5.Text = newtonAlgorithm.argument5;
            numberOfXsText.Text = newtonAlgorithm.numberOfXs.ToString();
            numberOfIterationsText.Text = newtonAlgorithm.numberOfIterations.ToString();
        }
        
        private void ChangeStart_Click(object sender, RoutedEventArgs e)
        {
            string equals = " = ";
            newtonAlgorithm.numberOfXs = int.Parse(numberOfXsText.Text);
            newtonAlgorithm.numberOfIterations = int.Parse(numberOfIterationsText.Text);
            functionText.Text = newtonAlgorithm.MakeFunction() + equals;
            double[] functionArguments = newtonAlgorithm.FunctionArguments();
            string function = textfun.Text;
            string realTextFun = newtonAlgorithm.MakeFunction() + equals + function;

            if (newtonAlgorithm.numberOfXs >= 1)
                functionArguments[0] = double.Parse(argumentText1.Text);
            if (newtonAlgorithm.numberOfXs >= 2)
                functionArguments[1] = double.Parse(argumentText2.Text);
            if (newtonAlgorithm.numberOfXs >= 3)
                functionArguments[2] = double.Parse(argumentText3.Text);
            if (newtonAlgorithm.numberOfXs >= 4)
                functionArguments[3] = double.Parse(argumentText4.Text);
            if (newtonAlgorithm.numberOfXs >= 5)
                functionArguments[4] = double.Parse(argumentText5.Text);

            newtonAlgorithm.LoadArguments(functionArguments);
            newtonAlgorithm.ParseFunction(realTextFun);

            newtonAlgorithm.ResetAlgotithm();
            newtonAlgorithm.ResetTable(myTable);
            newtonAlgorithm.ListIterationsValues(functionArguments, myTable);
            criterionText.Text = newtonAlgorithm.criterionText;
            optimalXText.Text = newtonAlgorithm.optimalX;
            optimalFunctionText.Text = newtonAlgorithm.optimalFunction;
            newtonAlgorithm.wasFunctionChanged = true;
        }

        private void PlotButton_Click(object sender, RoutedEventArgs e)
        {
            string equals = " = ";
            string function = textfun.Text;
            string realTextFun = newtonAlgorithm.MakeFunction() + equals + function;
            double[] X = new double[2];
            double[] Y = new double[2];
            double step;
            X[0] = double.Parse(startAxisX.Text);
            X[1] = double.Parse(endAxisX.Text);
            Y[0] = double.Parse(startAxisY.Text);
            Y[1] = double.Parse(endAxisY.Text);
            step = double.Parse(stepText.Text);

            if (newtonAlgorithm.wasFunctionChanged == true)
            {
                newtonAlgorithm.ParseFunction(realTextFun);
                newtonAlgorithm.PrepareSaveFile(X, Y, step);
                newtonAlgorithm.wasFunctionChanged = false;
            }
            newtonAlgorithm.MakePlot(X, Y);
        }

        private void ContourButton_Click(object sender, RoutedEventArgs e)
        {
            string equals = " = ";
            string function = textfun.Text;
            string realTextFun = newtonAlgorithm.MakeFunction() + equals + function;
            newtonAlgorithm.ParseFunction(realTextFun);
            double[] X = new double[2];
            double[] Y = new double[2];
            double step;
            X[0] = double.Parse(startAxisX.Text);
            X[1] = double.Parse(endAxisX.Text);
            Y[0] = double.Parse(startAxisY.Text);
            Y[1] = double.Parse(endAxisY.Text);
            step = double.Parse(stepText.Text);

            if (newtonAlgorithm.wasFunctionChanged == true)
            {
                newtonAlgorithm.ParseFunction(realTextFun);
                newtonAlgorithm.PrepareSaveFile(X, Y, step);
                newtonAlgorithm.wasFunctionChanged = false;
            }
            newtonAlgorithm.MakeContour(X, Y);
        }
    }
}