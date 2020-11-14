using SpurGearMechanismCreator.Calculations;
using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Documents;
using System.Windows.Media;
using System.Windows.Media.Animation;
using System.Windows.Shapes;

namespace SpurGearMechanismCreator
{
    public partial class MainWindow : Window
    {
        public Point OriginPoint { get; set; }

        public MainWindow()
        {
            InitializeComponent();

            OriginPoint = new Point(0, 0);
            DataContext = this;
            CollapsibleColumn.ElementStyle = this.FindResource("CollapsibleColumnStyle") as Style;
        }

        private void OnClick(object sender, RoutedEventArgs e)
        {
            CalculationsResultsData Data = DimensionCalculations.Calculate(
                double.Parse(ModuleTextBox.Text),
                int.Parse(Z1TextBox.Text),
                int.Parse(Z2TextBox.Text),
                double.Parse(X1TextBox.Text),
                double.Parse(X2TextBox.Text));

            DataTable.ItemsSource = TableCalculations.GetTableData(Data);
            TranslateTransformObject.X = Data.ActionPosition.X;
            TranslateTransformObject.Y = Data.ActionPosition.Y;

            GearCanvas.Children.Clear();
            foreach (UIElement Element in Data.GearGeometry)
			{
                GearCanvas.Children.Add(Element);
            }
        }

		private void OnSliderValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
		{
            if (ScaleTransformation != null)
			{
                ScaleTransformation.ScaleX =  e.NewValue;
                ScaleTransformation.ScaleY = -e.NewValue;
            }
		}

		private void OnStartButtonClick(object sender, RoutedEventArgs e)
		{
            /*
            GearCanvas.BeginStoryboard(aa);

            Storyboard storyboard = new Storyboard();
            storyboard.RepeatBehavior = RepeatBehavior.Forever;
            storyboard.Duration = new Duration(TimeSpan.FromSeconds(10.0));

            DoubleAnimation rotateAnimation = new DoubleAnimation()
            {
                From = 0,
                To = 360,
                Duration = storyboard.Duration
            };
            
            Storyboard.SetTarget(rotateAnimation, GearCanvas.Children[]);
            Storyboard.SetTargetProperty(rotateAnimation, new PropertyPath("(UIElement.RenderTransform).(RotateTransform.Angle)"));

            storyboard.Children.Add(rotateAnimation);

            Resources.Add("Storyboard", storyboard);
            */
        }

		private void OnStopButtonClick(object sender, RoutedEventArgs e)
		{

		}
	}
}
