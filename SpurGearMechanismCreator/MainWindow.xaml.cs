using Ookii.Dialogs.Wpf;
using SpurGearMechanismCreator.Calculations;
using System;
using System.Globalization;
using System.IO;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Animation;

namespace SpurGearMechanismCreator
{
	public partial class MainWindow : Window
    {
        public Point OriginPoint { get; set; }
        public Point GearPosition { get; set; }
        public Point PinionPosition { get; set; }

        public double GearRatio { get; set; }

        public PointCollection PinionPointData { get; set; }
        public PointCollection GearPointData { get; set; }

        public MainWindow()
        {
            InitializeComponent();

			var CustomCulture = (CultureInfo) Thread.CurrentThread.CurrentCulture.Clone();
            CustomCulture.NumberFormat.NumberDecimalSeparator = ".";
            Thread.CurrentThread.CurrentCulture = CustomCulture;

            OriginPoint = new Point(0, 0);
            DataContext = this;
            CollapsibleColumn.ElementStyle = this.FindResource("CollapsibleColumnStyle") as Style;
        }

        private void SetScale(double scaleX, double scaleY, double centerX, double centerY)
        {
            scaleY = -scaleY; //Horizontal flip
            RenderTranformationMatrix.Matrix = new Matrix(scaleX, 0, 0, scaleY, 0, 0);

            foreach (UIElement Element in GearCanvas.Children)
            {
                Canvas.SetLeft(Element, (GearCanvas.ActualWidth - centerX * scaleX * 2) / 2 / scaleX);
                Canvas.SetTop(Element, (GearCanvas.ActualHeight - centerY * scaleY * 2) / 2 / scaleY);
            }
        }

        private void OnClick(object sender, RoutedEventArgs e)
        {
            var Data = DimensionCalculations.Calculate(
                double.Parse(ModuleTextBox.Text),
                int.Parse(Z1TextBox.Text),
                int.Parse(Z2TextBox.Text),
                double.Parse(X1TextBox.Text),
                double.Parse(X2TextBox.Text));

            DataTable.ItemsSource = TableDataVisualisation.GetTableData(Data);
            OriginPoint = new Point(Data.ActionPosition.X, Data.ActionPosition.Y);
            GearPosition = Data.GearPosition;
            PinionPosition = Data.PinionPosition;

            GearRatio = Data.MechanismData.TransmissionRatio;

            PinionPointData = Data.PinionPoints;
            GearPointData = Data.GearPoints;

            GearCanvas.Children.Clear();
            foreach (var Element in Data.MechanismGeometry)
            {
                Canvas.SetLeft(Element, (GearCanvas.ActualWidth - Data.ActionPosition.X * 2) / 2);
                Canvas.SetTop(Element, (GearCanvas.ActualHeight - Data.ActionPosition.Y * 2) / 2);
                GearCanvas.Children.Add(Element);
            }

            SetScale(3, 3, OriginPoint.X, OriginPoint.Y);
        }

        private void OnSliderValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (RenderTranformationMatrix != null)
            {
                SetScale(e.NewValue, e.NewValue, OriginPoint.X, OriginPoint.Y);
            }
        }

        private void OnStartButtonClick(object sender, RoutedEventArgs e)
        {
            var Time = 60;

            var InitializeTransformPinion = new RotateTransform { CenterX = PinionPosition.X, CenterY = PinionPosition.Y };
            GearCanvas.Children[^2].RenderTransform = InitializeTransformPinion;

            var InitializeTransformGear = new RotateTransform { CenterX = GearPosition.X, CenterY = GearPosition.Y };
            GearCanvas.Children[^1].RenderTransform = InitializeTransformGear;

            var PinionAnimation = new DoubleAnimation(0, 360, new Duration(TimeSpan.FromSeconds(Time)));
            GearCanvas.Children[^2]
                .RenderTransform.BeginAnimation(RotateTransform.AngleProperty, PinionAnimation);

            var GearAnimation = new DoubleAnimation(360 / GearRatio, 0, new Duration(TimeSpan.FromSeconds(Time)));
            GearCanvas.Children[^1]
                .RenderTransform.BeginAnimation(RotateTransform.AngleProperty, GearAnimation);
        }

        private void OnStopButtonClick(object sender, RoutedEventArgs e)
        {
            GearCanvas.Children[^2]
                .RenderTransform.BeginAnimation(RotateTransform.AngleProperty, null);

            GearCanvas.Children[^1]
                .RenderTransform.BeginAnimation(RotateTransform.AngleProperty, null);
        }

        private void OnExportPinionClick(object sender, RoutedEventArgs e)
        {
            if (PinionPointData == null)
                return;

			var dialog = new VistaSaveFileDialog
			{
				Filter = "Text files (*.txt)|*.txt|All files (*.*)|*.*",
				DefaultExt = "txt",
				OverwritePrompt = true
			};

			if ((bool) dialog.ShowDialog(this))
			{
                File.WriteAllText(dialog.FileName, 
                    ExportData.GenerateTxtData(PinionPointData));
            }
        }

        private void OnExportGearClick(object sender, RoutedEventArgs e)
        {
            if (GearPointData == null)
                return;

			var dialog = new VistaSaveFileDialog
			{
				Filter = "Text files (*.txt)|*.txt|All files (*.*)|*.*",
				DefaultExt = "txt",
				OverwritePrompt = true
			};

			if ((bool) dialog.ShowDialog(this))
            {
                File.WriteAllText(dialog.FileName,
                    ExportData.GenerateTxtData(GearPointData));
            }
        }
    }
}
