using SpurGearMechanismCreator.Calculations;
using System;
using System.Globalization;
using System.IO;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Forms;
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

        public PointCollection? PinionPointData { get; set; }
        public PointCollection? GearPointData { get; set; }

        public MainWindow()
        {
            InitializeComponent();

            CultureInfo CustomCulture = (CultureInfo)Thread.CurrentThread.CurrentCulture.Clone();
            CustomCulture.NumberFormat.NumberDecimalSeparator = ".";
            Thread.CurrentThread.CurrentCulture = CustomCulture;

            OriginPoint = new Point(0, 0);
            DataContext = this;
            CollapsibleColumn.ElementStyle = FindResource("CollapsibleColumnStyle") as Style;
        }

        private void SetScale(double scaleX, double scaleY, double centerX, double centerY)
        {
            scaleY = -scaleY; //Horizontal flip
            RenderTranformationMatrix.Matrix = new Matrix(scaleX, 0, 0, scaleY, 0, 0);

            foreach (UIElement Element in GearCanvas.Children)
            {
                Canvas.SetLeft(Element, (GearCanvas.ActualWidth - (centerX * scaleX * 2)) / 2 / scaleX);
                Canvas.SetTop(Element, (GearCanvas.ActualHeight - (centerY * scaleY * 2)) / 2 / scaleY);
            }
        }

        private void OnClick(object sender, RoutedEventArgs e)
        {
            CalculationsResultsData Data = DimensionCalculations.Calculate(
                double.Parse(ModuleTextBox.Text, CultureInfo.CurrentCulture),
                int.Parse(Z1TextBox.Text, CultureInfo.CurrentCulture),
                int.Parse(Z2TextBox.Text, CultureInfo.CurrentCulture),
                double.Parse(X1TextBox.Text, CultureInfo.CurrentCulture),
                double.Parse(X2TextBox.Text, CultureInfo.CurrentCulture));

            if (Data is null || Data.MechanismData is null || Data.GearData is null || Data.PinionData is null || Data.MechanismGeometry is null)
            {
                return;
            }

            DataTable.ItemsSource = TableDataVisualisation.GetTableData(Data);
            OriginPoint = new Point(Data.ActionPosition.X, Data.ActionPosition.Y);
            GearPosition = Data.GearPosition;
            PinionPosition = Data.PinionPosition;

            GearRatio = Data.MechanismData.TransmissionRatio;

            PinionPointData = Data.PinionPoints;
            GearPointData = Data.GearPoints;

            GearCanvas.Children.Clear();
            foreach (UIElement Element in Data.MechanismGeometry)
            {
                Canvas.SetLeft(Element, (GearCanvas.ActualWidth - (Data.ActionPosition.X * 2)) / 2);
                Canvas.SetTop(Element, (GearCanvas.ActualHeight - (Data.ActionPosition.Y * 2)) / 2);
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
            int Time = 60;

            RotateTransform InitializeTransformPinion = new() { CenterX = PinionPosition.X, CenterY = PinionPosition.Y };
            GearCanvas.Children[^2].RenderTransform = InitializeTransformPinion;

            RotateTransform InitializeTransformGear = new() { CenterX = GearPosition.X, CenterY = GearPosition.Y };
            GearCanvas.Children[^1].RenderTransform = InitializeTransformGear;

            DoubleAnimation PinionAnimation = new(0, 360, new Duration(TimeSpan.FromSeconds(Time)));
            GearCanvas.Children[^2]
                .RenderTransform.BeginAnimation(RotateTransform.AngleProperty, PinionAnimation);

            DoubleAnimation GearAnimation = new(360 / GearRatio, 0, new Duration(TimeSpan.FromSeconds(Time)));
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
            if (PinionPointData is null || this is null)
            {
                return;
            }

            SaveFileDialog dialog = new()
            {
                Filter = "Text files (*.txt)|*.txt|All files (*.*)|*.*",
                DefaultExt = "txt",
                OverwritePrompt = true
            };

            if (dialog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                File.WriteAllText(dialog.FileName,
                    ExportData.GenerateTxtData(PinionPointData));
            }
        }

        private void OnExportGearClick(object sender, RoutedEventArgs e)
        {
            if (GearPointData == null)
            {
                return;
            }

            SaveFileDialog dialog = new()
            {
                Filter = "Text files (*.txt)|*.txt|All files (*.*)|*.*",
                DefaultExt = "txt",
                OverwritePrompt = true
            };

            if (dialog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                File.WriteAllText(dialog.FileName,
                    ExportData.GenerateTxtData(GearPointData));
            }
        }
    }
}
