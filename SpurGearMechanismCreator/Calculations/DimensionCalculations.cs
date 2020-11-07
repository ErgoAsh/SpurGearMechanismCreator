using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Security.Cryptography.Xml;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace SpurGearMechanismCreator
{
    public class CalculationsResultsData
    {
        public GearCharacteristicsData PinionData { get; set; }
        public GearCharacteristicsData GearData { get; set; }
        public GearMechanismData MechanismData { get; set; }
        public List<UIElement> PinionGeometry { get; set; }
        public List<UIElement> GearGeometry { get; set; }
        public Point PinionPosition { get; set; }
        public Point GearPosition { get; set; }
        public Point ActionPosition { get; set; }
    }

    public class GearCharacteristicsData
	{
        public int NumberOfTeeths { get; set; }
        public double ShiftCoefficient { get; set; }
        public double ReferencePitchDiameter { get; set; }
        public double OperatingPitchDiameter { get; set; }
        public double Addendum { get; set; }
        public double Dedendum { get; set; }
        public double DedendumDiameter { get; set; }
        public double AddendumDiameter { get; set; }
        public double BaseCircleDiameter { get; set; }
        public double OperatingClearance { get; set; }
        public double ThicknessReference { get; set; } //not needed?
        public double ThicknessOperating { get; set; } 
        public double ThinknessBase { get; set; }
        public double ThicknessTip { get; set; }
    }

    public class GearMechanismData
	{
        public double Module { get; set; }
        public double PressureAngle { get; set; }
        public double OperatingPressureAngle { get; set; }
        public double Clearance { get; set; }
        public double CenterDistance { get; set; }
        public double CenterDistanceCoefficient { get; set; }
        public double TransmissionRatio { get; set; }
        public double ContactRatio { get; set; }
        public double Pitch { get; set; }
		public double FilletRadius { get; internal set; }
	}

    public static class DimensionCalculations
    {
        public static double Radians(double angle)
        {
            return (Math.PI / 180) * angle;
        }

        public static double Degrees(double radians)
        {
            return (180 / Math.PI) * radians;
        }

        public static double Involute(double AngleRad)
        {
            return Math.Tan(AngleRad) - AngleRad;
        }

        public static double InverseInvolute(double Involute) {
            double Angle1 = 0;
            double Angle2 = 0;

            while(true)
            {
                Angle2 = Angle1;
                Angle1 = Math.Atan(Angle1 + Involute);

                double Diff = Math.Abs(Angle1 - Angle2);
                if (Diff < Math.Pow(10, -10))
                {
                    break;
                }
            }

            return Angle1;
        }

        public static Point Cartesian(double Rho, double phi)
		{
            return new Point
            {
                X = Rho * Math.Cos(phi),
                Y = Rho * Math.Sin(phi)
            };
		}
        /*
        public static Point Polar(double X, double Y, Point Center)
		{

		}
        */
        public static double Distance(Point A, Point B)
		{
            return Math.Sqrt(Math.Pow(B.X - A.X, 2) + Math.Pow(B.Y - A.Y, 2));
		}

        public static Point RotatePointAAroundB(Point A, Point B, double Angle)
        {
            return new Point{
                X = Math.Cos(Angle) * (A.X - B.X) - Math.Sin(Angle) * (A.Y - B.Y) + B.X,
                Y = Math.Sin(Angle) * (A.X - B.X) + Math.Cos(Angle) * (A.Y - B.Y) + B.Y
            };
        }

        public static PointCollection InvolutePoints(double StartAngle, double dTheta, double RMin, double RMax, Point Center, bool ReverseDirection)
		{
            PointCollection Coll = new PointCollection();
            double Theta = 0;
            while (true)
			{
                int neg = ReverseDirection ? -1 : 1;
                Point BasePoint = new Point
                {
                    X =       RMin * (Math.Cos(Theta) + Theta * Math.Sin(Theta)),
                    Y = neg * RMin * (Math.Sin(Theta) - Theta * Math.Cos(Theta))
                };

                Point ActualPoint = RotatePointAAroundB(BasePoint, Center, StartAngle);

                Theta += dTheta;
                if (Distance(ActualPoint, Center) > RMax)
				{
                    //Coll.Add(Cartesian(RMax, StartAngle + Theta)); //Small error for small dTheta
                    return Coll;
                } 
                else
				{
                    Coll.Add(ActualPoint);
				}
			}
		}

        public static CalculationsResultsData Calculate(double m, int z1, int z2, double x1, double x2)
        {
            double dTheta = 0.001;

            double ha = 1; //Coefficient of a tool hight
            double c_star = 0.25f; //Coefficient of radial spacing

            double alpha = Radians(20);

            double i = (double) z2 / z1;

            double inv_alpha_prime = 2 * Math.Tan(alpha) * (x1 + x2) / (z1 + z2) + Involute(alpha);
            double alpha_prime = InverseInvolute(inv_alpha_prime);

            double y = (z1 + z2) / 2 * ((Math.Cos(alpha) / Math.Cos(alpha_prime)) - 1);
            double a = ((z1 + z2) / 2 + y) * m;

            //Pitch
            double p = Math.PI * m;

            // Pitch circle
            double d1 = z1 * m;
            double d2 = z2 * m;

            // Base circle
            double d_b1 = d1 * Math.Cos(alpha);
            double d_b2 = d2 * Math.Cos(alpha);

            // Working pitch diameter
            double d_prime1 = d_b1 / Math.Cos(alpha_prime);
            double d_prime2 = d_b2 / Math.Cos(alpha_prime);

            // Dedendum circle
            double d_f1 = m * (z1 - 2 * ha - 2 * c_star + 2 * x1);
            double d_f2 = m * (z2 - 2 * ha - 2 * c_star + 2 * x2);

            // Addendum
            double h_a1 = (1 + y - x1) * m;
            double h_a2 = (1 + y - x2) * m;

            // Addendum circle
            double d_a1 = d1 + 2 * h_a1;
            double d_a2 = d2 + 2 * h_a1;

            // Overlap coefficient
            double epsilon = (Math.Sqrt(Math.Pow(d_a1 / 2, 2) - Math.Pow(d_b1 / 2, 2))
                + Math.Sqrt(Math.Pow(d_a2 / 2, 2) - Math.Pow(d_b2 / 2, 2))
                - a * Math.Sin(alpha_prime))/(Math.PI * m * Math.Cos(alpha));

            // Thickness (arc length) of tooth at base circle
            double s_1_angle = Radians(2 * (90 / z1 + (360 * x1 * Math.Tan(alpha)) / (Math.PI * z1)));

            // Tooth depth
            double h = (2.25 + y - (x1 + x2)) * m;

            double rho = 0.38 * m;

            PointCollection Dedendum1Points = new PointCollection();
            PointCollection Base1Points = new PointCollection();
            PointCollection RefPitch1Points = new PointCollection();
            PointCollection WorkPitch1Points = new PointCollection();
            PointCollection Addendum1Points = new PointCollection();

            for (double Theta = 0; Theta < 2 * Math.PI; Theta += dTheta)
            {
                Point dedendumPoint = new Point(
                    (d_f1 / 2) * Math.Cos(Theta),
                    (d_f1 / 2) * Math.Sin(Theta));
                Dedendum1Points.Add(dedendumPoint);

                Point basePoint = new Point(
                    (d_b1 / 2) * Math.Cos(Theta),
                    (d_b1 / 2) * Math.Sin(Theta));
                Base1Points.Add(basePoint);

                Point refPitchPoint = new Point(
                    (d1 / 2) * Math.Cos(Theta),
                    (d1 / 2) * Math.Sin(Theta));
                RefPitch1Points.Add(refPitchPoint);

                Point workPitchPoint = new Point(
                    (d_prime1 / 2) * Math.Cos(Theta),
                    (d_prime1 / 2) * Math.Sin(Theta));
                WorkPitch1Points.Add(workPitchPoint);

                Point addendumPoint = new Point(
                    (d_a1 / 2) * Math.Cos(Theta),
                    (d_a1 / 2) * Math.Sin(Theta));
                Addendum1Points.Add(addendumPoint);
            }

            PointCollection Dedendum2Points = new PointCollection();
            PointCollection Base2Points = new PointCollection();
            PointCollection RefPitch2Points = new PointCollection();
            PointCollection WorkPitch2Points = new PointCollection();
            PointCollection Addendum2Points = new PointCollection();

            for (double Theta = 0; Theta < 2 * Math.PI; Theta += dTheta)
            {
                Point dedendumPoint = new Point(
                    a + (d_f2 / 2) * Math.Cos(Theta),
                        (d_f2 / 2) * Math.Sin(Theta));
                Dedendum2Points.Add(dedendumPoint);

                Point basePoint = new Point(
                    a + (d_b2 / 2) * Math.Cos(Theta),
                        (d_b2 / 2) * Math.Sin(Theta));
                Base2Points.Add(basePoint);

                Point refPitchPoint = new Point(
                    a + (d2 / 2) * Math.Cos(Theta),
                        (d2 / 2) * Math.Sin(Theta));
                RefPitch2Points.Add(refPitchPoint);

                Point workPitchPoint = new Point(
                    a + (d_prime2 / 2) * Math.Cos(Theta),
                        (d_prime2 / 2) * Math.Sin(Theta));
                WorkPitch2Points.Add(workPitchPoint);

                Point addendumPoint = new Point(
                    a + (d_a2 / 2) * Math.Cos(Theta),
                        (d_a2 / 2) * Math.Sin(Theta));
                Addendum2Points.Add(addendumPoint);
            }

            List<UIElement> GearElements = new List<UIElement>();

			Polygon Dedendum1Polygon = new Polygon
			{
				Stroke = Brushes.Black,
				StrokeThickness = 0.5,
				Points = Dedendum1Points
            };
			GearElements.Add(Dedendum1Polygon);

			Polygon Base1Polygon = new Polygon
			{
				Stroke = Brushes.Black,
				StrokeThickness = 0.5,
				Points = Base1Points
            };
			GearElements.Add(Base1Polygon);

			Polygon RefPitch1Polygon = new Polygon
			{
				Stroke = Brushes.Black,
				StrokeThickness = 1,
                StrokeDashArray = DoubleCollection.Parse("3,1"),
                Points = RefPitch1Points
            };
			GearElements.Add(RefPitch1Polygon);

            Polygon WorkPitch1Polygon = new Polygon
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Points = WorkPitch1Points
            };
            GearElements.Add(WorkPitch1Polygon);

            Polygon Addendum1Polygon = new Polygon
			{
				Stroke = Brushes.Black,
				StrokeThickness = 0.5,
				Points = Addendum1Points
			};
			GearElements.Add(Addendum1Polygon);

            Polygon Dedendum2Polygon = new Polygon
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Points = Dedendum2Points
            };
            GearElements.Add(Dedendum2Polygon);

            Polygon Base2Polygon = new Polygon
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Points = Base2Points
            };
            GearElements.Add(Base2Polygon);

            Polygon RefPitch2Polygon = new Polygon
            {
                Stroke = Brushes.Black,
                StrokeThickness = 1,
                StrokeDashArray = DoubleCollection.Parse("3,1"),
                Points = RefPitch2Points
            };
            GearElements.Add(RefPitch2Polygon);

            Polygon WorkPitch2Polygon = new Polygon
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Points = WorkPitch2Points
            };
            GearElements.Add(WorkPitch2Polygon);

            Polygon Addendum2Polygon = new Polygon
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Points = Addendum2Points
            };
            GearElements.Add(Addendum2Polygon);

            double[] ThetaRange = Generate.LinearRange(0, p / (d1 / 2), 2 * Math.PI);
            foreach (double th in ThetaRange)
			{
                PointCollection Inv = InvolutePoints(th, dTheta, d_b1 / 2, d_a1 / 2, new Point(0, 0), false);
                Polyline InvoluteLine = new Polyline
                {
                    Stroke = Brushes.Red,
                    StrokeThickness = 1,
                    Points = Inv
                };
                GearElements.Add(InvoluteLine);
            }

            double[] Theta2Range = Generate.LinearRange(0 + s_1_angle, p / (d1 / 2), 2 * Math.PI + s_1_angle);
            foreach (double th in Theta2Range)
            {
                PointCollection Inv = InvolutePoints(th, dTheta, d_b1 / 2, d_a1 / 2, new Point(0, 0), true);
                Polyline InvoluteLine = new Polyline
                {
                    Stroke = Brushes.Red,
                    StrokeThickness = 1,
                    Points = Inv
                };
                GearElements.Add(InvoluteLine);
            }

            GearCharacteristicsData Pinion = new GearCharacteristicsData
            {
                NumberOfTeeths = z1,
                ShiftCoefficient = x1,
                ReferencePitchDiameter = d1,
                OperatingPitchDiameter = d_prime1,
                DedendumDiameter = d_f1,
                AddendumDiameter = d_a1,
                BaseCircleDiameter = d_b1,
                OperatingClearance = 0, //TODO
                ThicknessReference = 0,
                ThicknessOperating = 0,
                ThinknessBase = 0,
                ThicknessTip = 0,
            };

            GearCharacteristicsData Gear = new GearCharacteristicsData
            {
                NumberOfTeeths = z2,
                ShiftCoefficient = x2,
                ReferencePitchDiameter = d2,
                OperatingPitchDiameter = d_prime2,
                DedendumDiameter = d_f2,
                AddendumDiameter = d_a2,
                BaseCircleDiameter = d_b2,
                OperatingClearance = 0, //TODO
                ThicknessReference = 0,
                ThicknessOperating = 0,
                ThinknessBase = 0,
                ThicknessTip = 0,
            };

            GearMechanismData MechanismData = new GearMechanismData
            {
                Module = m,
                PressureAngle = 20,
                OperatingPressureAngle = Degrees(alpha_prime),
                Clearance = 0, //TODO
                CenterDistance = a,
                CenterDistanceCoefficient = y,
                TransmissionRatio = i,
                ContactRatio = epsilon,
                Pitch = p,
                FilletRadius = rho
            };

            CalculationsResultsData Result = new CalculationsResultsData
            {
                GearData = Gear,
                PinionData = Pinion,
                MechanismData = MechanismData,
                GearGeometry = GearElements,
                ActionPosition = new Point(-Pinion.OperatingPitchDiameter, 0)
            };

            return Result;
        }
    }
}
