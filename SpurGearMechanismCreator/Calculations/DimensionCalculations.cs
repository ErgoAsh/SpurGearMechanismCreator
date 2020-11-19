using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Security.Cryptography.Xml;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace SpurGearMechanismCreator
{
    public enum CurveType
    {
        Dedendum,
        RisingInvolute,
        ReturningInvolute,
        Addendum
    }

    public class PolarPoint { 
        public double Rho { get; set; }
        public double Theta { get; set; }

        public PolarPoint(double Rho = 0, double Theta = 0)
        {
            this.Rho = Rho;
            this.Theta = Theta;
        }
    }

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
		public PointCollection PinionPoints { get; set; }
        public PointCollection GearPoints { get; set; }
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
        public double ThicknessReference { get; set; }
        public double ThicknessOperating { get; set; } 
        public double ThinknessBase { get; set; }
        public double ThicknessTip { get; set; }
    }

    public class GearMechanismData
	{
        public double Module { get; set; }
        public double PressureAngle { get; set; }
        public double OperatingPressureAngle { get; set; }
        public double CenterDistance { get; set; }
        public double CenterDistanceCoefficient { get; set; }
        public double TransmissionRatio { get; set; }
        public double ContactRatio { get; set; }
        public double Pitch { get; set; }
		public double FilletRadius { get; set; }
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

            while(true)
            {
                double Angle2 = Angle1;
                Angle1 = Math.Atan(Angle1 + Involute);

                double Diff = Math.Abs(Angle1 - Angle2);
                if (Diff < Math.Pow(10, -10))
                {
                    break;
                }
            }

            return Angle1;
        }

        public static Point Cartesian(double Rho, double Phi)
		{
            return new Point
            {
                X = Rho * Math.Cos(Phi),
                Y = Rho * Math.Sin(Phi)
            };
		}

        public static PolarPoint Polar(Point P)
        {
            return new PolarPoint
            {
                Rho = Math.Sqrt(Math.Pow(P.X, 2) + Math.Pow(P.Y, 2)),
                Theta = Math.Atan(P.Y / P.X)
            };
        }

        public static double Distance(Point A, Point B)
		{
            return Math.Sqrt(Math.Pow(B.X - A.X, 2) + Math.Pow(B.Y - A.Y, 2));
		}

        public static Point TranslatePoint(Point P, double XOffset, double YOffset)
        {
            return new Point(P.X + XOffset, P.Y + YOffset);
        }

        public static Point RotatePointAAroundB(Point A, Point B, double Angle)
        {
            return new Point{
                X = Math.Cos(Angle) * (A.X - B.X) - Math.Sin(Angle) * (A.Y - B.Y) + B.X,
                Y = Math.Sin(Angle) * (A.X - B.X) + Math.Cos(Angle) * (A.Y - B.Y) + B.Y
            };
        }

        public static UIElement[] GenerateGearCirclesGeometry(Point Center, 
            double DedendumDiameter, double BaseDiameter, double ReferencePitchDiameter,
            double WorkingPitchDiameter, double AddendumDiameter)
        {
            var Elements = new List<UIElement>();

            EllipseGeometry DedendumGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = DedendumDiameter / 2,
                RadiusY = DedendumDiameter / 2
            };

            EllipseGeometry BaseGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = BaseDiameter / 2,
                RadiusY = BaseDiameter / 2
            };

            EllipseGeometry RefPitchGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = ReferencePitchDiameter / 2,
                RadiusY = ReferencePitchDiameter / 2
            };

            EllipseGeometry WorkPitchGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = WorkingPitchDiameter / 2,
                RadiusY = WorkingPitchDiameter / 2
            };

            EllipseGeometry AddendumGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = AddendumDiameter / 2,
                RadiusY = AddendumDiameter / 2
            };

            Path DedendumPath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Data = DedendumGeometry
            };
            Elements.Add(DedendumPath);

            Path BasePath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                StrokeDashArray = DoubleCollection.Parse("1,1"),
                Data = BaseGeometry
            };
            Elements.Add(BasePath);

            Path RefPitchPath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                StrokeDashArray = DoubleCollection.Parse("3,1"),
                Data = RefPitchGeometry
            };
            Elements.Add(RefPitchPath);

            Path WorkPitchPath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 1,
                Data = WorkPitchGeometry
            };
            Elements.Add(WorkPitchPath);

            Path AddendumPath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Data = AddendumGeometry
            };
            Elements.Add(AddendumPath);

            return Elements.ToArray();
        }

        public static SortedDictionary<double, CurveType> GenerateAngleData(
            double dTheta, double Teeths, double InvoluteAngle, 
            double ToothSpacingAngle, double TipAngle, double StartAngleOffset)
		{
            var GearAngleData = new SortedDictionary<double, CurveType>();

            double[] RisingInvolute = Generate.LinearRange(StartAngleOffset, dTheta, InvoluteAngle + StartAngleOffset);
            for (int j = 0; j < Teeths; j++)
            {
                GearAngleData.Add(RisingInvolute.Select(n => n + (j * ToothSpacingAngle)).First(), CurveType.RisingInvolute);
                GearAngleData.Add(RisingInvolute.Select(n => n + (j * ToothSpacingAngle)).Last(), CurveType.RisingInvolute);
            }

            double[] Tip = Generate.LinearRange(StartAngleOffset + InvoluteAngle, dTheta, StartAngleOffset + InvoluteAngle + TipAngle);
            for (int j = 0; j < Teeths; j++)
            {
                foreach (var Item in Tip.Select(n => n + (j * ToothSpacingAngle)))
                {
                    if (!GearAngleData.ContainsKey(Item))
                    {
                        GearAngleData.Add(Item, CurveType.Addendum);
                    }
                }
            }

            double[] ReturningInvolute = Generate.LinearRange(StartAngleOffset + InvoluteAngle + TipAngle, dTheta, StartAngleOffset + 2 * InvoluteAngle + TipAngle);
            for (int j = 0; j < Teeths; j++)
            {
                GearAngleData.Add(ReturningInvolute.Select(n => n + (j * ToothSpacingAngle)).First(), CurveType.ReturningInvolute);
                GearAngleData.Add(ReturningInvolute.Select(n => n + (j * ToothSpacingAngle)).Last(), CurveType.ReturningInvolute);
            }

            double[] Dwell = Generate.LinearRange(StartAngleOffset + 2 * InvoluteAngle + TipAngle, dTheta, StartAngleOffset + ToothSpacingAngle);
            for (int j = 0; j < Teeths; j++)
            {
                foreach (var Item in Dwell.Select(n => n + (j * ToothSpacingAngle)))
                {
                    if (!GearAngleData.ContainsKey(Item))
                    {
                        GearAngleData.Add(Item, CurveType.Dedendum);
                    }
                }
            }

            return GearAngleData;
        }

        public static Point[] GenerateInvoluteProfile(double dTheta, double BaseRadius, double InvoluteAngle, bool IsDirectionInverted)
		{
            double Alpha = InverseInvolute(InvoluteAngle);
            var NegateDirection = IsDirectionInverted ? -1 : 1;

            var List = new List<Point>();
            for (double i = 0; i < Alpha; i += dTheta)
			{
                var InvolutePoint = new Point
                {
                    X =                   BaseRadius * (Math.Cos(i) + i * Math.Sin(i)),
                    Y = NegateDirection * BaseRadius * (Math.Sin(i) - i * Math.Cos(i))
                };
                List.Add(InvolutePoint);
            }
            return List.ToArray();
        }

        public static PointCollection GenerateGearProfile(double dTheta,
            double BaseRadius, double DedendumRadius, double AddendumRadius,
            SortedDictionary<double, CurveType> AngleCollection, Point Center)
		{
            var Result = new PointCollection();

            double? InvoluteMaxAngle = null;
            for (int i = 1; i < AngleCollection.Count; i++)
            {
                var Data = AngleCollection.ElementAt(i);
                var Theta = Data.Key;
                switch (Data.Value)
                {
                    case CurveType.Dedendum:
                        Result.Add(TranslatePoint(new Point
                        {
                            X = DedendumRadius * Math.Cos(Theta),
                            Y = DedendumRadius * Math.Sin(Theta),
                        }, Center.X, Center.Y));
                        break;

                    case CurveType.RisingInvolute:
                        if (AngleCollection.ElementAt(i - 1).Value != CurveType.RisingInvolute)
                            continue;
                        
                        var RisingInvAlpha = AngleCollection.ElementAt(i).Key - AngleCollection.ElementAt(i - 1).Key;
                        var RisingInvolute = GenerateInvoluteProfile(dTheta, BaseRadius, RisingInvAlpha, false);

                        foreach (var Item in RisingInvolute)
						{
                            var RaisingClampedPoint = Item;

                            if (DedendumRadius > BaseRadius)
                            {
                                PolarPoint RaisingPolarPoint = Polar(Item);
                                if (RaisingPolarPoint.Rho < DedendumRadius)
                                    RaisingClampedPoint = Cartesian(DedendumRadius, RaisingPolarPoint.Theta);
                            }

                            var RisingTranslatedPoint = TranslatePoint(RaisingClampedPoint, Center.X, Center.Y);
                            var RisingRotatedPoint = RotatePointAAroundB(RisingTranslatedPoint, Center,
                                    AngleCollection.ElementAt(i - 1).Key);

                            Result.Add(RisingRotatedPoint);
                        }
                        break;

					case CurveType.ReturningInvolute:
                        if (AngleCollection.ElementAt(i - 1).Value != CurveType.ReturningInvolute)
                            continue;

                        var ReturningInvAlpha = AngleCollection.ElementAt(i).Key - AngleCollection.ElementAt(i - 1).Key;
                        var ReturningInvolute = GenerateInvoluteProfile(dTheta, BaseRadius, ReturningInvAlpha, true);

                        foreach (var Item in ReturningInvolute.Reverse())
                        {
                            var ReturningClampedPoint = Item;
                            
                            if (DedendumRadius > BaseRadius)
							{
                                PolarPoint ReturningPolarPoint = Polar(Item);
                                if (ReturningPolarPoint.Rho < DedendumRadius)
                                    ReturningClampedPoint = Cartesian(DedendumRadius, ReturningPolarPoint.Theta);
                            }

                            var ReturningTranslatedPoint = TranslatePoint(ReturningClampedPoint, Center.X, Center.Y);
                            var ReturningRotatedPoint = RotatePointAAroundB(ReturningTranslatedPoint, Center,
                                    AngleCollection.ElementAt(i).Key);

                            Result.Add(ReturningRotatedPoint);
                        }
                        break;

                    case CurveType.Addendum:
                        if (InvoluteMaxAngle == null) InvoluteMaxAngle = Theta - AngleCollection.ElementAt(i - 1).Key;
                        
                        Result.Add(TranslatePoint(new Point
                        {
                            X = AddendumRadius * Math.Cos(Theta),
                            Y = AddendumRadius * Math.Sin(Theta),
                        }, Center.X, Center.Y));
                        break;
				}
			}
            return Result;
		}

        public static CalculationsResultsData Calculate(double m, int z1, int z2, double x1, double x2)
        {
            double dTheta = 0.01;
            double alpha = Radians(20);

            double i = (double) z2 / z1;

            double inv_alpha_prime = 2 * Math.Tan(alpha) * (x1 + x2) / (z1 + z2) + Involute(alpha);
            double alpha_prime = InverseInvolute(inv_alpha_prime);

            double y = (z1 + z2) / 2 * ((Math.Cos(alpha) / Math.Cos(alpha_prime)) - 1);
            double a = ((z1 + z2) / 2 + y) * m;

            // Pitch circle
            double d1 = z1 * m;
            double d2 = z2 * m;

            // Base circle
            double d_b1 = d1 * Math.Cos(alpha);
            double d_b2 = d2 * Math.Cos(alpha);

            // Working pitch diameter
            double d_prime1 = d_b1 / Math.Cos(alpha_prime);
            double d_prime2 = d_b2 / Math.Cos(alpha_prime);

            // Addendum
            double h_a1 = (1 + y - x1) * m;
            double h_a2 = (1 + y - x2) * m;

            // Addendum circle
            double d_a1 = d1 + 2 * h_a1;
            double d_a2 = d2 + 2 * h_a2;

            // Dedendum circle
            double h = (2.25 + y - (x1 + x2)) * m;
            double d_f1 = d_a1 - 2 * h;
            double d_f2 = d_a2 - 2 * h;

            // Overlap coefficient
            double epsilon = (Math.Sqrt(Math.Pow(d_a1 / 2, 2) - Math.Pow(d_b1 / 2, 2))
                + Math.Sqrt(Math.Pow(d_a2 / 2, 2) - Math.Pow(d_b2 / 2, 2))
                - a * Math.Sin(alpha_prime)) / (Math.PI * m * Math.Cos(alpha));

            //Pitch 
            double p1 = Math.PI * d1 / z1;
            double p2 = Math.PI * d2 / z2;
            double p = Math.PI * m;
            //double spacing_1 = p / (d1 / 2);

            // Arc length of tooth at the reference pitch circle
            double s_1 = m * (Math.PI / 2 + 2 * x1 * Math.Tan(alpha));
            double s_2 = m * (Math.PI / 2 + 2 * x2 * Math.Tan(alpha));

            // Arc length of tooth at the working pitch circle
            double sw_1 = d_prime1 * (s_1 / d_prime1 + Involute(alpha) - Involute(alpha_prime));
            double sw_2 = d_prime2 * (s_2 / d_prime2 + Involute(alpha) - Involute(alpha_prime));

            // Angle of tooth thickness at the base pitch circle
            double theta_b1 = 2 * (sw_1 / d1 + Involute(alpha_prime));
            double theta_b2 = 2 * (sw_2 / d2 + Involute(alpha_prime));

            // Arc length of tooth at the working pitch circle
            double sb_1 = d_b1 * theta_b1 / 2;
            double sb_2 = d_b2 * theta_b2 / 2;

            double alpha_a1 = Math.Acos(d1 / d_a1 * Math.Cos(alpha));
            double alpha_a2 = Math.Acos(d2 / d_a2 * Math.Cos(alpha));

            double tip_angle1 = theta_b1 - 2 * Involute(alpha_a1);
            double tip_angle2 = theta_b2 - 2 * Involute(alpha_a2);

            double sa_1 = d_a1 * Involute(alpha_a1) / 2;
            double sa_2 = d_a2 * Involute(alpha_a2) / 2;

            double rho = 0.38 * m;

            List<UIElement> GearElements = new List<UIElement>();
            GearElements.AddRange(GenerateGearCirclesGeometry(new Point(0, 0), d_f1, d_b1, d1, d_prime1, d_a1));
            GearElements.AddRange(GenerateGearCirclesGeometry(new Point(a, 0), d_f2, d_b2, d2, d_prime2, d_a2));
            
            var Data1 = GenerateAngleData(dTheta, z1, Involute(alpha_a1), 2 * Math.PI / z1, tip_angle1, Involute(alpha_prime));
            var Points1 = GenerateGearProfile(dTheta, d_b1 / 2, d_f1 / 2, d_a1 / 2, Data1, new Point(0, 0));
            Polygon InvoluteLine1 = new Polygon
            {
                Stroke = Brushes.DarkOrange,
                StrokeThickness = 0.75,
                Points = Points1
            };
            GearElements.Add(InvoluteLine1);

			var Data2 = GenerateAngleData(dTheta, z2, Involute(alpha_a2), 2 * Math.PI / z2, tip_angle2, (double) 1 / 2 * Math.PI - theta_b2 + Involute(alpha_prime));
            var Points2 = GenerateGearProfile(dTheta, d_b2 / 2, d_f2 / 2, d_a2 / 2, Data2, new Point(a, 0));
            Polygon InvoluteLine2 = new Polygon
            {
                Stroke = Brushes.Red,
                StrokeThickness = 0.75,
                Points = Points2
            };
            GearElements.Add(InvoluteLine2);

            GearCharacteristicsData Pinion = new GearCharacteristicsData
            {
                NumberOfTeeths = z1,
                ShiftCoefficient = x1,
                ReferencePitchDiameter = d1,
                OperatingPitchDiameter = d_prime1,
                DedendumDiameter = d_f1,
                AddendumDiameter = d_a1,
                BaseCircleDiameter = d_b1,
                ThicknessReference = s_1,
                ThicknessOperating = sw_1,
                ThinknessBase = sb_1,
                ThicknessTip = sa_1,
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
                ThicknessReference = s_2,
                ThicknessOperating = sw_2,
                ThinknessBase = sb_2,
                ThicknessTip = sa_2,
            };

            GearMechanismData MechanismData = new GearMechanismData
            {
                Module = m,
                PressureAngle = 20,
                OperatingPressureAngle = Degrees(alpha_prime),
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
                PinionPoints = Points1,
                GearPoints = Points2,
                ActionPosition = new Point(Pinion.OperatingPitchDiameter, 0), //what / 2 
                GearPosition = new Point(Pinion.OperatingPitchDiameter / 2 + Gear.OperatingPitchDiameter / 2, 0),
                PinionPosition = new Point(0, 0)
    };

            return Result;
        }
    }
}
