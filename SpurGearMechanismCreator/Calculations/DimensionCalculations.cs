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

        public static Point Cartesian(double Rho, double phi)
		{
            return new Point
            {
                X = Rho * Math.Cos(phi),
                Y = Rho * Math.Sin(phi)
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

        public static SortedDictionary<double, CurveType> GenerateAngleData(
            double dTheta, double Teeths, double InvoluteAngle, 
            double ToothSpacingAngle, double TipAngle)
		{
            var GearAngleData = new SortedDictionary<double, CurveType>();

            double[] RisingInvolute = Generate.LinearRange(0, dTheta, InvoluteAngle);
            for (int j = 0; j < Teeths; j++)
            {
                foreach (var Item in RisingInvolute.Select(n => n + (j * ToothSpacingAngle)))
                {
                    if (!GearAngleData.ContainsKey(Item))
                    {
                        GearAngleData.Add(Item, CurveType.RisingInvolute);
                    }
                }
            }

            double[] Tip = Generate.LinearRange(InvoluteAngle, dTheta, InvoluteAngle + TipAngle);
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

            double[] ReturningInvolute = Generate.LinearRange(InvoluteAngle + TipAngle, dTheta, 2 * InvoluteAngle + TipAngle);
            for (int j = 0; j < Teeths; j++)
            {
                foreach (var Item in ReturningInvolute.Select(n => n + (j * ToothSpacingAngle)))
                {
                    if (!GearAngleData.ContainsKey(Item))
                    {
                        GearAngleData.Add(Item, CurveType.ReturningInvolute);
                    }
                }
            }

            double[] Dwell = Generate.LinearRange(2 * InvoluteAngle + TipAngle, dTheta, ToothSpacingAngle);
            for (int j = 0; j < Teeths; j++)
            {
                foreach (var Item in Dwell.Select(n => n + (j * ToothSpacingAngle)))
                {
                    if (!GearAngleData.ContainsKey(Item))
                    {
                        GearAngleData.Add(Item, CurveType.ReturningInvolute);
                    }
                }
            }

            return GearAngleData;
        }

        public static PointCollection GenerateGearProfile(
            double BaseRadius, double DedendumRadius, double AddendumRadius, double StartAngleOffset,
            SortedDictionary<double, CurveType> AngleCollection, Point Center)
		{
            var Result = new PointCollection();
            //TODO StartAngleOffset
            double? InvoluteStartAngle = null;
            foreach (var Data in AngleCollection)
            {
                var Theta = Data.Key;
                switch (Data.Value)
                {
                    case CurveType.Dedendum:
                        InvoluteStartAngle = null;
                        Result.Add(TranslatePoint(new Point
                        {
                            X = DedendumRadius * Math.Cos(Data.Key),
                            Y = DedendumRadius * Math.Sin(Data.Key),
                        }, Center.X, Center.Y));
                        break;

                    case CurveType.RisingInvolute:
                        if (InvoluteStartAngle == null) InvoluteStartAngle = Theta;
                        var RisingThetaInv = Theta - InvoluteStartAngle ?? throw new ArgumentException("Invalid theta value");

                        var RisingInvolutePoint = new Point
                        {
                            X = BaseRadius * (Math.Cos(RisingThetaInv) + RisingThetaInv * Math.Sin(RisingThetaInv)),
                            Y = BaseRadius * (Math.Sin(RisingThetaInv) - RisingThetaInv * Math.Cos(RisingThetaInv))
                        };
                        var RisingTranslatedPoint = TranslatePoint(RisingInvolutePoint, Center.X, Center.Y);
                        var RisingRotatedPoint = RotatePointAAroundB(RisingTranslatedPoint, Center, 
                            InvoluteStartAngle ?? throw new ArgumentException("Invalid theta value"));

                        Result.Add(RisingRotatedPoint);
                        break;

					case CurveType.ReturningInvolute:
                        if (InvoluteStartAngle == null) InvoluteStartAngle = Theta;
                        var ReturningThetaInv = Theta - InvoluteStartAngle ?? throw new ArgumentException("Invalid theta value");

                        var ReturningInvolutePoint = new Point
                        {
                            X =  BaseRadius * (Math.Cos(ReturningThetaInv) + ReturningThetaInv * Math.Sin(ReturningThetaInv)),
                            Y = -BaseRadius * (Math.Sin(ReturningThetaInv) - ReturningThetaInv * Math.Cos(ReturningThetaInv))
                        };
                        var ReturningTranslatedPoint = TranslatePoint(ReturningInvolutePoint, Center.X, Center.Y);
                        var ReturningRotatedPoint = RotatePointAAroundB(ReturningTranslatedPoint, Center,
                                InvoluteStartAngle ?? throw new ArgumentException("Invalid theta value"));

                        Result.Add(ReturningRotatedPoint);
                        break;

                    case CurveType.Addendum:
                        InvoluteStartAngle = null;
                        Result.Add(TranslatePoint(new Point
                        {
                            X = AddendumRadius * Math.Cos(Data.Key),
                            Y = AddendumRadius * Math.Sin(Data.Key),
                        }, Center.X, Center.Y));
                        break;
				}
			}
            return Result;
		}

        public static UIElement[] GenerateGearCirclesGeometry(
            Point Center, double DedendumDiameter, double BaseDiameter, double ReferencePitchDiameter, 
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

        public static CalculationsResultsData Calculate(double m, int z1, int z2, double x1, double x2)
        {
            double dTheta = 0.001;
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

            //double d_f1 = m * (z1 - 2.5 + 2 * x1);
            //double d_f2 = m * (z2 - 2.5 + 2 * x2);

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

            double alpha_a1 = Math.Acos(d1 / d_a1 * Math.Cos(alpha)); //p1 / d1;
            double alpha_a2 = Math.Acos(d2 / d_a2 * Math.Cos(alpha)); //p2 / d2;

            double tip_angle1 = theta_b1 - 2 * Involute(alpha_a1);
            double tip_angle2 = theta_b2 - 2 * Involute(alpha_a2);

            double sa_1 = d_a1 * Involute(alpha_a1);
            double sa_2 = d_a2 * Involute(alpha_a2);

            double rho = 0.38 * m;

            List<UIElement> GearElements = new List<UIElement>();
            GearElements.AddRange(GenerateGearCirclesGeometry(new Point(0, 0), d_f1, d_b1, d1, d_prime1, d_a1));
            GearElements.AddRange(GenerateGearCirclesGeometry(new Point(a, 0), d_f2, d_b2, d2, d_prime2, d_a2));
            
            var Data1 = GenerateAngleData(dTheta, z1, Involute(alpha_a1), 2 * Math.PI / z1, tip_angle1);
            Polyline InvoluteLine1 = new Polyline
            {
                Stroke = Brushes.Red,
                StrokeThickness = 1,
                Points = GenerateGearProfile(d_b1 / 2, d_f1 / 2, d_a1 / 2, Involute(alpha_prime), Data1, new Point(0, 0))
            };
            GearElements.Add(InvoluteLine1);

            var Data2 = GenerateAngleData(dTheta, z1, Involute(alpha_a2), 2 * Math.PI / z1, tip_angle2);
            Polyline InvoluteLine2 = new Polyline
            {
                Stroke = Brushes.Red,
                StrokeThickness = 1,
                Points = GenerateGearProfile(d_b2 / 2, d_f2 / 2, d_a2 / 2, theta_b2 - Involute(alpha_prime), Data2, new Point(a, 0))
            };
            GearElements.Add(InvoluteLine2);


            /*
            double spacing_1 = p / (d1 / 2);
            double base_1 = Involute(alpha_prime);
            double[] ThetaRange_1 = Generate.LinearRange(base_1, spacing_1, 2 * Math.PI + base_1 - dTheta);
            foreach (double th in ThetaRange_1)
			{
                PointCollection Inv = InvolutePoints(th, dTheta, d_b1 / 2, d_f1 / 2, d_a1 / 2, new Point(0, 0), false);
                Polyline InvoluteLine = new Polyline
                {
                    Stroke = Brushes.Red,
                    StrokeThickness = 1,
                    Points = Inv
                };
                GearElements.Add(InvoluteLine);
            }

            double[] Theta2Range_1 = Generate.LinearRange(0 + theta_b1, spacing_1, 2 * Math.PI + theta_b1 - dTheta);
            foreach (double th in Theta2Range_1)
            {
                PointCollection Inv = InvolutePoints(th, dTheta, d_b1 / 2, d_f1 / 2, d_a1 / 2, new Point(0, 0), true);
                Polyline InvoluteLine = new Polyline
                {
                    Stroke = Brushes.Red,
                    StrokeThickness = 1,
                    Points = Inv
                };
                GearElements.Add(InvoluteLine);
            }

            double spacing_2 = p / (d2 / 2);
            double startPosition = theta_b2 - Involute(alpha_prime);
            double[] ThetaRange_2 = Generate.LinearRange(((double) 1/2) * Math.PI - startPosition, spacing_2, ((double) 5/2) * Math.PI - dTheta - startPosition);
            foreach (double th in ThetaRange_2)
            {
                PointCollection Inv = InvolutePoints(th, dTheta, d_b2 / 2, d_f2 / 2, d_a2 / 2, new Point(a, 0), false);
                Polyline InvoluteLine = new Polyline
                {
                    Stroke = Brushes.Red,
                    StrokeThickness = 1,
                    Points = Inv
                };
                GearElements.Add(InvoluteLine);
            }

            double[] Theta2Range_2 = Generate.LinearRange(((double) 1/2) * Math.PI + theta_b2 - startPosition, spacing_2, ((double) 5/2) * Math.PI + theta_b2 - dTheta - startPosition);
            foreach (double th in Theta2Range_2)
            {
                PointCollection Inv = InvolutePoints(th, dTheta, d_b2 / 2, d_f2 / 2, d_a2 / 2, new Point(a, 0), true);
                Polyline InvoluteLine = new Polyline
                {
                    Stroke = Brushes.Red,
                    StrokeThickness = 1,
                    Points = Inv
                };
                GearElements.Add(InvoluteLine);
            }
            */

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
                ActionPosition = new Point(-Pinion.OperatingPitchDiameter, 0)
            };

            return Result;
        }
    }
}
