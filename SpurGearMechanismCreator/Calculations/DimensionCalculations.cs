using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using System.Windows.Media;
using System.Windows.Shapes;

namespace SpurGearMechanismCreator.Calculations
{
    public enum CurveType
    {
        Dedendum,
        RisingInvolute,
        ReturningInvolute,
        Addendum
    }

    public class PolarPoint
    {
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
        public GearCharacteristicsData? PinionData { get; set; }
        public GearCharacteristicsData? GearData { get; set; }
        public GearMechanismData? MechanismData { get; set; }
        public List<UIElement>? MechanismGeometry { get; set; }
        public Point PinionPosition { get; set; }
        public Point GearPosition { get; set; }
        public Point ActionPosition { get; set; }
        public PointCollection? PinionPoints { get; set; }
        public PointCollection? GearPoints { get; set; }
    }

    public class GearCharacteristicsData
    {
        public int NumberOfTeeth { get; set; }
        public double ShiftCoefficient { get; set; }
        public double ReferencePitchDiameter { get; set; }
        public double OperatingPitchDiameter { get; set; }
        public double DedendumDiameter { get; set; }
        public double AddendumDiameter { get; set; }
        public double BaseCircleDiameter { get; set; }
        public double ThicknessReference { get; set; }
        public double ThicknessOperating { get; set; }
        public double ThicknessBase { get; set; }
        public double ThicknessTip { get; set; }
        public double AngleTip { get; set; }
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
        public static double Radians(double Angle)
        {
            return Math.PI / 180 * Angle;
        }

        public static double Degrees(double Radians)
        {
            return 180 / Math.PI * Radians;
        }

        public static double Involute(double AngleRad)
        {
            return Math.Tan(AngleRad) - AngleRad;
        }

        public static double InverseInvolute(double Involute)
        {
            double Angle1 = 0;

            while (true)
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

        public static Point TranslatePoint(Point P, double XOffset, double YOffset)
        {
            return new Point(P.X + XOffset, P.Y + YOffset);
        }

        public static Point RotatePointAAroundB(Point A, Point B, double Angle)
        {
            return new Point
            {
                X = (Math.Cos(Angle) * (A.X - B.X)) - (Math.Sin(Angle) * (A.Y - B.Y)) + B.X,
                Y = (Math.Sin(Angle) * (A.X - B.X)) + (Math.Cos(Angle) * (A.Y - B.Y)) + B.Y
            };
        }

        public static UIElement[] GenerateGearCirclesGeometry(Point Center,
            double DedendumDiameter, double BaseDiameter, double ReferencePitchDiameter,
            double WorkingPitchDiameter, double AddendumDiameter)
        {
            List<UIElement> Elements = new();

            EllipseGeometry DedendumGeometry = new()
            {
                Center = Center,
                RadiusX = DedendumDiameter / 2,
                RadiusY = DedendumDiameter / 2
            };

            EllipseGeometry BaseGeometry = new()
            {
                Center = Center,
                RadiusX = BaseDiameter / 2,
                RadiusY = BaseDiameter / 2
            };

            EllipseGeometry RefPitchGeometry = new()
            {
                Center = Center,
                RadiusX = ReferencePitchDiameter / 2,
                RadiusY = ReferencePitchDiameter / 2
            };

            EllipseGeometry WorkPitchGeometry = new()
            {
                Center = Center,
                RadiusX = WorkingPitchDiameter / 2,
                RadiusY = WorkingPitchDiameter / 2
            };

            EllipseGeometry AddendumGeometry = new()
            {
                Center = Center,
                RadiusX = AddendumDiameter / 2,
                RadiusY = AddendumDiameter / 2
            };

            Path DedendumPath = new()
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Data = DedendumGeometry
            };
            Elements.Add(DedendumPath);

            Path BasePath = new()
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                StrokeDashArray = DoubleCollection.Parse("1,1"),
                Data = BaseGeometry
            };
            Elements.Add(BasePath);

            Path RefPitchPath = new()
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                StrokeDashArray = DoubleCollection.Parse("3,1"),
                Data = RefPitchGeometry
            };
            Elements.Add(RefPitchPath);

            Path WorkPitchPath = new()
            {
                Stroke = Brushes.Black,
                StrokeThickness = 1,
                Data = WorkPitchGeometry
            };
            Elements.Add(WorkPitchPath);

            Path AddendumPath = new()
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Data = AddendumGeometry
            };
            Elements.Add(AddendumPath);

            return Elements.ToArray();
        }

        public static SortedDictionary<double, CurveType> GenerateAngleData(
            double dTheta, double Teeth, double InvoluteAngle,
            double ToothSpacingAngle, double TipAngle, double StartAngleOffset)
        {
            SortedDictionary<double, CurveType> GearAngleData = new();
            const double InvoluteOffset = 0.0001;

            for (int j = 0; j < Teeth; j++)
            {
                GearAngleData.Add(StartAngleOffset + (j * ToothSpacingAngle) + InvoluteOffset, CurveType.RisingInvolute);
                GearAngleData.Add(InvoluteAngle + StartAngleOffset + (j * ToothSpacingAngle) - InvoluteOffset, CurveType.RisingInvolute);
            }

            double[] Tip = Generate.LinearSpaced(5, StartAngleOffset + InvoluteAngle, StartAngleOffset + InvoluteAngle + TipAngle);
            for (int j = 0; j < Teeth; j++)
            {
                foreach (double Item in Tip.Select(n => n + (j * ToothSpacingAngle)))
                {
                    GearAngleData.Add(Item, CurveType.Addendum);
                }
            }

            for (int j = 0; j < Teeth; j++)
            {
                GearAngleData.Add(StartAngleOffset + InvoluteAngle + TipAngle + (j * ToothSpacingAngle) + InvoluteOffset, CurveType.ReturningInvolute);
                GearAngleData.Add(StartAngleOffset + (2 * InvoluteAngle) + TipAngle + (j * ToothSpacingAngle) - InvoluteOffset, CurveType.ReturningInvolute);
            }

            double[] Dwell = Generate.LinearSpaced(5, StartAngleOffset + (2 * InvoluteAngle) + TipAngle, StartAngleOffset + ToothSpacingAngle);
            for (int j = 0; j < Teeth; j++)
            {
                foreach (double Item in Dwell.Select(n => n + (j * ToothSpacingAngle)))
                {
                    GearAngleData.Add(Item, CurveType.Dedendum);
                }
            }


            return GearAngleData;
        }

        public static Point[] GenerateInvoluteProfile(double dTheta, double BaseRadius, double AddendumRadius,
            double InvoluteAngle, bool IsDirectionInverted)
        {
            double Alpha = InverseInvolute(InvoluteAngle);
            int NegateDirection = IsDirectionInverted ? -1 : 1;

            List<Point> List = new();
            for (double i = 0; i < Alpha * 1.5; i += dTheta / 2)
            {
                Point InvolutePoint = new()
                {
                    X = BaseRadius * (Math.Cos(i) + (i * Math.Sin(i))),
                    Y = NegateDirection * BaseRadius * (Math.Sin(i) - (i * Math.Cos(i)))
                };

                if (Polar(InvolutePoint).Rho <= AddendumRadius)
                {
                    List.Add(InvolutePoint);
                }
            }
            return List.ToArray();
        }

        public static PointCollection GenerateGearProfile(double dTheta,
            double BaseRadius, double DedendumRadius, double AddendumRadius,
            SortedDictionary<double, CurveType> AngleCollection, Point Center)
        {
            PointCollection Result = new();

            double? InvoluteMaxAngle = null;
            for (int i = 1; i < AngleCollection.Count; i++)
            {
                KeyValuePair<double, CurveType> Data = AngleCollection.ElementAt(i);
                double Theta = Data.Key;
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
                        {
                            continue;
                        }

                        double RisingInvAlpha = AngleCollection.ElementAt(i).Key - AngleCollection.ElementAt(i - 1).Key;
                        Point[] RisingInvolute = GenerateInvoluteProfile(dTheta, BaseRadius, AddendumRadius, RisingInvAlpha, false);

                        foreach (Point Item in RisingInvolute)
                        {
                            Point RaisingClampedPoint = Item;

                            if (DedendumRadius > BaseRadius)
                            {
                                PolarPoint RaisingPolarPoint = Polar(Item);
                                if (RaisingPolarPoint.Rho < DedendumRadius)
                                {
                                    RaisingClampedPoint = Cartesian(DedendumRadius, RaisingPolarPoint.Theta);
                                }
                            }

                            Point RisingTranslatedPoint = TranslatePoint(RaisingClampedPoint, Center.X, Center.Y);
                            Point RisingRotatedPoint = RotatePointAAroundB(RisingTranslatedPoint, Center,
                                    AngleCollection.ElementAt(i - 1).Key);

                            Result.Add(RisingRotatedPoint);
                        }
                        break;

                    case CurveType.ReturningInvolute:
                        if (AngleCollection.ElementAt(i - 1).Value != CurveType.ReturningInvolute)
                        {
                            continue;
                        }

                        double ReturningInvAlpha = AngleCollection.ElementAt(i).Key - AngleCollection.ElementAt(i - 1).Key;
                        Point[] ReturningInvolute = GenerateInvoluteProfile(dTheta, BaseRadius, AddendumRadius, ReturningInvAlpha, true);

                        foreach (Point Item in ReturningInvolute.Reverse())
                        {
                            Point ReturningClampedPoint = Item;

                            if (DedendumRadius > BaseRadius)
                            {
                                PolarPoint ReturningPolarPoint = Polar(Item);
                                if (ReturningPolarPoint.Rho < DedendumRadius)
                                {
                                    ReturningClampedPoint = Cartesian(DedendumRadius, ReturningPolarPoint.Theta);
                                }
                            }

                            Point ReturningTranslatedPoint = TranslatePoint(ReturningClampedPoint, Center.X, Center.Y);
                            Point ReturningRotatedPoint = RotatePointAAroundB(ReturningTranslatedPoint, Center,
                                    AngleCollection.ElementAt(i).Key);

                            Result.Add(ReturningRotatedPoint);
                        }
                        break;

                    case CurveType.Addendum:
                        InvoluteMaxAngle ??= Theta - AngleCollection.ElementAt(i - 1).Key;

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
            double dTheta = 0.1;
            double alpha = Radians(20);

            double i = (double)z2 / z1;

            double inv_alpha_prime = (2 * Math.Tan(alpha) * (x1 + x2) / (z1 + z2)) + Involute(alpha);
            double alpha_prime = InverseInvolute(inv_alpha_prime);

            double y = (double)(z1 + z2) / 2 * ((Math.Cos(alpha) / Math.Cos(alpha_prime)) - 1);
            double a = (((double)(z1 + z2) / 2) + y) * m;

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
            //double h_a1 = (1 + x1) * m;
            //double h_a2 = (1 + x2) * m;

            // Addendum circle
            double d_a1 = d1 + (2 * h_a1);
            double d_a2 = d2 + (2 * h_a2);

            // Dedendum circle
            double h = (2.25 + y - (x1 + x2)) * m;
            //double h = 2.25 * m;
            double d_f1 = d_a1 - (2 * h);
            double d_f2 = d_a2 - (2 * h);

            // Overlap coefficient
            double epsilon = (Math.Sqrt(Math.Pow(d_a1 / 2, 2) - Math.Pow(d_b1 / 2, 2))
                           + Math.Sqrt(Math.Pow(d_a2 / 2, 2) - Math.Pow(d_b2 / 2, 2))
                           - (a * Math.Sin(alpha_prime))) / (Math.PI * m * Math.Cos(alpha));

            //Pitch 
            double p1 = Math.PI * d1 / z1;
            double p2 = Math.PI * d2 / z2;
            double p = Math.PI * m;
            //double spacing_1 = p / (d1 / 2);

            // Arc length of tooth at the reference pitch circle
            double s_1 = m * ((Math.PI / 2) + (2 * x1 * Math.Tan(alpha)));
            double s_2 = m * ((Math.PI / 2) + (2 * x2 * Math.Tan(alpha)));

            // Arc length of tooth at the working pitch circle
            double sw_1 = d_prime1 * ((s_1 / d1) - Involute(alpha_prime) + Involute(alpha));
            double sw_2 = d_prime2 * ((s_2 / d2) - Involute(alpha_prime) + Involute(alpha));

            // Arc length of tooth at the base pitch circle
            double sb_1 = d_b1 * ((sw_1 / d_prime1) + Involute(alpha_prime));
            double sb_2 = d_b2 * ((sw_2 / d_prime2) + Involute(alpha_prime));

            // InverseInvolute angle of whole involute curve
            double alpha_a1 = Math.Acos(d1 / d_a1 * Math.Cos(alpha));
            double alpha_a2 = Math.Acos(d2 / d_a2 * Math.Cos(alpha));

            // Arc length of tooth at the base pitch circle
            double sa_1 = d_a1 * ((sb_1 / d_b1) - Involute(alpha_a1));
            double sa_2 = d_a2 * ((sb_2 / d_b2) - Involute(alpha_a2));

            double tip_angle1 = 2 * sa_1 / d_a1;
            double tip_angle2 = 2 * sa_2 / d_a2;

            double ang = 2 * s_1 / d1;
            double angw = 2 * sw_1 / d_prime1;
            double angb = 2 * sb_1 / d_b1;
            double anga = 2 * sa_1 / d_a1;

            double test = Math.Acos(d1 / d1 * Math.Cos(alpha));
            double testw = Math.Acos(d1 / d_prime1 * Math.Cos(alpha));
            double testb = Math.Acos(d1 / d_b1 * Math.Cos(alpha));
            double testa = Math.Acos(d1 / d_a1 * Math.Cos(alpha));

            double rho = 0.38 * m;

            List<UIElement> GearElements = new();
            GearElements.AddRange(GenerateGearCirclesGeometry(new Point(0, 0), d_f1, d_b1, d1, d_prime1, d_a1));
            GearElements.AddRange(GenerateGearCirclesGeometry(new Point(a, 0), d_f2, d_b2, d2, d_prime2, d_a2));

            SortedDictionary<double, CurveType> Data1 = GenerateAngleData(dTheta, z1, Involute(alpha_a1), 2 * Math.PI / z1, tip_angle1, Involute(alpha_prime));
            PointCollection Points1 = GenerateGearProfile(dTheta, d_b1 / 2, d_f1 / 2, d_a1 / 2, Data1, new Point(0, 0));
            Polygon InvoluteLine1 = new()
            {
                Stroke = Brushes.DarkOrange,
                StrokeThickness = 0.75,
                Points = Points1
            };
            GearElements.Add(InvoluteLine1);

            double Offset = ((double)1 / 2 * Math.PI) - (2 * sb_2 / d_b2) + Involute(alpha_prime);
            SortedDictionary<double, CurveType> Data2 = GenerateAngleData(dTheta, z2, Involute(alpha_a2), 2 * Math.PI / z2, tip_angle2, Offset);
            PointCollection Points2 = GenerateGearProfile(dTheta, d_b2 / 2, d_f2 / 2, d_a2 / 2, Data2, new Point(a, 0));
            Polygon InvoluteLine2 = new()
            {
                Stroke = Brushes.Red,
                StrokeThickness = 0.75,
                Points = Points2
            };
            GearElements.Add(InvoluteLine2);

            GearCharacteristicsData Pinion = new()
            {
                NumberOfTeeth = z1,
                ShiftCoefficient = x1,
                ReferencePitchDiameter = d1,
                OperatingPitchDiameter = d_prime1,
                DedendumDiameter = d_f1,
                AddendumDiameter = d_a1,
                BaseCircleDiameter = d_b1,
                ThicknessReference = s_1,
                ThicknessOperating = sw_1,
                ThicknessBase = sb_1,
                ThicknessTip = sa_1,
                AngleTip = alpha_a1
            };

            GearCharacteristicsData Gear = new()
            {
                NumberOfTeeth = z2,
                ShiftCoefficient = x2,
                ReferencePitchDiameter = d2,
                OperatingPitchDiameter = d_prime2,
                DedendumDiameter = d_f2,
                AddendumDiameter = d_a2,
                BaseCircleDiameter = d_b2,
                ThicknessReference = s_2,
                ThicknessOperating = sw_2,
                ThicknessBase = sb_2,
                ThicknessTip = sa_2,
                AngleTip = alpha_a2
            };

            GearMechanismData MechanismData = new()
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

            CalculationsResultsData Result = new()
            {
                GearData = Gear,
                PinionData = Pinion,
                MechanismData = MechanismData,
                MechanismGeometry = GearElements,
                PinionPoints = Points1,
                GearPoints = Points2,
                ActionPosition = new Point(Pinion.OperatingPitchDiameter / 2, 0),
                GearPosition = new Point((Pinion.OperatingPitchDiameter / 2) + (Gear.OperatingPitchDiameter / 2), 0),
                PinionPosition = new Point(0, 0)
            };

            return Result;
        }
    }
}
