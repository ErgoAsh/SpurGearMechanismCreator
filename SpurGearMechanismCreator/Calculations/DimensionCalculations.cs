using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using System.Windows.Media;
using System.Windows.Shapes;
using MathNet.Numerics;

namespace SpurGearMechanismCreator.Calculations
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
        public List<UIElement> MechanismGeometry { get; set; }
        public Point PinionPosition { get; set; }
        public Point GearPosition { get; set; }
        public Point ActionPosition { get; set; }
		public PointCollection PinionPoints { get; set; }
        public PointCollection GearPoints { get; set; }
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
            return (Math.PI / 180) * Angle;
        }

        public static double Degrees(double Radians)
        {
            return (180 / Math.PI) * Radians;
        }

        public static double Involute(double AngleRad)
        {
            return Math.Tan(AngleRad) - AngleRad;
        }

        public static double InverseInvolute(double Involute) {
            double Angle1 = 0;

            while(true)
            {
                var Angle2 = Angle1;
                Angle1 = Math.Atan(Angle1 + Involute);

                var Diff = Math.Abs(Angle1 - Angle2);
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

            var DedendumGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = DedendumDiameter / 2,
                RadiusY = DedendumDiameter / 2
            };

            var BaseGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = BaseDiameter / 2,
                RadiusY = BaseDiameter / 2
            };

            var RefPitchGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = ReferencePitchDiameter / 2,
                RadiusY = ReferencePitchDiameter / 2
            };

            var WorkPitchGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = WorkingPitchDiameter / 2,
                RadiusY = WorkingPitchDiameter / 2
            };

            var AddendumGeometry = new EllipseGeometry
            {
                Center = Center,
                RadiusX = AddendumDiameter / 2,
                RadiusY = AddendumDiameter / 2
            };

            var DedendumPath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                Data = DedendumGeometry
            };
            Elements.Add(DedendumPath);

            var BasePath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                StrokeDashArray = DoubleCollection.Parse("1,1"),
                Data = BaseGeometry
            };
            Elements.Add(BasePath);

            var RefPitchPath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 0.5,
                StrokeDashArray = DoubleCollection.Parse("3,1"),
                Data = RefPitchGeometry
            };
            Elements.Add(RefPitchPath);

            var WorkPitchPath = new Path
            {
                Stroke = Brushes.Black,
                StrokeThickness = 1,
                Data = WorkPitchGeometry
            };
            Elements.Add(WorkPitchPath);

            var AddendumPath = new Path
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
            var GearAngleData = new SortedDictionary<double, CurveType>();
            const double InvoluteOffset = 0.0001;

            for (var j = 0; j < Teeth; j++)
            {
                GearAngleData.Add(StartAngleOffset + j * ToothSpacingAngle + InvoluteOffset, CurveType.RisingInvolute);
                GearAngleData.Add(InvoluteAngle + StartAngleOffset + j * ToothSpacingAngle - InvoluteOffset, CurveType.RisingInvolute);
            }

            var Tip = Generate.LinearSpaced(5, StartAngleOffset + InvoluteAngle, StartAngleOffset + InvoluteAngle + TipAngle);
            for (var j = 0; j < Teeth; j++)
            {
                foreach (var Item in Tip.Select(n => n + (j * ToothSpacingAngle)))
                {
	                GearAngleData.Add(Item, CurveType.Addendum);
                }
            }

            for (var j = 0; j < Teeth; j++)
            {
                GearAngleData.Add(StartAngleOffset + InvoluteAngle + TipAngle + j * ToothSpacingAngle + InvoluteOffset, CurveType.ReturningInvolute);
                GearAngleData.Add(StartAngleOffset + 2 * InvoluteAngle + TipAngle + j * ToothSpacingAngle - InvoluteOffset, CurveType.ReturningInvolute);
            }

            var Dwell = Generate.LinearSpaced(5, StartAngleOffset + 2 * InvoluteAngle + TipAngle, StartAngleOffset + ToothSpacingAngle);
            for (var j = 0; j < Teeth; j++)
            {
                foreach (var Item in Dwell.Select(n => n + (j * ToothSpacingAngle)))
                {
	                GearAngleData.Add(Item, CurveType.Dedendum);
                }
            }


            return GearAngleData;
        }

        public static Point[] GenerateInvoluteProfile(double dTheta, double BaseRadius, double AddendumRadius,
            double InvoluteAngle, bool IsDirectionInverted)
		{
            var Alpha = InverseInvolute(InvoluteAngle);
            var NegateDirection = IsDirectionInverted ? -1 : 1;

            var List = new List<Point>();
            for (double i = 0; i < Alpha * 1.5; i += dTheta / 2)
			{
                var InvolutePoint = new Point
                {
                    X =                   BaseRadius * (Math.Cos(i) + i * Math.Sin(i)),
                    Y = NegateDirection * BaseRadius * (Math.Sin(i) - i * Math.Cos(i))
                };

                if (Polar(InvolutePoint).Rho <= AddendumRadius)
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
            for (var i = 1; i < AngleCollection.Count; i++)
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
                        var RisingInvolute = GenerateInvoluteProfile(dTheta, BaseRadius, AddendumRadius, RisingInvAlpha, false);

                        foreach (var Item in RisingInvolute)
						{
                            var RaisingClampedPoint = Item;

                            if (DedendumRadius > BaseRadius)
                            {
                                var RaisingPolarPoint = Polar(Item);
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
                        var ReturningInvolute = GenerateInvoluteProfile(dTheta, BaseRadius, AddendumRadius, ReturningInvAlpha, true);

                        foreach (var Item in ReturningInvolute.Reverse())
                        {
                            var ReturningClampedPoint = Item;
                            
                            if (DedendumRadius > BaseRadius)
							{
                                var ReturningPolarPoint = Polar(Item);
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
            var dTheta = 0.1;
            var alpha = Radians(20);

            var i = (double) z2 / z1;

            var inv_alpha_prime = 2 * Math.Tan(alpha) * (x1 + x2) / (z1 + z2) + Involute(alpha);
            var alpha_prime = InverseInvolute(inv_alpha_prime);

            var y = (double) (z1 + z2) / 2 * ((Math.Cos(alpha) / Math.Cos(alpha_prime)) - 1);
            var a = ((double) (z1 + z2) / 2 + y) * m;

            // Pitch circle
            var d1 = z1 * m;
            var d2 = z2 * m;

            // Base circle
            var d_b1 = d1 * Math.Cos(alpha);
            var d_b2 = d2 * Math.Cos(alpha);

            // Working pitch diameter
            var d_prime1 = d_b1 / Math.Cos(alpha_prime);
            var d_prime2 = d_b2 / Math.Cos(alpha_prime);

            // Addendum
            var h_a1 = (1 + y - x1) * m;
            var h_a2 = (1 + y - x2) * m;
            //double h_a1 = (1 + x1) * m;
            //double h_a2 = (1 + x2) * m;

            // Addendum circle
            var d_a1 = d1 + 2 * h_a1;
            var d_a2 = d2 + 2 * h_a2;

            // Dedendum circle
            var h = (2.25 + y - (x1 + x2)) * m;
            //double h = 2.25 * m;
            var d_f1 = d_a1 - 2 * h;
            var d_f2 = d_a2 - 2 * h;

            // Overlap coefficient
            var epsilon = (Math.Sqrt(Math.Pow(d_a1 / 2, 2) - Math.Pow(d_b1 / 2, 2))
                           + Math.Sqrt(Math.Pow(d_a2 / 2, 2) - Math.Pow(d_b2 / 2, 2))
                           - a * Math.Sin(alpha_prime)) / (Math.PI * m * Math.Cos(alpha));

            //Pitch 
            var p1 = Math.PI * d1 / z1;
            var p2 = Math.PI * d2 / z2;
            var p = Math.PI * m;
            //double spacing_1 = p / (d1 / 2);

            // Arc length of tooth at the reference pitch circle
            var s_1 = m * (Math.PI / 2 + 2 * x1 * Math.Tan(alpha));
            var s_2 = m * (Math.PI / 2 + 2 * x2 * Math.Tan(alpha));

            // Arc length of tooth at the working pitch circle
            var sw_1 = d_prime1 * (s_1 / d1 - Involute(alpha_prime) + Involute(alpha));
            var sw_2 = d_prime2 * (s_2 / d2 - Involute(alpha_prime) + Involute(alpha));

            // Arc length of tooth at the base pitch circle
            var sb_1 = d_b1 * (sw_1 / d_prime1 + Involute(alpha_prime));
            var sb_2 = d_b2 * (sw_2 / d_prime2 + Involute(alpha_prime));

            // InverseInvolute angle of whole involute curve
            var alpha_a1 = Math.Acos(d1 / d_a1 * Math.Cos(alpha));
            var alpha_a2 = Math.Acos(d2 / d_a2 * Math.Cos(alpha));

            // Arc length of tooth at the base pitch circle
            var sa_1 = d_a1 * (sb_1 / d_b1 - Involute(alpha_a1));
            var sa_2 = d_a2 * (sb_2 / d_b2 - Involute(alpha_a2));

            var tip_angle1 = 2 * sa_1 / d_a1;
            var tip_angle2 = 2 * sa_2 / d_a2;

            var ang = 2 * s_1 / d1;
            var angw = 2 * sw_1 / d_prime1;
            var angb = 2 * sb_1 / d_b1;
            var anga = 2 * sa_1 / d_a1;

            var test = Math.Acos(d1 / d1 * Math.Cos(alpha));
            var testw = Math.Acos(d1 / d_prime1 * Math.Cos(alpha));
            var testb = Math.Acos(d1 / d_b1 * Math.Cos(alpha));
            var testa = Math.Acos(d1 / d_a1 * Math.Cos(alpha));

            var rho = 0.38 * m;

            var GearElements = new List<UIElement>();
            GearElements.AddRange(GenerateGearCirclesGeometry(new Point(0, 0), d_f1, d_b1, d1, d_prime1, d_a1));
            GearElements.AddRange(GenerateGearCirclesGeometry(new Point(a, 0), d_f2, d_b2, d2, d_prime2, d_a2));
            
            var Data1 = GenerateAngleData(dTheta, z1, Involute(alpha_a1), 2 * Math.PI / z1, tip_angle1, Involute(alpha_prime));
            var Points1 = GenerateGearProfile(dTheta, d_b1 / 2, d_f1 / 2, d_a1 / 2, Data1, new Point(0, 0));
            var InvoluteLine1 = new Polygon
            {
                Stroke = Brushes.DarkOrange,
                StrokeThickness = 0.75,
                Points = Points1
            };
            GearElements.Add(InvoluteLine1);

            var Offset = (double) 1 / 2 * Math.PI - (2 * sb_2 / d_b2) + Involute(alpha_prime);
            var Data2 = GenerateAngleData(dTheta, z2, Involute(alpha_a2), 2 * Math.PI / z2, tip_angle2, Offset);
            var Points2 = GenerateGearProfile(dTheta, d_b2 / 2, d_f2 / 2, d_a2 / 2, Data2, new Point(a, 0));
            var InvoluteLine2 = new Polygon
            {
                Stroke = Brushes.Red,
                StrokeThickness = 0.75,
                Points = Points2
            };
            GearElements.Add(InvoluteLine2);

            var Pinion = new GearCharacteristicsData
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

            var Gear = new GearCharacteristicsData
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

            var MechanismData = new GearMechanismData
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

            var Result = new CalculationsResultsData
            {
                GearData = Gear,
                PinionData = Pinion,
                MechanismData = MechanismData,
                MechanismGeometry = GearElements,
                PinionPoints = Points1,
                GearPoints = Points2,
                ActionPosition = new Point(Pinion.OperatingPitchDiameter / 2, 0),
                GearPosition = new Point(Pinion.OperatingPitchDiameter / 2 + Gear.OperatingPitchDiameter / 2, 0),
                PinionPosition = new Point(0, 0)
            };

            return Result;
        }
    }
}
