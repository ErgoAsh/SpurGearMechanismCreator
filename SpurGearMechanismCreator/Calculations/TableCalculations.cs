using System;
using System.Collections.Generic;
using System.Text;
using System.Windows.Controls;

namespace SpurGearMechanismCreator.Calculations
{
	public class TableDataRow
	{
		public string Name { get; set; }
		public string Formula { get; set; }
		public double Value { get; set; }
		public double ValueSecondary { get; set; }
		public double SharedValue { get; set; }
		public bool AreValuesShared { get; set; }
	}

	public static class TableCalculations
	{
		public static List<TableDataRow> GetTableData(CalculationsResultsData Data)
		{
			List<TableDataRow> Result = new List<TableDataRow> {
				new TableDataRow { 
					Name = "Module",
					Formula = "m",
					Value = Data.MechanismData.Module,
					AreValuesShared = true
				},
				new TableDataRow {
					Name = "Number of teeth",
					Formula = "z",
					Value = Data.GearData.NumberOfTeeths,
					ValueSecondary = Data.PinionData.NumberOfTeeths,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Profile shift coefficients",
					Formula = "x",
					Value = Data.GearData.ShiftCoefficient,
					ValueSecondary = Data.PinionData.ShiftCoefficient,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Pressure angle",
					Formula = @"\alpha",
					Value = Data.MechanismData.PressureAngle,
					AreValuesShared = true
				},
				new TableDataRow {
					Name = "Operating pressure angle",
					Formula = @"\alpha' = inverseInvolute(2\tan{\alpha}\frac{x_1+x_2}{z_1+z_2}+inv{\alpha})",
					Value = Data.MechanismData.OperatingPressureAngle,
					AreValuesShared = true
				},
				new TableDataRow {
					Name = "Transmission ratio",
					Formula = @"i = \frac{z_2}{z_1}",
					Value = Data.MechanismData.TransmissionRatio,
					AreValuesShared = true
				},
				new TableDataRow {
					Name = "Center distance modification coefficient",
					Formula = @"y = \frac{z_1+z_2}{2}(\frac{\cos{\alpha}}{\cos{\alpha'}}-1)",
					Value = Data.MechanismData.CenterDistanceCoefficient,
					AreValuesShared = true
				},
				new TableDataRow {
					Name = "Center distance",
					Formula = @"a = (\frac{z_1+z_2}{2}+y)m",
					Value = Data.MechanismData.CenterDistance,
					AreValuesShared = true
				},
				new TableDataRow {
					Name = "Pitch circle diameter",
					Formula = "d = zm",
					Value = Data.GearData.ReferencePitchDiameter,
					ValueSecondary = Data.PinionData.ReferencePitchDiameter,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Base circle diameter",
					Formula = @"d_b = d\cos{\alpha}",
					Value = Data.GearData.BaseCircleDiameter,
					ValueSecondary = Data.PinionData.BaseCircleDiameter,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Operating pitch circle diameter",
					Formula = @"d' = \frac{d_b}{\cos{\alpha'}}",
					Value = Data.GearData.OperatingPitchDiameter,
					ValueSecondary = Data.PinionData.OperatingPitchDiameter,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Addendum circle diameter",
					Formula = @"d_a = d + 2h_a",
					Value = Data.GearData.AddendumDiameter,
					ValueSecondary = Data.PinionData.AddendumDiameter,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Dedendum circle diameter",
					Formula = @"d_f = d - 2h",
					Value = Data.GearData.DedendumDiameter,
					ValueSecondary = Data.PinionData.DedendumDiameter,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Pitch",
					Formula = @"p = \pi m",
					Value = Data.MechanismData.Pitch,
					AreValuesShared = true
				},
				new TableDataRow {
					Name = "Fillet radius",
					Formula = @"\rho = 0.38 m",
					Value = Data.MechanismData.FilletRadius,
					AreValuesShared = true
				},
				new TableDataRow {
					Name = "Tooth thickness at the reference pitch circle",
					Formula = "-",
					Value = Data.GearData.ThicknessReference,
					ValueSecondary = Data.PinionData.ThicknessReference,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Tooth thickness at the operating pitch circle",
					Formula = "-",
					Value = Data.GearData.ThicknessOperating,
					ValueSecondary = Data.PinionData.ThicknessOperating,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Tooth thickness at the addendum pitch circle",
					Formula = "-",
					Value = Data.GearData.ThicknessTip,
					ValueSecondary = Data.PinionData.ThicknessTip,
					AreValuesShared = false
				},
				new TableDataRow {
					Name = "Contact Ratio",
					Formula = @"\epsilon = \frac{\sqrt{\frac{d_{a1}}{2}^2 - \frac{d_{b1}}{2}^2} + \sqrt{\frac{d_{a2}}{2}^2 - \frac{d_{b2}}{2}^2} + a \sin{\alpha}}{\pi m \cos{\alpha}}",
					Value = Data.MechanismData.ContactRatio,
					AreValuesShared = true
				},
			};

			return Result;
		}
	}
}
