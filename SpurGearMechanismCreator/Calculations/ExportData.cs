using System;
using System.Collections.Generic;
using System.Text;
using System.Windows;
using System.Windows.Media;

namespace SpurGearMechanismCreator.Calculations
{
	public static class ExportData
	{
		public static string GenerateExportData(PointCollection Points)
		{
			for (int i = 0; i < 50; i++) 
				Points.RemoveAt(Points.Count - 1);

			var Builder = new StringBuilder();
			foreach (Point Item in Points)
			{
				Builder.Append(Item.X).Append(" ")
					   .Append(Item.Y).Append(" ")
					   .Append(0).Append(Environment.NewLine);
			}
			return Builder.ToString();
		}
	}
}
