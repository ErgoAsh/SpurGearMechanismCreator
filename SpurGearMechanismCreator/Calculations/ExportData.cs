using System;
using System.Collections.Generic;
using System.Text;
using System.Windows;
using System.Windows.Media;

namespace SpurGearMechanismCreator.Calculations
{
	public static class ExportData
	{
		public static string GenerateTxtData(PointCollection Points)
		{
			//Points.RemoveAt(Points.Count - 1);

			var Builder = new StringBuilder();
			foreach (var Item in Points)
			{
				Builder.Append(Item.X).Append(" ")
					   .Append(Item.Y).Append(" ")
					   .Append(0).Append(Environment.NewLine);
			}
			return Builder.ToString();
		}
	}
}
