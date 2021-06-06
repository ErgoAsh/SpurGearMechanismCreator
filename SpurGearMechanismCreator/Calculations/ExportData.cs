using System;
using System.Text;
using System.Windows.Media;

namespace SpurGearMechanismCreator.Calculations
{
    public static class ExportData
    {
        private const char WhiteSpaceCharConst = ' ';

        public static string GenerateTxtData(PointCollection Points)
        {
            //Points.RemoveAt(Points.Count - 1);

            StringBuilder Builder = new();
            foreach (System.Windows.Point Item in Points)
            {
                _ = Builder.Append(Item.X).Append(WhiteSpaceCharConst)
                       .Append(Item.Y).Append(WhiteSpaceCharConst)
                       .Append(0).Append(Environment.NewLine);
            }
            return Builder.ToString();
        }
    }
}
