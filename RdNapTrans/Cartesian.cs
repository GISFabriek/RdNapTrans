// ***********************************************************************
// Assembly         : RdNapTrans
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-02-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-02-2019
// ***********************************************************************
// C# PORT from https://github.com/PDOK/rdnaptrans-java
// ***********************************************************************
namespace RdNapTrans
{

    /// <summary>
    /// Wraps a Cartesian value with 3 dimensions
    /// </summary>
    public class Cartesian
	{

        /// <summary>
        /// Constructor for Cartesian.
        /// </summary>
        /// <param name="x">X coordinate.</param>
        /// <param name="y">Y coordinate.</param>
        /// <param name="z">Z coordinate.</param>
        public Cartesian(double x, double y, double z)
        {
            X = x;
			Y = y;
			Z = z;
		}

        /// <summary>
        /// Constructor for Cartesian.
        /// </summary>
        /// <param name="x">X coordinate.</param>
        /// <param name="y">Y coordinate.</param>
        public Cartesian(double x, double y) : this(x, y, 0)
		{
		}

        /// <summary>
        /// Gets the X Coordinate.
        /// </summary>
        /// <value>The X Coordinate.</value>
        public double X { get; }
        /// <summary>
        /// Gets the Y Coordinate.
        /// </summary>
        /// <value>The Y Coordinate.</value>
        public double Y { get; }
        /// <summary>
        /// Gets the Z Coordinate.
        /// </summary>
        /// <value>The Z Coordinate.</value>
        public double Z { get; }


        /// <summary>
        /// Creates a copy of the existing instance and adds a Z Coordinate to it.
        /// </summary>
        /// <param name="z">The Z Coordinate to be added.</param>
        /// <returns>Cartesian.</returns>
        public virtual Cartesian WithZ(double z)
		{
			return new Cartesian(X,Y, z);
		}

        /// <summary>
        /// Returns a <see cref="System.String" /> that represents this instance.
        /// </summary>
        /// <returns>A <see cref="System.String" /> that represents this instance.</returns>
        public override string ToString()
        {
            return $"X: {X}; Y: {Y}; Z {Z}";
        }
    }

}