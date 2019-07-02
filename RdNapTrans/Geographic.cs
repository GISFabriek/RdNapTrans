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
    /// Wraps a Geographic value with 3 dimensions
    /// </summary>
    public class Geographic
	{
        /// <summary>
        /// Initializes a new instance of the <see cref="T:RdNapTrans.Geographic" /> class.
        /// </summary>
        /// <param name="phi">The Latitude.</param>
        /// <param name="lambda">The Longitude.</param>
        /// <param name="h">The ellipsoidal height.</param>
        public Geographic(double phi, double lambda, double h)
		{
			Phi = phi;
			Lambda = lambda;
			H = h;
		}

        /// <summary>
        /// Initializes a new instance of the <see cref="T:RdNapTrans.Geographic" /> class.
        /// </summary>
        /// <param name="phi">The Latitude.</param>
        /// <param name="lambda">The Longitude.</param>
        public Geographic(double phi, double lambda) : this(phi, lambda, 0)
		{
		}

        /// <summary>
        /// Latitude in degrees
        /// </summary>
        /// <value>The Latitude.</value>
        public double Phi { get; }
        /// <summary>
        /// Longitude in degrees
        /// </summary>
        /// <value>The Longitude.</value>
        public double Lambda { get; }
        /// <summary>
        /// Ellipsoidal height
        /// </summary>
        /// <value>The ellipsoidal height.</value>
        public double H { get; }


        /// <summary>
        /// Creates a copy of the existing instance and adds an ellipsoidal height to it.
        /// </summary>
        /// <param name="h">The ellipsoidal height to be added.</param>
        /// <returns>Geographic.</returns>
        public virtual Geographic WithH(double h)
		{
			return new Geographic(Phi, Lambda, h);
		}
        /// <summary>
        /// Returns a <see cref="System.String" /> that represents this instance.
        /// </summary>
        /// <returns>A <see cref="System.String" /> that represents this instance.</returns>
        public override string ToString()
        {
            return $"Phi: {Phi}; Lambda: {Lambda}; H {H}";
        }
    }

}