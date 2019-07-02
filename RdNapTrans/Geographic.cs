namespace RdNapTrans
{
	/// <summary>
	/// <para>Geographic class.</para>
	/// 
	/// @author raymond
	/// @version $Id: $Id
	/// </summary>
	public class Geographic
	{

		/*
		 **    Phi      latitude in degrees
		 **    Lambda   longitude in degrees
		 **    H        ellipsoidal height
		*/

		public readonly double Phi;
		public readonly double Lambda;
		public readonly double H;

        /// <summary>
        /// <para>Constructor for Geographic.</para>
        /// </summary>
        /// <param name="phi">  Phi coordinate. </param>
        /// <param name="lambda"> Lambda coordinate. </param>
        /// <param name="h"> H coordinate. </param>
        public Geographic(double phi, double lambda, double h)
		{
			Phi = phi;
			Lambda = lambda;
			H = h;
		}

		/// <summary>
		/// <para>Constructor for Geographic.</para>
		/// </summary>
		/// <param name="phi"> Phi coordinate. </param>
		/// <param name="lambda"> Lambda coordinate. </param>
		public Geographic(double phi, double lambda) : this(phi, lambda, 0)
		{
		}

		/// <summary>
		/// <para>withH.</para>
		/// </summary>
		/// <param name="h"> H Coordinate. </param>
		/// <returns> a <seealso cref="Geographic"/> object. </returns>
		public virtual Geographic WithH(double h)
		{
			return new Geographic(Phi, Lambda, h);
		}
        public override string ToString()
        {
            return $"Phi: {Phi}; Lambda: {Lambda}; H {H}";
        }
    }

}