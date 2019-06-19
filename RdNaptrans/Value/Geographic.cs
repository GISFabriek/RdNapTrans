namespace RdNaptrans.Value
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
		 **    phi      latitude in degrees
		 **    lambda   longitude in degrees
		 **    h        ellipsoidal height
		*/

        /// <summary>
		/// <para>Constructor for Geographic.</para>
		/// </summary>
		/// <param name="phi"> a double. </param>
		/// <param name="lambda"> a double. </param>
		/// <param name="h"> a double. </param>
		public Geographic(double phi, double lambda, double h)
		{
			Phi = phi;
			Lambda = lambda;
			H = h;
		}

		/// <summary>
		/// <para>Constructor for Geographic.</para>
		/// </summary>
		/// <param name="phi"> a double. </param>
		/// <param name="lambda"> a double. </param>
		public Geographic(double phi, double lambda) : this(phi, lambda, 0)
		{
		}

        public double Phi { get; }

        public double Lambda { get; }

        public double H { get; }

        /// <summary>
		/// <para>withH.</para>
		/// </summary>
		/// <param name="h"> a double. </param>
		/// <returns> a <seealso cref="rdnaptrans.value.Geographic"/> object. </returns>
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