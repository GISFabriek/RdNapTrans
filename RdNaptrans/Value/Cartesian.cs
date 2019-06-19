namespace RdNaptrans.Value
{
	/// <summary>
	/// <para>Cartesian class.</para>
	/// 
	/// @author raymond
	/// @version $Id: $Id
	/// </summary>
	public class Cartesian
	{
		public readonly double X;
		public readonly double Y;
		public readonly double Z;

		/// <summary>
		/// <para>Constructor for Cartesian.</para>
		/// </summary>
		/// <param name="x"> a double. </param>
		/// <param name="y"> a double. </param>
		/// <param name="z"> a double. </param>
		public Cartesian(double x, double y, double z)
		{
			X = x;
			Y = y;
			Z = z;
		}

		/// <summary>
		/// <para>Constructor for Cartesian.</para>
		/// </summary>
		/// <param name="x"> a double. </param>
		/// <param name="y"> a double. </param>
		public Cartesian(double x, double y) : this(x, y, 0)
		{
		}

		/// <summary>
		/// <para>withZ.</para>
		/// </summary>
		/// <param name="z"> a double. </param>
		/// <returns> a <seealso cref="rdnaptrans.value.Cartesian"/> object. </returns>
		public virtual Cartesian WithZ(double z)
		{
			return new Cartesian(X,Y, z);
		}

        public override string ToString()
        {
            return $"X: {X}; Y: {Y}; Z {Z}";
        }
    }

}