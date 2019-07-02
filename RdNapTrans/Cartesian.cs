namespace RdNapTrans
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
        /// <param name="x"> X coordinate. </param>
        /// <param name="y"> Y coordinate. </param>
        /// <param name="z"> Z coordinate. </param>
        public Cartesian(double x, double y, double z)
		{
			X = x;
			Y = y;
			Z = z;
		}

		/// <summary>
		/// <para>Constructor for Cartesian.</para>
		/// </summary>
		/// <param name="x"> X coordinate. </param>
		/// <param name="y"> Y coordinate. </param>
		public Cartesian(double x, double y) : this(x, y, 0)
		{
		}

        /// <summary>
        /// <para>withZ.</para>
        /// </summary>
        /// <param name="z"> Z coordinate. </param>
        /// <returns> a <seealso cref="Cartesian"/> object. </returns>
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