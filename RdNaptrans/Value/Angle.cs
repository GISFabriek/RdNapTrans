namespace RdNaptrans.Value
{
	/// <summary>
	/// <para>Angle class.</para>
	/// 
	/// @author raymond
	/// @version $Id: $Id
	/// </summary>
	public class Angle
	{
        /// <summary>
		/// <para>Constructor for Angle.</para>
		/// </summary>
		/// <param name="degrees"> a double. </param>
		/// <param name="minutes"> a double. </param>
		/// <param name="seconds"> a double. </param>
		public Angle(double degrees, double minutes, double seconds)
		{
			Degrees = degrees;
			Minutes = minutes;
			Seconds = seconds;
		}

        public double Degrees { get; }

        public double Minutes { get; }

        public double Seconds { get; }
    }

}