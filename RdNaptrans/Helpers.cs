using System;
using System.IO;

namespace RdNaptrans
{
	
	using static Constants;
	using Angle = Value.Angle;
	using Cartesian = Value.Cartesian;
	using Geographic = Value.Geographic;
	using GrdFile = Value.GrdFile;

    /// <summary>
    ///   <para>Helpers class.</para>
    /// Converted and adapted from  <a href="https://github.com/PDOK/rdnaptrans-java">https://github.com/PDOK/rdnaptrans-java</a></summary>
    public class Helpers
	{
        /// <summary>Sine for angles in degrees.</summary>
        /// <param name="angleInDegrees">The angle.</param>
        /// <returns>System.Double.</returns>
        internal static double DegSin(double angleInDegrees)
		{
			return Math.Sin(angleInDegrees / 180.0 * Math.PI);
		}

        /// <summary>Cosine for angles in degrees.</summary>
        /// <param name="angleInDegrees">The angle.</param>
        /// <returns>System.Double.</returns>
        internal static double DegCos(double angleInDegrees)
		{
			return Math.Cos(angleInDegrees / 180.0 * Math.PI);
		}

        /// <summary>Tangent for angles in degrees.</summary>
        /// <param name="angleInDegrees">The angle.</param>
        /// <returns>System.Double.</returns>
        internal static double DegTan(double angleInDegrees)
		{
			return Math.Tan(angleInDegrees / 180.0 * Math.PI);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: deg_asin
		**    Description:   inverse sine for angles in degrees
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    a              double      in     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    asin(a)
		**--------------------------------------------------------------
		*/
		internal static double DegAsin(double a)
		{
			return (Math.Asin(a) * 180.0 / Math.PI);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: deg_atan
		**    Description:   inverse tangent for angles in degrees
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    a              double in     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    atan(a)
		**--------------------------------------------------------------
		*/
		internal static double DegAtan(double a)
		{
			return (Math.Atan(a) * 180.0 / Math.PI);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: Atanh
		**    Description:   inverse hyperbolic tangent
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    a              double      in     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    atanh(a)
		**--------------------------------------------------------------
		*/
        internal static double Atanh(double a)
		{
			return (0.5 * Math.Log((1.0 + a) / (1.0 - a)));
		}

		/*
		**--------------------------------------------------------------
		**    Function name: deg_min_sec2decimal
		**    Description:   converts from degrees, minutes and seconds to decimal degrees
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    deg            double      in     req     none
		**    min            double      in     req     none
		**    sec            double      in     req     none
		**    dec_deg        double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    All parameters are doubles, so one can also enter decimal minutes or degrees.
		**    Note: Nonsense input is accepted too.
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		internal static double DegreesMinutesSeconds2Decimal(Angle angle)
		{
			return (angle.Degrees + angle.Minutes / 60.0 + angle.Seconds / 3600.0);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: decimal2deg_min_sec
		**    Description:   converts from decimal degrees to degrees, minutes and seconds
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    dec_deg        double      in     req     none
		**    deg            int         out    -       none
		**    min            int         out    -       none
		**    sec            double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		internal static Angle Decimal2DegreesMinutesSeconds(double decDeg)
		{
			var deg = (int)(decDeg);
			var mathMin = (int)((decDeg - deg) * 60.0);
			var sec = ((decDeg - deg) * 60.0 - mathMin) * 60.0;

			return new Angle(deg, mathMin, sec);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: geographic2cartesian
		**    Description:   from geographic coordinates to cartesian coordinates
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    phi            double      in     req     none
		**    lambda         double      in     req     none
		**    h              double      in     req     none
		**    a              double      in     req     none
		**    inv_f          double      in     req     none
		**    x              double      out    -       none
		**    y              double      out    -       none
		**    z              double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    phi      latitude in degrees
		**    lambda   longitude in degrees
		**    h        ellipsoidal height
		**    a        half major axis of the ellisoid
		**    inv_f    inverse flattening of the ellipsoid
		**    x, y, z  output of cartesian coordinates
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		internal static Cartesian Geographic2Cartesian(Geographic g, double a, double invF)
		{
			/*
			**--------------------------------------------------------------
			**    Source: G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
			**--------------------------------------------------------------
			*/

			/*
			**--------------------------------------------------------------
			**    Explanation of the meaning of variables:
			**        f    flattening of the ellipsoid
			**        ee   first eccentricity squared (e squared in some notations)
			**        n    second (East West) principal radius of curvature (N in some notations)
			**--------------------------------------------------------------
			*/
			var f = 1.0 / invF;
			var ee = f * (2.0 - f);
			var n = a / Math.Sqrt(1.0 - ee * Math.Pow(DegSin(g.Phi),2));

			var x = (n + g.H) * DegCos(g.Phi) * DegCos(g.Lambda);
			var y = (n + g.H) * DegCos(g.Phi) * DegSin(g.Lambda);
			var z = (n * (1.0 - ee) + g.H) * DegSin(g.Phi);

			return new Cartesian(x, y, z);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: cartesian2geographic
		**    Description:   from cartesian coordinates to geographic coordinates
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    x              double      in     req     none
		**    y              double      in     req     none
		**    z              double      in     req     none
		**    a              double      in     req     none
		**    inv_f          double      in     req     none
		**    phi            double      out    -       none
		**    lambda         double      out    -       none
		**    h              double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    x, y, z  input of cartesian coordinates
		**    a        half major axis of the ellisoid
		**    inv_f    inverse flattening of the ellipsoid
		**    phi      output latitude in degrees
		**    lambda   output longitude in degrees
		**    h        output ellipsoidal height
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		internal static Geographic Cartesian2Geographic(Cartesian c, double a, double invF)
		{
			/*
			**--------------------------------------------------------------
			**    Source: G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
			**--------------------------------------------------------------
			*/

			/*
			**--------------------------------------------------------------
			**    Explanation of the meaning of variables:
			**        f    flattening of the ellipsoid
			**        ee   first eccentricity squared (e squared in some notations)
			**        rho  distance to minor axis
			**        n    second (East West) principal radius of curvature (N in some notations)
			**--------------------------------------------------------------
			*/
			var f = 1.0 / invF;
			var ee = f * (2.0 - f);
			var rho = Math.Sqrt(c.X * c.X + c.Y * c.Y);
			double n = 0;

			/*
			**--------------------------------------------------------------
			**    Iterative calculation of phi
			**--------------------------------------------------------------
			*/
			double phi = 0;
            double diff = 90;
			while (diff > DegPrecision)
			{
				var previous = phi;
				n = a / Math.Sqrt(1.0 - ee * Math.Pow(DegSin(phi),2));
				phi = DegAtan(c.Z / rho + n * ee * DegSin(phi) / rho);
				diff = Math.Abs(phi - previous);
			}

			/*
			**--------------------------------------------------------------
			**     Calculation of lambda and h
			**--------------------------------------------------------------
			*/
			var lambda = DegAtan(c.Y / c.X);
			var h = rho * DegCos(phi) + c.Z * DegSin(phi) - n * (1.0 - ee * Math.Pow(DegSin(phi),2));

			return new Geographic(phi, lambda, h);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: sim_trans
		**    Description:   3 dimensional similarity transformation (7 parameters) around another pivot point "a" than the origin
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    x_in           double      in     req     none
		**    y_in           double      in     req     none
		**    z_in           double      in     req     none
		**    tx             double      in     req     none
		**    ty             double      in     req     none
		**    tz             double      in     req     none
		**    alpha          double      in     req     none
		**    beta           double      in     req     none
		**    gamma          double      in     req     none
		**    delta          double      in     req     none
		**    xa             double      in     req     none
		**    ya             double      in     req     none
		**    za             double      in     req     none
		**    x_out          double      out    -       none
		**    y_out          double      out    -       none
		**    z_out          double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    x_in, y_in, z_in     input coordinates
		**    tx                   translation in direction of x axis
		**    ty                   translation in direction of y axis
		**    tz                   translation in direction of z axis
		**    alpha                rotation around x axis in radials
		**    beta                 rotation around y axis in radials
		**    gamma                rotation around z axis in radials
		**    delta                scale parameter (scale = 1 + delta)
		**    xa, ya, za           coordinates of pivot point a (in case of rotation around the center of the ellipsoid these parameters are zero)
		**    x_out, y_out, z_out  output coordinates
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		internal static Cartesian SimTrans(Cartesian input, Cartesian translate, double alpha, double beta, double gamma, double delta, Cartesian pivot)

		{
			/*
			**--------------------------------------------------------------
			**    Source: HTW, "Handleiding voor de Technische Werkzaamheden van het Kadaster". Apeldoorn: Kadaster, 1996.
			**--------------------------------------------------------------
			*/

			/*
			**--------------------------------------------------------------
			**    Calculate the elements of the rotation_matrix:
			**
			**    a b c
			**    d e f
			**    g h i
			**
			**--------------------------------------------------------------
			*/
			var a = Math.Cos(gamma) * Math.Cos(beta);
			var b = Math.Cos(gamma) * Math.Sin(beta) * Math.Sin(alpha) + Math.Sin(gamma) * Math.Cos(alpha);
			var c = -Math.Cos(gamma) * Math.Sin(beta) * Math.Cos(alpha) + Math.Sin(gamma) * Math.Sin(alpha);
			var d = -Math.Sin(gamma) * Math.Cos(beta);
			var e = -Math.Sin(gamma) * Math.Sin(beta) * Math.Sin(alpha) + Math.Cos(gamma) * Math.Cos(alpha);
			var f = Math.Sin(gamma) * Math.Sin(beta) * Math.Cos(alpha) + Math.Cos(gamma) * Math.Sin(alpha);
			var g = Math.Sin(beta);
			var h = -Math.Cos(beta) * Math.Sin(alpha);
			var i = Math.Cos(beta) * Math.Cos(alpha);

			/*
			**--------------------------------------------------------------
			**    Calculate the elements of the vector input_point:
			**    point_2 = input_point - pivot_point
			**--------------------------------------------------------------
			*/
			var x = input.X - pivot.X;
			var y = input.Y - pivot.Y;
			var z = input.Z - pivot.Z;

			/*
			**--------------------------------------------------------------
			**    Calculate the elements of the output vector:
			**    output_point = scale * rotation_matrix * point_2 + translation_vector + pivot_point
			**--------------------------------------------------------------
			*/
			var xOut = (1.0 + delta) * (a * x + b * y + c * z) + translate.X + pivot.X;
			var yOut = (1.0 + delta) * (d * x + e * y + f * z) + translate.Y + pivot.Y;
			var zOut = (1.0 + delta) * (g * x + h * y + i * z) + translate.Z + pivot.Z;

			return new Cartesian(xOut, yOut, zOut);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: rd_projection
		**    Description:   stereographic double projection
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    phi            double      in     req     none
		**    lambda         double      in     req     none
		**    x_rd           double      out    -       none
		**    y_rd           double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    phi         input Bessel latitude in degrees
		**    lambda      input Bessel longitude in degrees
		**    x_rd, rd_y  output RD coordinates
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		internal static Cartesian RdProjection(Geographic input)
		{
			/*
			**--------------------------------------------------------------
			**    Source: G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
			**            G. Strang van Hees, "Globale en lokale geodetische systemen". Delft: Nederlandse Commissie voor Geodesie (NCG), 1997.
			**--------------------------------------------------------------
			*/

			/*
			**--------------------------------------------------------------
			**    Explanation of the meaning of constants:
			**        f                         flattening of the ellipsoid
			**        ee                        first eccentricity squared (e squared in some notations)
			**        e                         first eccentricity
			**        eea                       second eccentricity squared (e' squared in some notations)
			**
			**        phi_amersfoort_sphere     latitude of projection base point Amersfoort on sphere in degrees
			**        lambda_amersfoort_sphere  longitude of projection base point Amersfoort on sphere in degrees
			**
			**        r1                        first (North South) principal radius of curvature in Amersfoort (M in some notations)
			**        r2                        second (East West) principal radius of curvature in Amersfoort (N in some notations)
			**        r_sphere                  radius of sphere
			**
			**        n                         constant of Gaussian projection n = 1.000475...
			**        q_amersfoort              isometric latitude of Amersfoort on ellipsiod
			**        w_amersfoort              isometric latitude of Amersfoort on sphere
			**        m                         constant of Gaussian projection m = 0.003773... (also named c in some notations)
			**--------------------------------------------------------------
			*/
            var f = 1 / InvFBessel;
            var ee = f * (2 - f);
            var e = Math.Sqrt(ee);
            var eea = ee / (1.0 - ee);
            var phiAmersfoortSphere = DegAtan(DegTan(PhiAmersfoortBessel) / Math.Sqrt(1 + eea * Math.Pow(DegCos(PhiAmersfoortBessel),2)));
            var lambdaAmersfoortSphere = LambdaAmersfoortBessel;
            var r1 = ABessel * (1 - ee) / Math.Pow(Math.Sqrt(1 - ee * Math.Pow(DegSin(PhiAmersfoortBessel),2)),3);
            var r2 = ABessel / Math.Sqrt(1.0 - ee * Math.Pow(DegSin(PhiAmersfoortBessel),2));
            var rSphere = Math.Sqrt(r1 * r2);
            var n = Math.Sqrt(1 + eea * Math.Pow(DegCos(PhiAmersfoortBessel),4));
			var qAmersfoort = Atanh(DegSin(PhiAmersfoortBessel)) - e * Atanh(e * DegSin(PhiAmersfoortBessel));
            var wAmersfoort = Math.Log(DegTan(45 + 0.5 * phiAmersfoortSphere));
            var m = wAmersfoort - n * qAmersfoort;

			/*
			**--------------------------------------------------------------
			**    Explanation of the meaning of variables:
			**        q                    isometric latitude on ellipsiod
			**        w                    isometric latitude on sphere
			**        phi_sphere           latitide on sphere in degrees
			**        delta_lambda_sphere  difference in longitude on sphere with Amersfoort in degrees
			**        psi                  distance angle from Amersfoort on sphere
			**        alpha                azimuth from Amersfoort
			**        r                    distance from Amersfoort in projection plane
			**--------------------------------------------------------------
			*/
			var q = Atanh(DegSin(input.Phi)) - e * Atanh(e * DegSin(input.Phi));
			var w = n * q + m;
			var phiSphere = 2 * DegAtan(Math.Exp(w)) - 90;
			var deltaLambdaSphere = n * (input.Lambda - lambdaAmersfoortSphere);
			var sinHalfPsiSquared = Math.Pow(DegSin(0.5 * (phiSphere - phiAmersfoortSphere)),2) + Math.Pow(DegSin(0.5 * deltaLambdaSphere),2) * DegCos(phiSphere) * DegCos(phiAmersfoortSphere);
			var sinHalfPsi = Math.Sqrt(sinHalfPsiSquared);
			var cosHalfPsi = Math.Sqrt(1 - sinHalfPsiSquared);
			var tanHalfPsi = sinHalfPsi / cosHalfPsi;
			var sinPsi = 2 * sinHalfPsi * cosHalfPsi;
			var cosPsi = 1 - 2 * sinHalfPsiSquared;
			var sinAlpha = DegSin(deltaLambdaSphere) * (DegCos(phiSphere) / sinPsi);
			var cosAlpha = (DegSin(phiSphere) - DegSin(phiAmersfoortSphere) * cosPsi) / (DegCos(phiAmersfoortSphere) * sinPsi);
			var r = 2 * ScaleRd * rSphere * tanHalfPsi;

			var xRd = r * sinAlpha + XAmersfoortRd;
			var yRd = r * cosAlpha + YAmersfoortRd;

			return new Cartesian(xRd, yRd);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: inv_rd_projection
		**    Description:   inverse stereographic double projection
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    x_rd           double      in     req     none
		**    y_rd           double      in     req     none
		**    phi            double      out    -       none
		**    lambda         double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    x_rd, rd_y  input RD coordinates
		**    phi         output Bessel latitude in degrees
		**    lambda      output Bessel longitude in degrees
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		internal static Geographic InvRdProjection(Cartesian input)
		{
			/*
			**--------------------------------------------------------------
			**    Source: G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
			**            G. Strang van Hees, "Globale en lokale geodetische systemen". Delft: Nederlandse Commissie voor Geodesie (NCG), 1997.
			**--------------------------------------------------------------
			*/

			/*
			**--------------------------------------------------------------
			**    Explanation of the meaning of constants:
			**        f                         flattening of the ellipsoid
			**        ee                        first eccentricity squared (e squared in some notations)
			**        e                         first eccentricity
			**        eea                       second eccentricity squared (e' squared in some notations)
			**
			**        phi_amersfoort_sphere     latitude of projection base point Amersfoort on sphere in degrees
			**
			**        r1                        first (North South) principal radius of curvature in Amersfoort (M in some notations)
			**        r2                        second (East West) principal radius of curvature in Amersfoort (N in some notations)
			**        r_sphere                  radius of sphere
			**
			**        n                         constant of Gaussian projection n = 1.000475...
			**        q_amersfoort              isometric latitude of Amersfoort on ellipsiod
			**        w_amersfoort              isometric latitude of Amersfoort on sphere
			**        m                         constant of Gaussian projection m = 0.003773... (also named c in some notations)
			**--------------------------------------------------------------
			*/
            var f = 1 / InvFBessel;
            var ee = f * (2 - f);
            var e = Math.Sqrt(ee);
            var eea = ee / (1.0 - ee);
            var phiAmersfoortSphere = DegAtan(DegTan(PhiAmersfoortBessel) / Math.Sqrt(1 + eea * Math.Pow(DegCos(PhiAmersfoortBessel),2)));
            var r1 = ABessel * (1 - ee) / Math.Pow(Math.Sqrt(1 - ee * Math.Pow(DegSin(PhiAmersfoortBessel),2)),3);
			var r2 = ABessel / Math.Sqrt(1.0 - ee * Math.Pow(DegSin(PhiAmersfoortBessel),2));
			var rSphere = Math.Sqrt(r1 * r2);
            var n = Math.Sqrt(1 + eea * Math.Pow(DegCos(PhiAmersfoortBessel),4));
            var qAmersfoort = Atanh(DegSin(PhiAmersfoortBessel)) - e * Atanh(e * DegSin(PhiAmersfoortBessel));
            var wAmersfoort = Math.Log(DegTan(45 + 0.5 * phiAmersfoortSphere));
            var m = wAmersfoort - n * qAmersfoort;

			/*
			**--------------------------------------------------------------
			**    Explanation of the meaning of variables:
			**        r                    distance from Amersfoort in projection plane
			**        alpha                azimuth from Amersfoort
			**        psi                  distance angle from Amersfoort on sphere in degrees
			**        phi_sphere           latitide on sphere in degrees
			**        delta_lambda_sphere  difference in longitude on sphere with Amersfoort in degrees
			**        w                    isometric latitude on sphere
			**        q                    isometric latitude on ellipsiod
			**--------------------------------------------------------------
			*/
			var r = Math.Sqrt(Math.Pow(input.X - XAmersfoortRd,2) + Math.Pow(input.Y - YAmersfoortRd,2));
			var sinAlpha = (input.X - XAmersfoortRd) / r;
			if (r < Precision)
			{
				sinAlpha = 0;
			}
			var cosAlpha = (input.Y - YAmersfoortRd) / r;
			if (r < Precision)
			{
				cosAlpha = 1;
			}
			var psi = 2 * DegAtan(r / (2 * ScaleRd * rSphere));
			var phiSphere = DegAsin(cosAlpha * DegCos(phiAmersfoortSphere) * DegSin(psi) + DegSin(phiAmersfoortSphere) * DegCos(psi));
			var deltaLambdaSphere = DegAsin((sinAlpha * DegSin(psi)) / DegCos(phiSphere));
            var lambda = deltaLambdaSphere / n + LambdaAmersfoortBessel;
            var w = Atanh(DegSin(phiSphere));
			var q = (w - m) / n;

			/*
			**--------------------------------------------------------------
			**    Iterative calculation of phi
			**--------------------------------------------------------------
			*/
			double phi = 0;
            double diff = 90;
			while (diff > DegPrecision)
			{
				var previous = phi;
				phi = 2 * DegAtan(Math.Exp(q + 0.5 * e * Math.Log((1 + e * DegSin(phi)) / (1 - e * DegSin(phi))))) - 90;
				diff = Math.Abs(phi - previous);
			}

			return new Geographic(phi, lambda);
		}

		/// <summary>
		/// <para>read_double.</para>
		/// </summary>
		/// <param name="stream"> a <seealso cref="System.IO.Stream"/> object. </param>
		/// <returns> a double. </returns>
	
        public static double ReadDouble(Stream stream)
        {
            var bytes = new byte[8];
			stream.Read(bytes, 0, bytes.Length);
            if (!BitConverter.IsLittleEndian)
            {
                Array.Reverse(bytes);
              
            }
            return BitConverter.ToDouble(bytes, 0);
        }

		/// <summary>
		/// <para>read_float.</para>
		/// </summary>
		/// <param name="bytes"> an array of byte. </param>
		/// <returns> a float. </returns>
		public static float ReadFloat(byte[] bytes)
		{
            if (!BitConverter.IsLittleEndian)
            {
                Array.Reverse(bytes);

            }
            return BitConverter.ToSingle(bytes, 0);
        }

		/// <summary>
		/// <para>read_float.</para>
		/// </summary>
		/// <param name="stream"> a <seealso cref="System.IO.Stream"/> object. </param>
		/// <returns> a float. </returns>
	
        public static float ReadFloat(Stream stream)
		{
			var bytes = new byte[4];
			stream.Read(bytes, 0, bytes.Length);
			return ReadFloat(bytes);
		}

		/// <summary>
		/// <para>read_short.</para>
		/// </summary>
		/// <param name="stream"> a <seealso cref="System.IO.Stream"/> object. </param>
		/// <returns> a short. </returns>
        public static short ReadShort(Stream stream)
		{
			var bytes = new byte[2];
			stream.Read(bytes, 0, bytes.Length);
            if (!BitConverter.IsLittleEndian)
            {
                Array.Reverse(bytes);

            }
            return BitConverter.ToInt16(bytes, 0);
        }

		/*
		**--------------------------------------------------------------
		**    Function name: rd_correction
		**    Description:   apply the modelled distortions in the RD coordinate system
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    x_pseudo_rd    double      in     req     none
		**    y_pseudo_rd    double      in     req     none
		**    x_rd           double      out    -       none
		**    y_rd           double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    x_pseudo_rd, y_pseudo_rd  input coordinates in undistorted pseudo RD
		**    x_rd, y_rd                output coordinates in real RD
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
        internal static Cartesian RdCorrection(Cartesian pseudo)
		{

			var dx = GrdFile.GridFileDx.grid_interpolation(pseudo.X, pseudo.Y);
			var dy = GrdFile.GridFileDy.grid_interpolation(pseudo.X, pseudo.Y);

			return new Cartesian(pseudo.X -(dx ?? 0), pseudo.Y - (dy ?? 0), pseudo.Z);
		}

		/*
		**--------------------------------------------------------------
		**    Function name: inv_rd_correction
		**    Description:   remove the modelled distortions in the RD coordinate system
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    x_rd           double      in     req     none
		**    y_rd           double      in     req     none
		**    x_pseudo_rd    double      out    -       none
		**    x_pseudo_rd    double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    x_rd, y_rd                input coordinates in real RD
		**    x_pseudo_rd, y_pseudo_rd  output coordinates in undistorted pseudo RD
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
        internal static Cartesian InvRdCorrection(Cartesian rd)
		{

			/*
			**--------------------------------------------------------------
			**    The grid values are formally in pseudo RD. For the interpolation below the RD values are used. The intoduced error is certainly smaller than 0.0001 m for the X2c.grd and Y2c.grd.
			**--------------------------------------------------------------
			*/
			var dx = GrdFile.GridFileDx.grid_interpolation(rd.X, rd.Y);
			var dy = GrdFile.GridFileDy.grid_interpolation(rd.X, rd.Y);
			return new Cartesian(rd.X + dx ?? 0, rd.Y + dy ??0, rd.Z);
		}

	}

}