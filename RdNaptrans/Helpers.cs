// ***********************************************************************
// Assembly         : RdNapTrans
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 06-19-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-02-2019
// ***********************************************************************
// C# PORT from https://github.com/PDOK/rdnaptrans-java
// ***********************************************************************
using System;
using System.IO;

namespace RdNapTrans
{
    using static Constants;

    /// <summary>
    /// Class Helpers.
    /// </summary>
    public class Helpers
	{
        /*
		**--------------------------------------------------------------
		**    Functions
		**--------------------------------------------------------------
		*/

        /*
		**--------------------------------------------------------------
		**    Function name: DegSin
		**    Description:   sine for angles cartesian degrees
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    angleInDegrees          double      cartesian     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    sin(angleInDegrees)
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Return the sine for angles in cartesian degrees.
        /// </summary>
        /// <param name="angleInDegrees">The angle in degrees.</param>
        /// <returns>System.Double.</returns>
        internal static double DegSin(double angleInDegrees)
		{
			return Math.Sin(angleInDegrees / 180.0 * Math.PI);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: DegCos
		**    Description:   cosine for angles cartesian degrees
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    angleInDegrees          double      cartesian     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    cos(angleInDegrees)
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Return the cosine for angles in cartesian degrees.
        /// </summary>
        /// <param name="angleInDegrees">The angle in degrees.</param>
        /// <returns>System.Double.</returns>
        internal static double DegCos(double angleInDegrees)
		{
			return Math.Cos(angleInDegrees / 180.0 * Math.PI);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: DegTan
		**    Description:   tangent for angles cartesian degrees
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    angleInDegrees          double      cartesian     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    tan(angleInDegrees)
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Return the tangent for angles in cartesian degrees.
        /// </summary>
        /// <param name="angleInDegrees">The angle in degrees.</param>
        /// <returns>System.Double.</returns>
        internal static double DegTan(double angleInDegrees)
		{
			return Math.Tan(angleInDegrees / 180.0 * Math.PI);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: DegAsin
		**    Description:   inverse sine for angles cartesian degrees
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    angleInDegrees              double      cartesian     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    asin(angleInDegrees)
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Return the inverse sine for angles in cartesian degrees.
        /// </summary>
        /// <param name="angleInDegrees">The angle in degrees.</param>
        /// <returns>System.Double.</returns>
        internal static double DegAsin(double angleInDegrees)
		{
			return (Math.Asin(angleInDegrees) * 180.0 / Math.PI);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: DegAtan
		**    Description:   inverse tangent for angles cartesian degrees
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    angleInDegrees              double cartesian     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    atan(angleInDegrees)
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Return the inverse tangent for angles in cartesian degrees.
        /// </summary>
        /// <param name="angleInDegrees">The angle in degrees.</param>
        /// <returns>System.Double.</returns>
        internal static double DegAtan(double angleInDegrees)
		{
			return (Math.Atan(angleInDegrees) * 180.0 / Math.PI);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: Atanh
		**    Description:   inverse hyperbolic tangent
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    angleInDegrees              double      cartesian     req     none
		**
		**    Additional explanation of the meaning of parameters
		**    none
		**
		**    Return value: (besides the standard return values)
		**    Atanh(angleInDegrees)
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Return the inverse hyperbolic tangent for angles in cartesian degrees.
        /// </summary>
        /// <param name="angleInDegrees">The angle in degrees.</param>
        /// <returns>System.Double.</returns>
        internal static double Atanh(double angleInDegrees)
		{
			return (0.5 * Math.Log((1.0 + angleInDegrees) / (1.0 - angleInDegrees)));
		}


        /*
		**--------------------------------------------------------------
		**    Function name: Geographic2Cartesian
		**    Description:   from geographic coordinates to cartesian coordinates
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    phi            double      cartesian     req     none
		**    lambda         double      cartesian     req     none
		**    h              double      cartesian     req     none
		**    angleInDegrees              double      cartesian     req     none
		**    inverseFlattening          double      cartesian     req     none
		**    x              double      out    -       none
		**    y              double      out    -       none
		**    z              double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    phi      latitude cartesian degrees
		**    lambda   longitude cartesian degrees
		**    h        ellipsoidal height
		**    angleInDegrees        half major axis of the ellisoid
		**    inverseFlattening    inverse flattening of the ellipsoid
		**    x, y, z  output of cartesian coordinates
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Converts from geographic coordinate to cartesian coordinate.
        /// </summary>
        /// <param name="geographic">The geographic coordinate.</param>
        /// <param name="angleInDegrees">The angle in degrees.</param>
        /// <param name="inverseFlattening">The inverse flattening.</param>
        /// <returns>Cartesian coordinate.</returns>
        internal static Cartesian Geographic2Cartesian(Geographic geographic, double angleInDegrees, double inverseFlattening)
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
			**        ee   first eccentricity squared (e squared cartesian some notations)
			**        n    second (East West) principal radius of curvature (N cartesian some notations)
			**--------------------------------------------------------------
			*/
			var f = 1.0 / inverseFlattening;
			var ee = f * (2.0 - f);
			var n = angleInDegrees / Math.Sqrt(1.0 - ee * Math.Pow(DegSin(geographic.Phi),2));

			var x = (n + geographic.H) * DegCos(geographic.Phi) * DegCos(geographic.Lambda);
			var y = (n + geographic.H) * DegCos(geographic.Phi) * DegSin(geographic.Lambda);
			var z = (n * (1.0 - ee) + geographic.H) * DegSin(geographic.Phi);

			return new Cartesian(x, y, z);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: Cartesian2Geographic
		**    Description:   from cartesian coordinates to geographic coordinates
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    x              double      cartesian     req     none
		**    y              double      cartesian     req     none
		**    z              double      cartesian     req     none
		**    angleInDegrees              double      cartesian     req     none
		**    inverseFlattening          double      cartesian     req     none
		**    phi            double      out    -       none
		**    lambda         double      out    -       none
		**    h              double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    x, y, z  input of cartesian coordinates
		**    angleInDegrees        half major axis of the ellisoid
		**    inverseFlattening    inverse flattening of the ellipsoid
		**    phi      output latitude cartesian degrees
		**    lambda   output longitude cartesian degrees
		**    h        output ellipsoidal height
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Converts from cartesian coordinate to geographic coordinate.
        /// </summary>
        /// <param name="cartesian">The cartesian coordinate.</param>
        /// <param name="angleInDegrees">The angle in degrees.</param>
        /// <param name="inverseFlattening">The inverse flattening.</param>
        /// <returns>Geographic coordinate.</returns>
        internal static Geographic Cartesian2Geographic(Cartesian cartesian, double angleInDegrees, double inverseFlattening)
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
			**        ee   first eccentricity squared (e squared cartesian some notations)
			**        rho  distance to minor axis
			**        n    second (East West) principal radius of curvature (N cartesian some notations)
			**--------------------------------------------------------------
			*/
			var f = 1.0 / inverseFlattening;
			var ee = f * (2.0 - f);
			var rho = Math.Sqrt(cartesian.X * cartesian.X + cartesian.Y * cartesian.Y);
			var n = 0D;

			/*
			**--------------------------------------------------------------
			**    Iterative calculation of phi
			**--------------------------------------------------------------
			*/
			var phi = 0D;
            var diff = 90.0;
			while (diff > DegPrecision)
			{
				var previous = phi;
				n = angleInDegrees / Math.Sqrt(1.0 - ee * Math.Pow(DegSin(phi),2));
				phi = DegAtan(cartesian.Z / rho + n * ee * DegSin(phi) / rho);
				diff = Math.Abs(phi - previous);
			}

			/*
			**--------------------------------------------------------------
			**     Calculation of lambda and h
			**--------------------------------------------------------------
			*/
			var lambda = DegAtan(cartesian.Y / cartesian.X);
			var h = rho * DegCos(phi) + cartesian.Z * DegSin(phi) - n * (1.0 - ee * Math.Pow(DegSin(phi),2));

			return new Geographic(phi, lambda, h);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: SimTrans
		**    Description:   3 dimensional similarity transformation (7 parameters) around another pivot point "a" than the origin
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    xIn           double      cartesian     req     none
		**    yIn           double      cartesian     req     none
		**    zIn           double      cartesian     req     none
		**    tx             double      cartesian     req     none
		**    ty             double      cartesian     req     none
		**    tz             double      cartesian     req     none
		**    alpha          double      cartesian     req     none
		**    beta           double      cartesian     req     none
		**    gamma          double      cartesian     req     none
		**    delta          double      cartesian     req     none
		**    xa             double      cartesian     req     none
		**    ya             double      cartesian     req     none
		**    za             double      cartesian     req     none
		**    xOut          double      out    -       none
		**    yOut          double      out    -       none
		**    zOut          double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    xIn, yIn, zIn     input coordinates
		**    tx                   translation cartesian direction of x axis
		**    ty                   translation cartesian direction of y axis
		**    tz                   translation cartesian direction of z axis
		**    alpha                rotation around x axis cartesian radials
		**    beta                 rotation around y axis cartesian radials
		**    gamma                rotation around z axis cartesian radials
		**    delta                scale parameter (scale = 1 + delta)
		**    xa, ya, za           coordinates of pivot point a (cartesian case of rotation around the center of the ellipsoid these parameters are zero)
		**    xOut, yOut, zOut  output coordinates
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Performs a 3 dimensional similarity transformation (7 parameters) around another pivot point "a" than the origin
        /// </summary>
        /// <param name="cartesian">The cartesian coordinate.</param>
        /// <param name="translate">The translation directions in x,y, and z directions.</param>
        /// <param name="alpha">The rotation around x axis in cartesian radials.</param>
        /// <param name="beta">The rotation around y axis in cartesian radials.</param>
        /// <param name="gamma">The rotation around z axis in cartesian radials.</param>
        /// <param name="delta">Scale parameter (scale = 1 + delta).</param>
        /// <param name="pivot">The coordinates of pivot point a (cartesian case of rotation around the center of the ellipsoid these parameters are zero).</param>
        /// <returns>Cartesian.</returns>
        internal static Cartesian SimTrans(Cartesian cartesian, Cartesian translate, double alpha, double beta, double gamma, double delta, Cartesian pivot)

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
			**    angleInDegrees b cartesian
			**    d e f
			**    geographic h i
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
			**    point2 = inputPoint - pivotPoint
			**--------------------------------------------------------------
			*/
			var x = cartesian.X - pivot.X;
			var y = cartesian.Y - pivot.Y;
			var z = cartesian.Z - pivot.Z;

			/*
			**--------------------------------------------------------------
			**    Calculate the elements of the output vector:
			**    outputPoint = scale * rotationMatrix * point2 + translationVector + pivotPoint
			**--------------------------------------------------------------
			*/
			var xOut = (1.0 + delta) * (a * x + b * y + c * z) + translate.X + pivot.X;
			var yOut = (1.0 + delta) * (d * x + e * y + f * z) + translate.Y + pivot.Y;
			var zOut = (1.0 + delta) * (g * x + h * y + i * z) + translate.Z + pivot.Z;

			return new Cartesian(xOut, yOut, zOut);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: RdProjection
		**    Description:   stereographic double projection
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    phi            double      cartesian     req     none
		**    lambda         double      cartesian     req     none
		**    xRd           double      out    -       none
		**    yRd           double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    phi         input Bessel latitude cartesian degrees
		**    lambda      input Bessel longitude cartesian degrees
		**    xRd, yRd output RD coordinates
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Projects Geographic coordinates (EPSG:4258) to RD_New coordinates(EPSG:28992).
        /// </summary>
        /// <param name="geographic">The geographic coordinates container.</param>
        /// <returns>Cartesian object containing the converted coordinates.</returns>
        internal static Cartesian RdProjection(Geographic geographic)
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
			**        ee                        first eccentricity squared (e squared cartesian some notations)
			**        e                         first eccentricity
			**        eea                       second eccentricity squared (e' squared cartesian some notations)
			**
			**        phiAmersfoortSphere     latitude of projection base point Amersfoort on sphere cartesian degrees
			**        lambdaAmersfoortSphere  longitude of projection base point Amersfoort on sphere cartesian degrees
			**
			**        r1                        first (North South) principal radius of curvature cartesian Amersfoort (M cartesian some notations)
			**        r2                        second (East West) principal radius of curvature cartesian Amersfoort (N cartesian some notations)
			**        rSphere                  radius of sphere
			**
			**        n                         constant of Gaussian projection n = 1.000475...
			**        qAmersfoort              isometric latitude of Amersfoort on ellipsiod
			**        wAmersfoort              isometric latitude of Amersfoort on sphere
			**        m                         constant of Gaussian projection m = 0.003773... (also named cartesian cartesian some notations)
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
			**        q                    isometric latitude on ellipsoid
			**        w                    isometric latitude on sphere
			**        phiSphere           latitude on sphere cartesian degrees
			**        deltaLambdaSphere  difference cartesian longitude on sphere with Amersfoort cartesian degrees
			**        psi                  distance angle from Amersfoort on sphere
			**        alpha                azimuth from Amersfoort
			**        r                    distance from Amersfoort cartesian projection plane
			**--------------------------------------------------------------
			*/
			var q = Atanh(DegSin(geographic.Phi)) - e * Atanh(e * DegSin(geographic.Phi));
			var w = n * q + m;
			var phiSphere = 2 * DegAtan(Math.Exp(w)) - 90;
			var deltaLambdaSphere = n * (geographic.Lambda - lambdaAmersfoortSphere);
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
		**    Function name: InverseRdProjection
		**    Description:   inverse stereographic double projection
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    xRd           double      cartesian     req     none
		**    yRd           double      cartesian     req     none
		**    phi            double      out    -       none
		**    lambda         double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    xRd, yRd  input RD coordinates
		**    phi         output Bessel latitude cartesian degrees
		**    lambda      output Bessel longitude cartesian degrees
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/

        /// <summary>
        /// Projects RD_New coordinates(EPSG:28992) to Geographic coordinates (EPSG:4258).
        /// </summary>
        /// <param name="cartesian">The Cartesian coordinates container.</param>
        /// <returns>Geographic object containing the converted coordinates.</returns>
        internal static Geographic InverseRdProjection(Cartesian cartesian)
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
			**        ee                        first eccentricity squared (e squared cartesian some notations)
			**        e                         first eccentricity
			**        eea                       second eccentricity squared (e' squared cartesian some notations)
			**
			**        phiAmersfoortSphere     latitude of projection base point Amersfoort on sphere cartesian degrees
			**
			**        r1                        first (North South) principal radius of curvature cartesian Amersfoort (M cartesian some notations)
			**        r2                        second (East West) principal radius of curvature cartesian Amersfoort (N cartesian some notations)
			**        rSphere                  radius of sphere
			**
			**        n                         constant of Gaussian projection n = 1.000475...
			**        qAmersfoort              isometric latitude of Amersfoort on ellipsiod
			**        wAmersfoort              isometric latitude of Amersfoort on sphere
			**        m                         constant of Gaussian projection m = 0.003773... (also named cartesian cartesian some notations)
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
			**        r                    distance from Amersfoort cartesian projection plane
			**        angleInDegrees                azimuth from Amersfoort
			**        psi                  distance angle from Amersfoort on sphere cartesian degrees
			**        phiSphere           latitide on sphere cartesian degrees
			**        deltaLambdaSphere  difference cartesian longitude on sphere with Amersfoort cartesian degrees
			**        w                    isometric latitude on sphere
			**        q                    isometric latitude on ellipsiod
			**--------------------------------------------------------------
			*/
			var r = Math.Sqrt(Math.Pow(cartesian.X - XAmersfoortRd,2) + Math.Pow(cartesian.Y - YAmersfoortRd,2));
			var sinAlpha = (cartesian.X - XAmersfoortRd) / r;
			if (r < Precision)
			{
				sinAlpha = 0;
			}
			var cosAlpha = (cartesian.Y - YAmersfoortRd) / r;
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
			var phi = 0.0;
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
        /// ReadDouble.
        /// </summary>
        /// <param name="stream"><seealso cref="System.IO.Stream" /> object.</param>
        /// <returns>angleInDegrees double.</returns>
        public static double ReadDouble(Stream stream)
		{
			byte[] bytes = new byte[8];
			stream.Read(bytes, 0, bytes.Length);
            var result = BitConverter.ToDouble(bytes, 0);
            return result;
        }

        /// <summary>
        /// ReadFloat.
        /// </summary>
        /// <param name="bytes">an array of sbyte.</param>
        /// <returns>angleInDegrees float.</returns>
        public static float ReadFloat(sbyte[] bytes)
		{
            var result = BitConverter.ToSingle((byte[])(object)bytes, 0);
            return result;
        }

        /// <summary>
        /// ReadShort.
        /// </summary>
        /// <param name="stream"><seealso cref="System.IO.Stream" /></param>
        /// <returns>angleInDegrees short.</returns>
        public static short ReadShort(Stream stream)
		{
			var bytes = new byte[2];
            stream.Read(bytes, 0, bytes.Length);
            var result = BitConverter.ToInt16(bytes, 0);
            return result;
        }

        /*
		**--------------------------------------------------------------
		**    Function name: RdCorrection
		**    Description:   apply the modeled distortions cartesian the RD coordinate system
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    xPseudoRd    double      cartesian     req     none
		**    yPseudoRd    double      cartesian     req     none
		**    xRd           double      out    -       none
		**    yRd           double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    xPseudoRd, yPseudoRd  input coordinates cartesian undistorted pseudo RD
		**    xRd, yRd                output coordinates cartesian real RD
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/

        /// <summary>
        /// Applies a correction to the RD_New coordinates, using correction grids.
        /// </summary>
        /// <param name="pseudo">The Rd_New coordinates.</param>
        /// <returns>Cartesian containing the corrected coordinates.</returns>
        internal static Cartesian RdCorrection(Cartesian pseudo)
		{

			var dx = GrdFile.GridFileDx.InterpolateGrid(pseudo.X, pseudo.Y);
			var dy = GrdFile.GridFileDy.InterpolateGrid(pseudo.X, pseudo.Y);

			return new Cartesian(pseudo.X - (dx ?? 0), pseudo.Y - (dy?? 0), pseudo.Z);
		}

        /*
		**--------------------------------------------------------------
		**    Function name: InverseRdCorrection
		**    Description:   remove the modeled distortions cartesian the RD coordinate system
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    xRd           double      cartesian     req     none
		**    yRd           double      cartesian     req     none
		**    xPseudoRd    double      out    -       none
		**    xPseudoRd    double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    xRd, yRd                input coordinates cartesian real RD
		**    xPseudoRd, yPseudoRd  output coordinates cartesian undistorted pseudo RD
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Inverts a correction to the RD_New coordinates, using correction grids.
        /// </summary>
        /// <param name="rd">The Rd_New coordinates.</param>
        /// <returns>Cartesian containing the inversely corrected coordinates.</returns>
        internal static Cartesian InverseRdCorrection(Cartesian rd)
		{

			/*
			**--------------------------------------------------------------
			**    The grid values are formally cartesian pseudo RD. For the interpolation below the RD values are used. The introduced error is certainly smaller than 0.0001 m for the X2c.grd and Y2c.grd.
			**--------------------------------------------------------------
			*/
			var dx = GrdFile.GridFileDx.InterpolateGrid(rd.X, rd.Y);
			var dy = GrdFile.GridFileDy.InterpolateGrid(rd.X, rd.Y);
			return new Cartesian(rd.X + (dx ?? 0), rd.Y + (dy ??0), rd.Z);
		}

	}

}