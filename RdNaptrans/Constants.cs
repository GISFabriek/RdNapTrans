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
namespace RdNapTrans
{

    /// <summary>
    /// Class Constants.
    /// </summary>
    public static class Constants
	{
        /*
		**--------------------------------------------------------------
		**    Static data declarations
		**
		**    Geographic NL-Bessel coordinates of Amersfoort (pivot point and projection base point)
		**        phi     latitude in decimal degrees
		**        lambda  longitude in decimal degrees
		**        h       ellipsoidal height in meters
		**    Source of constants:
		**        Hk.J. Heuvelink, "De stereografische kaartprojectie in hare toepassing bij de Rijksdriehoeksmeting". Delft: Rijkscommissie voor Graadmeting en Waterpassing, 1918.
		**        HTW, "Handleiding voor de Technische Werkzaamheden van het Kadaster". Apeldoorn: Kadaster, 1996.
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Constant <code>PhiAmersfoortBessel=52.0+ 9.0/60.0+22.178/3600.0</code>
        /// </summary>
        public const double PhiAmersfoortBessel = 52.0 + 9.0 / 60.0 + 22.178 / 3600.0;
        /// <summary>
        /// Constant <code>LambdaAmersfoortBessel=5.0+23.0/60.0+15.500/3600.0</code>
        /// </summary>
        public const double LambdaAmersfoortBessel = 5.0 + 23.0 / 60.0 + 15.500 / 3600.0;
        /// <summary>
        /// Constant <code>HAmersfoortBessel=0.0</code>
        /// </summary>
        public const double HAmersfoortBessel = 0.0;
        /*
		**--------------------------------------------------------------
		**    Continuation of static data declarations
		**    Parameters of ellipsoids Bessel1841 and GRS80
		**        a      half major axis in meters
		**        inv_f  inverse flattening
		**    Source of constants: HTW, "Handleiding voor de Technische Werkzaamheden van het Kadaster". Apeldoorn: Kadaster, 1996.
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Constant <code>ABessel=6377397.155</code>
        /// </summary>
        public const double ABessel = 6377397.155;
        /// <summary>
        /// Constant <code>InvFBessel=299.1528128</code>
        /// </summary>
        public const double InvFBessel = 299.1528128;
        /// <summary>
        /// Constant <code>AEtrs=6378137</code>
        /// </summary>
        public const double AEtrs = 6378137;
        /// <summary>
        /// Constant <code>InvFEtrs=298.257222101</code>
        /// </summary>
        public const double InvFEtrs = 298.257222101;
        /*
		**--------------------------------------------------------------
		**    Continuation of static data declarations
		**    Transformation parameters relative to pivot point Amersfoort. Note: Do NOT confuse with parameters for the center of the ellipsoid!
		**        tx     translation in direction of x axis in meters
		**        ty     translation in direction of y axis in meters
		**        tz     translation in direction of z axis in meters
		**        alpha  rotation around x axis in radials
		**        beta   rotation around y axis in radials
		**        gamma  rotation around z axis in radials
		**        delta  scale parameter (scale = 1 + delta)
		**    Source of constants: A. de Bruijne, J. van Buren, A. K\u0148sters and H. van der Marel, "De geodetische referentiestelsels van Nederland; Definitie en vastlegging van ETRS89, RD en NAP en hun onderlinge relatie". Delft: Nederlandse Commissie voor Geodesie (NCG), to be published in 2005.
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Constant <code>TxBesselEtrs=593.0248</code>
        /// </summary>
        public const double TxBesselEtrs = 593.0248;
        /// <summary>
        /// Constant <code>TyBesselEtrs=25.9984</code>
        /// </summary>
        public const double TyBesselEtrs = 25.9984;
        /// <summary>
        /// Constant <code>TzBesselEtrs=478.7459</code>
        /// </summary>
        public const double TzBesselEtrs = 478.7459;
        /// <summary>
        /// Constant <code>AlphaBesselEtrs=1.9342e-6</code>
        /// </summary>
        public const double AlphaBesselEtrs = 1.9342e-6;
        /// <summary>
        /// Constant <code>BetaBesselEtrs=-1.6677e-6</code>
        /// </summary>
        public const double BetaBesselEtrs = -1.6677e-6;
        /// <summary>
        /// Constant <code>GammaBesselEtrs=9.1019e-6</code>
        /// </summary>
        public const double GammaBesselEtrs = 9.1019e-6;
        /// <summary>
        /// Constant <code>DeltaBesselEtrs=4.0725e-6</code>
        /// </summary>
        public const double DeltaBesselEtrs = 4.0725e-6;

        /// <summary>
        /// Constant <code>TxEtrsBessel=-593.0248</code>
        /// </summary>
        public const double TxEtrsBessel = -593.0248;
        /// <summary>
        /// Constant <code>TyEtrsBessel=-25.9984</code>
        /// </summary>
        public const double TyEtrsBessel = -25.9984;
        /// <summary>
        /// Constant <code>TzEtrsBessel=-478.7459</code>
        /// </summary>
        public const double TzEtrsBessel = -478.7459;
        /// <summary>
        /// Constant <code>AlphaEtrsBessel=-1.9342e-6</code>
        /// </summary>
        public const double AlphaEtrsBessel = -1.9342e-6;
        /// <summary>
        /// Constant <code>BetaEtrsBessel=1.6677e-6</code>
        /// </summary>
        public const double BetaEtrsBessel = 1.6677e-6;
        /// <summary>
        /// Constant <code>GammaEtrsBessel=-9.1019e-6</code>
        /// </summary>
        public const double GammaEtrsBessel = -9.1019e-6;
        /// <summary>
        /// Constant <code>DELTA_ETRS_BESSEL=-4.0725e-6</code>
        /// </summary>
        public const double DeltaEtrsBessel = -4.0725e-6;
        /*
		**--------------------------------------------------------------
		**    Continuation of static data declarations
		**    Parameters of RD projection
		**        scale         scale factor (k in some notations)
		**                      this factor was first defined by Hk.J. Heuvelink as pow(10,-400e-7), nowadays we define it as exactly 0.9999079
		**        x_amersfoort  false Easting
		**        y_amersfoort  false Northing
		**    Source of constants:
		**        G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
		**        G. Strang van Hees, "Globale en lokale geodetische systemen". Delft: Nederlandse Commissie voor Geodesie (NCG), 1997.
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Constant <code>ScaleRd=0.9999079</code>
        /// </summary>
        public const double ScaleRd = 0.9999079;
        /// <summary>
        /// Constant <code>XAmersfoortRd=155000</code>
        /// </summary>
        public const double XAmersfoortRd = 155000;
        /// <summary>
        /// Constant <code>YAmersfoortRd=463000</code>
        /// </summary>
        public const double YAmersfoortRd = 463000;

        /*
		**--------------------------------------------------------------
		**    Continuation of static data declarations
		**    Precision parameters for iterations (respectively in meters and degrees)
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Constant <code>Precision=0.0001</code>
        /// </summary>
        public const double Precision = 0.0001;
        /// <summary>
        /// Constant <code>DegPrecision=PRECISION/40e6*360</code>
        /// </summary>
        public static readonly double DegPrecision = Precision / 40e6 * 360;
        /*
		**--------------------------------------------------------------
		**    Continuation of static data declarations
		**    Mean difference between NAP and ellipsoidal Bessel height. This is only used for getting from x, y in RD to phi, lambda in ETRS89.
		**--------------------------------------------------------------
		*/
        /// <summary>
        /// Constant <code>MeanGeoidHeightBessel=0.0</code>
        /// </summary>
        public const double MeanGeoidHeightBessel = 0.0;
	}

}