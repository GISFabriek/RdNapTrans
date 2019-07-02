// ***********************************************************************
// Assembly         : RdNapTrans
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-02-2019
//
// Last Modified By :  Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-02-2019
// ***********************************************************************
// C# PORT from https://github.com/PDOK/rdnaptrans-java
// ***********************************************************************
namespace RdNapTrans
{
    using static Constants;
    using static Helpers;

    /// <summary>
    /// Class Transformer.
    /// </summary>
    public class Transformer
    {
        /*
	**--------------------------------------------------------------
	**    RDNAPTRANS(TM)2008
	**
	**    Authors: Jochem Lesparre, Joop van Buren, Marc Crombaghs, Frank Dentz, Arnoud Pol, Sander Oude Elberink
	**             http://www.rdnap.nl
	**    Based on RDNAPTRANS2004
	**    Main changes:
	**    - 7 similarity transformation parameters
	**    - 0.0088 offset in the transformation between ellipsoidal height (h) and orthometric heights (NAP)
	**    - coordinates are computed also outside the validity regions of the grid files x2c.grd, y2c.grd and nlgeo04.grd
	**--------------------------------------------------------------
	*/

        /*
        **--------------------------------------------------------------
        **    Function name: Etrs2Rd
        **    Description:   convert ETRS89 coordinates to RD coordinates
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    phiEtrs       double      in     req     none
        **    lambdaEtrs    double      in     req     none
        **    hEtrs         double      in     req     none
        **    xRd           double      out    -       none
        **    yRd           double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    phiEtrs, lambdaEtrs, hEtrs  input ETRS89 coordinates
        **    xRd, yRd                     output RD coordinates
        **
        **    Return value: (besides the standard return values)
        **    none
        **--------------------------------------------------------------
        */
        /// <summary>
        /// Converts an ETRS89 (EPSG:4258) coordinate to an RD_New (EPSG:28992) coordinate.
        /// </summary>
        /// <param name="etrs">a <seealso cref="Geographic" /> object containing ETRS89 coordinates.</param>
        /// <returns>a <seealso cref="Cartesian" /> object containing RD coordinates.</returns>
        public static Cartesian Etrs2Rd(Geographic etrs)
        {
            /*
            **--------------------------------------------------------------
            **    Calculate the cartesian ETRS89 coordinates of the pivot point Amersfoort
            **--------------------------------------------------------------
            */
            var amersfoortBessel =
                Geographic2Cartesian(new Geographic(PhiAmersfoortBessel, LambdaAmersfoortBessel, HAmersfoortBessel),
                    ABessel, InvFBessel);
            var xAmersfoortEtrs = amersfoortBessel.X + TxBesselEtrs;
            var yAmersfoortEtrs = amersfoortBessel.Y + TyBesselEtrs;
            var zAmersfoortEtrs = amersfoortBessel.Z + TzBesselEtrs;

            /*
            **--------------------------------------------------------------
            **    Convert ETRS89 coordinates to RD coordinates
            **    (To convert from degrees, minutes and seconds use the function deg_min_sec2decimal() here)
            **--------------------------------------------------------------
            */
            var cartesianEtrs = Geographic2Cartesian(etrs, AEtrs, InvFEtrs);
            var cartesianBessel = SimTrans(cartesianEtrs, new Cartesian(TxEtrsBessel, TyEtrsBessel, TzEtrsBessel),
                AlphaEtrsBessel, BetaEtrsBessel, GammaEtrsBessel, DeltaEtrsBessel,
                new Cartesian(xAmersfoortEtrs, yAmersfoortEtrs, zAmersfoortEtrs));

            var geographicBessel = Cartesian2Geographic(cartesianBessel, ABessel, InvFBessel);

            var pseudoRd = RdProjection(geographicBessel);
            return RdCorrection(pseudoRd).WithZ(geographicBessel.H);
        }

        /*
        **--------------------------------------------------------------
        **    Function name: Rd2Etrs
        **    Description:   convert RD coordinates to ETRS89 coordinates
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    xRd           double      in     req     none
        **    yRd           double      in     req     none
        **    nap            double      in     req     none
        **    phiEtrs       double      out    -       none
        **    lambdaEtrs    double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    xRd, yRd, nap        input RD and NAP coordinates
        **    phiEtrs, lambdaEtrs  output ETRS89 coordinates
        **
        **    Return value: (besides the standard return values)
        **    none
        **--------------------------------------------------------------
        */
        /// <summary>
        /// Converts RD_New (EPSG:28992) coordinateS to ETRS89 (EPSG:4258) coordinateS.
        /// </summary>
        /// <param name="rd">a <seealso cref="Cartesian" /> object containing RD coordinates.</param>
        /// <returns>a <seealso cref="Geographic" /> object containing ETRS89 coordinates.</returns>
        public static Geographic Rd2Etrs(Cartesian rd)
        {
            /*
            **--------------------------------------------------------------
            **    Calculate the cartesian Bessel coordinates of the pivot point Amersfoort
            **--------------------------------------------------------------
            */
            var amersfoortBessel =
                Geographic2Cartesian(new Geographic(PhiAmersfoortBessel, LambdaAmersfoortBessel, HAmersfoortBessel),
                    ABessel, InvFBessel);

            /*
            **--------------------------------------------------------------
            **    Calculate approximated value of ellipsoidal Bessel height
            **    The error made by using a constant for de Bessel geoid height is max. circa 1 meter in the ellipsoidal height (for the NLGEO2004 geoid model). This intoduces an error in the phi, lambda position too, this error is nevertheless certainly smaller than 0.0001 m.
            **--------------------------------------------------------------
            */
            var hBessel = rd.Z + MeanGeoidHeightBessel;

            /*
            **--------------------------------------------------------------
            **    Convert RD coordinates to ETRS89 coordinates
            **--------------------------------------------------------------
            */
            var pseudoRd = InverseRdCorrection(rd);
            var etrsBessel = InverseRdProjection(pseudoRd);
            var cartesianBessel = Geographic2Cartesian(etrsBessel.WithH(hBessel), ABessel, InvFBessel);
            var cartesianEtrs = SimTrans(cartesianBessel, new Cartesian(TxBesselEtrs, TyBesselEtrs, TzBesselEtrs),
                AlphaBesselEtrs, BetaBesselEtrs, GammaBesselEtrs, DeltaBesselEtrs, amersfoortBessel);
            return Cartesian2Geographic(cartesianEtrs, AEtrs, InvFEtrs);
            /*
            **--------------------------------------------------------------
            **    To convert to degrees, minutes and seconds use the function decimal2deg_min_sec() here
            **--------------------------------------------------------------
            */
        }

        /*
        **--------------------------------------------------------------
        **    Function name: Etrs2Nap
        **    Description:   convert ellipsoidal ETRS89 height to NAP height
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    phi            double      in     req     none
        **    lambda         double      in     req     none
        **    h              double      in     req     none
        **    nap            double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    phi, lambda, h  input ETRS89 coordinates
        **    nap             output NAP height
        **
        **    Return value: (besides the standard return values) none
        **    on error (outside geoid grid) nap is not computed here
        **    instead in Etrs2RdNap nap=hBessel
        **--------------------------------------------------------------
        */
        /// <summary>
        /// Converts AN ellipsoidal ETRS89 (EPSG:28992) height to an NAP (EPSG:4258) height.
        /// </summary>
        /// <param name="etrs">a <seealso cref="Geographic" /> object containing ETRS89 coordinates.</param>
        /// <returns>a nullable double containing the NAP height (if all went well).</returns>

        public static double? Etrs2Nap(Geographic etrs)
        {
            /*
            **--------------------------------------------------------------
            **    Explanation of the meaning of variables:
            **        n  geoid height
            **    on error (outside geoid grid) nap is not computed
            **    instead in etrs2rdnap nap=hBessel
            **--------------------------------------------------------------
            */

            var n = GrdFile.GridFileGeoid.InterpolateGrid(etrs.Lambda, etrs.Phi);

            if (n.HasValue)
            {
                return etrs.H - n.Value + 0.0088;
            }

            return null;

        }

        /*
        **--------------------------------------------------------------
        **    Function name: Nap2Etrs
        **    Description:   convert NAP height to ellipsoidal ETRS89 height
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    phi            double      in     req     none
        **    lambda         double      in     req     none
        **    nap            double      in     req     none
        **    h              double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    phi, lambda  input ETRS89 position
        **    nap          input NAP height at position phi, lambda
        **    h            output ellipsoidal ETRS89 height
        **
        **    Return value: (besides the standard return values)
        **    none
        **    on error (outside geoid grid) h is not computed here
        **    instead in rdnap2etrs h=hEtrsSim (from similarity transformation)
        **--------------------------------------------------------------
        */
        /// <summary>
        ///Converts an NAP (EPSG:4258) height to ellipsoidal ETRS89 (EPSG:28992) height.
        /// </summary>
        /// <param name="phi">Latitude in degrees.</param>
        /// <param name="lambda">Longitude in degrees.</param>
        /// <param name="nap">NAP height at position phi, lambda.</param>
        /// <returns>a nullable double containing the ellipsoidal height (if all went well).</returns>
        public static double? Nap2Etrs(double phi, double lambda, double nap)
        {
            /*
            **--------------------------------------------------------------
            **    Explanation of the meaning of variables:
            **        n  geoid height
            **--------------------------------------------------------------
            */
            var n = GrdFile.GridFileGeoid.InterpolateGrid(lambda, phi);

            return nap + n - 0.0088;

        }

        /*
        **--------------------------------------------------------------
        **    Function name: Etrs2Rdnap
        **    Description:   convert ETRS89 coordinates to RD and NAP coordinates
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    phi            double      in     req     none
        **    lambda         double      in     req     none
        **    h              double      in     req     none
        **    xRd           double      out    -       none
        **    yRd           double      out    -       none
        **    nap            double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    phi, lambda, h   input ETRS89 coordinates
        **    xRd, yRd, nap  output RD and NAP coordinates
        **
        **    Return value: (besides the standard return values)
        **    none
        **--------------------------------------------------------------
        */
        /// <summary>
        /// Converts ETRS89 (EPSG:4258) coordinates to RD_New (EPSG:28992) coordinates (including a vertical NAP coordinate).
        /// </summary>
        /// <param name="etrs">a <seealso cref="Geographic" /> object containing ETRS89 coordinates.</param>
        /// <returns>a <seealso cref="Cartesian" /> object containing RD coordinates.</returns>
        public static Cartesian Etrs2Rdnap(Geographic etrs)
        {
            var rd = Etrs2Rd(etrs);
            var betterH = Etrs2Nap(etrs);
            if (betterH.HasValue)
            {
                return rd.WithZ(betterH.Value);
            }

            return rd;
        }

        /*
        **--------------------------------------------------------------
        **    Function name: Rdnap2Etrs
        **    Description:   convert RD and NAP coordinates to ETRS89 coordinates
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    xRd           double      in     req     none
        **    yRd           double      in     req     none
        **    nap            double      in     req     none
        **    phi            double      out    -       none
        **    lambda         double      out    -       none
        **    h              double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    xRd, yRd, nap  input RD and NAP coordinates
        **    phi, lambda, h   output ETRS89 coordinates
        **
        **    Return value: (besides the standard return values)
        **    none
        **--------------------------------------------------------------
        */
        /// <summary>
        /// Converts RD_New (EPSG:28992) coordinates to ETRS89 (EPSG:4258) coordinates (including a vertical h coordinate).
        /// </summary>
        /// <param name="rdnap">a <seealso cref="Cartesian" /> object containing RD coordinates.</param>
        /// <returns>a <seealso cref="Geographic" /> object containing ETRS89 coordinates.</returns>
        public static Geographic Rdnap2Etrs(Cartesian rdnap)
        {
            var etrs = Rd2Etrs(rdnap);
            var betterH = Nap2Etrs(etrs.Phi, etrs.Lambda, rdnap.Z);

            if (betterH.HasValue)
            {
                return etrs.WithH(betterH.Value);
            }

            return etrs;
        }
    }

}