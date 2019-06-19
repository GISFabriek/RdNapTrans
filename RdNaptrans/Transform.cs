namespace RdNaptrans
{
    using static Constants;
    using static Helpers;
    using Cartesian = Value.Cartesian;
    using Geographic = Value.Geographic;
    using GrdFile = Value.GrdFile;

  
    /// <summary>
    /// <para>Transform class.</para>
    /// 
    /// @author raymond
    /// @version $Id: $Id
    /// </summary>
    public class Transform
    {
        /* JAVA PORT
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
        **    Function name: etrs2rd
        **    Description:   convert ETRS89 coordinates to RD coordinates
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    phi_etrs       double      in     req     none
        **    lambda_etrs    double      in     req     none
        **    h_etrs         double      in     req     none
        **    x_rd           double      out    -       none
        **    y_rd           double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    phi_etrs, lambda_etrs, h_etrs  input ETRS89 coordinates
        **    x_rd, y_rd                     output RD coordinates
        **
        **    Return value: (besides the standard return values)
        **    none
        **--------------------------------------------------------------
        */
        /// <summary>
        /// <para>etrs2rd.</para>
        /// </summary>
        /// <param name="etrs"> a <seealso cref="rdnaptrans.value.Geographic"/> object. </param>
        /// <returns> a <seealso cref="rdnaptrans.value.Cartesian"/> object. </returns>
              public static Cartesian Etrs2Rd(Geographic etrs)
        {
            /*
            **--------------------------------------------------------------
            **    Calculate the cartesian ETRS89 coordinates of the pivot point Amersfoort
            **--------------------------------------------------------------
            */
            var amersfoortBessel =
                Geographic2Cartesian(
                    new Geographic(PhiAmersfoortBessel, LambdaAmersfoortBessel, HAmersfoortBessel), ABessel,
                    InvFBessel);
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
            var cartesianBessel = SimTrans(cartesianEtrs,
                new Cartesian(TxEtrsBessel, TyEtrsBessel, TzEtrsBessel), AlphaEtrsBessel, BetaEtrsBessel,
                GammaEtrsBessel, DeltaEtrsBessel,
                new Cartesian(xAmersfoortEtrs, yAmersfoortEtrs, zAmersfoortEtrs));

            var geographicBessel = Cartesian2Geographic(cartesianBessel, ABessel, InvFBessel);

            var pseudoRd = RdProjection(geographicBessel);
            return RdCorrection(pseudoRd).WithZ(geographicBessel.H);
        }

        /*
        **--------------------------------------------------------------
        **    Function name: rd2etrs
        **    Description:   convert RD coordinates to ETRS89 coordinates
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    x_rd           double      in     req     none
        **    y_rd           double      in     req     none
        **    nap            double      in     req     none
        **    phi_etrs       double      out    -       none
        **    lambda_etrs    double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    x_rd, y_rd, nap        input RD and NAP coordinates
        **    phi_etrs, lambda_etrs  output ETRS89 coordinates
        **
        **    Return value: (besides the standard return values)
        **    none
        **--------------------------------------------------------------
        */
        /// <summary>
        /// <para>rd2etrs.</para>
        /// </summary>
        /// <param name="rd"> a <seealso cref="rdnaptrans.value.Cartesian"/> object. </param>
        /// <returns> a <seealso cref="rdnaptrans.value.Geographic"/> object. </returns>

        public static Geographic Rd2Etrs(Cartesian rd)
        {
            /*
            **--------------------------------------------------------------
            **    Calculate the cartesian Bessel coordinates of the pivot point Amersfoort
            **--------------------------------------------------------------
            */
            var amersfoortBessel =
                Geographic2Cartesian(
                    new Geographic(PhiAmersfoortBessel, LambdaAmersfoortBessel, HAmersfoortBessel), ABessel,
                    InvFBessel);

            /*
            **--------------------------------------------------------------
            **    Calculate appoximated value of ellipsoidal Bessel height
            **    The error made by using a constant for de Bessel geoid height is max. circa 1 meter in the ellipsoidal height (for the NLGEO2004 geoid model). This intoduces an error in the phi, lambda position too, this error is nevertheless certainly smaller than 0.0001 m.
            **--------------------------------------------------------------
            */
            var hBessel = rd.Z + MeanGeoidHeightBessel;

            /*
            **--------------------------------------------------------------
            **    Convert RD coordinates to ETRS89 coordinates
            **--------------------------------------------------------------
            */
            var pseudoRd = InvRdCorrection(rd);
            var etrsBessel = InvRdProjection(pseudoRd);
            var cartesianBessel = Geographic2Cartesian(etrsBessel.WithH(hBessel), ABessel, InvFBessel);
            var cartesianEtrs = SimTrans(cartesianBessel,
                new Cartesian(TxBesselEtrs, TyBesselEtrs, TzBesselEtrs), AlphaBesselEtrs, BetaBesselEtrs,
                GammaBesselEtrs, DeltaBesselEtrs, amersfoortBessel);
            return Cartesian2Geographic(cartesianEtrs, AEtrs, InvFEtrs);
            /*
            **--------------------------------------------------------------
            **    To convert to degrees, minutes and seconds use the function decimal2deg_min_sec() here
            **--------------------------------------------------------------
            */
        }

        /*
        **--------------------------------------------------------------
        **    Function name: etrs2nap
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
        **    on error (outside geoid grid) nap is not compted here
        **    instead in etrs2rdnap nap=h_bessel
        **--------------------------------------------------------------
        */
        /// <summary>
        /// <para>etrs2nap.</para>
        /// </summary>
        /// <param name="etrs"> a <seealso cref="rdnaptrans.value.Geographic"/> object. </param>
        /// <returns> a double. </returns>
        
        public static double? Etrs2Nap(Geographic etrs)
        {
            /*
            **--------------------------------------------------------------
            **    Explanation of the meaning of variables:
            **        n  geoid height
            **    on error (outside geoid grid) nap is not compted
            **    instead in etrs2rdnap nap=h_bessel
            **--------------------------------------------------------------
            */

            var n = GrdFile.GridFileGeoid.grid_interpolation(etrs.Lambda, etrs.Phi);

            if (n.HasValue)
            {
                return etrs.H - n.Value + 0.0088;
            }

            return null;

        }

        /*
        **--------------------------------------------------------------
        **    Function name: nap2etrs
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
        **    on error (outside geoid grid) h is not compted here
        **    instead in rdnap2etrs h=h_etrs_sim (from similarity transformation)
        **--------------------------------------------------------------
        */
        /// <summary>
        /// <para>nap2etrs.</para>
        /// </summary>
        /// <param name="phi"> a double. </param>
        /// <param name="lambda"> a double. </param>
        /// <param name="nap"> a double. </param>
        /// <returns> a double. </returns>
     

        public static double? Nap2Etrs(double phi, double lambda, double nap)
        {
            /*
            **--------------------------------------------------------------
            **    Explanation of the meaning of variables:
            **        n  geoid height
            **--------------------------------------------------------------
            */
            var n = GrdFile.GridFileGeoid.grid_interpolation(lambda, phi);

            return nap + n - 0.0088;
        }

        /*
        **--------------------------------------------------------------
        **    Function name: etrs2rdnap
        **    Description:   convert ETRS89 coordinates to RD and NAP coordinates
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    phi            double      in     req     none
        **    lambda         double      in     req     none
        **    h              double      in     req     none
        **    x_rd           double      out    -       none
        **    y_rd           double      out    -       none
        **    nap            double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    phi, lambda, h   input ETRS89 coordinates
        **    x_rd, y_rd, nap  output RD and NAP coordinates
        **
        **    Return value: (besides the standard return values)
        **    none
        **--------------------------------------------------------------
        */
        /// <summary>
        /// <para>etrs2rdnap.</para>
        /// </summary>
        /// <param name="etrs"> a <seealso cref="rdnaptrans.value.Geographic"/> object. </param>
        /// <returns> a <seealso cref="rdnaptrans.value.Cartesian"/> object. </returns>
      
        public static Cartesian Etrs2Rdnap(Geographic etrs)
        {
            var rd = Etrs2Rd(etrs);
            var betterH = Etrs2Nap(etrs);
            return betterH.HasValue ? rd.WithZ(betterH.Value) : rd;
        }

        /*
        **--------------------------------------------------------------
        **    Function name: rdnap2etrs
        **    Description:   convert RD and NAP coordinates to ETRS89 coordinates
        **
        **    Parameter      Type        In/Out Req/Opt Default
        **    x_rd           double      in     req     none
        **    y_rd           double      in     req     none
        **    nap            double      in     req     none
        **    phi            double      out    -       none
        **    lambda         double      out    -       none
        **    h              double      out    -       none
        **
        **    Additional explanation of the meaning of parameters
        **    x_rd, y_rd, nap  input RD and NAP coordinates
        **    phi, lambda, h   output ETRS89 coordinates
        **
        **    Return value: (besides the standard return values)
        **    none
        **--------------------------------------------------------------
        */
        /// <summary>
        /// <para>rdnap2etrs.</para>
        /// </summary>
        /// <param name="rdnap"> a <seealso cref="rdnaptrans.value.Cartesian"/> object. </param>
        /// <returns> a <seealso cref="rdnaptrans.value.Geographic"/> object. </returns>
       
        public static Geographic Rdnap2Etrs(Cartesian rdnap)
        {
            var etrs = Rd2Etrs(rdnap);
            var betterH = Nap2Etrs(etrs.Phi, etrs.Lambda, rdnap.Z);

            return betterH.HasValue ? etrs.WithH(betterH.Value) : etrs;
        }

    }

}