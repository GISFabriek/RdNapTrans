using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;

namespace RdNaptrans.Value
{
    using static Constants;
    using static Helpers;

    /// <summary>
    /// <para>GrdFile class.</para>
    /// 
    /// @author raymond
    /// @version $Id: $Id
    /// </summary>
    public class GrdFile
	{

        /*
    **--------------------------------------------------------------
    **    Continuation of static data declarations
    **    Names of grd files
    **
    **    Grd files are binary grid files in the format of the program Surfer(R)
    **    The header contains information on the number of grid points, bounding box and extreme values.
    **
    **    RD-corrections in x and y
    **
    **          -8000 meters < RD Easting  (stepsize 1 km) < 301000 meters
    **         288000 meters < RD Northing (stepsize 1 km) < 630000 meters
    **
    **    Geoid model NLGEO2004
    **
    **        50.525   degrees < ETRS89 latitude  (stepsize 0.050000 degrees) < 53.675 degrees
    **         3.20833 degrees < ETRS89 longitude (stepsize 0.083333 degrees) <  7.45833 degrees
    **
    **        Alternative notation:
    **        50\u0248 31' 30" < ETRS89_latitude  (stepsize 0\u0248 3' 0") < 53\u0248 40' 30"
    **         3\u0248 12' 30" < ETRS89_longitude (stepsize 0\u0248 5' 0") <  7\u0248 27' 30"
    **
    **        The stepsizes correspond to about 5,5 km x 5,5 km in the Netherlands.
    **--------------------------------------------------------------
    */
        /// <summary>
        /// Constant <code>GridFileDx</code> </summary>
        public static readonly GrdFile GridFileDx = new GrdFile("RdNaptrans.resources.rdnaptrans.x2c.grd");
        /// <summary>
        /// Constant <code>GridFileDy</code> </summary>
        public static readonly GrdFile GridFileDy = new GrdFile("RdNaptrans.resources.rdnaptrans.y2c.grd");
        /// <summary>
        /// Constant <code>GridFileGeoid</code> </summary>
        public static readonly GrdFile GridFileGeoid = new GrdFile("RdNaptrans.resources.rdnaptrans.nlgeo04.grd");

		private readonly sbyte[] _grdInner;
		private readonly GrdFileHeader _header;

        /// <summary>
        /// <para>Constructor for GrdFile.</para>
        /// </summary>
        /// <param name="resourceName"> name of the grd resource</param>
        public GrdFile(string resourceName)
        {
            var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resourceName);

            _header = read_grd_file_header(stream);
            using (var ms = new MemoryStream())
            {
                if (stream != null)
                {
                    stream.CopyTo(ms);
                    _grdInner = (sbyte[])(object)ms.ToArray();
                }
            }
        }

		/*
		**--------------------------------------------------------------
		**    Function name: grid_interpolation
		**    Description:   grid interpolation using Overhauser splines
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    x              double      in     req     none
		**    y              double      in     req     none
		**    grd_file       string      in     req     none
		**    value          double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    x, y           coordinates of the point for which a interpolated value is desired
		**    grd_file       name of the grd file to be read
		**    record_value   output of the interpolated value
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		/// <summary>
		/// <para>grid_interpolation.</para>
		/// </summary>
		/// <param name="x"> a double. </param>
		/// <param name="y"> a double. </param>
		/// <returns> a double. </returns>
		/// <exception cref="IOException"> if any. </exception>
		public virtual double? grid_interpolation(double x, double y)
		{
			var recordNumber = new int[16];
			var recordValue = new float[16];
            var f = new double[4];
			var g = new double[4];
			var gfac = new double[16];

			/*
			**--------------------------------------------------------------
			**    Explanation of the meaning of variables:
			**    size_x     number of grid values in x direction (row)
			**    size_y     number of grid values in y direction (col)
			**    min_x      minimum of x
			**    max_x      maximum of x
			**    min_y      minimum of y
			**    max_y      maximum of x
			**    min_value  minimum value in grid (besides the error values)
			**    max_value  maximum value in grid (besides the error values)
			**--------------------------------------------------------------
			*/

			/*
			**--------------------------------------------------------------
			**    Check for location safely inside the bounding box of grid
			**--------------------------------------------------------------
			*/
			if (x <= _header.SafeMinX || x >= _header.SafeMaxX || y <= _header.SafeMinY || y >= _header.SafeMaxY)
			{
				return null;
			}

			/*
			**--------------------------------------------------------------
			**    The selected grid points are situated around point X like this:
			**
			**        12  13  14  15
			**
			**         8   9  10  11
			**               X
			**         4   5   6   7
			**
			**         0   1   2   3
			**
			**    ddx and ddy (in parts of the grid interval) are defined relative to grid point 9, respectively to the right and down.
			**--------------------------------------------------------------
			*/
			var ddx = (x - _header.MinX) / _header.StepSizeX - Math.Floor((x - _header.MinX) / _header.StepSizeX);
			var ddy = 1 - ((y - _header.MinY) / _header.StepSizeY - Math.Floor((y - _header.MinY) / _header.StepSizeY));

			/*
			**--------------------------------------------------------------
			**    Calculate the record numbers of the selected grid points
			**    The records are numbered from lower left corner to the uper right corner starting with 0:
			**
			**    size_x*(size_y-1) . . size_x*size_y-1
			**                   .                    .
			**                   .                    .
			**                   0 . . . . . . size_x-1
			**--------------------------------------------------------------
			*/
			recordNumber[5] = (int)((x - _header.MinX) / _header.StepSizeX + Math.Floor((y - _header.MinY) / _header.StepSizeY) * _header.SizeX);
			recordNumber[0] = recordNumber[5] - _header.SizeX - 1;
			recordNumber[1] = recordNumber[5] - _header.SizeX;
			recordNumber[2] = recordNumber[5] - _header.SizeX + 1;
			recordNumber[3] = recordNumber[5] - _header.SizeX + 2;
			recordNumber[4] = recordNumber[5] - 1;
			recordNumber[6] = recordNumber[5] + 1;
			recordNumber[7] = recordNumber[5] + 2;
			recordNumber[8] = recordNumber[5] + _header.SizeX - 1;
			recordNumber[9] = recordNumber[5] + _header.SizeX;
			recordNumber[10] = recordNumber[5] + _header.SizeX + 1;
			recordNumber[11] = recordNumber[5] + _header.SizeX + 2;
			recordNumber[12] = recordNumber[5] + 2 * _header.SizeX - 1;
			recordNumber[13] = recordNumber[5] + 2 * _header.SizeX;
			recordNumber[14] = recordNumber[5] + 2 * _header.SizeX + 1;
			recordNumber[15] = recordNumber[5] + 2 * _header.SizeX + 2;

			/*
			**--------------------------------------------------------------
			**    Read the record values of the selected grid point
			**    Outside the validity area the records have a very large value (circa 1.7e38).
			**--------------------------------------------------------------
			*/
            MemoizedReadGrdFileBody = new MemoizeAnonymousInnerClass(this);
			for (var i = 0; i < 16; i++)
			{
				recordValue[i] = MemoizedReadGrdFileBody.Get(recordNumber[i]);
                if (recordValue[i] > _header.MaxValue + Precision || recordValue[i] < _header.MinValue - Precision)
				{
					return null;
				}
			}

			/*
			**--------------------------------------------------------------
			**    Calculation of the multiplication factors
			**--------------------------------------------------------------
			*/
			f[0] = -0.5 * ddx + ddx * ddx - 0.5 * ddx * ddx * ddx;
			f[1] = 1.0 - 2.5 * ddx * ddx + 1.5 * ddx * ddx * ddx;
			f[2] = 0.5 * ddx + 2.0 * ddx * ddx - 1.5 * ddx * ddx * ddx;
			f[3] = -0.5 * ddx * ddx + 0.5 * ddx * ddx * ddx;
			g[0] = -0.5 * ddy + ddy * ddy - 0.5 * ddy * ddy * ddy;
			g[1] = 1.0 - 2.5 * ddy * ddy + 1.5 * ddy * ddy * ddy;
			g[2] = 0.5 * ddy + 2.0 * ddy * ddy - 1.5 * ddy * ddy * ddy;
			g[3] = -0.5 * ddy * ddy + 0.5 * ddy * ddy * ddy;

			gfac[12] = f[0] * g[0];
			gfac[8] = f[0] * g[1];
			gfac[4] = f[0] * g[2];
			gfac[0] = f[0] * g[3];
			gfac[13] = f[1] * g[0];
			gfac[9] = f[1] * g[1];
			gfac[5] = f[1] * g[2];
			gfac[1] = f[1] * g[3];
			gfac[14] = f[2] * g[0];
			gfac[10] = f[2] * g[1];
			gfac[6] = f[2] * g[2];
			gfac[2] = f[2] * g[3];
			gfac[15] = f[3] * g[0];
			gfac[11] = f[3] * g[1];
			gfac[7] = f[3] * g[2];
			gfac[3] = f[3] * g[3];

			/*
			**--------------------------------------------------------------
			**    Calculation of the interpolated value
			**    Applying the multiplication factors on the selected grid values
			**--------------------------------------------------------------
			*/
			var value = 0.0;
			for (var i = 0; i < 16; i = i + 1)
			{
				value += gfac[i] * recordValue[i];
			}

			return value;
		}

		/*
		**--------------------------------------------------------------
		**    Function name: read_grd_file_header
		**    Description:   reads the header of a grd file
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    filename       string      in     req     none
		**    size_x         short int   out    -       none
		**    size_y         short int   out    -       none
		**    min_x          double      out    -       none
		**    max_x          double      out    -       none
		**    min_y          double      out    -       none
		**    max_y          double      out    -       none
		**    min_value      double      out    -       none
		**    max_value      double      out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    filename   name of the to be read binary file
		**    size_x     number of grid values in x direction (row)
		**    size_y     number of grid values in y direction (col)
		**    min_x      minimum of x
		**    max_x      maximum of x
		**    min_y      minimum of y
		**    max_y      maximum of x
		**    min_value  minimum value in grid (besides the error values)
		**    max_value  maximum value in grid (besides the error values)
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		private static GrdFileHeader read_grd_file_header(Stream input)
		{

            var id = new byte[4];
			input.Read(id, 0, 4);
			var idString = System.Text.Encoding.UTF8.GetString(id, 0, 4);

            /*
			**--------------------------------------------------------------
			**    Checks
			**--------------------------------------------------------------
			*/

            if (!idString.Equals("DSBB"))
			{
				throw new FormatException("Not a valid grd file");
			}

			/*
			**--------------------------------------------------------------
			**    Read output parameters
			**--------------------------------------------------------------
			*/

			var sizeX = ReadShort(input);
			var sizeY = ReadShort(input);
			var minX = ReadDouble(input);
			var maxX = ReadDouble(input);
			var minY = ReadDouble(input);
			var maxY = ReadDouble(input);
			var minValue = ReadDouble(input);
			var maxValue = ReadDouble(input);

			return new GrdFileHeader(sizeX, sizeY, minX, maxX, minY, maxY, minValue, maxValue);
			
		}

		/*
		**--------------------------------------------------------------
		**    Function name: read_grd_file_body
		**    Description:   reads a value from a grd file
		**
		**    Parameter      Type        In/Out Req/Opt Default
		**    filename       string      in     req     none
		**    number         long int    in     req     none
		**    value          float       out    -       none
		**
		**    Additional explanation of the meaning of parameters
		**    filename       name of the grd file to be read
		**    record_number  number defining the position in the file
		**    record_value   output of the read value
		**
		**    Return value: (besides the standard return values)
		**    none
		**--------------------------------------------------------------
		*/
		private float read_grd_file_body(int recordNumber)
		{
			const int recordLength = 4;
			const int headerLength = 56;


            /*
			**--------------------------------------------------------------
			**    Read
			**    Grd files are binary grid files in the format of the program Surfer(R)
			**    The first "header_length" bytes are the header of the file
			**    The body of the file consists of records of "record_length" bytes
			**    The records have a "record_number", starting with 0,1,2,...
			**--------------------------------------------------------------
			*/

            var start = headerLength + recordNumber * recordLength;
            var record
                = new byte[recordLength];
            Array.Copy(_grdInner, start, record, 0, recordLength);
            return ReadFloat(record);
		}

        
		internal MemoizeBase<int, float> MemoizedReadGrdFileBody;

		private class MemoizeAnonymousInnerClass : MemoizeBase<int, float>
		{
            public MemoizeAnonymousInnerClass(GrdFile outerInstance) :base(outerInstance)
            {
            }
			public override float Calculate(int recordNumber)
            {
                return OuterInstance.read_grd_file_body(recordNumber);
			}
		}

		private class GrdFileHeader
		{

			/*
			**    Additional explanation of the meaning of parameters
			**    filename   name of the to be read binary file
			**    size_x     number of grid values in x direction (row)
			**    size_y     number of grid values in y direction (col)
			**    min_x      minimum of x
			**    max_x      maximum of x
			**    min_y      minimum of y
			**    max_y      maximum of x
			**    min_value  minimum value in grid (besides the error values)
			**    max_value  maximum value in grid (besides the error values)
			*/

            public GrdFileHeader(short sizeX, short sizeY, double minX, double maxX, double minY, double maxY, double minValue, double maxValue)
			{
				SizeX = sizeX;
				SizeY = sizeY;
				MinX = minX;
				MaxX = maxX;
				MinY = minY;
				MaxY = maxY;
				MinValue = minValue;
				MaxValue = maxValue;

				StepSizeX = (MaxX - MinX) / (SizeX - 1);
				StepSizeY = (MaxY - MinY) / (SizeY - 1);

				SafeMinX = MinX + StepSizeX;
				SafeMaxX = MaxX - StepSizeX;
				SafeMinY = MinY + StepSizeY;
				SafeMaxY = MaxY - StepSizeY;
			}

            public short SizeX { get; }

            private short SizeY { get; }

            public double MinX { get; }

            private double MaxX { get; }

            public double MinY { get; }

            private double MaxY { get; }

            public double MinValue { get; }

            public double MaxValue { get; }

            public double StepSizeX { get; }

            public double StepSizeY { get; }

            public double SafeMinX { get; }

            public double SafeMaxX { get; }

            public double SafeMinY { get; }

            public double SafeMaxY { get; }
        }

	}

    public abstract class MemoizeBase<TP, TV>
    {
        protected readonly GrdFile OuterInstance;
        protected IDictionary<TP, TV> Values = new Dictionary<TP, TV>();
        protected MemoizeBase(GrdFile outerInstance)
        {
            OuterInstance = outerInstance;
        }


        internal virtual TV Get(TP p)
        {
            if (!Values.ContainsKey(p))
            {
                Values[p] = Calculate(p);
            }
            return Values[p];
        }

        /// <summary>
        /// Will implement the calculations that are
        /// to be remembered thanks to this class
        /// (one calculation per distinct parameter) </summary>
        /// <param name="p"> Calc V for p </param>
        /// <returns> v based on p </returns>
        public abstract TV Calculate(TP p);
    }

}