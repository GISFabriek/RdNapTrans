using RdNapTrans;
using Xunit;

namespace RdNaptransUnitTestProject
{
    using System;
 
    using System.Collections.Generic;

    public class ConversionTest
    {
        private readonly List<(string name, Geographic geographic, Cartesian cartesian)> _testItems;
        public static double MaxDeltaRd = 0.001;
        public static double MaxDeltaAngle = 0.00000001;
        public static double MaxDeltaH = 0.001;

        public ConversionTest()
        {
            _testItems = new List<(string name, Geographic geographic, Cartesian cartesian)>();
            var item = ("Texel", new Geographic(53.160753042, 4.824761912, 42.8614), new Cartesian(117380.1200, 575040.3400, 1.0000));
            _testItems.Add(item);
            item = ("Noord-Groningen", new Geographic(53.419482050, 6.776726674, 42.3586),
                new Cartesian(247380.5600, 604580.7800, 2.0000));
            _testItems.Add(item);
            item = ("Amersfoort", new Geographic(52.155172897, 5.387203657, 43.2551),
                new Cartesian(155000.0000, 463000.0000, 0.0000));
            _testItems.Add(item);
            item = ("Amersfoort 100m", new Geographic(52.155172910, 5.387203658, 143.2551),
                new Cartesian(155000.0000, 463000.0000, 100.0000));
            _testItems.Add(item);
            item = ("Zeeuws-Vlaanderen", new Geographic(51.368607152, 3.397588595, 47.4024),
                new Cartesian(16460.9100, 377380.2300, 3.0000));
            _testItems.Add(item);
            item = ("Zuid-Limburg", new Geographic(50.792584916, 5.773795548, 245.9478),
                new Cartesian(182260.4500, 311480.6700, 200.0000));
            _testItems.Add(item);
            item = ("Maasvlakte", new Geographic(51.947393898, 4.072887101, 47.5968),
                new Cartesian(64640.8900, 440700.0101, 4.0000));
            _testItems.Add(item);
            item = ("outside", new Geographic(48.843030210, 8.723260235, 52.0289),
                new Cartesian(400000.2300, 100000.4500, 5.0000));
            _testItems.Add(item);
            item = ("no_rd&geoid", new Geographic(50.687420392, 4.608971813, 51.6108),
                new Cartesian(100000.6700, 300000.8900, 6.0000));
            _testItems.Add(item);
            item = ("no_geoid", new Geographic(51.136825197, 4.601375361, 50.9672),
                new Cartesian(100000.6700, 350000.8900, 6.0000));
            _testItems.Add(item);
            item = ("no_rd", new Geographic(52.482440839, 4.268403889, 49.9436),
                new Cartesian(79000.0100, 500000.2300, 7.0000));
            _testItems.Add(item);
            item = ("edge_rd", new Geographic(51.003976532, 3.891247830, 52.7427),
                new Cartesian(50000.4500, 335999.6700, 8.0000));
            _testItems.Add(item);
            
        }

        public bool IsWithinRange(double first, double second, double tolerance)
        {
            var difference = Math.Abs(first - second);
            return difference <= tolerance;
        }

        [Fact]
        public void TestEtrs2RdNap()
        {
            foreach (var item in _testItems)
            {
                var result = Transformer.Etrs2Rdnap(item.geographic);
                Assert.True(IsWithinRange(result.X, item.cartesian.X, MaxDeltaRd));
                Assert.True(IsWithinRange(result.Y, item.cartesian.Y, MaxDeltaRd));
                Assert.True(IsWithinRange(result.Z, item.cartesian.Z, MaxDeltaH));
            }
        }

        [Fact]
        public void TestRdNap2Etrs()
        {
            foreach (var item in _testItems)
            {
                var result = Transformer.Rdnap2Etrs(item.cartesian);
                Assert.True(IsWithinRange(result.Lambda, item.geographic.Lambda, MaxDeltaAngle));
                Assert.True(IsWithinRange(result.Phi, item.geographic.Phi, MaxDeltaAngle));
                Assert.True(IsWithinRange(result.H, item.geographic.H, MaxDeltaH));
            }
        }
    }
}
