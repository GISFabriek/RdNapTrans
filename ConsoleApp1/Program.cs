using System;
using RdNapTrans;

namespace ConsoleApp1
{ 

    class Program
    {
        static void Main()
        {
            var input = new Geographic(53.160753042, 4.824761912, 42.8614); // Texel
            Console.WriteLine($"Input: " + input);
            var output = Transformer.Etrs2Rdnap(input);
            Console.WriteLine($"Output: " + output);
            Console.ReadKey();
        }
    }
}
