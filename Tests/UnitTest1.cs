using NUnit.Framework;
using WindowsFormsApp10;
namespace Tests
{
    public class Tests
    {
        [SetUp]
        public void Setup()
        {
        }

        [Test]
        public void Test1()
        {
            SystemSolution solution = new SystemSolution(1,2,3,4,5,6,7);
            double
                x = 0,
                y = 0,
                dx = 0,
                dy = 0,
                h = 0.1;
            double[] arrayofsolutions = new double[4];
            for (double tau = 0; tau < 200; tau += h)
            {
                arrayofsolutions = solution.RungeKutta(x, y, dy, tau, h, dx);
                for (int i = 0; i < arrayofsolutions.Length; i++)
                {
                    x = arrayofsolutions[i];
                }
                Assert.Pass();
            }
        }
    }
}