using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace WindowsFormsApp10
{
    public class SystemSolution
    {
        private double _eta;
        private double _ro;
        private double _alpha;
        private double _mu0;
        private double _fi;
        private double _mu;
        private double _X10;
        private double _q;
        private double _alpha0;
        private int _score;
        private int _count;
        private double h; 
        private const double e_tol = 1e-4;
        private double _tau;

        public int Score { get { return _score; } set { this._score = value; } }
        public double Time { get { return _tau;} set { this._tau = value;} }
        public double DeltaTime { get { return h; } set { this.h = value; } }
        public int Count { get { return _count; } set { this._count = value; } }
        public double X0 { get { return _X10; } set { this._X10 = value; } }
        public SystemSolution(double eta, double ro, double alpha, double mu0, double fi, double q, double X10 , double mu, double alpha0)
        {
            _eta = eta;
            _ro = ro;
            _alpha = alpha;
            _mu0 = mu0;
            _fi = fi;
            _q = q;
            _X10 = X10;
            _mu = mu;
            _alpha0 = alpha0;
        }
        // 
        private double[] InitializeStep(double x, double y, double z1, double z , double _time, double _deltatime)
        {
            double[] k = new double[4];
            k[0] = _deltatime * getYfunction(z);
            k[1] = _deltatime * getXfunction(z1);
            k[2] = _deltatime * getZfunction(x, y, _time);
            k[3] = _deltatime * getZ1function(x, y, z1);
            return k;
        }
        public double[] RungeKuttaMerson(double x, double y, double z1, double z)
        {
            double[] Solution = new double[4] { y, x, z, z1 };
            double[] EstimateError = new double[4];
            double max = 0;
            double[] k1 = InitializeStep(x, y, z1, z, Time, DeltaTime / 3.0);
            double[] k2 = InitializeStep(x +  k1[1], y + k1[0], z1 + k1[3], z + k1[2], Time + (DeltaTime / 3.0), DeltaTime / 3.0);
            double[] k3 = InitializeStep(x + 0.5 * (k1[1] + k2[1]), y + 0.5 * (k1[0] + k2[0]), z1 + 0.5 * (k1[3] + k2[3]), z + 0.5 * (k1[2] + k2[2]), Time + (DeltaTime / 3.0), DeltaTime);
            double[] k4 = InitializeStep(x + 0.375 * (k1[1] + k3[1]), y + 0.375 * (k1[0] + k3[0]), z1 + 0.375 * (k1[3] + k3[3]), z + 0.375 * (k1[2] + k3[2]), Time + (DeltaTime / 2.0), (4.0 * DeltaTime) / 3.0);
            for (int i = 0; i < k1.Length; i++) k4[i] += k1[i];
            double[] k5 = InitializeStep(x + 1.5 * (-k3[1] + k4[1]), y + 1.5 * (-k3[0] + 4 * k4[0]), z1 + 1.5 * (-k3[3] + k4[3]), z + 1.5 * (-k3[2] + k4[2]), Time + DeltaTime, DeltaTime / 3.0);

            for (int i = 0; i < Solution.Length; i++)
            {
                Solution[i] += (1.0 / 2.0) * (k4[i] + k5[i]);
                EstimateError[i] = (1.0 / 10.0) * (2 * k4[i] - 3 * k3[i] - k5[i]);
            }
            double[] Tmp = new double[4] { y, x, z, z1 };
            max = EstimateError.Select(t => Math.Abs(t)).Prepend(0.0).Max();
            if(max > e_tol)
            {
                DeltaTime /= 2;
                return Tmp;
            }
            else 
            {
                Time += DeltaTime;
                if(max <= e_tol / 32.0)
                {
                    DeltaTime *= 2;
                }
                return Solution;
            }
        }

        public int getscore(double x, double y)
        {
            if ((Score == -1 || Score == 3) && y > x + _eta) return Score = -1;
            else if (Score == -1 && y <= x + _eta) return Score = 0;
            else if ((Score == 0 || Score == 1) && Math.Pow(_alpha0, 2) * (x - y + _eta) < _mu0 * _fi - _ro + _q) return Score = 1;
            else if ((Score == 1 || Score == 2) && Math.Pow(_alpha0, 2) * (x - y + _eta) >= _mu0 * _fi - _ro + _q) return Score = 2;
            else if (Score == 2 && Math.Pow(_alpha0, 2) * (x - y + _eta) < _mu0 * _fi - _ro + _q) return Score = 3;
            else return Score;
        }
        public double getZfunction(double x, double y, double tau)
        {
            if (Score == -1) return -_ro - (_mu * Math.Cos(tau));
            else if (Score == 1) return -_ro - (_mu * Math.Cos(tau)) - Math.Pow(_alpha0, 2) * (y - X0 - _eta);
            else if (Score == 2) return -_ro - (_mu * Math.Cos(tau)) - Math.Pow(_alpha0, 2) * (y - x - _eta);
            else if (Score == 4) return -_ro;
            else return 0;

        }
        public double getZ1function(double x, double y, double z1)
        {
            if (Score == 2) return -Math.Pow(_alpha,2) * (x - y + _eta) - _mu0 * _fi * Convert.ToDouble(Math.Sign(z1)) - _ro + _q;
            else return 0;
        }
        public double getYfunction(double z)
        {
            return z;
        }
        public double getXfunction(double z1)
        {
            if(Score == 2) return z1;
            else return 0;
        }
    }
    
}
