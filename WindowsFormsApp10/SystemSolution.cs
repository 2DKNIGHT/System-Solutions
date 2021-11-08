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

        public int Score { get { return _score; } set { this._score = value; } }

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
        public double[] RungeKutta(double x, double y, double z, double tau, double h, double z1)
        {
            double[] k1, k2, k3, k4;
            k1 = new double[4];
            k2 = new double[4];
            k3 = new double[4];
            k4 = new double[4];
            double
                _x = x,
                _y = y,
                _z = z,
                _z1 = z1,
                _tau = tau;
            double[] _RungeKutta = new double[4] { z, z1, y, x};
            getscore(x, y);
            if(Score == 1)
            {
                _x = _X10;
            }
            for (int i = 0; i < k1.Length; i++)
            {
                if (i < 3 && i != 0)
                {
                    _x = x + k4[i - 1] / 2;
                    _y = y + k3[i - 1] / 2;
                    _z = z + k1[i - 1] / 2;
                    _z1 = z1 + k2[i - 1] / 2;
                    _tau = tau + h / 2;
                }
                if (i == 3)
                {
                    _x = x + k4[i - 1];
                    _y = y + k3[i - 1];
                    _z = z + k1[i - 1];
                    _z1 = z1 + k2[i - 1];
                    _tau = tau + h;
                }
                k2[i] = h * getXfunction(_z1);
                k1[i] = h * getYfunction(_z);
                k4[i] = h * getZ1function(_x, _y, _z1);
                k3[i] = h * getZfunction(_x, _y, _tau);
            }
            _RungeKutta[0] += (k1[0] + k1[1] * 2 + k1[2] * 2 + k1[3]) / 6; // z
            _RungeKutta[1] += (k2[0] + k2[1] * 2 + k2[2] * 2 + k2[3]) / 6; // z1
            _RungeKutta[2] += (k3[0] + k3[1] * 2 + k3[2] * 2 + k3[3]) / 6; // y
            _RungeKutta[3] += (k4[0] + k4[1] * 2 + k4[2] * 2 + k4[3]) / 6; // x
            if (k3[0] == 0)
            {
                _RungeKutta[1] = z1;
                _RungeKutta[3] = x;
            }

            return _RungeKutta;
        }
        private void getscore(double x, double y)
        {
            Score = -1;
            if (y > x + _eta && (Score == -1 || Score == 0))
            {
                Score = 0;
            }
            if ((Score != 2 && Score != -1) && (y == x + _eta || (_alpha * _alpha) * (y - _X10 - _eta) < _mu0 * _fi - _ro + _q))
            {
                Score = 1;
            }
            if (Score != 0 && ((_alpha * _alpha) * (y - _X10 - _eta) >= _mu0 * _fi - _ro + _q))
            {
                Score = 2;
            }
            else Score = -1;
        }
        
        public double getZfunction(double x, double y, double tau)
        {
            if (Score  == 0)
            {
                return -_ro - (_mu * Math.Cos(tau));
            }
            if (Score == 1)
            {
                return -_ro - (_mu * Math.Cos(tau)) - (_alpha0 * _alpha0) * (y - _X10 - _eta);
            }
            if (Score == 2)
            {
                return -_ro - (_mu * Math.Cos(tau)) - (_alpha0 * _alpha0) * (y - x - _eta);
            }
            return 0;
        }
        public double getZ1function(double x, double y, double z1)
        {
            if (Score == 2)
            {
                return -(_alpha * _alpha) * (x - y + _eta) - _mu0 * _fi * Math.Sign(z1) - _ro + _q;
            }
            return 0;
        }
        public double getYfunction(double z)
        {
            return z;
        }
        public double getXfunction(double z1)
        {
            return z1;
        }
    }
    
}
