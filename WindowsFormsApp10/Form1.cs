﻿using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace WindowsFormsApp10
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }
        private void button1_Click(object sender, EventArgs e)
        {
            this.chart1.Series.Clear();
            double
            nu = double.Parse(textBox1.Text),
            ro = double.Parse(textBox2.Text),
            alpha = double.Parse(textBox3.Text),
            mu0 = double.Parse(textBox4.Text),
            fi = double.Parse(textBox5.Text),
            q = double.Parse(textBox6.Text),
            X10 = double.Parse(textBox7.Text),
            h = double.Parse(textBox8.Text),
            dy = double.Parse(textBox9.Text),
            dx = double.Parse(textBox10.Text),
            y = double.Parse(textBox11.Text),
            x = double.Parse(textBox12.Text),
            mu = double.Parse(textBox13.Text),
            alpha0 = double.Parse(textBox14.Text),
            tau = 0;
            double[] arrayofsolutions = new double[4] {y,x,dy,dx};
            string[] seriesname = new string[arrayofsolutions.Length];
            for(int k = 0; k < seriesname.Length; k++)
            {
                seriesname[k] = k.ToString();
                this.chart1.Series.Add(seriesname[k]);
                this.chart1.Series[seriesname[k]].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
                this.chart1.Series[seriesname[k]].BorderWidth = 2;
            }
            for(int k = 0; k < this.chart1.Series.Count; k++)
            {
                this.chart1.Series[k].Points.Clear();
            }
            for (int k = 0; k < this.chart1.Series.Count;k++)
            {
                this.chart1.Series[k].Points.AddXY(tau, arrayofsolutions[k]);
            }
            SystemSolution solutions = new SystemSolution(nu,ro,alpha,mu0,fi,q,X10,mu,alpha0); // y,x,z,z1
            for (tau = h; tau < 200; tau += h)
            {
                arrayofsolutions = solutions.RungeKutta(x, y, dy, tau, h, dx);
                x = arrayofsolutions[1];
                y = arrayofsolutions[0];
                dx = arrayofsolutions[3];
                dy = arrayofsolutions[2];
                this.chart1.Series[3].Points.AddXY(tau, dx);
                this.chart1.Series[2].Points.AddXY(tau, dy);
            }
            

        }


    }
}
