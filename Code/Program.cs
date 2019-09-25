using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CFDProject2
{
    class Program
    {
        //Reference Type for each Node
        struct Node
        {
            public double x, y, T, dx, dy, aE, aP, aN, aW, aS, Sc, Sp, b;

        }
        //General Dimensions:
        static double Width = 1.2f, Height = 0.9f;
        static double xS1 = 1f, yS1 = 0.7f, WS1 = 0.1f, HS1 = 0.1f, xS2 = 0.2f, yS2 = 0.15f, WS2 = 0.3f, HS2 = 0.6f, xS3 = 1.05f, yS3 = 0.2f, WS3 = 0.15f, HS3 = 0.3f;
        //Grid Size Scale ==> for Grid Study
        static double GridScale = 1f;
        //Total Number of Nodes
        static int Ni = 120, Nj = 90;
        //Conductivity
        static double K = 2;
        //h,Emissivity,Sigma
        static double h = 60, e = 0.6f, sigma = 5.67E-8f;
        //Temperatures
        static double T1 = 400, T2 = 500, Tinf = 300, InitAssumpTemp = 350, InitTwall = 350;
        static Node[,] Nodes;
        static double[,] OldTemps;
        static double[] OldTwall;
        static double  Res = 100, Epsilon = 0.00000001f;
        //Relaxation Coefficient
        static double W = 1f;
       // static double QD, QR, QS;
        static int Counter = 0;
        static double S = 4800;
        static double dx, dy;
        static void Main(string[] args)
        {
            Grid();
            AssignTemps();
            AssignSource();
            AssignCoeffs();
            SOLVE();

            //North Side
            for (int i = 1; i < Ni + 1; i++)
            {
                Nodes[i, Nj + 1].T = Nodes[i, Nj].T;
            }

            //Left Side
            for (int j = 0; j < Nj + 2; j++)
            {
                Nodes[0, j].T = Nodes[1, j].T;
            }
            Export();
            CalculateQ();
        }

        static void Grid()
        {
            //Apply The Scale
            Ni = (int)(Ni / GridScale);
            Nj = (int)(Nj / GridScale);

            dx = Width / Ni;
            dy = Height / Nj;
            //Initialization
            Nodes = new Node[Ni + 2, Nj + 2];
            OldTemps = new double[Ni + 2, Nj + 2];
            OldTwall = new double[Ni + 2];

            //Assign dx and dy to Every Node
            for (int i = 0; i < Ni + 2; i++)
            {
                for (int j = 0; j < Nj + 2; j++)
                {
                    Nodes[i, j].dx = dx;
                    Nodes[i, j].dy = dy;
                }
            }

            //Assign x and y Position to each Node
            for (int i = 0; i < Ni + 2; i++)
            {

                for (int j = 0; j < Nj + 2; j++)
                {
                    //X Position
                    //set 0 to the Boundary Node
                    if (i == 0) Nodes[i, j].x = 0;
                    //The First Node is dx/2 away from the Boundary
                    else if (i == 1) Nodes[i, j].x = Nodes[i, j].dx / 2;
                    //Other Nodes Distance dx from each other
                    else if (i == Ni + 1) Nodes[i, j].x = Nodes[i - 1, j].x + Nodes[i, j].dx / 2;
                    else Nodes[i, j].x = Nodes[i - 1, j].x + Nodes[i, j].dx;

                    //Same for Y Position
                    if (j == 0) Nodes[i, j].y = 0;
                    else if (j == 1) Nodes[i, j].y = Nodes[i, j].dy / 2;
                    else if (j == Nj + 1) Nodes[i, j].y = Nodes[i, j - 1].y + Nodes[i, j].dy / 2;
                    else Nodes[i, j].y = Nodes[i, j - 1].y + Nodes[i, j].dy;
                }
            }

        }



        static void AssignCoeffs()
        {
            for (int i = 1; i < Ni + 1; i++)
            {
                for (int j = 1; j < Nj + 1; j++)
                {
                    //aE
                    Nodes[i, j].aE = K / (Nodes[i + 1, j].x - Nodes[i, j].x) * Nodes[i, j].dy;
                    //aW
                    if (i != 1)
                        Nodes[i, j].aW = K / (Nodes[i, j].x - Nodes[i - 1, j].x) * Nodes[i, j].dy;
                    // at Western Node ==> Assign Insulate Boundary Condition
                    else
                        Nodes[i, j].aW = 0;
                    //aN
                    if (j != Nj)
                        Nodes[i, j].aN = K / (Nodes[i, j + 1].y - Nodes[i, j].y) * Nodes[i, j].dx;
                    // at Northen Node ==> Assign Insulate Boundary Condition
                    else
                        Nodes[i, j].aN = 0;

                    //aS,aP,b
                    Nodes[i, j].aS = K / (Nodes[i, j].y - Nodes[i, j - 1].y) * Nodes[i, j].dx;
                    Nodes[i, j].aP = Nodes[i, j].aN + Nodes[i, j].aS + Nodes[i, j].aW + Nodes[i, j].aE - Nodes[i, j].Sp * Nodes[i, j].dx * Nodes[i, j].dy;
                    Nodes[i, j].b = Nodes[i, j].Sc * Nodes[i, j].dx * Nodes[i, j].dy;
                }
            }
        }

        static void AssignTemps()
        {
            for (int i = 0; i < Ni + 2; i++)
            {
                for (int j = 0; j < Nj + 2; j++)
                {
                    if (i == Ni + 1) Nodes[i, j].T = T2;
                    else if (j == 0) Nodes[i, 0].T = InitTwall;
                    else
                        Nodes[i, j].T = InitAssumpTemp;
                }

            }
        }

        static void AssignSource()
        {
            for (int i = 0; i < Ni + 2; i++)
            {
                for (int j = 0; j < Nj + 2; j++)
                {
                    //Virtual Sources for Fixing The Hole Temperature at T1 and T2
                    if (Nodes[i, j].x >= xS2 && Nodes[i, j].x <= xS2 + WS2 && Nodes[i, j].y >= yS2 && Nodes[i, j].y <= yS2 + HS2)
                    {
                        Nodes[i, j].Sc = T1 * (float)Math.Pow(10, 30);
                        Nodes[i, j].Sp = -(float)Math.Pow(10, 30);
                    }
                    else if (Nodes[i, j].x >= xS3 && Nodes[i, j].x <= xS3 + WS3 && Nodes[i, j].y >= yS3 && Nodes[i, j].y <= yS3 + HS3)
                    {
                        Nodes[i, j].Sc = T2 * (float)Math.Pow(10, 30);
                        Nodes[i, j].Sp = -(float)Math.Pow(10, 30);
                    }
                    else if (Nodes[i, j].x >= xS1 && Nodes[i, j].x <= xS1 + WS1 && Nodes[i, j].y >= yS1 && Nodes[i, j].y <= yS1 + HS1)
                    {
                        Nodes[i, j].Sc = S;
                        Nodes[i, j].Sp = 0;
                    }
                    else
                    {
                        Nodes[i, j].Sc = 0;
                        Nodes[i, j].Sp = 0;
                    }
                }
            }
        }


        static void SOLVE()
        {
            //TDMA Variables Initialization
            double[] a, b, c, d;
            a = new double[Ni + 1];
            b = new double[Ni + 1];
            c = new double[Ni + 1];
            d = new double[Ni + 1];
            //Count the Number of Iteration
            Counter = 0;
            while (Res > (Epsilon))
            {
                Counter++;
                //Line By Line Loop
                for (int j = 1; j < Nj + 1; j++)
                {
                    //Set the Proper Coefficient for Each Node Equation in order to Solve by TDMA Algorithm
                    for (int i = 1; i < Ni + 1; i++)
                    {
                        if (i == Ni) a[i] = 0;
                        else
                            a[i] = (-1) * Nodes[i, j].aE;
                        if (i == 1) b[i] = 0;
                        else
                            b[i] = (-1) * Nodes[i, j].aW;
                        d[i] = Nodes[i, j].aP;
                        if (i == 1) c[i] = Nodes[i, j].aN * Nodes[i, j + 1].T + Nodes[i, j].aS * Nodes[i, j - 1].T + Nodes[i, j].aW * Nodes[i - 1, j].T + Nodes[i, j].b;
                        else if (i == Ni) c[i] = Nodes[i, j].aN * Nodes[i, j + 1].T + Nodes[i, j].aS * Nodes[i, j - 1].T + Nodes[i, j].aE * Nodes[i + 1, j].T + Nodes[i, j].b;
                        else c[i] = Nodes[i, j].aN * Nodes[i, j + 1].T + Nodes[i, j].aS * Nodes[i, j - 1].T + Nodes[i, j].b;
                    }
                    //TDMA First Loop
                    for (int k = 2; k < Ni + 1; k++)
                    {
                        d[k] -= b[k] * a[k - 1] / d[k - 1];
                        c[k] -= b[k] * c[k - 1] / d[k - 1];
                    }
                    Nodes[Ni, j].T = c[Ni] / d[Ni];
                    //TDMA Second Loop
                    for (int k = Ni - 1; k > 0; k--)
                    {
                        Nodes[k, j].T = (c[k] - a[k] * Nodes[k + 1, j].T) / d[k];
                    }

                }
                //South Wall Temperature is Unknown,Some Iteration is Required to Correct the Wall Temperature
                for (int k = 0; k < 7; k++)
                {
                    for (int i = 1; i < Ni + 1; i++)
                    {
                        OldTwall[i] = Nodes[i, 0].T;
                        Nodes[i, 0].T = ((K * Nodes[i, 1].T / (Nodes[i, 1].y - Nodes[i, 0].y)) + h * Tinf + e * sigma * (float)Math.Pow(Tinf, 4.0) + 3 * e * sigma * (float)Math.Pow(Nodes[i, 0].T, 4)) / ((K / (Nodes[i, 1].y - Nodes[i, 0].y)) + h + 4 * e * sigma * (float)Math.Pow(Nodes[i, 0].T, 3));
                    }
                }

                CalculateResidual();
                Console.WriteLine(Res);
                for (int i = 1; i < Ni + 1; i++)
                {
                    for (int j = 1; j < Nj + 1; j++)
                    {
                        //Apply Relaxation
                        Nodes[i, j].T += (W - 1) * (Nodes[i, j].T - OldTemps[i, j]);
                        //Save The Old Temps for Residual Calculation
                        OldTemps[i, j] = Nodes[i, j].T;

                    }
                }
            }


        }

        static void CalculateResidual()
        {
            Res = 0;
            for (int i = 1; i < Ni + 1; i++)
            {
                for (int j = 1; j < Nj + 1; j++)
                {
                    Res += Math.Abs(Nodes[i, j].T - OldTemps[i, j]);
                }
            }
            Res /= Ni * Nj;
        }
        
        //Export Tecplot Format Data
        static void Export()
        {
            string[] Data = new string[(Ni + 2) * (Nj + 2) + 2];
            Data[0] = "VARIABLES = X, Y, T";
            Data[1] = string.Format("ZONE I = {0} , J = {1}", Nj + 2, Ni + 2);
            for (int i = 0; i < Ni + 2; i++)
            {
                for (int j = 0; j < Nj + 2; j++)
                {
                    Data[i * (Nj + 2) + j + 2] = Nodes[i, j].x + "\t" + Nodes[i, j].y + "\t" + Nodes[i, j].T;
                }
            }
            File.WriteAllLines(@"F:\Result.plt", Data);

        }
            
        static void CalculateQ()
        {
            double QD = 0, QDL = 0, QDR = 0;

            //Index Range of Source Area & Holes
            //====================
            int iS21 = (int)(xS2 / dx + 1 + 0.005);
            int iS22 = (int)((xS2 + WS2) / dx + 0.005);
            int jS21 = (int)(yS2 / dy + 1 + 0.005);
            int jS22 = (int)((yS2 + HS2) / dy + 0.005);
            //======================
            int iS11 = (int)(xS1 / dx + 1 + 0.005);
            int iS12 = (int)((xS1 + WS1) / dx + 0.005);
            int jS11 = (int)(yS1 / dy + 1 + 0.005);
            int jS12 = (int)((yS1 + HS1) / dy + 0.005);
            //======================
            int iS31 = (int)(xS3 / dx + 1 + 0.005);
            int iS32 = (int)((xS3 + WS3) / dx + 0.005);
            int jS31 = (int)(yS3 / dy + 1 + 0.005);
            int jS32 = (int)((yS3 + HS3) / dy + 0.005);
            //======================

            //Down
            for (int i = 1; i < Ni+1; i++)
            {
                QD += (h * (Nodes[i, 0].T - Tinf) + (e * sigma * (float)(Math.Pow(Nodes[i, 0].T, 4) - Math.Pow(Tinf, 4)))) * dx;
            }
            //Corners
            QDL += (h * (Nodes[0, 0].T - Tinf) + (e * sigma * (float)(Math.Pow(Nodes[0, 0].T, 4) - Math.Pow(Tinf, 4)))) * (dx / 2);
            QDR += (h * (Nodes[Ni+1, 0].T - Tinf) + (e * sigma * (float)(Math.Pow(Nodes[Ni+1, 0].T, 4) - Math.Pow(Tinf, 4)))) * (dx / 2);

            //Right
            double QR = 0, QRD = 0, QRU = 0;
            for (int j = 1; j < Nj+1; j++)
            {
                if(!(j> jS31 && j< jS32))
                QR += (K * (Nodes[Ni, j].T - Nodes[Ni+1, j].T) / (dx / 2)) * dy;
            }
            //Corners
            QRD += (K * (Nodes[Ni, 0].T - Nodes[Ni+1, 0].T) / (dx / 2)) * (dy / 2);
            QRU += (K * (Nodes[Ni, Nj].T - Nodes[Ni+1, Nj].T) / (dy)) * (dx / 2);

            //S2
            double QS2R=0, QS2L=0, QS2U=0, QS2D=0;
            //S2 Up and Down
            for (int i = iS21+1; i < iS22; i++)
            {
                QS2D += (K * (Nodes[i,jS21-1].T - Nodes[i,jS21].T) / dy) * dx;
                QS2U += (K * (Nodes[i, jS22+1].T - Nodes[i, jS22].T) / dy) * dx;
            }
            //Left And Right Corners
            var QS2DL = (K * (Nodes[iS21, jS21 - 1].T - Nodes[iS21, jS21].T) / dy) * (dx / 2);
            var QS2DR = (K * (Nodes[iS22, jS21 - 1].T - Nodes[iS22, jS21].T) / dy) * (dx / 2);

            var QS2UL = (K * (Nodes[iS21, jS22+1].T - Nodes[iS21, jS22].T) / dy) * (dx / 2);
            var QS2UR = (K * (Nodes[iS22, jS22+1].T - Nodes[iS22, jS22].T) / dy) * (dx / 2);

            //S2 Left And Right
            for (int j = jS21+1; j < jS22; j++)
            {
                QS2L += (K * (Nodes[iS21 - 1, j].T - Nodes[iS21, j].T) / dx) * dy;
                QS2R += (K * (Nodes[iS22 + 1, j].T - Nodes[iS22, j].T) / dx) * dy;
            }
            //Up and Down Corners
            var QS2LD = (K * (Nodes[iS21 - 1, jS21].T - Nodes[iS21, jS21].T) / dx) * (dy / 2);
            var QS2LU = (K * (Nodes[iS21 - 1, jS22].T - Nodes[iS21, jS22].T) / dx) * (dy / 2);

            var QS2RD = (K * (Nodes[iS22 + 1, jS21].T - Nodes[iS22, jS21].T) / dx) * (dy / 2);
            var QS2RU = (K * (Nodes[iS22 + 1, jS22].T - Nodes[iS22, jS22].T) / dx) * (dy / 2);

            double QS3D=0,QS3U = 0;
            //S3 Up and Down
            for (int i = iS31 + 1; i < iS32; i++)
            {
                QS3D += (K * (Nodes[i, jS31 - 1].T - Nodes[i, jS31].T) / dy) * dx;
                QS3U += (K * (Nodes[i, jS32 + 1].T - Nodes[i, jS32].T) / dy) * dx;
            }
            //Left And Right Corners
            var QS3DL = (K * (Nodes[iS31, jS31 - 1].T - Nodes[iS31, jS31].T) / dy) * (dx / 2);
            var QS3DR = (K * (Nodes[iS32, jS31 - 1].T - Nodes[iS32, jS31].T) / dy) * (dx / 2);

            var QS3UL = (K * (Nodes[iS31, jS32 + 1].T - Nodes[iS31, jS32].T) / dy) * (dx / 2);
            var QS3UR = (K * (Nodes[iS32, jS32 + 1].T - Nodes[iS32, jS32].T) / dy) * (dx / 2);

            double QS3L = 0;
            //S3 Left And Right
            for (int j = jS31 + 1; j < jS32; j++)
            {
                QS3L += (K * (Nodes[iS31 - 1, j].T - Nodes[iS31, j].T) / dx) * dy;
            }
            //Up and Down Corners
            var QS3LD = (K * (Nodes[iS31 - 1, jS31].T - Nodes[iS31, jS31].T) / dx) * (dy / 2);
            var QS3LU = (K * (Nodes[iS31 - 1, jS32].T - Nodes[iS31, jS32].T) / dx) * (dy / 2);

            double qs = 0;
            for (int i = iS11; i <= iS12; i++)
            {
                for (int j = jS11; j <= jS12; j++)
                {
                    qs += Nodes[i, j].Sc * Nodes[i, j].dx * Nodes[i, j].dy;
                }
            }

            var Q =qs + QD + QDL + QDR + QR + QRD + QRU + QS2D + QS2DL + QS2DR + QS2L + QS2LD + QS2LU + QS2R + QS2RD + QS2RU + QS2U + QS2UL + QS2UR + QS3D + QS3DL + QS3DR + QS3L + QS3LD + QS3LU + QS3U + QS3UL + QS3UR;
            Console.WriteLine(Q);
            Console.ReadLine();
        }
    }

}
