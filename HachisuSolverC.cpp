using namespace std;
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include<string>
#include<iomanip>
#include <math.h>
#define pi 3.141592653589793


//fix params
double rho_pc = 1;
//rho_ec = 0.66965
//rho_ec = 0.669642 #is close to no atmosphere
double rho_ec = 0.66964175;
double TOL = 0.0001;
double magnetOn = 0;
double protonRatio=0.95;
double electronRatio=1;

//integration params
int LMAX = 50;
double me = 0.1;
double h = 10;
double C = 2;
double c = 2 * 137 / h;
double mu0 = 4 * pi * h / (c * 137);
double sigma = pow((3.0 / pi), (2.0 / 3)) * (5 * h * h / (40 * me * 4 * pi));
//grid will be a square, must be an odd number
int gridsize = 129;
//need to change this manually. First two entries are the grid size, and last entry should be half of LMAX
double FGrid[129][129][50];

string Mstr = "MagneticGrid1C.txt";
string Estr = "ElectricGrid1C.txt";
string Gstr = "GravityGrid1C.txt";
string Hpstr = "HpGrid1C.txt";
string Hestr = "HeGrid1C.txt";
string elecstr = "iteration1electronC.txt";
string protstr = "iteration1protonC.txt";
string Ecurvestr = "ElectronCurve1.txt";
string Pcurvestr = "ProtonCurve1.txt";


double f(double t, double ue, double up, double sp) {
    if (t == 0) {
        return (C - 1) * pow(up,(3.0 / 2)) / (sigma * me) - (C + me) *pow( ue ,(3.0 / 2)) / (sigma * me);
    }

    else if (up == 0) {
        return 0;
    }
    else {
        return -2 / t * sp + (C - 1) *pow( up ,(3.0 / 2)) / (sigma * me) - (C + me) * pow(ue,(3.0 / 2)) / (sigma * me);
    }
}

double g(double t, double ue, double up,double  se) {
    if (t == 0) {
        return -(C + me) * pow(up, (3.0 / 2)) / (sigma)+(C - me*me) * pow(ue, (3.0 / 2)) / (sigma);
    }
    else if (up == 0) {
        return 0;
    }
    else {
        return -2 / t * se - (C + me) * pow(up, (3.0 / 2)) / (sigma)+(C - me*me) * pow(ue, (3.0 / 2)) / (sigma);
    }
}

double F(double r1, double r2, int n) {
    if (r1 < r2) {
        return pow(r1, (n + 2)) / pow(r2, (n + 1));
    }
    else if (r2 < r1)
    {
        return pow(r2,n) / pow(r1,(n - 1));
    }
    else {
        return r1;
    }
}


int main() {

    // To save time later, we precompute the legendre polynomial grid, and some other grids
    double thetaStep = pi /(2* (gridsize - 1));
    double LegendreGrid[LMAX/2][gridsize];
    //is this numbered correctly?
    for (int k = 0; k < LMAX/2; k++) {
        for (int n = 0;n < gridsize; n++) {
            LegendreGrid[k][n] = legendre(2*k, cos(thetaStep * n)) * sin(thetaStep * n);
        }
    }
    double SinGrid[gridsize];
    for (int n = 0;n < gridsize; n++) {
        SinGrid[n] = sin(thetaStep * n);
    }


    //Step 1: solve nonrotating problem.
    double t = 0;
    double ue = rho_ec;
    double up = rho_pc;
    double se = 0;
    double sp = 0;
    int flag = 0;
    int flagp = 0;
    int flage = 0;
    int i = 0;
    double step = 0.0000001;

    double nsp;
    double nse;
    double nue;
    double nup;
    double radius_e = 0;
    double radius_p = 0;

    ofstream logfile;
    logfile.open("outputC.txt");
    do {

        if (sp > 0.01) {
            cout << "sp" << '\n';
            flag = 1;
            break;
        }

        if (se > 0.01) {
            cout << "se" << '\n';
            flag = 1;
            break;
        }
        if (ue + step * se < 0) {
            if (flage == 0) {
                flage = 1;
                radius_e = t;}
            ue = 0;
            se = 0;
        }
        if (up + step * sp < 0) {
            if (flagp == 0) {
                flagp = 1;
                radius_p = t;}
            up = 0;
            sp = 0;
        }

        nsp = sp + step * f(t, ue, up, sp);
        nse = se + step * g(t, ue, up, se);
        nup = up + step * sp;
        nue = ue + step * se;

        sp = nsp;
        up = nup;
        se = nse;
        ue = nue;
        if (logfile.is_open() && remainder(i,100000)==0)
            logfile << t << " " << pow(up, 3.0 / 2) << " " << pow(ue, 3.0 / 2) << "\n";
        t = t + step;
        i += 1;
    } while ((up > 0 || ue > 0) && i < 5000000000);
    logfile.close();
    logfile.clear();
    
    if ((flage == 0 or flagp == 0)) {
        flag = 1;
    }
    cout.precision(8);
    cout  << "initial electron radius: "  <<radius_e << "\n";
    cout  << "initial proton radius: "  << radius_p << "\n";


    //if flag=1, the program did not conclude succesfully, and these may be bad values
    if (flag == 1) {
        cout << "not a good central pair" << "\n";
        exit(EXIT_FAILURE);
    }

    //Step 2
    double rmax = 32.0 / 31 * max(radius_p, radius_e);

    //closest entry to actual radius
    double closestStartingeR = radius_e / rmax * (gridsize - 1);
    double closestStartingpR = radius_p / rmax * (gridsize - 1);
    cout << "closest starting point on grid for electron radius: "<<closestStartingeR << "\n";
    cout << "closest starting point on grid for proton radius: " << closestStartingpR << "\n";

    //put the starting values into the grid
    //so instead, we will iterate again and choose the rmax / (10 * gridsize) as the step size

    double ProtonGrid[gridsize][gridsize];
    double ElectronGrid[gridsize][gridsize];
    for (int j = 0;j < gridsize; j++) {
        for (int i = 0;i < gridsize; i++) {
            ProtonGrid[j][i] =0;
            ElectronGrid[j][i] = 0;
        }
    }
    t = 0;
    ue = rho_ec;
    up = rho_pc;
    se = 0;
    sp = 0;
    int j = 0;

    step = rmax / (1000000 * (gridsize-1));

    do {

        if (ue + step * se < 0) {
            ue = 0;
            se = 0;
        }
        if (up + step * sp < 0) {
            up = 0;
            sp = 0;
        }

        nsp = sp + step * f(t, ue, up, sp);
        nse = se + step * g(t, ue, up, se);
        nup = up + step * sp;
        nue = ue + step * se;

        sp = nsp;
        up = nup;
        se = nse;
        ue = nue;

        t = t + step;
        if (remainder(j, 1000000) == 0) {
            for (int i = 0;i < gridsize; i++) {
                ProtonGrid[j / 1000000][i] = pow(up, (3.0 / 2));
                ElectronGrid[j / 1000000][i] = pow(ue, (3.0 / 2));
            }
        }
        j += 1;
        

    } while ((up > 0 || ue > 0) && j < 1000000 * gridsize);


    double rStep = rmax / (gridsize - 1);
    cout << rStep << "\n";

    /*
    logfile.open(Pcurvestr);
    //cout << "printing curve" << "\n";
    for (int i = 0; i < gridsize; i++) {
        for (int j = 0; j < gridsize; j++) {
            if (ProtonGrid[j][i] == 0) {
                logfile << j * rStep * cos(pi / 2 - i * thetaStep) << " " << j * rStep * cos(i * thetaStep) << "\n";
                break;
            }
        }
    }
    logfile.close();
    logfile.open(Ecurvestr);
    //cout << "printing curve" << "\n";
    for (int i = 0; i < gridsize; i++) {
        for (int j = 0; j < gridsize; j++) {
            if (ElectronGrid[j][i] == 0) {
                logfile << j * rStep * cos(pi / 2 - i * thetaStep) << " " << j * rStep * cos(i * thetaStep) << "\n";
                break;
            }
        }
    }
    */
    logfile.close();
    //step 2: compute the potentials
    double ElectricGrid[gridsize][gridsize];
    double GravityGrid[gridsize][gridsize];
    double MagneticGrid[gridsize][gridsize];
    
    //compute potentials. since we have to integrate over the grid multiple times, should do all potentials at once
    // seems like we need to provide an initial guess for omega
    //double OmegaGrid[gridsize];
    //for (int i = 0;i < gridsize; i++) {
    //    OmegaGrid[i] = omega * rStep * i;
    //}
   
    for (int m = 0; m < gridsize; m++) {
        for (int j = 0; j < gridsize; j++) {
            for (int k = 0; k < LMAX/2; k++) {
                FGrid[m][j][k] = F(rStep * m, rStep * j, 2*k);
            }
        }
    }
    
    double g = 2.0/45* pi /(2* (gridsize - 1));
    double factor = (2.0 / 45) * rStep;

    double D_rEk;
    double D_rMk;
    double D_rGk;
    double D_tEm0;
    double D_tMm0;
    double D_tGm0;
    int m;
    int n;

    for (int j = 0; j < gridsize; j++) {
        for (int i = 0; i < gridsize; i++) {
            ElectricGrid[j][i] = 0;
            MagneticGrid[j][i] = 0;
            GravityGrid[j][i] = 0;
        }
    }
    
    int upper = gridsize - 1;
    for (int j=0; j <gridsize ; j++) {
        //cout <<j<<'\n';
        // j is radial coordinate
        for (int i = 0; i < gridsize; i++) {
            // i is theta coordinate
            for (int k = 0; k < LMAX/2; k++) {
                D_rEk = 0;
                D_rMk = 0;
                D_rGk = 0;
                m = 0;
                do {
                    //SHOULD THESE BE RAISED POWERS? WHAT IS IN PROTONGRID? no those grids account for the power.
                    D_tEm0 = 7*LegendreGrid[k][0] * (ProtonGrid[m][0] - ElectronGrid[m][0]) - 7*LegendreGrid[k][upper] * (ProtonGrid[m][upper] - ElectronGrid[m][upper]);
                    //D_tMm0 = 7*LegendreGrid[k][0] * (ProtonGrid[m][0] - ElectronGrid[m][0]) * OmegaGrid[m] * SinGrid[0] - 7*LegendreGrid[k][upper] * (ProtonGrid[m][upper] - ElectronGrid[m][upper]) * OmegaGrid[m] * SinGrid[upper];
                    D_tGm0 = 7*LegendreGrid[k][0] * (ProtonGrid[m][0] + me * ElectronGrid[m][0]) - 7*LegendreGrid[k][upper] * (ProtonGrid[m][upper] + me * ElectronGrid[m][upper]);

                    n = 1;
                    do {
                        //this loop is Simpson's along the theta direction. n here is the (first) i in the Hachisu paper
                        D_tEm0 += 32 * LegendreGrid[k][n] * (ProtonGrid[m][n] - ElectronGrid[m][n]) + 12 * LegendreGrid[k][n + 1] * (ProtonGrid[m][n + 1] - ElectronGrid[m][n + 1])+ 32 * LegendreGrid[k][n + 2] * (ProtonGrid[m][n + 2] - ElectronGrid[m][n + 2])+ 14 * LegendreGrid[k][n + 3] * (ProtonGrid[m][n + 3] - ElectronGrid[m][n + 3]);
                        // for the magnetic potential, since we have axisymmetry, we assume that the y coordinate of Rperp is zero, and therefore we have only the y coordinate in the cross product
                        //D_tMm0 += 32 * LegendreGrid[k][n] * (ProtonGrid[m][n] - ElectronGrid[m][n]) * OmegaGrid[m] * SinGrid[n] + 12 * LegendreGrid[k][n + 1] * (ProtonGrid[m][n + 1] - ElectronGrid[m][n + 1]) * OmegaGrid[m] * SinGrid[n + 1]+ 32 * LegendreGrid[k][n + 2] * (ProtonGrid[m][n + 2] - ElectronGrid[m][n + 2]) * OmegaGrid[m] * SinGrid[n + 2]+ 14 * LegendreGrid[k][n + 3] * (ProtonGrid[m][n + 3] - ElectronGrid[m][n + 3]) * OmegaGrid[m] * SinGrid[n + 3];
                        D_tGm0 += 32 * LegendreGrid[k][n] * (ProtonGrid[m][n] + me * ElectronGrid[m][n]) + 12 * LegendreGrid[k][n + 1] * (ProtonGrid[m][n + 1] + me * ElectronGrid[m][n + 1])+ 32 * LegendreGrid[k][n + 2] * (ProtonGrid[m][n + 2] + me * ElectronGrid[m][n + 2])+ 14 * LegendreGrid[k][n + 3] * (ProtonGrid[m][n + 3] + me * ElectronGrid[m][n + 3]);
                        n += 4;
                        //cout << D_tMm0 << " ";
                   
                    }while (n < upper);
                    //Simpson's in r direction. m here is k in the hachisu paper
                    if (m == 0 || m == upper){
                        D_rEk += 7*FGrid[m][j][k] * D_tEm0;
                        //D_rMk += 7*FGrid[m][j][k] * D_tMm0;
                        D_rGk += 7*FGrid[m][j][k] * D_tGm0;
                    }
                    else if (remainder(m, 4) == 0) {
                        //cout << 0 <<" "<<FGrid[m][j][k]<<" ";
                        D_rEk += 14 * FGrid[m][j][k] * D_tEm0;
                        //D_rMk += 14 * FGrid[m][j][k] * D_tMm0;
                        D_rGk += 14 * FGrid[m][j][k] * D_tGm0;
                    }
                    else if (remainder(m, 4) == 2) {
                        //cout << 0 <<" "<<FGrid[m][j][k]<<" ";
                        D_rEk += 12 * FGrid[m][j][k] * D_tEm0;
                        //D_rMk += 12 * FGrid[m][j][k] * D_tMm0;
                        D_rGk += 12 * FGrid[m][j][k] * D_tGm0;
                    } else if (remainder(m,2) ==1 || remainder(m,2)==-1) {
                        //cout << 1 << " " << FGrid[m][j][k] << " ";
                        D_rEk += 32 * FGrid[m][j][k] * D_tEm0;
                        //D_rMk += 32 * FGrid[m][j][k] * D_tMm0;
                        D_rGk += 32 * FGrid[m][j][k] * D_tGm0;
                    }
                    //cout <<D_rMk<< " " <<m<< "\n";
                    m += 1;
                }while (m < gridsize);
                ElectricGrid[j][i] += factor * g * 4 * pi * C * D_rEk * legendre(2*k,cos(thetaStep * i));
                //MagneticGrid[j][i] += factor * g  * mu0 * D_rMk * legendre(2*k,cos(thetaStep * i));
                GravityGrid[j][i] += factor * g * 4 * pi * D_rGk * legendre(2*k,cos(thetaStep * i));
            }
        }
    }

    //step 4: compute lambdap and lambdae
    // what are we claiming the proton and electron radius are ? pRad and eRad should be integers
    int pRad =(31.0/32)*(gridsize-1);
    int eRad = round(pRad*electronRatio);
    int paxisRatio = round(pRad * protonRatio);
    cout << "new electron radius: " << rStep * eRad << "\n";
    cout << "new proton radius: " << rStep * pRad << "\n";
    

    // now compute lambdap, lambdae, and omega
    double Ulambdap = ElectricGrid[paxisRatio][0] - GravityGrid[paxisRatio][0];
    double omega = pow(abs(2.0/pow((rStep * pRad), 2) *(ElectricGrid[pRad][upper] - GravityGrid[pRad][upper]-Ulambdap)), 0.5);
    double Ulambdae = -ElectricGrid[eRad][upper] - me * GravityGrid[eRad][upper] - magnetOn * MagneticGrid[eRad][upper] * omega * rStep * eRad - me * pow(omega, 2) * (1.0 / 2) * pow((rStep * eRad), 2);

    
    double lambdap = 0;
    double lambdae = 0;

    cout << Ulambdap<< "\n";
    cout <<omega << "\n";
    //step 5: Compute H's

    double HpGrid[gridsize][gridsize];
    double HeGrid[gridsize][gridsize];
    double UHpGrid[gridsize][gridsize];
    double UHeGrid[gridsize][gridsize];

    for (int j = 0;j < gridsize; j++) {
        for (int i = 0;i < gridsize; i++) {
            HpGrid[j][i] = 0;
            HeGrid[j][i] = 0;
            // THERE WAS AN ERROR HERE WITH ADDING THE ROTATION AND MULTIPLY BY SIN(THETA). THIS MIGHT BE APPEARING IN OTHER PLACES AS WELL.
            UHpGrid[j][i] = -ElectricGrid[j][i] + GravityGrid[j][i] - magnetOn * MagneticGrid[j][i] * omega * rStep * j*sin(thetaStep*i) + pow(omega,2) * (1.0 / 2) * pow((rStep * j) * sin(thetaStep * i),2) + Ulambdap;
            UHeGrid[j][i] = ElectricGrid[j][i] + me * GravityGrid[j][i] + magnetOn * MagneticGrid[j][i] * omega * rStep * j*sin(thetaStep*i)+ me * pow(omega,2) * (1.0 / 2) * pow((rStep * j) * sin(thetaStep * i),2) + Ulambdae;
        }
    }

    //step 6: Repeat
    int count = 1;
    bool cond = true;
   
    logfile.open(protstr);
    for (int j = 0;j < gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            logfile << ProtonGrid[j][i] <<" " ;
        }
        logfile << "\n";
        
    }
    logfile.close();
    logfile.clear();

    logfile.open(elecstr);
    for (int j = 0;j < gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            logfile << ElectronGrid[j][i] << " ";
        }
        logfile << "\n";

    }
    logfile.close();
    logfile.clear();

    logfile.open(Hpstr);
    for (int j = 0;j < gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            logfile << UHpGrid[j][i] << " ";
        }
        logfile << "\n";

    }
    logfile.close();
    logfile.clear();

    logfile.open(Hestr);
    for (int j = 0;j < gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            logfile << UHeGrid[j][i] << " ";
        }
        logfile << "\n";

    }
    logfile.close();
    logfile.clear();

    logfile.open(Estr);
    for (int j = 0;j < gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            logfile << ElectricGrid[j][i] << " ";
        }
        logfile << "\n";

    }
    logfile.close();
    logfile.clear();

    logfile.open(Mstr);
    for (int j = 0;j < gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            logfile << MagneticGrid[j][i] << " ";
        }
        logfile << "\n";

    }
    logfile.close();
    logfile.clear();

    logfile.open(Gstr);
    for (int j = 0;j < gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            logfile << GravityGrid[j][i] << " ";
        }
        logfile << "\n";

    }
    logfile.close();
    logfile.clear();

    double Hpmax;
    double Hpdiffmax;
    double Hemax;
    double Hediffmax;
    double lpdiff;
    double lediff;

    char Charcount;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // start the iteration here
    do {
        count += 1;
        Charcount = '0' + count;
        protstr[9] = Charcount;
        elecstr[9] = Charcount;
        Mstr[12] = Charcount;
        Estr[12] = Charcount;
        Gstr[11] = Charcount;
        Hpstr[6] = Charcount;
        Hestr[6] = Charcount;
        Ecurvestr[13] = Charcount;
        Pcurvestr[11] = Charcount;
        cout << "iteration: " << count << "\n";
        //need to transfer updated values to old values

        if (UHpGrid[0][0] <=0) {
            cout << "Something has gone terribly wrong.\n";
            break;
        }

        for (int j = 0;j < gridsize; j++) {
            for (int i = 0;i < gridsize; i++) {
                HpGrid[j][i] = UHpGrid[j][i];
                HeGrid[j][i] = UHeGrid[j][i];
                if (UHpGrid[j][i] > 0){
                    ProtonGrid[j][i] = pow((UHpGrid[j][i] / UHpGrid[0][0]) ,(3.0 / 2));}
                else {
                    ProtonGrid[j][i] = 0;
                }
                if (UHeGrid[j][i] > 0){

                    ElectronGrid[j][i] = pow((me * UHeGrid[j][i] / UHpGrid[0][0]),(3.0 / 2));
                }
                else {
                    ElectronGrid[j][i] = 0;
                }
            }
        }

        
        logfile.open(protstr);
        for (int j = 0;j < gridsize; j++) {
            logfile << j << ": ";
            for (int i = 0; i < gridsize; i++) {
                logfile << ProtonGrid[j][i] << " ";
            }
            logfile << "\n";

        }
        logfile.close();
        logfile.clear();

        logfile.open(elecstr);
        for (int j = 0;j < gridsize; j++) {
            logfile << j << ": ";
            for (int i = 0; i < gridsize; i++) {
                logfile << ElectronGrid[j][i] << " ";
            }
            logfile << "\n";

        }
        logfile.close();
        logfile.clear();

        for (int j = 0; j < gridsize; j++) {
            for (int i = 0; i < gridsize; i++) {
                ElectricGrid[j][i] = 0;
                MagneticGrid[j][i] = 0;
                GravityGrid[j][i] = 0;
            }
        }

        /*
        logfile.open(Pcurvestr);
        //cout << "printing curve" << "\n";
        for (int i = 0; i < gridsize; i++) {
            for (int j = 0; j < gridsize; j++) {
                if (ProtonGrid[j][i] == 0) {
                    logfile << j * rStep * cos(pi / 2 - i * thetaStep) << " " << j * rStep * cos(i * thetaStep) << "\n";
                    break;
                }
            }
        }
        logfile.close();
        logfile.open(Ecurvestr);
        //cout << "printing curve" << "\n";
        for (int i = 0; i < gridsize; i++) {
            for (int j = 0; j < gridsize; j++) {
                if (ElectronGrid[j][i] == 0) {
                    logfile << j * rStep * cos(pi / 2 - i * thetaStep) << " " << j * rStep * cos(i * thetaStep) << "\n";
                    break;
                }
            }
        }
        logfile.close();
        */

        for (int j = 0; j < gridsize; j++) {
            //cout << j << '\n';
            // j is radial coordinate
            for (int i = 0; i < gridsize; i++) {
                // i is theta coordinate
                for (int k = 0; k < LMAX/2; k++) {
                    D_rEk = 0;
                    D_rMk = 0;
                    D_rGk = 0;
                    m = 0;

                    //note the we go from 0 to 256 while Hachisu goes from 1 to 257
                    do {
                        D_tEm0 = 7 * LegendreGrid[k][0] * (ProtonGrid[m][0] - ElectronGrid[m][0]) - 7 * LegendreGrid[k][upper] * (ProtonGrid[m][upper] - ElectronGrid[m][upper]);
                        //D_tMm0 = 7 * LegendreGrid[k][0] * (ProtonGrid[m][0] - ElectronGrid[m][0]) * OmegaGrid[m] * SinGrid[0] - 7 * LegendreGrid[k][upper] * (ProtonGrid[m][upper] - ElectronGrid[m][upper]) * OmegaGrid[m] * SinGrid[upper];
                        D_tGm0 = 7 * LegendreGrid[k][0] * (ProtonGrid[m][0] + me * ElectronGrid[m][0]) - 7 * LegendreGrid[k][upper] * (ProtonGrid[m][upper] + me * ElectronGrid[m][upper]);

                        n = 1;
                        do {
                            //this may have some double counting for the last term.

                            //this loop is Simpson's along the theta direction. n here is the (first) i in the Hachisu paper
                            D_tEm0 += 32 * LegendreGrid[k][n] * (ProtonGrid[m][n] - ElectronGrid[m][n]) + 12 * LegendreGrid[k][n + 1] * (ProtonGrid[m][n + 1] - ElectronGrid[m][n + 1]) + 32 * LegendreGrid[k][n + 2] * (ProtonGrid[m][n + 2] - ElectronGrid[m][n + 2]) + 14 * LegendreGrid[k][n + 3] * (ProtonGrid[m][n + 3] - ElectronGrid[m][n + 3]);
                            // for the magnetic potential, since we have axisymmetry, we assume that the y coordinate of Rperp is zero, and therefore we have only the y coordinate in the cross product
                            //D_tMm0 += 32 * LegendreGrid[k][n] * (ProtonGrid[m][n] - ElectronGrid[m][n]) * OmegaGrid[m] * SinGrid[n] + 12 * LegendreGrid[k][n + 1] * (ProtonGrid[m][n + 1] - ElectronGrid[m][n + 1]) * OmegaGrid[m] * SinGrid[n + 1] + 32 * LegendreGrid[k][n + 2] * (ProtonGrid[m][n + 2] - ElectronGrid[m][n + 2]) * OmegaGrid[m] * SinGrid[n + 2] + 14 * LegendreGrid[k][n + 3] * (ProtonGrid[m][n + 3] - ElectronGrid[m][n + 3]) * OmegaGrid[m] * SinGrid[n + 3];
                            D_tGm0 += 32 * LegendreGrid[k][n] * (ProtonGrid[m][n] + me * ElectronGrid[m][n]) + 12 * LegendreGrid[k][n + 1] * (ProtonGrid[m][n + 1] + me * ElectronGrid[m][n + 1]) + 32 * LegendreGrid[k][n + 2] * (ProtonGrid[m][n + 2] + me * ElectronGrid[m][n + 2]) + 14 * LegendreGrid[k][n + 3] * (ProtonGrid[m][n + 3] + me * ElectronGrid[m][n + 3]);
                            n += 4;
                            //cout << D_tMm0 << " ";

                        } while (n < upper);
                        //Simpson's in r direction. m here is k in the hachisu paper
                        if (m == 0 || m == upper) {
                            D_rEk += 7 * FGrid[m][j][k] * D_tEm0;
                            //D_rMk += 7 * FGrid[m][j][k] * D_tMm0;
                            D_rGk += 7 * FGrid[m][j][k] * D_tGm0;
                        }
                        else if (remainder(m, 4) == 0) {
                            //cout << 0 <<" "<<FGrid[m][j][k]<<" ";
                            D_rEk += 14 * FGrid[m][j][k] * D_tEm0;
                            //D_rMk += 14 * FGrid[m][j][k] * D_tMm0;
                            D_rGk += 14 * FGrid[m][j][k] * D_tGm0;
                        }
                        else if (remainder(m, 4) == 2) {
                            //cout << 0 <<" "<<FGrid[m][j][k]<<" ";
                            D_rEk += 12 * FGrid[m][j][k] * D_tEm0;
                            //D_rMk += 12 * FGrid[m][j][k] * D_tMm0;
                            D_rGk += 12 * FGrid[m][j][k] * D_tGm0;
                        }
                        else if (remainder(m, 2) == 1 || remainder(m, 2) == -1) {
                            //cout << 1 << " " << FGrid[m][j][k] << " ";
                            D_rEk += 32 * FGrid[m][j][k] * D_tEm0;
                            //D_rMk += 32 * FGrid[m][j][k] * D_tMm0;
                            D_rGk += 32 * FGrid[m][j][k] * D_tGm0;
                        }

                        m += 1;
                    }while (m < gridsize);

                    ElectricGrid[j][i] += factor * g * 4 * pi * C * D_rEk * legendre(2*k, cos(thetaStep * i));
                    MagneticGrid[j][i] += factor * g  * mu0 * D_rMk * legendre(2*k, cos(thetaStep * i));
                    GravityGrid[j][i] += factor * g * 4 * pi * D_rGk * legendre(2*k, cos(thetaStep * i));
                }
            }
        }

            

        logfile.open(Estr);
        for (int j = 0;j < gridsize; j++) {
            logfile << j << ": ";
            for (int i = 0; i < gridsize; i++) {
                logfile << ElectricGrid[j][i] << " ";
            }
            logfile << "\n";

        }
        logfile.close();
        logfile.clear();

        logfile.open(Mstr);
        for (int j = 0;j < gridsize; j++) {
            logfile << j << ": ";
            for (int i = 0; i < gridsize; i++) {
                logfile << MagneticGrid[j][i] << " ";
            }
            logfile << "\n";

        }
        logfile.close();
        logfile.clear();

        logfile.open(Gstr);
        for (int j = 0;j < gridsize; j++) {
            logfile << j << ": ";
            for (int i = 0; i < gridsize; i++) {
                logfile << GravityGrid[j][i] << " ";
            }
            logfile << "\n";

        }
        logfile.close();
        logfile.clear();
        //step 4: compute lambdap and lambdae

        // now compute lambdap, lambdae, and omega

        lambdap = Ulambdap;
        lambdae = Ulambdae;
        Ulambdap = ElectricGrid[paxisRatio][0] - GravityGrid[paxisRatio][0];
        omega = pow(abs(2.0 / pow((rStep * pRad), 2) * (ElectricGrid[pRad][upper] - GravityGrid[pRad][upper] - Ulambdap)), 0.5);
        Ulambdae = -ElectricGrid[eRad][upper] - me * GravityGrid[eRad][upper] - magnetOn * MagneticGrid[eRad][upper] * omega * rStep * eRad - me * pow(omega, 2) * (1.0 / 2) * pow((rStep * eRad), 2);
       

        cout << Ulambdap << "\n";
        cout << omega << "\n";

        for (int j = 0;j < gridsize; j++) {
            for (int i = 0;i < gridsize; i++) {
                UHpGrid[j][i] = -ElectricGrid[j][i] + GravityGrid[j][i] - magnetOn * MagneticGrid[j][i] *  omega * rStep * j* sin(thetaStep * i) + pow(omega,2) * (1.0 / 2) * pow((rStep * j) * sin(thetaStep * i),2) + Ulambdap;
                UHeGrid[j][i] = ElectricGrid[j][i] + me * GravityGrid[j][i] + magnetOn * MagneticGrid[j][i]  * omega * rStep * j* sin(thetaStep * i) + me * pow(omega,2) * (1.0 / 2) * pow((rStep * j) * sin(thetaStep * i),2) + Ulambdae;
            }
        }

        logfile.open(Hpstr);
        for (int j = 0;j < gridsize; j++) {
            logfile << j << ": ";
            for (int i = 0; i < gridsize; i++) {
                logfile << UHpGrid[j][i] << " ";
            }
            logfile << "\n";

        }
        logfile.close();
        logfile.clear();

        logfile.open(Hestr);
        for (int j = 0;j < gridsize; j++) {
            logfile << j << ": ";
            for (int i = 0; i < gridsize; i++) {
                logfile << UHeGrid[j][i] << " ";
            }
            logfile << "\n";

        }
        logfile.close();
        logfile.clear();

            //now we need to compute cond
        lpdiff = abs((Ulambdap - lambdap) / Ulambdap);
        lediff = abs((Ulambdae - lambdae) / Ulambdae);
        Hpmax = 0;
        Hpdiffmax = 0;
        Hemax = 0;
        Hediffmax = 0;
        for (int j = 0;j < gridsize; j++) {
            for (int i = 0;i < gridsize; i++) {
                if (abs(UHpGrid[j][i] - HpGrid[j][i]) > Hpdiffmax) {
                    Hpdiffmax = abs(UHpGrid[j][i] - HpGrid[j][i]);
                }
                if (abs(UHpGrid[j][i]) > Hpmax) {
                    Hpmax = abs(UHpGrid[j][i]);
                }
                if (abs(UHeGrid[j][i] - HeGrid[j][i]) > Hediffmax) {
                    Hediffmax = abs(UHeGrid[j][i] - HeGrid[j][i]);
                }
                if (abs(UHeGrid[j][i]) > Hemax) {
                    Hemax = abs(UHeGrid[j][i]);
                }
            }
        }
        cout << lpdiff << " " << lediff << " " << Hpdiffmax / Hpmax << " " << Hediffmax / Hemax << "\n";
        if (lpdiff < TOL && lediff < TOL && Hpdiffmax / Hpmax < TOL && Hediffmax / Hemax < TOL) {
            cond = false;
            cout << "success!\n";
        }
    }while (cond && count <5);
}