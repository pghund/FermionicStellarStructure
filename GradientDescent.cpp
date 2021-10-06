using namespace std;
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include<string>
#include<iomanip>
#define pi 3.141592653589793


//fix params
double rho_pc = 1;
double DescentParameter = 0.001;
//rho_ec = 0.66965
//rho_ec = 0.669642 #is close to no atmosphere
double rho_ec = 0.66964175;
double omega = 0.1;
double TOL = 0.001;
double magnetOn = 0;
double targeteradius = 9.243;
double targetpradius = 9.243;

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
double FGrid[129][129][25];

string Mstr = "MagneticGrid1C.txt";
string Estr = "ElectricGrid1C.txt";
string Gstr = "GravityGrid1C.txt";
string Hpstr = "HpGrid1C.txt";
string Hestr = "HeGrid1C.txt";
string elecstr = "iteration1electronC.txt";
string protstr = "iteration1protonC.txt";
string Ecurvestr = "ElectronCurve1.txt";
string Pcurvestr = "ProtonCurve1.txt";
string gradstr = "Gradient1.txt";
string updatestr = "updated1.txt";




double f(double t, double ue, double up, double sp) {
    if (t == 0) {
        return (C - 1) * pow(up, (3.0 / 2)) / (sigma * me) - (C + me) * pow(ue, (3.0 / 2)) / (sigma * me);
    }

    else if (up == 0) {
        return 0;
    }
    else {
        return -2 / t * sp + (C - 1) * pow(up, (3.0 / 2)) / (sigma * me) - (C + me) * pow(ue, (3.0 / 2)) / (sigma * me);
    }
}

double g(double t, double ue, double up, double  se) {
    if (t == 0) {
        return -(C + me) * pow(up, (3.0 / 2)) / (sigma)+(C - me * me) * pow(ue, (3.0 / 2)) / (sigma);
    }
    else if (up == 0) {
        return 0;
    }
    else {
        return -2 / t * se - (C + me) * pow(up, (3.0 / 2)) / (sigma)+(C - me * me) * pow(ue, (3.0 / 2)) / (sigma);
    }
}

double F(double r1, double r2, int n) {
    if (r1 < r2) {
        return pow(r1, (n + 2)) / pow(r2, (n + 1));
    }
    else if (r2 < r1)
    {
        return pow(r2, n) / pow(r1, (n - 1));
    }
    else {
        return r1;
    }
}


int main() {

    // To save time later, we precompute the legendre polynomial grid, and some other grids
    double thetaStep = pi / (2 * (gridsize - 1));
    double LegendreGrid[LMAX / 2][gridsize];
    double LGrid[LMAX / 2][gridsize];
    //is this numbered correctly?
    for (int k = 0; k < LMAX / 2; k++) {
        for (int n = 0;n < gridsize; n++) {
            LegendreGrid[k][n] = legendre(2 * k, cos(thetaStep * n)) * sin(thetaStep * n);
            LGrid[k][n] = legendre(2 * k, cos(thetaStep * n));
        }
    }
    double SinGrid[gridsize];
    for (int n = 0;n < gridsize; n++) {
        SinGrid[n] = sin(thetaStep * n);
    }

    double MGrid[gridsize];
    for (int n = 0;n < gridsize; n++) {
        if (n == 0 || n == gridsize-1) {
            MGrid[n] = 7;
        }
        else if (remainder(n, 4) == 0) {
            MGrid[n] = 14;
        }
        else if (remainder(n, 4) == 2) {
            MGrid[n] = 12;
        }
        else if (remainder(n, 2) == 1 || remainder(n, 2) == -1) {
            MGrid[n] = 32;
        }
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
                radius_e = t;
            }
            ue = 0;
            se = 0;
        }
        if (up + step * sp < 0) {
            if (flagp == 0) {
                flagp = 1;
                radius_p = t;
            }
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
        if (logfile.is_open() && remainder(i, 100000) == 0)
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
    cout << "initial electron radius: " << radius_e << "\n";
    cout << "initial proton radius: " << radius_p << "\n";


    //if flag=1, the program did not conclude succesfully, and these may be bad values
    if (flag == 1) {
        cout << "not a good central pair" << "\n";
        exit(EXIT_FAILURE);
    }

    //Step 2
    double rmax = 16.0 / 15 * max(radius_p, radius_e);

    //closest entry to actual radius
    double closestStartingeR = radius_e / rmax * (gridsize - 1);
    double closestStartingpR = radius_p / rmax * (gridsize - 1);
    cout << "closest starting point on grid for electron radius: " << closestStartingeR << "\n";
    cout << "closest starting point on grid for proton radius: " << closestStartingpR << "\n";

    //put the starting values into the grid
    //so instead, we will iterate again and choose the rmax / (10 * gridsize) as the step size

    double ProtonGrid[gridsize][gridsize];
    double ElectronGrid[gridsize][gridsize];
    double UProtonGrid[gridsize][gridsize];
    double UElectronGrid[gridsize][gridsize];
    for (int j = 0;j < gridsize; j++) {
        for (int i = 0;i < gridsize; i++) {
            ProtonGrid[j][i] = 0;
            ElectronGrid[j][i] = 0;
        }
    }
    t = 0;
    ue = rho_ec;
    up = rho_pc;
    se = 0;
    sp = 0;
    int j = 0;

    step = rmax / (1000000 * (gridsize - 1));

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
    //cout << rStep << "\n";
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
    //step 2: compute the potentials
    double ElectricGrid[gridsize][gridsize];
    double GravityGrid[gridsize][gridsize];
    double MagneticGrid[gridsize][gridsize];

    //compute potentials. since we have to integrate over the grid multiple times, should do all potentials at once
    double OmegaGrid[gridsize];
    for (int i = 0;i < gridsize; i++) {
        OmegaGrid[i] = omega * rStep * i;
    }

    for (int m = 0; m < gridsize; m++) {
        for (int j = 0; j < gridsize; j++) {
            for (int k = 0; k < LMAX / 2; k++) {
                FGrid[m][j][k] = F(rStep * m, rStep * j, 2 * k);
            }
        }
    }

    double g = 2.0 / 45 * pi / (2 * (gridsize - 1));
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
    for (int j = 0; j < gridsize; j++) {
        //cout <<j<<'\n';
        // j is radial coordinate
        for (int i = 0; i < gridsize; i++) {
            // i is theta coordinate
            for (int k = 0; k < LMAX / 2; k++) {
                D_rEk = 0;
                D_rMk = 0;
                D_rGk = 0;
                m = 0;
                for (int m = 0; m < gridsize; m++) {
                    D_tEm0 = 0;
                    D_tGm0 = 0;
                    //D_tMm0=0;

                    for (int n = 0; n < gridsize; n++) {
                        D_tEm0 += MGrid[n] * LegendreGrid[k][n] * (ProtonGrid[m][n] - ElectronGrid[m][n]);
                        D_tGm0 += MGrid[n] * LegendreGrid[k][n] * (ProtonGrid[m][n] + me * ElectronGrid[m][n]);
                        //D_tMm0 += MGrid[n] * LegendreGrid[k][n] * (ProtonGrid[m][n] - ElectronGrid[m][n]) * OmegaGrid[m] * SinGrid[n];
                    }
                    D_rEk += MGrid[m]* FGrid[m][j][k] * D_tEm0;
                    //D_rMk += MGrid[m] * FGrid[m][j][k] * D_tMm0;
                    D_rGk += MGrid[m]* FGrid[m][j][k] * D_tGm0;
                }
                ElectricGrid[j][i] += factor * g * 4 * pi * C * D_rEk * LGrid[k][i];
                //MagneticGrid[j][i] += factor * g  * mu0 * D_rMk * LGrid[k][i];
                GravityGrid[j][i] += factor * g * 4 * pi * D_rGk * LGrid[k][i];
            }
        }
    }

    //step 4: compute lambdap and lambdae
    // what are we claiming the proton and electron radius are ?
    int pRad = floor(targetpradius / rmax * (gridsize - 1));
    int eRad = floor(targeteradius / rmax * (gridsize - 1));
    cout << "new electron radius: " << rStep * eRad << "\n";
    cout << "new proton radius: " << rStep * pRad << "\n";


    double lambdap = ElectricGrid[pRad][upper] - GravityGrid[pRad][upper] + magnetOn * MagneticGrid[pRad][upper] * omega * rStep * pRad - pow(omega, 2) * (1.0 / 2) * pow((rStep * pRad), 2);
    double lambdae = -ElectricGrid[eRad][upper] - me * GravityGrid[eRad][upper] - magnetOn * MagneticGrid[eRad][upper] * omega * rStep * eRad - me * pow(omega, 2) * (1.0 / 2) * pow((rStep * eRad), 2);


    //cout << Ulambdap<< "\n";
    //cout <<Ulambdae << "\n";
    //step 5: Compute H's

    double HpGrid[gridsize][gridsize];
    double HeGrid[gridsize][gridsize];

    for (int j = 0;j < gridsize; j++) {
        for (int i = 0;i < gridsize; i++) {
            // THERE WAS AN ERROR HERE WITH ADDING THE ROTATION AND MULTIPLY BY SIN(THETA). THIS MIGHT BE APPEARING IN OTHER PLACES AS WELL.
            HpGrid[j][i] = -ElectricGrid[j][i] + GravityGrid[j][i] - magnetOn * MagneticGrid[j][i] * omega * rStep * j * sin(thetaStep * i) + pow(omega, 2) * (1.0 / 2) * pow((rStep * j) * sin(thetaStep * i), 2) + lambdap;
            HeGrid[j][i] = ElectricGrid[j][i] + me * GravityGrid[j][i] + magnetOn * MagneticGrid[j][i] * omega * rStep * j * sin(thetaStep * i) + me * pow(omega, 2) * (1.0 / 2) * pow((rStep * j) * sin(thetaStep * i), 2) + lambdae;
        }
    }

 
    for (int j = 0;j < gridsize; j++) {
        for (int i = 0;i < gridsize; i++) {
            if (HpGrid[j][i] > 0) {
                UProtonGrid[j][i] = pow((HpGrid[j][i] / HpGrid[0][0]), (3.0 / 2));
            }
            else {
                UProtonGrid[j][i] = 0;
            }
            if (HeGrid[j][i] > 0) {

                UElectronGrid[j][i] = pow((me * HeGrid[j][i] / HpGrid[0][0]), (3.0 / 2));
            }
            else {
                UElectronGrid[j][i] = 0;
            }
        }
    }
    double ProtonByProton;
    double ElectronByProton;
    double ProtonByElectron;
    double ElectronByElectron;
    double innersum;
    double sumP=0;
    double sumE = 0;
    double radial;
    double factor2 = factor * g * 4 * pi;
    double Gradient[2 * gridsize][gridsize];
    double factor3 = 2*rStep * sin(thetaStep / 2);
    //Now we need to compute the gradient.
    for (int alpha = 0; alpha < gridsize; alpha++) {
        for (int beta = 0; beta < gridsize; beta++) {
            for (int j = 0; j < gridsize; j++) {
                radial = j * rStep;
                for (int i = 0; i < gridsize; i++) {
                    innersum = 0;
                    for (int k = 0; k < LMAX / 2; k++) {
                        innersum+=(LGrid[k][i] * FGrid[alpha][j][k] - LGrid[k][upper] * FGrid[alpha][pRad][k])* LegendreGrid[k][beta];
                    }
                    ProtonByProton = innersum * factor2 * MGrid[alpha] * MGrid[beta] * (-C + 1);
                    ElectronByProton= innersum * factor2 * MGrid[alpha] * MGrid[beta] * (C + me);
                    ProtonByElectron = innersum * factor2 * MGrid[alpha] * MGrid[beta] * (C+me);
                    ElectronByElectron = innersum * factor2 * MGrid[alpha] * MGrid[beta] * (-C + me*me);

                    if (i == beta && j == alpha) {
                        sumP += 2 *(radial*factor3)* (ProtonGrid[j][i] - UProtonGrid[j][i]) * (1 - ProtonByProton) + 2 * (ElectronGrid[j][i] - UElectronGrid[j][i]) * (-ElectronByProton);
                        sumE += 2 * (radial * factor3) *(ProtonGrid[j][i] - UProtonGrid[j][i]) * (-ProtonByElectron) + 2 * (ElectronGrid[j][i] - UElectronGrid[j][i]) * (1 - ElectronByElectron);
                    }
                    else {
                        sumP += 2 * (radial * factor3) *(ProtonGrid[j][i] - UProtonGrid[j][i]) * (-ProtonByProton) + 2 * (ElectronGrid[j][i] - UElectronGrid[j][i]) * (-ElectronByProton);
                        sumE += 2 * (radial * factor3) *(ProtonGrid[j][i] - UProtonGrid[j][i]) * (-ProtonByElectron) + 2 * (ElectronGrid[j][i] - UElectronGrid[j][i]) * (-ElectronByElectron);
                    }

                }
                
            }
            //this takes care of the proton terms in the gradient
            Gradient[alpha][beta] = sumP;
            sumP = 0;
            Gradient[alpha + gridsize-1][beta] = sumE;
            sumE = 0;

        }
    }
    double max = 0;
    logfile.open(gradstr);
    for (int j = 0; j < 2 * gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            if (abs(Gradient[j][i]) > max) {
                max = abs(Gradient[j][i]);
            }
            logfile << Gradient[j][i] << " ";
        }
        logfile << "\n";
    }
    logfile.close();
    logfile.clear();
    cout << max << "\n";

 
    double error = 0;
    for (int j = 0; j < gridsize; j++) {
        for (int i = 0; i < gridsize; i++) {
            error += j*rStep*factor3*(pow(ProtonGrid[j][i] - UProtonGrid[j][i], 2) + pow(ElectronGrid[j][i] - UElectronGrid[j][i], 2));
        }
    }
    cout << "error: " << error << "\n";

    logfile.open(updatestr);
    for (int j = 0;j < gridsize; j++) {
        logfile << j << ": ";
        for (int i = 0; i < gridsize; i++) {
            logfile << UElectronGrid[j][i] << " ";
        }
        logfile << "\n";

    }
    logfile.close();
    logfile.clear();


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

    

    /*
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
    */
   
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
    

    //step 6: Repeat
    int count = 1;

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
        updatestr[7] = Charcount;
        cout << "iteration: " << count << "\n";
        //update values
        
        for (int i = 0;i < gridsize; i++) {
            flage = 0;
            flagp = 0;
            for (int j = 0; j < gridsize; j++) {
                if (flage == 0){
                    ElectronGrid[j][i] -= DescentParameter * Gradient[j + gridsize][i];
                    if (ElectronGrid[j][i] <= 0) {
                        ElectronGrid[j][i] = 0;
                        flage = 1;
                        eRad = j;
                    }
                }
                else {
                    ElectronGrid[j][i] = 0;
                }
                if (flagp == 0) {
                    ProtonGrid[j][i] -= DescentParameter * Gradient[j][i];
                    if (ProtonGrid[j][i] < 0) {
                        ProtonGrid[j][i] = 0;
                        flagp = 1;
                        pRad = j;
                    }
                }
                else {
                    ProtonGrid[j][i] = 0;
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
            //cout <<j<<'\n';
            // j is radial coordinate
            for (int i = 0; i < gridsize; i++) {
                // i is theta coordinate
                for (int k = 0; k < LMAX / 2; k++) {
                    D_rEk = 0;
                    D_rMk = 0;
                    D_rGk = 0;
                    m = 0;
                    for (int m = 0; m < gridsize; m++) {
                        D_tEm0 = 0;
                        D_tGm0 = 0;
                        //D_tMm0=0;

                        for (int n = 0; n < gridsize; n++) {
                            D_tEm0 += MGrid[n] * LegendreGrid[k][n] * (ProtonGrid[m][n] - ElectronGrid[m][n]);
                            D_tGm0 += MGrid[n] * LegendreGrid[k][n] * (ProtonGrid[m][n] + me * ElectronGrid[m][n]);
                            //D_tMm0 += MGrid[n] * LegendreGrid[k][n] * (ProtonGrid[m][n] - ElectronGrid[m][n]) * OmegaGrid[m] * SinGrid[n];
                        }
                        D_rEk += MGrid[m] * FGrid[m][j][k] * D_tEm0;
                        //D_rMk += MGrid[m] * FGrid[m][j][k] * D_tMm0;
                        D_rGk += MGrid[m] * FGrid[m][j][k] * D_tGm0;
                    }
                    ElectricGrid[j][i] += factor * g * 4 * pi * C * D_rEk * LGrid[k][i];
                    //MagneticGrid[j][i] += factor * g  * mu0 * D_rMk * LGrid[k][i];
                    GravityGrid[j][i] += factor * g * 4 * pi * D_rGk * LGrid[k][i];
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
        

        /*
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
        */
        
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

        lambdap = ElectricGrid[pRad][upper] - GravityGrid[pRad][upper] + magnetOn * MagneticGrid[pRad][upper] * omega * rStep * pRad  - pow(omega, 2) * (1.0 / 2) * pow((rStep * pRad), 2);
        lambdae = -ElectricGrid[eRad][upper] - me * GravityGrid[eRad][upper] - magnetOn * MagneticGrid[eRad][upper] * omega * rStep * eRad  - me * pow(omega, 2) * (1.0 / 2) * pow((rStep * eRad), 2);

        for (int j = 0;j < gridsize; j++) {
            for (int i = 0;i < gridsize; i++) {
                HpGrid[j][i] = -ElectricGrid[j][i] + GravityGrid[j][i] - magnetOn * MagneticGrid[j][i] * omega * rStep * j * sin(thetaStep * i) + pow(omega, 2) * (1.0 / 2) * pow((rStep * j) * sin(thetaStep * i), 2) + lambdap;
                HeGrid[j][i] = ElectricGrid[j][i] + me * GravityGrid[j][i] + magnetOn * MagneticGrid[j][i] * omega * rStep * j * sin(thetaStep * i) + me * pow(omega, 2) * (1.0 / 2) * pow((rStep * j) * sin(thetaStep * i), 2) + lambdae;
            }
        }


        for (int j = 0;j < gridsize; j++) {
            for (int i = 0;i < gridsize; i++) {
                if (HpGrid[j][i] > 0) {
                    UProtonGrid[j][i] = pow((HpGrid[j][i] / HpGrid[0][0]), (3.0 / 2));
                }
                else {
                    UProtonGrid[j][i] = 0;
                }
                if (HeGrid[j][i] > 0) {

                    UElectronGrid[j][i] = pow((me * HeGrid[j][i] / HpGrid[0][0]), (3.0 / 2));
                }
                else {
                    UElectronGrid[j][i] = 0;
                }
            }
        }

        logfile.open(updatestr);
        for (int j = 0;j < gridsize; j++) {
            logfile << j << ": ";
            for (int i = 0; i < gridsize; i++) {
                logfile << UElectronGrid[j][i] << " ";
            }
            logfile << "\n";

        }
        logfile.close();
        logfile.clear();

        for (int alpha = 0; alpha < gridsize; alpha++) {
            for (int beta = 0; beta < gridsize; beta++) {
                for (int j = 0; j < gridsize; j++) {
                    radial = j * rStep;
                    for (int i = 0; i < gridsize; i++) {
                        innersum = 0;
                        for (int k = 0; k < LMAX / 2; k++) {
                            innersum += (LGrid[k][i] * FGrid[alpha][j][k] - LGrid[k][upper] * FGrid[alpha][pRad][k]) * LegendreGrid[k][beta];
                        }
                        ProtonByProton = innersum * factor2 * MGrid[alpha] * MGrid[beta] * (-C + 1);
                        ElectronByProton = innersum * factor2 * MGrid[alpha] * MGrid[beta] * (C + me);
                        ProtonByElectron = innersum * factor2 * MGrid[alpha] * MGrid[beta] * (C + me);
                        ElectronByElectron = innersum * factor2 * MGrid[alpha] * MGrid[beta] * (-C + me * me);

                        if (i == beta && j == alpha) {
                            sumP += 2 * (radial * factor3) * (ProtonGrid[j][i] - UProtonGrid[j][i]) * (1 - ProtonByProton) + 2 * (ElectronGrid[j][i] - UElectronGrid[j][i]) * (-ElectronByProton);
                            sumE += 2 * (radial * factor3) * (ProtonGrid[j][i] - UProtonGrid[j][i]) * (-ProtonByElectron) + 2 * (ElectronGrid[j][i] - UElectronGrid[j][i]) * (1 - ElectronByElectron);
                        }
                        else {
                            sumP += 2 * (radial * factor3) * (ProtonGrid[j][i] - UProtonGrid[j][i]) * (-ProtonByProton) + 2 * (ElectronGrid[j][i] - UElectronGrid[j][i]) * (-ElectronByProton);
                            sumE += 2 * (radial * factor3) * (ProtonGrid[j][i] - UProtonGrid[j][i]) * (-ProtonByElectron) + 2 * (ElectronGrid[j][i] - UElectronGrid[j][i]) * (-ElectronByElectron);
                        }

                    }

                }
                //this takes care of the proton terms in the gradient
                Gradient[alpha][beta] = sumP;
                sumP = 0;
                Gradient[alpha + gridsize - 1][beta] = sumE;
                sumE = 0;

            }
        }


        max = 0;
        for (int j = 0; j < 2 * gridsize; j++) {
            for (int i = 0; i < gridsize; i++) {
                if (abs(Gradient[j][i]) > max) {
                    max = abs(Gradient[j][i]);
                }
            }
        }
        cout << max << "\n";


        error = 0;
        for (int j = 0; j < gridsize; j++) {
            for (int i = 0; i < gridsize; i++) {
                error += j * rStep * factor3 * (pow(ProtonGrid[j][i] - UProtonGrid[j][i], 2) + pow(ElectronGrid[j][i] - UElectronGrid[j][i], 2));
            }
        }
        cout << "error: " << error << "\n";
        
    } while (error>TOL && count < 1);
}