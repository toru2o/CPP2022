#pragma warning( disable:4996 )
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
//    Program ising
// 2021.11.29 random_device を使って
//=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

#include "cip.h"
#include <opencv2/opencv.hpp>
#include<iostream>
#include "Eigen/Core"
#include <random>
//#include "Eigen/Dense"



using namespace cv;
using namespace std;
using namespace Eigen;

constexpr double PI = 3.14159265358979323846;
constexpr int width = 1000;
constexpr int height = 1000;



//---- definition of fundermental constants
//---- physical constants
const double J = 1.0;           // bonding constant
const double H = 0.0;           // external magnetic field
const int SIZE = 16;            // size of lattice

//---- experimental settings
// settings about temperature
const double T_INIT = 0.25;          // initial temperature
const double T_FINA = 5.00;          // final   temperature
const double T_STEP = 0.25;          // increments

// settings about sweep the lattice
const int RELAX = 400;           // steps till relaxation
const int INTERVAL = 10;            // steps per one measure

// settings about number of data to measure fluctuations of ...
const int SAMPLES = 40;            // energy, magnetization
const int GROUPS = 20;            // susceptibility, specific heat

// file names of output data
const char* fnames[4] = {
  "energy.dat", "magnet.dat", "suscep.dat", "spheat.dat"
};
// titles of each quantities
const char* titles[4] = {
  "Energy", "Magnetization", "Susceptibility", "Specific heat"
};
// file pointers
static FILE* fptrs[4];

//---- index indicators
const int ENERGY = 0;
const int MAGNET = 1;
const int SUSCEP = 2;
const int SPHEAT = 3;

//---- graphics settings
const int    WIN_WIDTH = 256;
const int    WIN_HEIGHT = 256;

random_device rd;
mt19937 mt(rd());
//浮動小数点の一様分布 
uniform_real_distribution<double> rand_dist(0.0, 1.0);


//---- declaration and definition of class Statistic
class Statistic
{
    int n;                // current total of data
    double sum1, sum2;    // current sum of data and its square
public:
    inline Statistic(void) {
        n = 0;
        sum1 = sum2 = 0.0;
    }
    //-- adds a new datum
    inline void Add(double datum) {
        n++;
        sum1 += datum;
        sum2 += datum * datum;
    }
    //-- current mean value
    inline double Mean(void) {
        return sum1 / n;
    }
    //-- current variance value
    inline double Variance(void) {
        return (sum2 - sqr(sum1) / n) / (n - 1);
    }
    //-- current deviation value
    inline double Deviation(void) {
        return sqrt(Variance());
    }
    //-- current deviation of mean value, namely error value
    inline double Error(void) {
        return sqrt(Variance() / n);
    }
};


//---- declaration and difinition of class Spin
class Spin
{
    int spin;      // the value of spin, it takes +1 or -1
public:
    //-- sets this spin randomly +1 or -1
    void RandSet(void) {
        //if (Urand() > 0.5) spin = +1; else spin = -1;
        if (rand_dist(mt) > 0.5) spin = +1; else spin = -1;
    }
    //-- represents the value of spin when this object is called
    operator int() {
        return spin;
    }
    //-- transits this spin by given probability
    void Transit(int s_sum, double prob) {
        //if (Urand() < prob) spin = +1; else spin = -1;
        if (rand_dist(mt) < prob) spin = +1; else spin = -1;
    }
};


//---- declaration of class Lattice
class Lattice
{
    Spin spin[SIZE][SIZE];    // spins in this lattice
    double prob[5];           // table of transition probability
public:
    Lattice(double T);
    void Sweep(int interval);
    void Measure(Statistic& Energy, Statistic& Magnet);
    void Draw(cv::Mat& image);
private:
    inline void Neighbor(int& im, int& ip, int i);
};

//---- makes transition propability table and let spins randomly turn
Lattice::Lattice(double T)
{
    for (int s_sum = -4; s_sum <= +4; s_sum += 2) {
        double eplus = exp((J * s_sum + H) / T);
        prob[s_sum / 2 + 2] = eplus / (eplus + 1.0 / eplus);
    }
    for (int ix = 0; ix < SIZE; ix++) {
        for (int iy = 0; iy < SIZE; iy++) {
            spin[ix][iy].RandSet();
        }
    }
}


//---- gets neighbors' indices of given index
inline void Lattice::Neighbor(int& im, int& ip, int i)
{
    if (i == 0) im = SIZE - 1; else im = i - 1;
    if (i == SIZE - 1) ip = 0; else ip = i + 1;
}

//---- sweeps this lattice [interval] times
void Lattice::Sweep(int interval)
{
    int ix, iy, ixm, ixp, iym, iyp;
    int s_sum;

    while (interval--) {
        for (ix = 0; ix < SIZE; ix++) {
            Neighbor(ixm, ixp, ix);
            for (iy = 0; iy < SIZE; iy++) {
                Neighbor(iym, iyp, iy);

                s_sum = spin[ixp][iy] + spin[ixm][iy]
                    + spin[ix][iyp] + spin[ix][iym];

                spin[ix][iy].Transit(s_sum, prob[s_sum / 2 + 2]);
            }
        }
    }
}

//---- evaluates energy and magnetization of this lattice
void Lattice::Measure(Statistic& Energy, Statistic& Magnet)
{
    int ix, iy, ixm, ixp, iym, iyp;
    int s_sum = 0, ss_sum = 0;

    for (ix = 0; ix < SIZE; ix++) {
        Neighbor(ixm, ixp, ix);
        for (iy = 0; iy < SIZE; iy++) {
            Neighbor(iym, iyp, iy);

            s_sum += spin[ix][iy];
            ss_sum += spin[ix][iy] * (spin[ixm][iy] + spin[ix][iym]);
        }
    }

    const double energy = -J * ss_sum - H * s_sum;
    const double magnet = (double)s_sum;

    Energy.Add(energy);
    Magnet.Add(magnet);
}

//---- draws spins as colored square
void Lattice::Draw(cv::Mat& image)
{
    const int dX = WIN_WIDTH / SIZE, dY = WIN_HEIGHT / SIZE;
    double mp = 2.0;

    for (int ix = 0; ix < SIZE; ix++) {
        for (int iy = 0; iy < SIZE; iy++) {
            int B, G, R;
            if (spin[ix][iy] == +1) {
                B = 0;
                G = 255;
                R = 0;
            }
            else {
                B = 255;
                G = 0;
                R = 0;
            }
            
            int x1, y1, x2, y2;
            x1 = (int)(ix * dX * mp);
            y1 = (int)(iy * dY * mp);
            x2 = (int)((ix * dX + dX - 1) * mp);
            y2 = (int)((iy * dY + dY - 1) * mp);
            rectangle(image, Point(x1, y1), Point(x2, y2), Scalar(B, G, R), -1);
        }
    }
    cv::imshow("ising", image);
    cv::waitKey(100);
}


//---- opens four files to save data
void OpenFiles(void)
{
    for (int i = 0; i < 4; i++) {
        fptrs[i] = fopen(fnames[i], "wt");
        if (fptrs[i] == NULL) {
            perror("Error at fopen.");
            exit(1);
        }
    }
}

//---- prints and saves each quantities in the four files
void Output(double T, Statistic quantity[])
{
    printf("\nTemperature    = % f\n", T);

    for (int i = 0; i < 4; i++) {
        printf("%-14s = % f +/- %f\n", titles[i], quantity[i].Mean(), quantity[i].Error());
        fprintf(fptrs[i], "%f % f %f\n", T, quantity[i].Mean(), quantity[i].Error());
    }
}

//---- closes the files
void CloseFiles(void)
{
    for (int i = 0; i < 4; i++) {
        fclose(fptrs[i]);
    }
}

//---- executes experiment at given temperature
void Experiment(cv::Mat& image, double T)
{
    Lattice   latt(T);
    Statistic quantity[4];

    latt.Sweep(RELAX);     // makes the lattice equilibrium

    for (int group = 0; group < GROUPS; group++) {
        Statistic Energy0, Magnet0;

        // measures fluctuations of energy and magnetization
        for (int sample = 0; sample < SAMPLES; sample++) {
            latt.Sweep(INTERVAL);
            latt.Measure(Energy0, Magnet0);
        }
        // calculates other data per spin statistically
        quantity[ENERGY].Add(Energy0.Mean() / (SIZE * SIZE));
        quantity[MAGNET].Add(Magnet0.Mean() / (SIZE * SIZE));
        quantity[SUSCEP].Add(Magnet0.Variance() / T / (SIZE * SIZE));
        quantity[SPHEAT].Add(Energy0.Variance() / sqr(T) / (SIZE * SIZE));

        latt.Draw(image);
    }

    Output(T, quantity); // saves data in files
}

void gnuplot() {
    FILE* gp = _popen("C:/PROGRA~1/gnuplot/bin/gnuplot.exe -persist", "w");
    if (gp == NULL) return;
    fputs("load 'D:/VCPP/VC++2022/計算物理のためのCC++言語入門/ising/ising/gChart_2D.txt'\n", gp);
    fflush(gp);
    _pclose(gp);
}


int main(void)
{
    //mt19937 mt(rd());
    cv::Mat image(height, width, CV_8UC3);
    OpenFiles();

    int cnt = 0;
    for (double T = T_INIT; T <= T_FINA; T += T_STEP) {  // temperature loop
        Experiment(image, T);
        cnt++;
        if (cnt > 30) break;

    }

    CloseFiles();
    gnuplot();

    return 0;
}