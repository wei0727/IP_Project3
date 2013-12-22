#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <string>
#include <complex>
#include <climits>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;

void setTable() ;

bool mCmp(int a, int b) ;

Mat centering(Mat &m) ;
Mat zeroPadding(Mat &m, int type=0) ;

Mat_<complex<double>> convertToCoplex(Mat &m) ;
vector< vector< complex<double>>> matToVector(Mat &m) ;
Mat vectorToMat(vector<vector<complex<double>>> &v) ;
Mat vectorToMat_powerSpec(vector<vector<complex<double>>> &v) ;
Mat vectorToMat_enhanced(vector<vector<complex<double>>> &v) ;
Mat vectorToMat_real(vector<vector<complex<double>>> &v) ;

vector< complex<double>> FFT_1D(vector< complex<double>> &m) ;
vector< complex<double>> IFFT_1D(vector< complex<double>> &m) ;

vector<vector<complex<double>>> FFT_2D(vector<vector<complex<double>>> &m) ;
vector<vector<complex<double>>> IFFT_2D(vector<vector<complex<double>>> &m) ;

bitset<16> reverseBit(bitset<16> &b) ;
vector< complex<double>> rearrange(vector< complex<double>> &m) ;

//difference (uchar)
Mat difMat(Mat &m1, Mat &m2) ;
//Multiply
//m1->img(uchar), m2->filter(double)
Mat mulMat(Mat &m1, Mat &m2) ;
vector<vector<complex<double>>> mulVec(vector<vector<complex<double>>> &v1, vector<vector<complex<double>>> &v2) ;

//low pass gaussian
vector<vector<complex<double>>> BLPF(vector<vector<complex<double>>> &m, double d0=1) ;
vector<vector<complex<double>>> LPF_Gaussian(vector<vector<complex<double>>> &m, double d0=1) ;
//hight pass gaussian
vector<vector<complex<double>>> HPF_Gaussian(vector<vector<complex<double>>> &m, double d0=1) ;
//hight frequency emphasis filter(gausian)
vector<vector<complex<double>>> HFEF_Gaussian(vector<vector<complex<double>>> &m, double d0=1) ;

vector<vector<complex<double>>> Butterworth_reject(vector<vector<complex<double>>> &m, double d0=1, double w=64, double n=1) ;
vector<vector<complex<double>>> Butterworth_pass(vector<vector<complex<double>>> &m, double d0=1, double w=64, double n=1) ;
vector<vector<complex<double>>> Butterworth_reject_byimg(vector<vector<complex<double>>> &m, string name) ;
vector<vector<complex<double>>> Butterworth_pass_byimg(vector<vector<complex<double>>> &m, string name) ;

Mat reverseMat(Mat &m) ;


