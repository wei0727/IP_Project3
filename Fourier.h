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
Mat vectorToMat_enhanced(vector<vector<complex<double>>> &v) ;

vector< complex<double>> FFT_1D(vector< complex<double>> &m) ;
Mat_<complex<double>> FFT_1D(Mat_<complex<double>> &m) ;
vector< complex<double>> IFFT_1D(vector< complex<double>> &m) ;

vector<vector<complex<double>>> FFT_2D(vector<vector<complex<double>>> &m) ;
vector<vector<complex<double>>> IFFT_2D(vector<vector<complex<double>>> &m) ;

bitset<16> reverseBit(bitset<16> &b) ;
vector< complex<double>> rearrange(vector< complex<double>> &m) ;
//Mat_<complex<double>> rearrange(Mat_<complex<double>> &m) ;

Mat enhanceSpectrum(Mat &m) ;


