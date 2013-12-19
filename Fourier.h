#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <string>
#include <complex>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;

bool mCmp(int a, int b) ;

Mat centering(Mat &m) ;
Mat zeroPadding(Mat &m, int type=0) ;

Mat_<complex<double>> convertToCoplex(Mat &m) ;
vector< vector< complex<double>>> matToVector(Mat &m) ;

vector< complex<double>> FFT_1D(vector< complex<double>> &m) ;
Mat_<complex<double>> FFT_1D(Mat_<complex<double>> &m) ;

vector<vector<complex<double>>> FFT_2D(vector<vector<complex<double>>> &m) ;

bitset<32> reverseBit(bitset<32> &b) ;
vector< complex<double>> rearrange(vector< complex<double>> &m) ;
//Mat_<complex<double>> rearrange(Mat_<complex<double>> &m) ;
