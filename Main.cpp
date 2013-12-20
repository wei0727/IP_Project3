#include <ctime>
#include "Fourier.h"


void test(){
	clock_t t1, t2 ;
	//cout << (t2-t1)/(double)CLOCKS_PER_SEC << endl ;
	Mat t = (Mat_<double>(8, 8) << 0, 0, 0, 0, 0, 0, 0, 0,
								  0, 0, 70, 80, 90, 0, 0, 0,
								  0, 0, 90, 100, 110, 0, 0, 0,
								  0, 0, 110, 120, 130, 0, 0, 0,
								  0, 0, 130, 140, 150, 0, 0, 0,
								  0, 0, 0, 0, 0, 0, 0, 0,
								  0, 0, 0, 0, 0, 0, 0, 0,
								  0, 0, 0, 0, 0, 0, 0, 0) ;
	//IplImage *iplimg = cvLoadImage("rectangle.tif", CV_LOAD_IMAGE_GRAYSCALE) ;
	//Mat img(iplimg) ;
	//Mat m = centering(img) ;
	//vector<vector<complex<double>>> vImg = matToVector(m) ;
	//t1 = clock() ;
	//setTable() ;
	//t2 = clock() ;
	//cout << (t2-t1)/(double)CLOCKS_PER_SEC << endl ;
	//t1 = clock() ;
	//vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	//t2 = clock() ;
	//cout << (t2-t1)/(double)CLOCKS_PER_SEC << endl ;
	vector<complex<double>> v(4) ;
	v[0] = complex<double>(0, 0) ;
	v[1] = complex<double>(20, 0) ;
	v[2] = complex<double>(10, 0) ;
	v[3] = complex<double>(30, 0) ;
	v = FFT_1D(v) ;
	cout << v[0] << endl ;
	cout << v[1] << endl ;
	cout << v[2] << endl ;
	cout << v[3] << endl ;
}

void FFT_Spectrum(string imgName="Fig0424(a)(rectangle).tif"){
	IplImage *iplimg = cvLoadImage(imgName.c_str(), CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat img(iplimg) ;
	Mat m = centering(img) ;
	vector<vector<complex<double>>> vImg = matToVector(m) ;
	vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	//Mat spectrum = vectorToMat(vFFT) ;
	Mat spectrum_enhanced = vectorToMat_enhanced(vFFT) ;
	cvShowImage("enhanced spectrum", &IplImage(spectrum_enhanced)) ;
	cvWaitKey(0) ;
}

void FFT_HFEF_HIST(string imgName="Fig0459(a)(orig_chest_xray).tif"){
	IplImage *iplimg = cvLoadImage(imgName.c_str(), CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat img(iplimg) ;
	Mat m = centering(img) ;
	m = zeroPadding(m, 1) ;
	//m = zeroPadding(m, 0) ;
	vector<vector<complex<double>>> vImg = matToVector(m) ;
	vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	vector<vector<complex<double>>> vHFFT = vFFT ;
	vFFT = HFEF_Gaussian(vFFT, 40) ;
	vHFFT = HPF_Gaussian(vHFFT, 40) ;
	vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT) ;
	vector<vector<complex<double>>> vIHFFT = IFFT_2D(vHFFT) ;
	Mat ifImg = vectorToMat(vIFFT) ;
	Mat hpfImg = vectorToMat(vIHFFT) ;
	Mat hpfResult = hpfImg(cvRect(0, 0, img.cols, img.rows)) ;
	Mat result = ifImg(cvRect(0, 0, img.cols, img.rows)) ;
	Mat HFEFMat ;
	result.copyTo(HFEFMat) ;
	equalizeHist(result, result) ;
	equalizeHist(img, img) ;
	cvShowImage("b", &IplImage(hpfResult)) ;
	cvShowImage("c", &IplImage(HFEFMat)) ;
	cvShowImage("d", &IplImage(result)) ;
	cvShowImage("Only histogram equalization", &IplImage(img)) ;
	cvWaitKey(0) ;
}

int main(){
	//test() ;
	setTable() ;
	//FFT_Spectrum() ;
	FFT_HFEF_HIST() ;
	//IplImage *iplimg = cvLoadImage("Fig0333(a)(test_pattern_blurring_orig).tif", CV_LOAD_IMAGE_GRAYSCALE) ;
	//Mat img(iplimg) ;
	//Mat m = centering(img) ;
	//m = zeroPadding(m, 1) ;
	//vector<vector<complex<double>>> vImg = matToVector(m) ;
	//vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	//Mat spectrum = vectorToMat(vFFT) ;
	////Mat spectrum_enhanced = vectorToMat_enhanced(vFFT) ;
	//cvShowImage("enhanced spectrum", &IplImage(spectrum)) ;
	//cvWaitKey(0) ;
	return 0 ;
}