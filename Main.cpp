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

int main(){
	//test() ;
	setTable() ;
	IplImage *iplimg = cvLoadImage("Fig0459(a)(orig_chest_xray).tif", CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat img(iplimg) ;
	//cvShowImage("input", iplimg) ;
	Mat m = centering(img) ;
	m = zeroPadding(m, 0) ;
	m = zeroPadding(m, 1) ;
	vector<vector<complex<double>>> vImg = matToVector(m) ;
	vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	//vFFT = LPF_Gaussian(vFFT, 10) ;
	vFFT = HFEF_Gaussian(vFFT, 40) ;
	//vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT) ;
	Mat fImg = vectorToMat(vFFT) ;
	//Mat fImg = vectorToMat_enhanced(vFFT) ;
	//fImg = LPF_Gaussian(fImg, 10) ;
	//fImg = HPF_Gaussian(fImg, 20) ;
	//fImg = HFEF_Gaussian(fImg, 40) ;
	//fImg = enhanceSpectrum(fImg) ;
	//vector<vector<complex<double>>> vIFFT = FFT_2D(vFFT) ;
	vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT) ;
	Mat ifImg = vectorToMat(vIFFT) ;
	//Mat reconstruct = ifImg(cvRect(ifImg.cols-img.cols-1, ifImg.rows-img.rows-1, img.cols, img.rows)) ;
	//reconstruct = reverseMat(reconstruct) ;
	//equalizeHist(reconstruct, reconstruct) ;
	IplImage *iplfImg = cvCloneImage(&IplImage(fImg)) ;
	IplImage *iplifImg = cvCloneImage(&IplImage(ifImg)) ;
	cvSaveImage("fftSpectrum.jpg", iplfImg) ;
	cvSaveImage("ifftSpectrum.jpg", iplifImg) ;
	cvShowImage("output", iplfImg) ;
	cvShowImage("output2", iplifImg) ;
	cvWaitKey(0) ;
	return 0 ;
}