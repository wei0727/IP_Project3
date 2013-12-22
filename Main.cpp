 #include <ctime>
#include "Fourier.h"


void test(){
	clock_t t1, t2 ;
	//Mat t = (Mat_<double>(4, 4) << 0, 0, 0, 0,
	//							   1, 2, 3, 4,
	//							   2, 3, 4, 5,
	//							   3, 4, 5, 6) ;
	//vector<vector<complex<double>>> vImg = matToVector(t) ;
	//vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	//vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT) ;
	//cout << fixed ;
	//for(int r=0; r<vFFT.size(); r++){
	//	for(int c=0; c<vFFT[0].size(); c++)
	//		cout << vIFFT[r][c] << "\t" ;
	//	cout << endl ;
	//}
	//vector<complex<double>> v ;
	//v.push_back(complex<double>(1, 0)) ;
	//v.push_back(complex<double>(2, 0)) ;
	//v.push_back(complex<double>(3, 0)) ;
	//v.push_back(complex<double>(4, 0)) ;
	//v = FFT_1D(v) ;
	//v[0] /= 4 ;
	//v[1] /= 4 ;
	//v[2] /= 4 ;
	//v[3] /= 4 ;
	//cout << v[0] << endl ;
	//cout << v[1] << endl ;
	//cout << v[2] << endl ;
	//cout << v[3] << endl ;
	//
	//IplImage *iplimg = cvLoadImage("Fig0516(a)(applo17_boulder_noisy).tif", CV_LOAD_IMAGE_GRAYSCALE) ;
	//Mat img(iplimg) ;
	//Mat m = centering(img) ;
	//m = zeroPadding(m, 0) ;
	//m = zeroPadding(m, 1) ;
	//vector<vector<complex<double>>> vImg = matToVector(m) ;
	//vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	//vector<vector<complex<double>>> vFFT_resize = specResize(vFFT, img.rows, img.cols) ;
	//vector<vector<complex<double>>> vFFT_resize2 = specResize(vFFT_resize, m.rows, m.cols) ;
	//vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT_resize2) ;
	//Mat result = vectorToMat(vIFFT) ;
	//cvShowImage("d", &IplImage(result)) ;
	//cvWaitKey(0) ;
}

void FFT_Spectrum(string imgName="Fig0424(a)(rectangle).tif"){
	IplImage *iplimg = cvLoadImage(imgName.c_str(), CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat img(iplimg) ;
	Mat m = centering(img) ;
	vector<vector<complex<double>>> vImg = matToVector(m) ;
	vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	Mat spectrum = vectorToMat(vFFT) ;
	//Mat pSpec = vectorToMat_powerSpec(vFFT) ;
	Mat spectrum_enhanced = vectorToMat_enhanced(vFFT) ;
	cvShowImage("enhanced spectrum", &IplImage(spectrum_enhanced)) ;
	cvShowImage("spectrum", &IplImage(spectrum)) ;
	cvSaveImage("1_Spectrum.jpg", &IplImage(spectrum)) ;
	cvSaveImage("1_EnhancedSpectrum.jpg", &IplImage(spectrum_enhanced)) ;
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
	Mat spec = vectorToMat(vFFT) ;
	vector<vector<complex<double>>> vHFFT = vFFT ;
	vFFT = HFEF_Gaussian(vFFT, 40) ;
	vHFFT = HPF_Gaussian(vHFFT, 40) ;
	vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT) ;
	vector<vector<complex<double>>> vIHFFT = IFFT_2D(vHFFT) ;
	Mat ifImg = vectorToMat_real(vIFFT) ;
	Mat hpfImg = vectorToMat_real(vIHFFT) ;
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
	//cvShowImage("spec", &IplImage(spec)) ;
	cvSaveImage("2_b.jpg", &IplImage(hpfResult)) ;
	cvSaveImage("2_c.jpg", &IplImage(HFEFMat)) ;
	cvSaveImage("2_d.jpg", &IplImage(result)) ;
	cvSaveImage("2_equalization.jpg", &IplImage(img)) ;
	cvWaitKey(0) ;
}

void FFT_ButterworthReject(string name="Fig0516(a)(applo17_boulder_noisy).tif"){
	IplImage *iplimg = cvLoadImage("Fig0516(a)(applo17_boulder_noisy).tif", CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat img(iplimg) ;
	Mat m = centering(img) ;
	m = zeroPadding(m, 0) ;
	m = zeroPadding(m, 1) ;
	vector<vector<complex<double>>> vImg = matToVector(m) ;
	vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	vector<vector<complex<double>>> vFFT_noise = FFT_2D(vImg) ;
	Mat spectrum = vectorToMat(vFFT) ;
	Mat specResize(img.rows, img.cols, img.type()) ;
	resize(spectrum, specResize, cvSize(img.cols, img.rows)) ;
	//uchar *p ;
	//int cr = spectrum.rows/2 ;
	//int cc = spectrum.cols/2 ;
	//for(int r=0; r<m.rows; r++){
	//	p = spectrum.ptr<uchar>(r) ;
	//	for(int c=0; c<m.cols; c++){
	//		if(p[c] > 50)
	//			cout << r-cr << ", " << c-cc << ": " << sqrt((r-cr)*(r-cr) + (c-cc)*(c-cc)) << endl ;
	//	}
	//}

	//vFFT = Butterworth_reject(vFFT, 766, 8, 5) ;
	//vFFT = Butterworth_reject(vFFT, 676, 8, 5) ;
	//vFFT = Butterworth_reject(vFFT, 575, 8, 5) ;
	//vFFT_noise = Butterworth_pass(vFFT_noise, 766, 8, 5) ;
	//vFFT_noise = Butterworth_pass(vFFT_noise, 676, 8, 5) ;
	//vFFT_noise = Butterworth_pass(vFFT_noise, 575, 8, 5) ;

	vFFT = Butterworth_reject_byimg(vFFT, "Fig0516(c)(BW_banreject_order4).tif") ;
	vFFT_noise = Butterworth_pass_byimg(vFFT_noise, "Fig0516(c)(BW_banreject_order4).tif") ;

	Mat spectrum_filtered = vectorToMat(vFFT) ;
	Mat specFilterdResize(img.rows, img.cols, img.type()) ;
	resize(spectrum_filtered, specFilterdResize, cvSize(img.cols, img.rows)) ;
	vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT) ;
	vector<vector<complex<double>>> vIFFT_noise = IFFT_2D(vFFT_noise) ;
	//vector<vector<complex<double>>> vIFFT = FFT_2D(vFFT) ;
	Mat ifImg = vectorToMat_real(vIFFT) ;
	Mat result = ifImg(cvRect(0, 0, img.cols, img.rows)) ;
	Mat noise = vectorToMat(vIFFT_noise) ;
	Mat result_noise = noise(cvRect(0, 0, img.cols, img.rows)) ;
	cvShowImage("b", &IplImage(specResize)) ;
	cvShowImage("spectrum_filtered", &IplImage(specFilterdResize)) ;
	cvShowImage("d", &IplImage(result)) ;
	cvShowImage("noise", &IplImage(result_noise)) ;
	cvSaveImage("3_butter.jpg", &IplImage(result)) ;
	cvSaveImage("3_spectrum.jpg", &IplImage(specResize)) ;
	cvSaveImage("3_spectrum_filtered.jpg", &IplImage(specFilterdResize)) ;
	cvSaveImage("3_noise.jpg", &IplImage(result_noise)) ;
	cvWaitKey(0) ;
}

int main(){
	setTable() ;
	//test() ;
	//FFT_Spectrum() ;
	//FFT_HFEF_HIST() ;
	FFT_ButterworthReject() ;
	return 0 ;
}