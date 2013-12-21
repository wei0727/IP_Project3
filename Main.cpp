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
	vector<vector<complex<double>>> vImg = matToVector(t) ;
	vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
	vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT) ;
	vector<vector<complex<double>>> vt = vImg ;
	for(int r=0; r<8; r++){
		for(int c=0; c<8; c++){
			double d = vIFFT[r][c].real() ;
			d *= ((r+c)&1==1? -1: 1) ;
			d *= 64 ;
			//cout << abs(vIFFT[r][c]) << "  " ;
			cout << vImg[r][c].real() << "\t" ;
		}
		cout << endl ;
	}
	cout << endl ;
	for(int r=0; r<8; r++){
		for(int c=0; c<8; c++){
			double d = vIFFT[r][c].real() ;
			d *= ((r+c)&1==1? -1: 1) ;
			//d *= 64 ;
			//d = vIFFT[r][c].real() ;
			d = abs(vIFFT[r][c]) ;
			vt[7-r][7-c] = d ;
			if(d<0.01)
				d = 0 ;
			//cout << abs(vIFFT[r][c]) << "  " ;
			cout << d << "\t" ;
		}
		cout << endl ;
	}
	cout << endl ;
	for(int r=0; r<8; r++){
		for(int c=0; c<8; c++){
			double d = vt[r][c].real() ; 
			d = abs(vt[r][c]) ;
			if(d<0.01)
				d = 0 ;
			//cout << abs(vIFFT[r][c]) << "  " ;
			cout << d << "\t" ;
		}
		cout << endl ;
	}
	
	
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
	//Mat ifImg = vectorToMat_real(vIFFT) ;
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

void FFT_ButterworthReject(string name="Fig0516(a)(applo17_boulder_noisy).tif"){
	IplImage *iplimg = cvLoadImage("Fig0516(a)(applo17_boulder_noisy).tif", CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat img(iplimg) ;
	Mat m = centering(img) ;
	m = zeroPadding(m, 0) ;
	m = zeroPadding(m, 1) ;
	vector<vector<complex<double>>> vImg = matToVector(m) ;
	vector<vector<complex<double>>> vFFT = FFT_2D(vImg) ;
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
	vFFT = Butterworth_reject(vFFT, 766, 8, 5) ;
	vFFT = Butterworth_reject(vFFT, 676, 8, 5) ;
	vFFT = Butterworth_reject(vFFT, 575, 8, 5) ;
	Mat spectrum_filtered = vectorToMat(vFFT) ;
	Mat specFilterdResize(img.rows, img.cols, img.type()) ;
	resize(spectrum_filtered, specFilterdResize, cvSize(img.cols, img.rows)) ;
	vector<vector<complex<double>>> vIFFT = IFFT_2D(vFFT) ;
	//vector<vector<complex<double>>> vIFFT = FFT_2D(vFFT) ;
	Mat ifImg = vectorToMat(vIFFT) ;
	Mat result = ifImg(cvRect(0, 0, img.cols, img.rows)) ;
	cvShowImage("b", &IplImage(specResize)) ;
	cvShowImage("spectrum_filtered", &IplImage(specFilterdResize)) ;
	cvShowImage("d", &IplImage(result)) ;
	//cvSaveImage("butter.jpg", &IplImage(result)) ;
	//cvSaveImage("spectrum.jpg", &IplImage(specResize)) ;
	cvWaitKey(0) ;
}

int main(){
	setTable() ;
	test() ;
	//FFT_Spectrum() ;
	//FFT_HFEF_HIST() ;
	//FFT_ButterworthReject() ;
	return 0 ;
}