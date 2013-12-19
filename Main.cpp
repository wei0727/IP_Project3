#include "Fourier.h"

void test(){
	/*Mat t = (Mat_<double>(8, 8) << 0, 0, 0, 0, 0, 0, 0, 0,
								  0, 0, 70, 80, 90, 0, 0, 0,
								  0, 0, 90, 100, 110, 0, 0, 0,
								  0, 0, 110, 120, 130, 0, 0, 0,
								  0, 0, 130, 140, 150, 0, 0, 0,
								  0, 0, 0, 0, 0, 0, 0, 0,
								  0, 0, 0, 0, 0, 0, 0, 0,
								  0, 0, 0, 0, 0, 0, 0, 0) ;*/
	Mat t = (Mat_<double>(1, 8) << 0, 0, 2, 3, 4, 0, 0, 0) ;
	vector< vector< complex<double>>> v = matToVector(t) ;
	for(int i=0; i<8; i++)
		cout << v[0][i] << endl ;
	//for(int r=0; r<8; r++){
	//	for(int c=0; c<8; c++){
	//		cout << v[r][c].real() << ", " ;
	//	}
	//	cout << endl ;
	//}
	vector<complex<double>> f = FFT_1D(v[0]) ;
	for(int i=0; i<8; i++)
		cout << f[i] << endl ;
	//vector< vector< complex<double>>> f = FFT_2D(v) ;
	//for(int r=0; r<8; r++){
	//	for(int c=0; c<8; c++){
	//		cout << f[r][c].real() << ", " ;
	//	}
	//	cout << endl ;
	//}
}

int main(){
	//Mat test = (Mat_<int>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1) ;
	//Mat r = test.row(0) ;
	//cout << r << endl ;
	//Mat c = test.col(0) ;
	//cout << c << endl ;
	test() ;
	Mat img = imread("Fig0516(a)(applo17_boulder_noisy).tif", CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat img2 = centering(img) ;
	img2 = zeroPadding(img, 1) ;
	imshow("img", img) ;
	imshow("img2", img2) ;
	cvWaitKey(0) ;
	return 0 ;
}