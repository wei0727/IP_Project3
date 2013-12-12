#include "Fourier.h"

int main(){
	//Mat test = (Mat_<int>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1) ;
	//Mat r = test.row(0) ;
	//cout << r << endl ;
	Mat img = imread("Fig0516(a)(applo17_boulder_noisy).tif", CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat img2 = centering(img) ;
	img2 = zeroPadding(img) ;
	imshow("img", img) ;
	imshow("img2", img2) ;
	cvWaitKey(0) ;
	return 0 ;
}