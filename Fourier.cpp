#include "Fourier.h"

Mat centering(Mat &m){
	Mat img ;
	m.convertTo(img, CV_64FC1) ;
	double *p ;
	for(int r=0; r<m.rows; r++){
		p = img.ptr<double>(r) ;
		for(int c=0; c<m.cols; c++){
			//r+c odd->-1 even->1
			//img.at<unsigned char>(r, c) *= ((r+c)&1==1? -1: 1) ;
			//cout << (double)p[c] << ", " ;
			p[c] *= ((r+c)&1==1? -1: 1) ;
			//cout << (double)p[c] << endl ;
		}
	}
	return img ;
}

Mat zeroPadding(Mat &m){
	Mat img = Mat::zeros(m.rows*2, m.cols*2, m.type()) ;
	m.copyTo(img(cvRect(0, 0, m.cols, m.rows))) ;
	return img ;
}