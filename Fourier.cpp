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

Mat zeroPadding(Mat &m, int type){
	Mat img ;
	//col*2 row*2
	if(type==0){
		img = Mat::zeros(m.rows*2, m.cols*2, m.type()) ;
	}
	//padding to 2^k
	else{
		int num = 1 ;
		while(num < m.rows || num <m.cols){
			num = num<<1 ;
		}
		img = Mat::zeros(num, num, m.type()) ;
	}
	m.copyTo(img(cvRect(0, 0, m.cols, m.rows))) ;
	return img ;
}

vector< vector< complex<double>>> matToVector(Mat &m){
	vector< vector< complex<double>>> v(m.rows, vector< complex<double>>(m.cols)) ;
	double *p ;
	for(int r=0; r<m.rows; r++){
		p = m.ptr<double>(r) ;
		for(int c=0; c<m.cols; c++){
			v[r][c] = complex<double>(p[c], 0) ;
			//v[r][c] = 0 ;
		}
	}
	return v ;
}

bitset<32> reverseBit(bitset<32> &b){
	string bs = b.to_string() ;
	reverse(bs.begin(), bs.end()) ;
	return bitset<32>(bs) ;
}

vector<complex<double>> FFT_1D(vector<complex<double>> &m){
	vector<complex<double>> v = rearrange(m) ;
	vector<complex<double>> f = v ;
	for(int n=2; n!=v.size(); n=n<<1){
		vector<complex<double>> tmp(v.size()) ;
		int nGroup = v.size()/n ;
		int k = n >> 1 ;
		for(int index=0; index<nGroup; index++){
			for(int j=0; j<k; j++){
				int u = index*n + j ;
				double angle = CV_PI*j/k ;
				complex<double> w(cos(angle), -sin(angle)) ;
				tmp[u] = f[u] + f[u+k]*w ;
				tmp[u+k] = f[u] - f[u+k]*w ;
			}
		}
		f = tmp ;
	}
	return f ;
}

vector<complex<double>> rearrange(vector<complex<double>> &m){
	vector<complex<double>> v(m.size()) ;
	for(int i=0; i<m.size(); i++){
		bitset<32> b(i) ;
		b = reverseBit(b) ;
		int j = b.to_ulong() ;
		v[j] = m[i] ;
	}
	return v ;
}