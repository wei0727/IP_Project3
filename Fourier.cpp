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

Mat_<complex<double>> convertToCoplex(Mat &m){
	Mat_<complex<double>> cm(m.rows, m.cols) ;
	double *p ;
	complex<double> *cp ;
	for(int r=0; r<m.rows; r++){
		p = m.ptr<double>(r) ;
		cp = cm.ptr<complex<double>>(r) ;
		for(int c=0; c<m.cols; c++){
			cp[c] = complex<double>(p[c], 0) ;
			//v[r][c] = complex<double>(p[c], 0) ;
			//v[r][c] = 0 ;
		}
	}
	return cm ;
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

//Mat_<complex<double>> FFT_1D(Mat_<complex<double>> &m){
//	Mat_<complex<double>> v = rearrange(m) ;
//	Mat_<complex<double>> f(v) ;
//	for(int n=2; n!=v.cols; n=n<<1){
//		Mat_<complex<double>> tmp(1, f.cols) ;
//		complex<double> *tmpPtr = tmp.ptr<complex<double>>(0) ;
//		complex<double> *fPtr = f.ptr<complex<double>>(0) ;
//		int nGroup = v.cols/n ;
//		int k = n >> 1 ;
//		for(int index=0; index<nGroup; index++){
//			for(int j=0; j<k; j++){
//				int u = index*n + j ;
//				double angle = CV_PI*j/k ;
//				complex<double> w(cos(angle), -sin(angle)) ;
//				tmpPtr[u] = fPtr[u] + fPtr[u+k]*w ;
//				tmpPtr[u+k] = fPtr[u] - fPtr[u+k]*w ;
//			}
//		}
//		//f = tmp ;
//		tmp.copyTo(f) ;
//	}
//	return f ;
//}

vector<vector<complex<double>>> FFT_2D(vector<vector<complex<double>>> &m){
	vector<vector<complex<double>>> f(m) ;
	//rows
	for(int r=0; r<m.size(); r++){
		f[r] = FFT_1D(m[r]) ;
	}
	//cols
	for(int c=0; c<m[0].size(); c++){
		vector<complex<double>> col(m.size()) ;
		//get col vector
		for(int i=0; i<m.size(); i++)
			col[i] = f[i][c] ;
		col = FFT_1D(col) ;
		//set col vector
		for(int i=0; i<m.size(); i++)
			f[i][c] = col[i] ;
	}
	return f ;
}

bool mCmp(int a, int b){
	bitset<32> ba(a) ;
	bitset<32> bb(b) ;
	ba = reverseBit(ba) ;
	bb = reverseBit(bb) ;
	return ba.to_ulong() < bb.to_ulong() ;
}

vector<complex<double>> rearrange(vector<complex<double>> &m){
	vector< complex<double>> v(m.size()) ;
	vector<int> tmp(m.size()) ;
	for(int i=0; i<m.size(); i++){
		tmp[i] = i ;
	}
	sort(tmp.begin(), tmp.end(), mCmp) ;
	for(int i=0; i<m.size(); i++){
		v[i] = m[tmp[i]] ;
	}
	return v ;
}

//Mat_<complex<double>> rearrange(Mat_<complex<double>> &m){
//	Mat_<complex<double>> v(1, m.cols) ;
//	complex<double> *p = v.ptr<complex<double>>(0) ;
//	complex<double> *mp = m.ptr<complex<double>>(0) ;
//	for(int i=0; i<m.cols; i++){
//		bitset<32> b(i) ;
//		b = reverseBit(b) ;
//		int j = b.to_ulong() ;
//		p[j] = mp[i] ;
//	}
//	return v ;
//}