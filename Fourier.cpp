#include "Fourier.h"

vector<unsigned int> bitReversalTable ;

void setTable(){
	int n = 1<<16 ;
	bitReversalTable = vector<unsigned int>(n) ;
	for(unsigned int i=0; i<n; i++){
		bitset<16> b(i) ;
		b = reverseBit(b) ;
		bitReversalTable[i] = b.to_ulong() ;
		//cout << i <<endl ;
	}
	//sort(bitReversalTable.begin(), bitReversalTable.end(), mCmp) ;
}

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
	vector<vector<complex<double>>> v(m.rows, vector< complex<double>>(m.cols)) ;
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

Mat vectorToMat(vector<vector<complex<double>>> &v){
	Mat m = Mat::zeros(v[0].size(), v.size(), CV_8UC1) ;
	uchar *p ;
	double pmin=abs(v[0][0]), pmax=abs(v[0][0]) ;
	for(int r=0; r<m.rows; r++){
		//p = m.ptr<double>(r) ;
		for(int c=0; c<m.cols; c++){
			//p[c] = abs(v[r][c]) ;
			double p = abs(v[r][c]) ;
			pmin = min(p, pmin) ;
			pmax = max(p, pmax) ;
		}
	}
	//cout << pmin << pmax << endl ;
	for(int r=0; r<m.rows; r++){
		p = m.ptr<uchar>(r) ;
		for(int c=0; c<m.cols; c++){
			p[c] = ((abs(v[r][c])-pmin)/(pmax-pmin))*255 ;
		}
	}
	return m ;
}

Mat vectorToMat_enhanced(vector<vector<complex<double>>> &v){
	Mat m = Mat::zeros(v[0].size(), v.size(), CV_8UC1) ;
	uchar *p ;
	double pmin=abs(v[0][0]), pmax=abs(v[0][0]) ;
	for(int r=0; r<m.rows; r++){
		//p = m.ptr<double>(r) ;
		for(int c=0; c<m.cols; c++){
			//p[c] = abs(v[r][c]) ;
			double p = log(1+abs(v[r][c])) ;
			pmin = min(p, pmin) ;
			pmax = max(p, pmax) ;
		}
	}
	//cout << pmin << pmax << endl ;
	for(int r=0; r<m.rows; r++){
		p = m.ptr<uchar>(r) ;
		for(int c=0; c<m.cols; c++){
			p[c] = (log(1+abs(v[r][c]))-pmin)/(pmax-pmin)*255 ;
			//p[c] = ((abs(v[r][c])-pmin)/(pmax-pmin))*255 ;
			//pmin = min(p[c], pmin) ;
			//pmax = max(p[c], pmax) ;
		}
	}
	return m ;
}

bitset<16> reverseBit(bitset<16> &b){
	string bs = b.to_string() ;
	reverse(bs.begin(), bs.end()) ;
	return bitset<16>(bs) ;
}

vector<complex<double>> FFT_1D(vector<complex<double>> &m){
	//vector<complex<double>> v = rearrange(m) ;
	vector<complex<double>> f(m) ;
	vector<complex<double>> tmp(f.size()) ;
	for(int n=2; n<=f.size(); n=n<<1){
		//vector<complex<double>> tmp(v.size()) ;
		int nGroup = f.size()/n ;
		int k = n >> 1 ;
		for(int index=0; index<nGroup; index++){
			for(int j=0; j<k; j++){
				int u = index*n + j ;
				double angle = CV_PI*j/k ;
				complex<double> w(cos(angle), -sin(angle)) ;
				tmp[u] = f[u] + f[u+k]*w ;
				tmp[u+k] = f[u] - f[u+k]*w ;
				//cout << u << ", " << u+k << endl ;
			}
		}
		f = tmp ;
	}
	return f ;
}

vector< complex<double>> IFFT_1D(vector< complex<double>> &m){
	vector<complex<double>> f = m ;
	for(int n=2; n<=f.size(); n=n<<1){
		vector<complex<double>> tmp(f.size()) ;
		int nGroup = f.size()/n ;
		int k = n >> 1 ;
		for(int index=0; index<nGroup; index++){
			for(int j=0; j<k; j++){
				int u = index*n + j ;
				double angle = CV_PI*j/k ;
				complex<double> w(cos(angle), -sin(angle)) ;
				tmp[u] = f[u] + f[u+k]*w ;
				tmp[u+k] = f[u] - f[u+k]*w ;
				//cout << u << ", " << u+k << endl ;
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
	//rearrange
	vector<int> tmp(m.size()) ;
	for(int i=0; i<m.size(); i++){
		tmp[i] = i ;
	}
	sort(tmp.begin(), tmp.end(), mCmp) ;
	for(int i=0; i<m.size(); i++){
		for(int r=0; r<m.size(); r++)
			f[r][i] = m[r][tmp[i]] ;
	}
	//FFT
	for(int r=0; r<m.size(); r++){
		f[r] = FFT_1D(f[r]) ;
	}
	vector<vector<complex<double>>> fr(f) ;
	//cols
	//rearrange
	for(int i=0; i<m.size(); i++){
		for(int c=0; c<m[0].size(); c++)
			f[i][c] = fr[tmp[i]][c] ;
	}
	//FFT
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

vector<vector<complex<double>>> IFFT_2D(vector<vector<complex<double>>> &m){
	vector<vector<complex<double>>> f(m) ;
	//rows
	for(int r=0; r<m.size(); r++){
		f[r] = IFFT_1D(m[r]) ;
	}
	//cols
	for(int c=0; c<m[0].size(); c++){
		vector<complex<double>> col(m.size()) ;
		//get col vector
		for(int i=0; i<m.size(); i++)
			col[i] = f[i][c] ;
		col = IFFT_1D(col) ;
		//set col vector
		for(int i=0; i<m.size(); i++)
			f[i][c] = col[i] ;
	}
	return f ;
}

bool mCmp(int a, int b){
	//bitset<16> ba(a) ;
	//bitset<16> bb(b) ;
	//ba = reverseBit(ba) ;
	//bb = reverseBit(bb) ;
	//return ba.to_ulong() < bb.to_ulong() ;
	return bitReversalTable[a] > bitReversalTable[b] ;
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

Mat enhanceSpectrum(Mat &m){
	Mat e = Mat::zeros(m.rows, m.cols, m.type()) ;
	uchar *p ;
	uchar *mp ;
	for(int r=0; r<m.rows; r++){
		p = e.ptr<uchar>(r) ;
		mp = m.ptr<uchar>(r) ;
		for(int c=0; c<m.cols; c++){
			p[c] = log(1+mp[c]) ;
			if(mp[c] != 0)
				cout << (int)mp[c] << ", " << (int)p[c] << endl ;
		}
	}
	return e ;
}

Mat mulMat(Mat &m1, Mat &m2){
	Mat img = Mat::zeros(m1.rows, m1.cols, m1.type()) ;
	if(m1.rows!=m2.rows && m1.cols!=m2.cols)
		return img ;
	uchar *p, *p1 ;
	double *p2 ;
	for(int r=0; r<m1.rows; r++){
		p = img.ptr<uchar>(r) ;
		p1 = m1.ptr<uchar>(r) ;
		p2 = m2.ptr<double>(r) ;
		for(int c=0; c<m1.cols; c++){
			p[c] = (p1[c]/255 * p2[c])*255 ;
			//cout << (double)p[c] << endl ;
		}
	}
	return img ;
}

Mat difMat(Mat &m1, Mat &m2){
	Mat img = Mat::zeros(m1.rows, m1.cols, m1.type()) ;
	uchar *p, *p1, *p2 ;
	for(int r=0; r<m1.rows; r++){
		p = img.ptr<uchar>(r) ;
		p1 = m1.ptr<uchar>(r) ;
		p2 = m2.ptr<uchar>(r) ;
		for(int c=0; c<m1.cols; c++){
			p[c] =abs(p1[c]-p2[c]) ;
			//cout << (double)p[c] << endl ;
		}
	}
	return img ;
}

Mat HFEF_Gausian(Mat &m, double d0){
	static double k1=0.4, k2=0.8 ;
	Mat f = Mat::zeros(m.rows, m.cols, CV_64FC1) ;
	double *p ;
	double sigma = 2*d0*d0 ;
	for(int r=0; r<f.rows; r++){
		p = f.ptr<double>(r) ;
		for(int c=0; c<f.cols; c++){
			double duv = sqrt(r*r + c*c) ;
			double g = exp(-duv/sigma) ;
			double h = k1 + k2*(1-g) ;
			p[c] = h ;
		}
	}
	return mulMat(m, f) ;
}

Mat reverseMat(Mat &m){
	Mat img = Mat::zeros(m.rows, m.cols, m.type()) ;
	uchar *p1, *p2 ;
	for(int r=0; r<m.rows; r++){
		p1 = img.ptr<uchar>(r) ;
		p2 = m.ptr<uchar>(m.rows-r-1) ;
		for(int c=0; c<m.cols; c++){
			p1[c] = p2[m.cols-c-1] ;
		}
	}
	return img ;
}