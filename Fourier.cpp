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

vector< vector< complex<double>>> matToVector(Mat &m){
	Mat md ;
	m.convertTo(md, CV_64FC1) ;
	vector<vector<complex<double>>> v(m.rows, vector< complex<double>>(m.cols)) ;
	double *p ;
	for(int r=0; r<m.rows; r++){
		p = md.ptr<double>(r) ;
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
	//cout << pmin << "\t" << pmax << endl ;
	for(int r=0; r<m.rows; r++){
		p = m.ptr<uchar>(r) ;
		for(int c=0; c<m.cols; c++){
			//p[c] = abs(v[r][c]) ;
			p[c] = ((abs(v[r][c])-pmin)/(pmax-pmin))*255 ;
			//p[c] = abs(v[r][c])-pmin ;
		}
	}
	return m ;
}

Mat vectorToMat_powerSpec(vector<vector<complex<double>>> &v){
	Mat m = Mat::zeros(v[0].size(), v.size(), CV_8UC1) ;
	uchar *p ;
	double pmin=norm(v[0][0]), pmax=norm(v[0][0]) ;
	for(int r=0; r<m.rows; r++){
		//p = m.ptr<double>(r) ;
		for(int c=0; c<m.cols; c++){
			//p[c] = abs(v[r][c]) ;
			double p = norm(v[r][c]) ;
			pmin = min(p, pmin) ;
			pmax = max(p, pmax) ;
		}
	}
	//cout << pmin << pmax << endl ;
	for(int r=0; r<m.rows; r++){
		p = m.ptr<uchar>(r) ;
		for(int c=0; c<m.cols; c++){
			p[c] = ((norm(v[r][c])-pmin)/(pmax-pmin))*255 ;
			//p[c] = abs(v[r][c])-pmin ;
		}
	}
	return m ;
}

Mat vectorToMat_real(vector<vector<complex<double>>> &v){
	Mat m = Mat::zeros(v[0].size(), v.size(), CV_8UC1) ;
	uchar *p ;
	double pmin=v[0][0].real(), pmax=v[0][0].real() ;
	for(int r=0; r<m.rows; r++){
		//p = m.ptr<double>(r) ;
		for(int c=0; c<m.cols; c++){
			//p[c] = abs(v[r][c]) ;
			double p = v[r][c].real()*((r+c)&1==1? -1: 1) ;
			pmin = min(p, pmin) ;
			pmax = max(p, pmax) ;
		}
	}
	//cout << pmin << pmax << endl ;
	for(int r=0; r<m.rows; r++){
		p = m.ptr<uchar>(r) ;
		for(int c=0; c<m.cols; c++){
			//p[c] = ((v[r][c].real()-pmin)/(pmax-pmin)) ;
			p[c] = ((v[r][c].real()*((r+c)&1==1? -1: 1)-pmin)/(pmax-pmin))*255 ;
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
	//cout << pmin << "\t" << pmax << endl ;
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
	//vector<complex<double>> tmp(f.size()) ;
	for(int n=2; n<=f.size(); n=n<<1){
		//vector<complex<double>> tmp(v.size()) ;
		int nGroup = f.size()/n ;
		int k = n >> 1 ;
		for(int index=0; index<nGroup; index++){
			for(int j=0; j<k; j++){
				int u = index*n + j ;
				double angle = CV_PI*j/k ;
				complex<double> w(cos(angle), -sin(angle)) ;
				complex<double> even = f[u] ;
				complex<double> odd = f[u+k]*w ;
				f[u] = even+odd ;
				f[u+k] = even-odd ;
				//tmp[u] = f[u] + f[u+k]*w ;
				//tmp[u+k] = f[u] - f[u+k]*w ;
				//cout << u << ", " << u+k << endl ;
			}
		}
		//f = tmp ;
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

vector<vector<complex<double>>> FFT_2D(vector<vector<complex<double>>> &m){
	vector<vector<complex<double>>> f(m) ;
	//rows
	//rearrange
	vector<int> tmp(m.size()) ;
	for(int i=0; i<m.size(); i++){
		tmp[i] = i ;
	}
	sort(tmp.begin(), tmp.end(), mCmp) ;
	for(int i=0; i<m[0].size(); i++){
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
		for(int c=0; c<m[0].size(); c++){
			f[i][c] = fr[tmp[i]][c] ;
		}
	}
	//FFT
	for(int c=0; c<m[0].size(); c++){
		vector<complex<double>> col(m.size()) ;
		//get col vector
		for(int i=0; i<m.size(); i++){
			col[i] = f[i][c] ;
		}
		col = FFT_1D(col) ;
		//set col vector
		for(int i=0; i<m.size(); i++){
			f[i][c] = col[i] ;
		}
	}
	return f ;
}

vector<vector<complex<double>>> IFFT_2D(vector<vector<complex<double>>> &m){
	vector<vector<complex<double>>> v = m ;
	for(int r=0; r<v.size(); r++){
		for(int c=0; c<v[0].size(); c++){
			v[r][c] = conj(v[r][c]) ;
		}
	}
	v = FFT_2D(v) ;
	double mn = v.size()*v[0].size() ;
	for(int r=0; r<v.size(); r++){
		for(int c=0; c<v[0].size(); c++){
			v[r][c] /= mn ;
		}
	}
	return v ;
}

bool mCmp(int a, int b){
	return bitReversalTable[a] < bitReversalTable[b] ;
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
			//p[c] = (p1[c]/255 * p2[c])*255 ;
			p[c] = p1[c] * p2[c] ;
			//cout << (double)p[c] << endl ;
		}
	}
	return img ;
}

vector<vector<complex<double>>> mulVec(vector<vector<complex<double>>> &v1, vector<vector<complex<double>>> &v2){
	vector<vector<complex<double>>> v = v1 ;
	for(int r=0; r<v1.size(); r++){
		for(int c=0; c<v1[0].size(); c++){
			v[r][c] = v1[r][c]*v2[r][c] ;
		}
	}
	return v ;
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

vector<vector<complex<double>>> LPF_Gaussian(vector<vector<complex<double>>> &m, double d0){
	vector<vector<complex<double>>> f(m.size(), vector<complex<double>>(m[0].size(), complex<double>(0, 0))) ;
	double sigma = 2*d0*d0 ;
	int cr = m.size()/2 ;
	int cc = m[0].size()/2 ;
	for(int r=0; r<f.size(); r++){
		for(int c=0; c<f[0].size(); c++){
			double duv = sqrt((r-cr)*(r-cr) + (c-cc)*(c-cc)) ;
			double g = exp(-duv/sigma) ;
			//double h = 1-g ;
			f[r][c] = g ;
		}
	}
	return mulVec(m, f) ;
}

vector<vector<complex<double>>> HPF_Gaussian(vector<vector<complex<double>>> &m, double d0){
	vector<vector<complex<double>>> f(m.size(), vector<complex<double>>(m[0].size(), complex<double>(0, 0))) ;
	double sigma = 2*d0*d0 ;
	int cr = m.size()/2 ;
	int cc = m[0].size()/2 ;
	for(int r=0; r<f.size(); r++){
		for(int c=0; c<f[0].size(); c++){
			double duv = ((r-cr)*(r-cr) + (c-cc)*(c-cc)) ;
			double g = exp(-duv/sigma) ;
			double h = 1-g ;
			f[r][c] = h ;
			//f[r][c] = complex<double>(k1 + k2*h, 0) ;
		}
	}
	//Mat mf = vectorToMat(f) ;
	//cvShowImage("high pass", &IplImage(mf)) ;
	return mulVec(m, f) ;
}

vector<vector<complex<double>>> HFEF_Gaussian(vector<vector<complex<double>>> &m, double d0){
	static double k1=0.5, k2=0.75 ;
	vector<vector<complex<double>>> f(m.size(), vector<complex<double>>(m[0].size(), complex<double>(0, 0))) ;
	double sigma = 2*d0*d0 ;
	int cr = m.size()/2 ;
	int cc = m[0].size()/2 ;
	for(int r=0; r<f.size(); r++){
		for(int c=0; c<f[0].size(); c++){
			double duv = ((r-cr)*(r-cr) + (c-cc)*(c-cc)) ;
			double g = exp(-duv/sigma) ;
			double h = 1-g ;
			f[r][c] = k1 + k2*h ;
			//f[r][c] = complex<double>(k1 + k2*h, 0) ;
		}
	}
	return mulVec(m, f) ;
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

vector<vector<complex<double>>> Butterworth_reject(vector<vector<complex<double>>> &m, double d0, double w, double n){
	vector<vector<complex<double>>> f(m.size(), vector<complex<double>>(m[0].size(), complex<double>(1, 0))) ;
	double d0square = d0*d0 ;
	int cr = m.size()/2 ;
	int cc = m[0].size()/2 ;
	for(int r=0; r<f.size(); r++){
		for(int c=0; c<f[0].size(); c++){
			double duv = sqrt((r-cr)*(r-cr) + (c-cc)*(c-cc)) ;
			if(duv>=(d0-w/2) && duv<=(d0+w/2)){
				double base = duv*w/(duv*duv-d0square) ;
				f[r][c] = 1/(1+pow(base, 2*n)) ;
			}
			//f[r][c] = complex<double>(k1 + k2*h, 0) ;
		}
	}
	//Mat mf = vectorToMat(f) ;
	//cvShowImage("butterworth reject", &IplImage(mf)) ;
	return mulVec(m, f) ;
}

vector<vector<complex<double>>> Butterworth_reject_byimg(vector<vector<complex<double>>> &m, string name){
	IplImage *img = cvLoadImage(name.c_str(), CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat mImg(img) ;
	Mat resize_m ;
	resize(mImg, resize_m, cvSize(m[0].size(), m.size())) ;
	vector<vector<complex<double>>> f(m.size(), vector<complex<double>>(m[0].size(), complex<double>())) ;
	uchar *p ;
	for(int r=0; r<f.size(); r++){
		p = resize_m.ptr<uchar>(r) ;
		for(int c=0; c<f[0].size(); c++){
			f[r][c] = complex<double>(p[c]/255.0, 0) ;
		}
	}
	//Mat mf = vectorToMat(f) ;
	//cvShowImage("butterworth reject", &IplImage(mf)) ;
	return mulVec(m, f) ;
}

vector<vector<complex<double>>> Butterworth_pass_byimg(vector<vector<complex<double>>> &m, string name){
	IplImage *img = cvLoadImage(name.c_str(), CV_LOAD_IMAGE_GRAYSCALE) ;
	Mat mImg(img) ;
	Mat resize_m ;
	resize(mImg, resize_m, cvSize(m[0].size(), m.size())) ;
	vector<vector<complex<double>>> f(m.size(), vector<complex<double>>(m[0].size(), complex<double>())) ;
	uchar *p ;
	for(int r=0; r<f.size(); r++){
		p = resize_m.ptr<uchar>(r) ;
		for(int c=0; c<f[0].size(); c++){
			f[r][c] = complex<double>(1-p[c]/255.0, 0) ;
		}
	}
	return mulVec(m, f) ;
}

vector<vector<complex<double>>> Butterworth_pass(vector<vector<complex<double>>> &m, double d0, double w, double n){
	vector<vector<complex<double>>> f(m.size(), vector<complex<double>>(m[0].size(), complex<double>(0, 0))) ;
	double d0square = d0*d0 ;
	int cr = m.size()/2 ;
	int cc = m[0].size()/2 ;
	for(int r=0; r<f.size(); r++){
		for(int c=0; c<f[0].size(); c++){
			double duv = sqrt((r-cr)*(r-cr) + (c-cc)*(c-cc)) ;
			if(duv>=(d0-w/2) && duv<=(d0+w/2)){
				double base = duv*w/(duv*duv-d0square) ;
				f[r][c] = 1 - 1/(1+pow(base, 2*n)) ;
			}
			//f[r][c] = complex<double>(k1 + k2*h, 0) ;
		}
	}
	//Mat mf = vectorToMat(f) ;
	//cvShowImage("butterworth reject", &IplImage(mf)) ;
	return mulVec(m, f) ;
}
