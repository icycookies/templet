#include<cstdio>
#include<algorithm>
#include<cstring>
#include<cmath>

//#define Integer
//#define Double
#define Module
//choose one

using namespace std;
typedef long long ll;
const int mod=998244353;
const double eps=1e-8;
class Int{
	public:
		Int();
		Int(int x);
		
		friend Int operator + (const Int &a,const Int &b);
		friend Int operator + (const Int &a,int b);
		friend Int operator + (int a,const Int &b);
		friend Int operator - (const Int &a);
		friend Int operator - (const Int &a,const Int &b);
		friend Int operator - (const Int &a,int b);
		friend Int operator - (int a,const Int &b);
		friend Int operator * (const Int &a,const Int &b);
		friend Int operator * (const Int &a,int b);
		friend Int operator * (int a,const Int &b);
		friend Int operator / (const Int &a,const Int &b);
		friend Int operator / (const Int &a,int b);
		friend Int operator / (int a,const Int &b);
		friend bool operator == (const Int &a,const Int &b);
		friend bool operator == (const Int &a,int b);
		friend bool operator == (int a,const Int &b);
		friend bool operator != (const Int &a,const Int &b);
		friend bool operator != (const Int &a,int b);
		friend bool operator != (int a,const Int &b);
		
		Int & operator += (const Int &b);
		Int & operator += (int b);
		Int & operator -= (const Int &b);
		Int & operator -= (int b);
		Int & operator *= (const Int &b);
		Int & operator *= (int b);
		Int & operator /= (const Int &a);
		Int & operator /= (int b);
		
		Int inv() const;
		Int pow(int x) const;
		Int pow(const Int &x) const;
		bool is_prime() const;
		
		void print();
	private:
		int w;
};
Int::Int(){w=0;}
Int::Int(int x){w=x;}
Int operator + (const Int &a,const Int &b){
	Int ref=a;
	ref.w+=b.w;
	if (ref.w>mod)ref.w-=mod;
	return ref;
}
Int operator + (const Int &a,int b){Int ref(b);return a+ref;}
Int operator + (int a,const Int &b){Int ref(a);return b+ref;}
Int operator - (const Int &a){Int ref=a;if (ref.w)ref.w=mod-ref.w;return ref;}
Int operator - (const Int &a,const Int &b){return a+(-b);}
Int operator - (const Int &a,int b){Int ref(b);return a-ref;}
Int operator - (int a,const Int &b){Int ref(a);return ref-b;}
Int operator * (const Int &a,const Int &b){
	Int ref=a;
	ref.w=(ll)ref.w*b.w%mod;
	return ref;
}
Int operator * (const Int &a,int b){Int ref(b);return a*ref;}
Int operator * (int a,const Int &b){Int ref(a);return b*ref;}
Int operator / (const Int &a,const Int &b){return a*b.inv();}
Int operator / (const Int &a,int b){Int ref(b);return a/ref;}
Int operator / (int a,const Int &b){Int ref(a);return ref/b;}
bool operator == (const Int &a,const Int &b){if (a.w==b.w)return true;else return false;}
bool operator == (const Int &a,int b){Int ref(b);return a==ref;}
bool operator == (int a,const Int &b){Int ref(a);return b==ref;}
bool operator != (const Int &a,const Int &b){return !(a==b);}
bool operator != (const Int &a,int b){return !(a==b);}
bool operator != (int a,const Int &b){return !(a==b);}
Int&Int::operator += (const Int &b){*this=*this+b;return *this;}
Int&Int::operator += (int b){*this=*this+b;return *this;}
Int&Int::operator -= (const Int &b){*this=*this-b;return *this;}
Int&Int::operator -= (int b){*this=*this-b;return *this;}
Int&Int::operator *= (const Int &b){*this=*this*b;return *this;}
Int&Int::operator *= (int b){*this=*this*b;return *this;}
Int&Int::operator /= (const Int &b){*this=*this/b;return *this;}
Int&Int::operator /= (int b){*this=*this/b;return *this;}
Int Int::inv() const{
	Int ref(w);
	return ref.pow(mod-2); 
}
Int Int::pow(const Int &x) const{
	int n=x.w;
	Int base(w),ans(1);
	while (n){
		if (n&1)ans*=base;
		base*=base;
		n>>=1;
	}
	return ans;
}
Int Int::pow(int x) const{
	Int ref(w),power(x);
	return ref.pow(power);
}
bool Int::is_prime() const{
	for (int i=2;i*i<=w;i++)
		if (!(w%i))return false;
	return true;
}
void Int::print(){printf("%d",w);}
template <class T> class vector{
	public:
		vector();
		vector(int x);
		void reset(int x);
		void init(T *x);
		
		template <class S> friend vector<S> operator + (const vector<S> &a,const vector<S> &b);
		template <class S> friend vector<S> operator - (const vector<S> &a,const vector<S> &b);
		template <class S> friend vector<S> operator - (const vector<S> &a);
		template <class S> friend vector<S> operator * (const vector<S> &a,S x);
		template <class S> friend vector<S> operator * (S x,const vector<S> &a);
		template <class S> friend vector<S> operator / (const vector<S> &a,S x);
		template <class S> friend S operator * (const vector<S> &a,const vector<S> &b);
		template <class S> friend bool operator == (const vector<S> &a,const vector<S> &b);
		template <class S> friend bool operator != (const vector<S> &a,const vector<S> &b);
		
		vector<T> & operator += (const vector<T> &b);
		vector<T> & operator -= (const vector<T> &b);
		vector<T> & operator *= (T x);
		vector<T> & operator /= (T x);
		
		double length();
		template <class S> friend vector<S> dim3_product(const vector<S> &a,const vector<S> &b);
		template <class S> friend S dim2_product(const vector<S> &a,const vector<S> &b);
		
		void print();
		void print(int x);
	private:
		static const int N=1000;
		int n;
		T a[N];
};
template <class T> vector<T>::vector(){n=0;for (int i=0;i<N;i++)a[i]=0;}
template <class T> vector<T>::vector(int x){n=x;for (int i=0;i<n;i++)a[i]=0;}
template <class T> void vector<T>::reset(int x){vector<T> ref(x);*this=ref;}
template <class T> void vector<T>::init(T *x){for (int i=0;i<n;i++)a[i]=*(x+i);}
template <class T> vector<T> operator + (const vector<T> &a,const vector<T> &b){
	vector<T> ret=a;
	if (a.n==b.n)for (int i=0;i<ret.n;i++)ret.a[i]+=b.a[i];else puts("The length of a and b is not equal!");
	return ret;
}
template <class T> vector<T> operator - (const vector<T> &a){
	vector<T> ret=a;
	for (int i=0;i<ret.n;i++)ret.a[i]=-ret.a[i];
	return ret;
}
template <class T> vector<T> operator - (const vector<T> &a,const vector<T> &b){
	return a+(-b);
}
template <class T> vector<T> operator * (const vector<T> &a,T x){
	vector<T> ret=a;
	for (int i=0;i<ret.n;i++)ret.a[i]*=x;
	return ret;
}
template <class T> vector<T> operator * (T x,const vector<T> &a){
	vector<T> ret=a;
	return ret*x;
}
template <class T> vector<T> operator / (const vector<T> &a,T x){
	vector<T> ret=a;
	if (x)for (int i=0;i<ret.n;i++)ret.a[i]/=x;else puts("Divided by 0!");
	return ret;
}
template <class T> double operator * (const vector<T> &a,const vector<T> &b){
	double ans=0;
	if (a.n==b.n)for (int i=0;i<a.n;i++)ans+=a.a[i]*b.a[i];else puts("The length of a and b is not equal!");
	return ans;
}
template <class T> bool operator == (const vector<T> &a,const vector<T> &b){
	if (a.n!=b.n)return false;
	for (int i=0;i<a.n;i++)if (a.a[i]!=b.a[i])return false;
	return true;
}
template <class T> bool operator != (const vector<T> &a,const vector<T> &b){
	return !(a==b);
}
template <class T> vector<T>&vector<T>::operator +=(const vector<T> &b){*this=*this+b;return *this;}
template <class T> vector<T>&vector<T>::operator -=(const vector<T> &b){*this=*this-b;return *this;}
template <class T> vector<T>&vector<T>::operator *=(T x){*this=*this*x;return *this;}
template <class T> vector<T>&vector<T>::operator /=(T x){*this=*this/x;return *this;}
template <class T> vector<T> dim3_product(const vector<T> &a,const vector<T> &b){
	vector<T> ret=a;
	if (a.n==3&&b.n==3){
		ret.a[0]=a.a[1]*b.a[2]-a.a[2]*b.a[1];
		ret.a[1]=a.a[2]*b.a[0]-a.a[0]*b.a[2];
		ret.a[2]=a.a[0]*b.a[1]-a.a[1]*b.a[0];
	}
	return ret;
}
template <class T> T dim2_product(const vector<T> &a,const vector<T> &b){
	if (a.n==2&&b.n==2)return a.a[0]*b.a[1]-a.a[1]*b.a[0];
		else return 0;
}
template <class T> double vector<T>::length(){
	double ans=0;
	for (int i=0;i<n;i++)ans+=a[i]*a[i];
	return sqrt(ans);
}
template <class T> void vector<T>::print(){
	#ifdef Integer
		printf("[%d",a[0]);
		for (int i=1;i<n;i++)printf(",%d",a[i]);
	#endif
	#ifdef Double
		printf("[%.2lf",a[0]);
		for (int i=1;i<n;i++)printf(",%.2lf",a[i]);
	#endif
	#ifdef Module
		printf("[");a[0].print();
		for (int i=1;i<n;i++){printf(",");a[i].print();}
	#endif
	puts("]");
}
template <class T> void vector<T>::print(int x){
	#ifdef Integer
		printf("%d\n",a[x]);
	#endif
	#ifdef Double
		printf("%.2lf\n",a[x]);
	#endif
	#ifdef Module
		a[x].print();puts("");
	#endif
}
class complex{
	public:
		complex();
		complex(int X);
		complex(double X);
		complex(int X,int Y);
		complex(double X,double Y);
		
		friend complex operator + (const complex &a,const complex &b);
		friend complex operator + (const complex &a,int x);
		friend complex operator + (const complex &a,double x);
		friend complex operator + (int x,const complex &a);
		friend complex operator + (double x,const complex &a);
		friend complex operator - (const complex &a,const complex &b);
		friend complex operator - (const complex &a,int x);
		friend complex operator - (const complex &a,double x);
		friend complex operator - (int x,const complex &a);
		friend complex operator - (double x,const complex &a);
		friend complex operator - (const complex &a);
		friend complex operator * (const complex &a,const complex &b);
		friend complex operator * (const complex &a,int x);
		friend complex operator * (const complex &a,double x);
		friend complex operator * (int x,const complex &a);
		friend complex operator * (double x,const complex &a);
		friend complex operator / (const complex &a,const complex &b);
		friend complex operator / (const complex &a,int x);
		friend complex operator / (const complex &a,double x);
		friend complex operator / (int x,const complex &a);
		friend complex operator / (double x,const complex &a);
		friend bool operator == (const complex &a,const complex &b);
		friend bool operator != (const complex &a,const complex &b);
		
		complex & operator += (const complex &b);
		complex & operator += (int b);
		complex & operator += (double b);
		complex & operator -= (const complex &b);
		complex & operator -= (int b);
		complex & operator -= (double b);
		complex & operator *= (const complex &b);
		complex & operator *= (int b);
		complex & operator *= (double b);
		complex & operator /= (const complex &b);
		complex & operator /= (int b);
		complex & operator /= (double b);
		
		double length() const;
		complex conj() const;
		
		void print();
		void print_real();
		void print_imag();
	private:
		double x,y;
};
complex::complex(){x=y=0;}
complex::complex(int X){x=X;y=0;}
complex::complex(double X){x=X;y=0;}
complex::complex(int X,int Y){x=X;y=Y;}
complex::complex(double X,double Y){x=X;y=Y;}
complex operator + (const complex &a,const complex &b){
	complex ret=a;
	ret.x+=b.x;ret.y+=b.y;
	return ret;
}
complex operator + (const complex &a,int x){
	complex ret(x);
	return ret+a;
}
complex operator + (const complex &a,double x){
	complex ret(x);
	return ret+a;
}
complex operator + (int x,const complex &a){return a+x;}
complex operator + (double x,const complex &a){return a+x;}
complex operator - (const complex &a){
	complex ret=a;
	ret.x=-ret.x;ret.y=-ret.y;
	return ret;
}
complex operator - (const complex &a,const complex &b){return a+(-b);}
complex operator - (const complex &a,int x){
	complex ret(x);
	return a-ret;
}
complex operator - (const complex &a,double x){
	complex ret(x);
	return a-ret;
}
complex operator - (int x,const complex &a){return -(a-x);}
complex operator - (double x,const complex &a){return -(a-x);}
complex operator * (const complex &a,const complex &b){
	complex ret;
	ret.x=a.x*b.x-a.y*b.y;
	ret.y=a.x*b.y+a.y*b.x;
	return ret;
}
complex operator * (const complex &a,int x){
	complex ret(x);
	return a*ret;
}
complex operator * (const complex &a,double x){
	complex ret(x);
	return a*ret;
}
complex operator * (int x,const complex &a){return a*x;}
complex operator * (double x,const complex &a){return a*x;}
complex operator / (const complex &a,const complex &b){
	double x=b.length();x=x*x;
	return a*b.conj()/x;
}
complex operator / (const complex &a,int x){
	complex ret=a;
	ret.x/=x;ret.y/=x;
	return ret;
}
complex operator / (const complex &a,double x){
	complex ret=a;
	ret.x/=x;ret.y/=x;
	return ret;
}
complex operator / (int x,const complex &a){
	complex ret(x);
	return ret/a;
}
complex operator / (double x,const complex &a){
	complex ret(x);
	return ret/a;
}
bool operator == (const complex &a,const complex &b){
	if (a.x==b.x&&a.y==b.y)return true;else return false;
}
bool operator != (const complex &a,const complex &b){return !(a==b);}
complex&complex::operator += (const complex &b){*this=*this+b;return *this;}
complex&complex::operator += (int b){*this=*this+b;return *this;}
complex&complex::operator += (double b){*this=*this+b;return *this;}
complex&complex::operator -= (const complex &b){*this=*this-b;return *this;}
complex&complex::operator -= (int b){*this=*this-b;return *this;}
complex&complex::operator -= (double b){*this=*this-b;return *this;}
complex&complex::operator *= (const complex &b){*this=*this*b;return *this;}
complex&complex::operator *= (int b){*this=*this*b;return *this;}
complex&complex::operator *= (double b){*this=*this*b;return *this;}
complex&complex::operator /= (const complex &b){*this=*this/b;return *this;}
complex&complex::operator /= (int b){*this=*this/b;return *this;}
complex&complex::operator /= (double b){*this=*this/b;return *this;}
double complex::length() const{return sqrt(x*x+y*y);}
complex complex::conj() const{
	complex ret=*this;
	ret.y=-ret.y;
	return ret;
}
void complex::print(){
	if (x){
		printf("%.2lf",x);
		if (y<-1e-8)printf("%.2lfi\n",y);
			else if (y>1e-8)printf("+%.2lfi\n",y);
	}else if (y<-1e-8||y>1e-8)printf("%.2lfi\n",y);
		else puts("0");
}
void complex::print_real(){printf("%.2lf\n",x);}
void complex::print_imag(){printf("%.2lf\n",y);}
template <class T> class matrix{
	public:
		matrix();
		matrix(int x,int y);
		void reset(int x,int y);
		void set_unitary();
		void init(T **x);
		
		template <class S> friend matrix<S> operator + (const matrix<S> &a,const matrix<S> &b);
		template <class S> friend matrix<S> operator - (const matrix<S> &a,const matrix<S> &b);
		template <class S> friend matrix<S> operator * (const matrix<S> &a,S x);
		template <class S> friend matrix<S> operator * (S x,const matrix<S> &a);
		template <class S> friend matrix<S> operator * (const matrix<S> &a,const matrix<S> &b);
		template <class S> friend bool operator == (const matrix<S> &a,const matrix<S> &b);
		template <class S> friend bool operator != (const matrix<S> &a,const matrix<S> &b);
		
		matrix<T> & operator += (const matrix<T> &b);
		matrix<T> & operator -= (const matrix<T> &b);
		matrix<T> & operator *= (const matrix<T> &b);
		matrix<T> & operator *= (const T x);
		
		template <class S> friend int eliminate(matrix<S> &a,matrix<S> &b);
		matrix<T> trans() const;
		matrix<T> rref() const;
		matrix<T> inv() const;
		T det() const;
		int rank() const;
		void print();
	private:
		static const int N=400;
		int n,m;
		T a[N][N];
};
template <class T> matrix<T>::matrix(){n=m=0;for (int i=0;i<N;i++)for (int j=0;j<N;j++)a[i][j]=0;}
template <class T> matrix<T>::matrix(int x,int y){n=x;m=y;for (int i=0;i<n;i++)for (int j=0;j<m;j++)a[i][j]=0;}
template <class T> void matrix<T>::reset(int x,int y){matrix<T> ref(x,y);*this=ref;}
template <class T> void matrix<T>::set_unitary(){
	if (n!=m){puts("This is a non-square matrix!");return;}
	for (int i=0;i<n;i++)
		for (int j=0;j<m;j++)
			if (i==j)a[i][j]=1;else a[i][j]=0;
}
template <class T> void matrix<T>::init(T **x){for (int i=0;i<n;i++)for (int j=0;j<m;j++)a[i][j]=*((T*)x+m*i+j);}
template <class T> matrix<T> operator + (const matrix<T> &a,const matrix<T> &b){
	matrix<T> ref=a;
	if (a.n==b.n&&a.m==b.m)
		for (int i=0;i<ref.n;i++)
			for (int j=0;j<ref.m;j++)ref.a[i][j]+=b.a[i][j];
	return ref;
}
template <class T> matrix<T> operator - (const matrix<T> &a,const matrix<T> &b){
	matrix<T> ref=a;
	if (a.n==b.n&&a.m==b.m)
		for (int i=0;i<ref.n;i++)
			for (int j=0;j<ref.m;j++){
				ref.a[i][j]-=b.a[i][j];
				if (ref.a[i][j]<0)ref.a[i][j]+=ref.mod;
			}
	return ref;
}
template <class T> matrix<T> operator * (const matrix<T> &a,T x){
	matrix<T> ref=a;
	for (int i=0;i<ref.n;i++)
		for (int j=0;j<ref.m;j++)ref.a[i][j]=ref.a[i][j]*x;
	return ref;
}
template <class T> matrix<T> operator * (T x,const matrix<T> &a){return a*x;}
template <class T> matrix<T> operator * (const matrix<T> &a,const matrix<T> &b){
	if (a.m==b.n){
		matrix<T> ref(a.n,b.m);
		for (int i=0;i<a.n;i++)
			for (int j=0;j<a.m;j++)
				for (int k=0;k<b.m;k++)
					ref.a[i][k]=ref.a[i][k]+a.a[i][j]*b.a[j][k];
		return ref;
	}
}
template <class T> bool operator == (const matrix<T> &a,const matrix<T> &b){
	if (a.n==b.n&&a.m==b.m){
		for (int i=0;i<a.n;i++)
			for (int j=0;j<a.m;j++)if (a.a[i][j]!=b.a[i][j])return false;
		return true;
	}
	return false;
}
template <class T> bool operator != (const matrix<T> &a,const matrix<T> &b){
	return !(a==b);
}
template <class T> matrix<T>&matrix<T>::operator += (const matrix<T> &b){*this=*this+b;return *this;}
template <class T> matrix<T>&matrix<T>::operator -= (const matrix<T> &b){*this=*this-b;return *this;}
template <class T> matrix<T>&matrix<T>::operator *= (T b){*this=*this*b;return *this;}
template <class T> matrix<T>&matrix<T>::operator *= (const matrix<T> &b){*this=b**this;return *this;}
template <class T> void matrix<T>::print(){
	for (int i=0;i<n;i++){
		for (int j=0;j<m;j++){
			#ifdef Integer
				printf("%d ",a[i][j]);
			#endif
			#ifdef Double
				printf("%.4lf ",a[i][j]);
			#endif
			#ifdef Module
				a[i][j].print();printf(" ");
			#endif
		}
		puts("");
	}
}
template <class T> matrix<T> matrix<T>::trans() const{
	matrix<T> ref(m,n);
	for (int i=0;i<m;i++)
		for (int j=0;j<n;j++)ref.a[i][j]=a[j][i];
	return ref;
}
template <class T> matrix<T> matrix<T>::rref() const{
	matrix<T> a=*this,b(n,m);
	eliminate(a,b);
	return a;
}
template <class T> matrix<T> matrix<T>::inv() const{
	if (n!=m){puts("This is a non-square matrix!");return *this;}
	matrix<T> a=*this,b(n,m);
	b.set_unitary();
	if (abs(eliminate(a,b))==n){
		for (int i=n-1;i>=0;i--){
			for (int j=0;j<n;j++)b.a[i][j]/=a.a[i][i];
			for (int j=0;j<i;j++)
				for (int k=0;k<n;k++)b.a[j][k]-=a.a[j][i]*b.a[i][k];
		}
		return b;
	}else{
		puts("The inverse does not exist!");
		return *this;
	}
}
template <class T> T matrix<T>::det() const{
	#ifdef Module
		Int zero(0); 
	#else
		T zero=0;
	#endif
	if (n!=m){puts("This is a non-square matrix!");return zero;}
	matrix<T> a=*this,b(n,m);
	int tmp=eliminate(a,b);
	if (abs(tmp)==n){
		#ifdef Module
			Int ans(1);
		#else
			T ans=1;
		#endif
		for (int i=0;i<n;i++)ans*=a.a[i][i];
		if (tmp<0)ans=-ans;
		return ans;
	}else return zero;
}
template <class T> int matrix<T>::rank() const{
	matrix<T> a=*this,b(n,m);
	return abs(eliminate(a,b));
}
template <class T> int eliminate(matrix<T> &a,matrix<T> &b){
	int row=0,n=a.n,m=a.m;
	int i,j,k,flag=1;
	for (i=0;i<m&&row<n;i++){
		for (j=row;j<n;j++){
			#ifdef Integer
				if (a.a[j][i])break;
			#endif
			#ifdef Double
				if (fabs(a.a[j][i])>eps)break;
			#endif
			#ifdef Module
				if (a.a[j][i]!=0)break;
			#endif
		}
		if (j<n){
			if (j!=row){for (k=0;k<m;k++)swap(a.a[row][k],a.a[j][k]),swap(b.a[row][k],b.a[j][k]);flag=-flag;}
			for (j=row+1;j<n;j++){
				T ratio=a.a[j][i]/a.a[row][i];
				for (k=0;k<m;k++)a.a[j][k]-=a.a[row][k]*ratio,b.a[j][k]-=b.a[row][k]*ratio;
			}
			row++;
		}
	}
	return row*flag;
}
int main(){
	puts("This is linear algebra test!");
	/*vector<Int> a(3),b(3);
	Int arr[3]={Int(1),Int(0),Int(0)},Arr[3]={Int(0),Int(1),Int(0)};
	a.init((Int *)arr);
	b.init((Int *)Arr);
	vector<Int> c=dim3_product(a,b);
	c.print();
	c+=b;
	c.print();*/
	matrix<Int> a(2,2);
	Int arr[2][2]={{Int(0),Int(1)},{Int(1),Int(0)}};
	a.init((Int **)arr);
	a.det().print();
}
