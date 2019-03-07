#ifndef _XVEC
#define _XVEC
#include <iostream>
#include <cmath>
using namespace std;

// Xvec contains vector calculations, header only

class Xvec{
	double xv, yv, zv;
	friend inline ostream & operator <<( ostream & os, Xvec const & v );
	friend inline void cross_product(const Xvec &v1, const Xvec &v2, Xvec &v3);
 public:

	Xvec(double xi=0.0, double yi=0.0, double zi=0.0) : xv(xi), yv(yi), zv(zi){}

	Xvec(const Xvec &from){xv=from.xv; yv=from.yv; zv=from.zv;}
	Xvec(const double *x){xv=*x; yv=x[1]; zv=x[2];}

	inline double operator[](int inx) const{ return (inx==0?xv:inx==1?yv:zv); }
	inline Xvec operator-() const{ return Xvec(-xv, -yv, -zv);}
	inline Xvec operator-(const Xvec &v) const{ return Xvec(xv-v.xv, yv-v.yv, zv-v.zv); }
	inline Xvec operator+(const Xvec &v) const{ return Xvec(xv+v.xv, yv+v.yv, zv+v.zv); }
	inline double operator*(Xvec &v) const{ return xv*v.xv + yv*v.yv + zv*v.zv; }
	inline Xvec operator/(double v) const{ return Xvec(xv/v, yv/v, zv/v); }
	inline Xvec operator*(double v) const{ return Xvec(xv*v, yv*v, zv*v); }

	inline Xvec &operator=(const Xvec& from){
		xv=from.xv; yv=from.yv; zv=from.zv; return(*this);
	}
	inline Xvec &operator+=(const Xvec& from){
		xv+=from.xv; yv+=from.yv; zv+=from.zv; return(*this);
	}
	inline Xvec &operator-=(const Xvec& from){
		xv-=from.xv; yv-=from.yv; zv-=from.zv; return(*this);
	}

    inline Xvec cross_product(const Xvec &v) {
        double x_v = yv * v.zv - zv * v.yv;
        double y_v = zv * v.xv - xv * v.zv;
        double z_v = xv * v.yv - yv * v.xv;

        return Xvec(x_v,y_v,z_v);
    }

	inline double normalize(){
		double r=sqrt(xv*xv + yv*yv + zv*zv), len=r;
		if(len < 1.0e-8) len=1.0e-8;
		xv/=len; yv/=len; zv/=len;
		return r;
	}
    inline double angle(Xvec &v){
        double val = xv*v.xv + yv*v.yv + zv*v.zv;
        return acos(val/this->normalize()/v.normalize());
    };

    inline double angle_cos(Xvec &v){
        double val = xv*v.xv + yv*v.yv + zv*v.zv;
        return val/this->normalize()/v.normalize();
    }
    inline char* c_str(){
        char * tmpstring;
        sprintf(tmpstring, "%8.4f %8.4f %8.4f", xv,yv,zv);
        return tmpstring;
    }
};
inline ostream & operator <<( ostream & os, Xvec const & v ){
	os<<v.xv<<' '<<v.yv<<' '<<v.zv; return os;
}
inline void cross_product(const Xvec &v1, const Xvec &v2, Xvec &v3){
	v3.xv = v1.yv*v2.zv - v1.zv*v2.yv;
	v3.yv = v1.zv*v2.xv - v1.xv*v2.zv;
	v3.zv = v1.xv*v2.yv - v1.yv*v2.xv;
}


inline static double torsion_xvec(const Xvec &x1, const Xvec &x2, const Xvec &x3, const Xvec &x4){

	Xvec b1 = x2 - x1;
    Xvec b2 = x3 - x2;
    Xvec b3 = x4 - x3;

    Xvec n1, n2, m1;
    n1 = b1.cross_product(b2);
    n2 = b2.cross_product(b3);
    m1 = n1.cross_product(b2);

    // cross_product(b1, b2, n1);
    // cross_product(b2, b3, n2);
    // cross_product(n1, b2, m1);
    n1.normalize();
    n2.normalize();
    m1.normalize();

    double x = n1*n2;
    double y = m1*n2;

    return atan2(y,x) * 180 /3.1415926;

}

#endif
