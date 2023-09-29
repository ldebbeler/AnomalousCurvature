# pragma once

// double3 class
struct double3 {
	//double x[2];
	double x;
	double y;
    double z;
	double3() { }
	double3(double xx, double yy, double zz) : x{xx}, y{yy}, z{zz}{ }
	double3(const double3& o) : x {o.x}, y{o.y}, z{o.z} { }
	~double3() { }
	double3& operator=(const double3& o)  { x = o.x; y = o.y; z=o.z; return *this; }
	double3& operator+=(const double3& o) { x += o.x; y += o.y; z += o.z; return *this; }
	double3& operator-=(const double3& o)	{ x -= o.x; y -= o.y; z -= o.z; return *this; }
    // are these operators necessary?
	double& operator[](int i) 
	{ return (i==0)?x:y; }
	double const& operator[](int i) const
	{ return (i==0)?x:y; }

	operator int() const = delete;
};
static inline double3 operator*(const double& a, const double3& v) 
{ return double3(v[0]*a,v[1]*a,v[2]*a); }
static inline double3 operator*(const double3& v, const double& a) 
{ return double3(v[0]*a,v[1]*a,v[2]*a); }
static inline double3 operator/(const double& a, const double3& v) 
{ return double3(v[0]/a,v[1]/a,v[2]/a); }
static inline double3 operator/(const double3& v, const double& a) 
{ return double3(v[0]/a,v[1]/a,v[2]/a); }

static inline double3 operator+(const double3& a, const double3& b)
{ 
	double3 ret = a;
	ret += b;
	return ret;
}

static inline double3 operator-(const double3& a, const double3& b)
{
	double3 ret = a;
	ret -= b;
	return ret;
}

static inline double3 operator-(const double3& b)
{ return operator-(double3(0.,0.,0.),b); }

/* scalar product */
static inline double operator*(const double3& a, const double3& b)
{ return a.x*b.x + a.y*b.y + a.z*b.z; }
