# pragma once

// double2 class
struct double2 {
	//double x[2];
	double x;
	double y;
	double2() { }
	double2(double xx, double yy) : x{xx}, y{yy}{ }
	double2(const double2& o) : x {o.x}, y{o.y} { }
	~double2() { }
	double2& operator=(const double2& o)  { x = o.x; y = o.y; return *this; }
	double2& operator+=(const double2& o) { x += o.x; y += o.y; return *this; }
	double2& operator-=(const double2& o)	{ x -= o.x; y -= o.y; return *this; }
	double& operator[](int i) 
	{ return (i==0)?x:y; }
	double const& operator[](int i) const
	{ return (i==0)?x:y; }

	operator int() const = delete;
};
static inline double2 operator*(const double& a, const double2& v) 
{ return double2(v[0]*a,v[1]*a); }
static inline double2 operator*(const double2& v, const double& a) 
{ return double2(v[0]*a,v[1]*a); }
static inline double2 operator/(const double& a, const double2& v) 
{ return double2(v[0]/a,v[1]/a); }
static inline double2 operator/(const double2& v, const double& a) 
{ return double2(v[0]/a,v[1]/a); }

static inline double2 operator+(const double2& a, const double2& b)
{ 
	double2 ret = a;
	ret += b;
	return ret;
}

static inline double2 operator-(const double2& a, const double2& b)
{
	double2 ret = a;
	ret -= b;
	return ret;
}

static inline double2 operator-(const double2& b)
{ return operator-(double2(0.,0.),b); }

/* scalar product */
static inline double operator*(const double2& a, const double2& b)
{ return a.x*b.x + a.y*b.y; }
