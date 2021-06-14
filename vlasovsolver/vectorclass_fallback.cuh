#include "../vlasovsolver_cuda/cuda_header.cuh"
template <typename T>
class Vec4Simple
{
  //T *_p;
  public:
    T val[4] __attribute__((aligned(32)));
    CUDA_HOSTDEV Vec4Simple();
    //CUDA_HOSTDEV Vec4Simple(T *p);
    CUDA_HOSTDEV Vec4Simple(T x);
    CUDA_HOSTDEV Vec4Simple(T a,T b,T c,T d);
    CUDA_HOSTDEV Vec4Simple(Vec4Simple const &x);
    CUDA_HOSTDEV Vec4Simple<T> & load(T const * p);
    CUDA_HOSTDEV Vec4Simple<T> & load_a(T const * p);
    CUDA_HOSTDEV Vec4Simple<T> & insert(int i,T const &x);
    CUDA_HOSTDEV void store(T * p) const;
    CUDA_HOSTDEV void store_a(T * p) const;
    CUDA_HOSTDEV Vec4Simple<T> & operator = (Vec4Simple<T> const & r);
    CUDA_HOSTDEV T operator [](int i) const;
    CUDA_HOSTDEV Vec4Simple<T> operator++ (int);
    static CUDA_HOSTDEV T getSquare(T b);
};
static CUDA_HOSTDEV void no_subnormals(){};
/*
static CUDA_HOSTDEV void no_subnormals();
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> abs(const Vec4Simple<T> &l);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> sqrt(const Vec4Simple<T> &l);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator + (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator + (const S &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator + (const Vec4Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator - (const Vec4Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator - (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator - (const S &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator - (const Vec4Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator * (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator * (const Vec4Simple<T> &l, const S &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator * (const S &l,const Vec4Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator / (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator / (const Vec4Simple<T> &l, const S &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator / (const S &l, const Vec4Simple<T> &r );
template <class T>
static CUDA_HOSTDEV inline  Vec4Simple<T> & operator += (Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline  Vec4Simple<T> & operator += (Vec4Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> & operator -= (Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> & operator -= (Vec4Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator || (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator && (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator == (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator == (const Vec4Simple<T> &l, const S& r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator != (const Vec4Simple<T> &l, const S& r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator ! (const Vec4Simple<T> &l);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator > (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator > (const Vec4Simple<T> &l, const S r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator > (const S l,const Vec4Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator >= (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator >= (const Vec4Simple<T> &l, const S r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator >= (const S l,const Vec4Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator < (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator < (const Vec4Simple<T> &l,const S &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator < (const S l, const Vec4Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator <= (const Vec4Simple<T> &l, const Vec4Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator <= (const Vec4Simple<T> &l,const S &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator <= (const S l, const Vec4Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> min(Vec4Simple<T> const & l, Vec4Simple<T> const & r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> min(S const & l, Vec4Simple<T> const & r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> max(Vec4Simple<T> const & l, Vec4Simple<T> const & r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> max(Vec4Simple<T> const & l, S const & r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> max(S const & l, Vec4Simple<T> const & r);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> select(Vec4Simple<bool> const & a, Vec4Simple<T> const & b, Vec4Simple<T> const & c);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> select(Vec4Simple<bool> const & a, S const & b, Vec4Simple<T> const & c);
template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> select(Vec4Simple<bool> const & a, Vec4Simple<T> const & b, S const & c);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> select(Vec4Simple<bool> const & a, T const & b, T const & c);
template <class T>
static CUDA_HOSTDEV inline bool horizontal_or(Vec4Simple<T> const & a);
template <class T>
static CUDA_HOSTDEV inline bool horizontal_and(Vec4Simple<T> const & a);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<int> truncate_to_int(Vec4Simple<T> const & a);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<double> to_double(Vec4Simple<T> const & a);
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<float> to_float(Vec4Simple<T> const & a);
*/

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> abs(const Vec4Simple<T> &l)
{
   return Vec4Simple<T>
   (
      fabs(l.val[0]),
      fabs(l.val[1]),
      fabs(l.val[2]),
      fabs(l.val[3])
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> sqrt(const Vec4Simple<T> &l){
   return Vec4Simple<T>(
      sqrt(l.val[0]),
      sqrt(l.val[1]),
      sqrt(l.val[2]),
      sqrt(l.val[3])
   );
}



template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator + (const Vec4Simple<T> &l, const Vec4Simple<T> &r){
   return Vec4Simple<T>(
      l.val[0]+r.val[0],
      l.val[1]+r.val[1],
      l.val[2]+r.val[2],
      l.val[3]+r.val[3]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator + (const S &l, const Vec4Simple<T> &r){
   return Vec4Simple<T>(
      l+r.val[0],
      l+r.val[1],
      l+r.val[2],
      l+r.val[3]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator + (const Vec4Simple<T> &l, const S &r){
   return Vec4Simple<T>(
      l.val[0]+r,
      l.val[1]+r,
      l.val[2]+r,
      l.val[3]+r
   );
}
template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator - (const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      -r.val[0],
      -r.val[1],
      -r.val[2],
      -r.val[3]
   );
}




template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator - (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      l.val[0]-r.val[0],
      l.val[1]-r.val[1],
      l.val[2]-r.val[2],
      l.val[3]-r.val[3]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator - (const S &l, const Vec4Simple<T> &r){
   return Vec4Simple<T>(
      l-r.val[0],
      l-r.val[1],
      l-r.val[2],
      l-r.val[3]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator - (const Vec4Simple<T> &l, const S &r){
   return Vec4Simple<T>(
      l.val[0]-r,
      l.val[1]-r,
      l.val[2]-r,
      l.val[3]-r
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator * (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      l.val[0]*r.val[0],
      l.val[1]*r.val[1],
      l.val[2]*r.val[2],
      l.val[3]*r.val[3]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator * (const Vec4Simple<T> &l, const S &r)
{
   return Vec4Simple<T>(
      l.val[0]*r,
      l.val[1]*r,
      l.val[2]*r,
      l.val[3]*r
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator * (const S &l,const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      l*r.val[0],
      l*r.val[1],
      l*r.val[2],
      l*r.val[3]
   );
}



template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> operator / (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      l.val[0]/r.val[0],
      l.val[1]/r.val[1],
      l.val[2]/r.val[2],
      l.val[3]/r.val[3]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator / (const Vec4Simple<T> &l, const S &r)
{
   return Vec4Simple<T>(
      l.val[0]/r,
      l.val[1]/r,
      l.val[2]/r,
      l.val[3]/r
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> operator / (const S &l, const Vec4Simple<T> &r )
{
   return Vec4Simple<T>(
      l/r.val[0],
      l/r.val[1],
      l/r.val[2],
      l/r.val[3]
   );
}

template <class T>
static CUDA_HOSTDEV inline  Vec4Simple<T> & operator += (Vec4Simple<T> &l, const Vec4Simple<T> &r){
   l=l+r;
   return l;
}

template <class T, class S>
static CUDA_HOSTDEV inline  Vec4Simple<T> & operator += (Vec4Simple<T> &l, const S &r){
   l = l+r;
   return l;
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> & operator -= (Vec4Simple<T> &l, const Vec4Simple<T> &r){
   l=l-r;
   return l;
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> & operator -= (Vec4Simple<T> &l, const S &r){
   l = l - r;
   return l;
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator || (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] || r.val[0],
      l.val[1] || r.val[1],
      l.val[2] || r.val[2],
      l.val[3] || r.val[3]
   );
}



template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator && (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] && r.val[0],
      l.val[1] && r.val[1],
      l.val[2] && r.val[2],
      l.val[3] && r.val[3]
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator == (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] == r.val[0],
      l.val[1] == r.val[1],
      l.val[2] == r.val[2],
      l.val[3] == r.val[3]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator == (const Vec4Simple<T> &l, const S& r)
{
   return Vec4Simple<bool>(
      l.val[0] == r,
      l.val[1] == r,
      l.val[2] == r,
      l.val[3] == r
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator != (const Vec4Simple<T> &l, const S& r)
{
   return Vec4Simple<bool>(
      l.val[0] != r,
      l.val[1] != r,
      l.val[2] != r,
      l.val[3] != r
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator ! (const Vec4Simple<T> &l)
{
   return Vec4Simple<bool>(
      !l.val[0],
      !l.val[1],
      !l.val[2],
      !l.val[3]
   );
}


template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator > (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] > r.val[0],
      l.val[1] > r.val[1],
      l.val[2] > r.val[2],
      l.val[3] > r.val[3]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator > (const Vec4Simple<T> &l, const S r)
{
   return Vec4Simple<bool>(
      l.val[0] > r,
      l.val[1] > r,
      l.val[2] > r,
      l.val[3] > r
   );
}



template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator > (const S l,const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l > r.val[0],
      l > r.val[1],
      l > r.val[2],
      l > r.val[3]
   );
}


template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator >= (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] >= r.val[0],
      l.val[1] >= r.val[1],
      l.val[2] >= r.val[2],
      l.val[3] >= r.val[3]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator >= (const Vec4Simple<T> &l, const S r)
{
   return Vec4Simple<bool>(
      l.val[0] >= r,
      l.val[1] >= r,
      l.val[2] >= r,
      l.val[3] >= r
   );
}



template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator >= (const S l,const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l >= r.val[0],
      l >= r.val[1],
      l >= r.val[2],
      l >= r.val[3]
   );
}



template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator < (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] < r.val[0],
      l.val[1] < r.val[1],
      l.val[2] < r.val[2],
      l.val[3] < r.val[3]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator < (const Vec4Simple<T> &l,const S &r)
{
   return Vec4Simple<bool>(
      l.val[0] < r,
      l.val[1] < r,
      l.val[2] < r,
      l.val[3] < r
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator < (const S l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l < r.val[0],
      l < r.val[1],
      l < r.val[2],
      l < r.val[3]
   );
}



template <class T>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator <= (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] <= r.val[0],
      l.val[1] <= r.val[1],
      l.val[2] <= r.val[2],
      l.val[3] <= r.val[3]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator <= (const Vec4Simple<T> &l,const S &r)
{
   return Vec4Simple<bool>(
      l.val[0] <= r,
      l.val[1] <= r,
      l.val[2] <= r,
      l.val[3] <= r
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<bool> operator <= (const S l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l <= r.val[0],
      l <= r.val[1],
      l <= r.val[2],
      l <= r.val[3]
   );
}




template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> min(Vec4Simple<T> const & l, Vec4Simple<T> const & r){
   return Vec4Simple<T>(
      l.val[0]<r.val[0]?l.val[0]:r.val[0],
      l.val[1]<r.val[1]?l.val[1]:r.val[1],
      l.val[2]<r.val[2]?l.val[2]:r.val[2],
      l.val[3]<r.val[3]?l.val[3]:r.val[3]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> min(S const & l, Vec4Simple<T> const & r){
   return Vec4Simple<T>(
      l <r.val[0] ? l:r.val[0],
      l <r.val[1] ? l:r.val[1],
      l <r.val[2] ? l:r.val[2],
      l <r.val[3] ? l:r.val[3]
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> max(Vec4Simple<T> const & l, Vec4Simple<T> const & r){
   return Vec4Simple<T>(
      l.val[0]>r.val[0]?l.val[0]:r.val[0],
      l.val[1]>r.val[1]?l.val[1]:r.val[1],
      l.val[2]>r.val[2]?l.val[2]:r.val[2],
      l.val[3]>r.val[3]?l.val[3]:r.val[3]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> max(Vec4Simple<T> const & l, S const & r){
   return Vec4Simple<T>(
      l.val[0] > r ? l.val[0] : r,
      l.val[1] > r ? l.val[1] : r,
      l.val[2] > r ? l.val[2] : r,
      l.val[3] > r ? l.val[3] : r
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> max(S const & l, Vec4Simple<T> const & r){
  return Vec4Simple<T>(
     r.val[0] > l ? r.val[0] : l,
     r.val[1] > l ? r.val[1] : l,
     r.val[2] > l ? r.val[2] : l,
     r.val[3] > l ? r.val[3] : l
  );
}



template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> select(Vec4Simple<bool> const & a, Vec4Simple<T> const & b, Vec4Simple<T> const & c){
   return Vec4Simple<T>(
      a.val[0] ? b.val[0] : c.val[0],
      a.val[1] ? b.val[1] : c.val[1],
      a.val[2] ? b.val[2] : c.val[2],
      a.val[3] ? b.val[3] : c.val[3]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> select(Vec4Simple<bool> const & a, S const & b, Vec4Simple<T> const & c){
   return Vec4Simple<T>(
      a.val[0] ? b : c.val[0],
      a.val[1] ? b : c.val[1],
      a.val[2] ? b : c.val[2],
      a.val[3] ? b : c.val[3]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec4Simple<T> select(Vec4Simple<bool> const & a, Vec4Simple<T> const & b, S const & c){
   return Vec4Simple<T>(
      a.val[0] ? b.val[0] : c,
      a.val[1] ? b.val[1] : c,
      a.val[2] ? b.val[2] : c,
      a.val[3] ? b.val[3] : c
   );
}


template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> select(Vec4Simple<bool> const & a, T const & b, T const & c){
   return Vec4Simple<T>(
      a.val[0] ? b : c,
      a.val[1] ? b : c,
      a.val[2] ? b : c,
      a.val[3] ? b : c
   );
}

template <class T>
static CUDA_HOSTDEV inline bool horizontal_or(Vec4Simple<T> const & a){
  return a.val[0] || a.val[1] || a.val[2] || a.val[3];
}


template <class T>
static CUDA_HOSTDEV inline bool horizontal_and(Vec4Simple<T> const & a){
  return a.val[0] && a.val[1] && a.val[2] && a.val[3];
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<int> truncate_to_int(Vec4Simple<T> const & a){
  return Vec4Simple<int>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<double> to_double(Vec4Simple<T> const & a){
  return Vec4Simple<double>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<float> to_float(Vec4Simple<T> const & a){
  return Vec4Simple<float>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

/*
template <typename T>
static CUDA_HOSTDEV Vec4Simple<T> truncate_to_int(Vec4Simple<T> const & a);
template <typename T>
static CUDA_HOSTDEV Vec4Simple<T> to_double(Vec4Simple<T> const & a);
*/
template <typename T>
class Vec8Simple
{
  //T *_p;
  public:
    T val[8] __attribute__((aligned(32)));
    CUDA_HOSTDEV Vec8Simple();
    CUDA_HOSTDEV Vec8Simple(T x);
    CUDA_HOSTDEV Vec8Simple(T a,T b,T c,T d, T e,T f,T g,T h);
    CUDA_HOSTDEV Vec8Simple(Vec8Simple const &x);
    CUDA_HOSTDEV Vec8Simple<T> & load(T const * p);
    CUDA_HOSTDEV Vec8Simple<T> & load_a(T const * p);
    CUDA_HOSTDEV Vec8Simple<T> & insert(int i,T const &x);
    CUDA_HOSTDEV void store(T * p) const;
    CUDA_HOSTDEV void store_a(T * p) const;
    CUDA_HOSTDEV Vec8Simple<T> & operator = (Vec8Simple<T> const & r);
    CUDA_HOSTDEV T operator [](int i) const;
    CUDA_HOSTDEV Vec8Simple<T> operator++ (int);
    static CUDA_HOSTDEV T getCube(T b);
};
/*
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> abs(const Vec8Simple<T> &l);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> sqrt(const Vec8Simple<T> &l);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator + (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator + (const S &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator + (const Vec8Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator - (const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator - (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator - (const S &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator - (const Vec8Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator * (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator * (const Vec8Simple<T> &l, const S &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator * (const S &l,const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator / (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator / (const Vec8Simple<T> &l, const S &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator / (const S &l, const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline  Vec8Simple<T> & operator += (Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline  Vec8Simple<T> & operator += (Vec8Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> & operator -= (Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> & operator -= (Vec8Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator || (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator && (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator == (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator == (const Vec8Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator != (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator != (const Vec8Simple<T> &l, const S &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator ! (const Vec8Simple<T> &l);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator > (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator > (const Vec8Simple<T> &l, const S r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator > (const S l,const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator >= (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator >= (const Vec8Simple<T> &l, const S r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator >= (const S l,const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator < (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator < (const Vec8Simple<T> &l,const S &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator < (const S l, const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator <= (const Vec8Simple<T> &l, const Vec8Simple<T> &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator <= (const Vec8Simple<T> &l,const S &r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator <= (const S l, const Vec8Simple<T> &r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> min(Vec8Simple<T> const & l, Vec8Simple<T> const & r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> min(S const & l, Vec8Simple<T> const & r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> max(Vec8Simple<T> const & l, Vec8Simple<T> const & r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> max(Vec8Simple<T> const & l, S const & r);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> max(S const & l, Vec8Simple<T> const & r);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> select(Vec8Simple<bool> const & a, Vec8Simple<T> const & b, Vec8Simple<T> const & c);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> select(Vec8Simple<bool> const & a, S const & b, Vec8Simple<T> const & c);
template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> select(Vec8Simple<bool> const & a, Vec8Simple<T> const & b, S const & c);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> select(Vec8Simple<bool> const & a, T const & b, T const & c);
template <class T>
static CUDA_HOSTDEV inline bool horizontal_or(Vec8Simple<T> const & a);
template <class T>
static CUDA_HOSTDEV inline bool horizontal_and(Vec8Simple<T> const & a);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<int> truncate_to_int(Vec8Simple<T> const & a);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<double> to_double(Vec8Simple<T> const & a);
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<float> to_float(Vec8Simple<T> const & a){
   return Vec8Simple<float>(a.val[0], a.val[1], a.val[2], a.val[3],
                            a.val[4], a.val[5], a.val[6], a.val[7]);
}
*/
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> abs(const Vec8Simple<T> &l)
{
   return Vec8Simple<T>
   (
      fabs(l.val[0]),
      fabs(l.val[1]),
      fabs(l.val[2]),
      fabs(l.val[3]),
      fabs(l.val[4]),
      fabs(l.val[5]),
      fabs(l.val[6]),
      fabs(l.val[7])
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> sqrt(const Vec8Simple<T> &l){
   return Vec8Simple<T>(
      sqrt(l.val[0]),
      sqrt(l.val[1]),
      sqrt(l.val[2]),
      sqrt(l.val[3]),
      sqrt(l.val[4]),
      sqrt(l.val[5]),
      sqrt(l.val[6]),
      sqrt(l.val[7])
   );
}



template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator + (const Vec8Simple<T> &l, const Vec8Simple<T> &r){
   return Vec8Simple<T>(
      l.val[0]+r.val[0],
      l.val[1]+r.val[1],
      l.val[2]+r.val[2],
      l.val[3]+r.val[3],
      l.val[4]+r.val[4],
      l.val[5]+r.val[5],
      l.val[6]+r.val[6],
      l.val[7]+r.val[7]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator + (const S &l, const Vec8Simple<T> &r){
   return Vec8Simple<T>(
      l+r.val[0],
      l+r.val[1],
      l+r.val[2],
      l+r.val[3],
      l+r.val[4],
      l+r.val[5],
      l+r.val[6],
      l+r.val[7]

   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator + (const Vec8Simple<T> &l, const S &r){
   return Vec8Simple<T>(
      l.val[0]+r,
      l.val[1]+r,
      l.val[2]+r,
      l.val[3]+r,
      l.val[4]+r,
      l.val[5]+r,
      l.val[6]+r,
      l.val[7]+r

   );
}
template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator - (const Vec8Simple<T> &r)
{
   return Vec8Simple<T>(
      -r.val[0],
      -r.val[1],
      -r.val[2],
      -r.val[3],
      -r.val[4],
      -r.val[5],
      -r.val[6],
      -r.val[7]

   );
}




template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator - (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<T>(
      l.val[0]-r.val[0],
      l.val[1]-r.val[1],
      l.val[2]-r.val[2],
      l.val[3]-r.val[3],
      l.val[4]-r.val[4],
      l.val[5]-r.val[5],
      l.val[6]-r.val[6],
      l.val[7]-r.val[7]

   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator - (const S &l, const Vec8Simple<T> &r){
   return Vec8Simple<T>(
      l-r.val[0],
      l-r.val[1],
      l-r.val[2],
      l-r.val[3],
      l-r.val[4],
      l-r.val[5],
      l-r.val[6],
      l-r.val[7]
      );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator - (const Vec8Simple<T> &l, const S &r){
   return Vec8Simple<T>(
      l.val[0]-r,
      l.val[1]-r,
      l.val[2]-r,
      l.val[3]-r,
      l.val[4]-r,
      l.val[5]-r,
      l.val[6]-r,
      l.val[7]-r

   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator * (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<T>(
      l.val[0]*r.val[0],
      l.val[1]*r.val[1],
      l.val[2]*r.val[2],
      l.val[3]*r.val[3],
      l.val[4]*r.val[4],
      l.val[5]*r.val[5],
      l.val[6]*r.val[6],
      l.val[7]*r.val[7]

   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator * (const Vec8Simple<T> &l, const S &r)
{
   return Vec8Simple<T>
   (
      l.val[0]*r,
      l.val[1]*r,
      l.val[2]*r,
      l.val[3]*r,
      l.val[4]*r,
      l.val[5]*r,
      l.val[6]*r,
      l.val[7]*r
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator * (const S &l,const Vec8Simple<T> &r)
{
   return Vec8Simple<T>(
      r.val[0]*l,
      r.val[1]*l,
      r.val[2]*l,
      r.val[3]*l,
      r.val[4]*l,
      r.val[5]*l,
      r.val[6]*l,
      r.val[7]*l

   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> operator / (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<T>(
      l.val[0]/r.val[0],
      l.val[1]/r.val[1],
      l.val[2]/r.val[2],
      l.val[3]/r.val[3],
      l.val[4]/r.val[4],
      l.val[5]/r.val[5],
      l.val[6]/r.val[6],
      l.val[7]/r.val[7]

   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator / (const Vec8Simple<T> &l, const S &r)
{
   return Vec8Simple<T>(
      l.val[0]/r,
      l.val[1]/r,
      l.val[2]/r,
      l.val[3]/r,
      l.val[4]/r,
      l.val[5]/r,
      l.val[6]/r,
      l.val[7]/r

   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> operator / (const S &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<T>(
      l/r.val[0],
      l/r.val[1],
      l/r.val[2],
      l/r.val[3],
      l/r.val[4],
      l/r.val[5],
      l/r.val[6],
      l/r.val[7]
   );
}

template <class T>
static CUDA_HOSTDEV inline  Vec8Simple<T> & operator += (Vec8Simple<T> &l, const Vec8Simple<T> &r){
   l=l+r;
   return l;
}

template <class T, class S>
static CUDA_HOSTDEV inline  Vec8Simple<T> & operator += (Vec8Simple<T> &l, const S &r){
   l = l+r;
   return l;
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> & operator -= (Vec8Simple<T> &l, const Vec8Simple<T> &r){
   l=l-r;
   return l;
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> & operator -= (Vec8Simple<T> &l, const S &r){
   l = l - r;
   return l;
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator || (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l.val[0] || r.val[0],
      l.val[1] || r.val[1],
      l.val[2] || r.val[2],
      l.val[3] || r.val[3],
      l.val[4] || r.val[4],
      l.val[5] || r.val[5],
      l.val[6] || r.val[6],
      l.val[7] || r.val[7]

   );
}



template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator && (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l.val[0] && r.val[0],
      l.val[1] && r.val[1],
      l.val[2] && r.val[2],
      l.val[3] && r.val[3],
      l.val[4] && r.val[4],
      l.val[5] && r.val[5],
      l.val[6] && r.val[6],
      l.val[7] && r.val[7]

   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator == (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l.val[0] == r.val[0],
      l.val[1] == r.val[1],
      l.val[2] == r.val[2],
      l.val[3] == r.val[3],
      l.val[4] == r.val[4],
      l.val[5] == r.val[5],
      l.val[6] == r.val[6],
      l.val[7] == r.val[7]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator == (const Vec8Simple<T> &l, const S &r)
{
   return Vec8Simple<bool>(
      l.val[0] == r,
      l.val[1] == r,
      l.val[2] == r,
      l.val[3] == r,
      l.val[4] == r,
      l.val[5] == r,
      l.val[6] == r,
      l.val[7] == r
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator != (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l.val[0] != r.val[0],
      l.val[1] != r.val[1],
      l.val[2] != r.val[2],
      l.val[3] != r.val[3],
      l.val[4] != r.val[4],
      l.val[5] != r.val[5],
      l.val[6] != r.val[6],
      l.val[7] != r.val[7]
   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator != (const Vec8Simple<T> &l, const S &r)
{
   return Vec8Simple<bool>(
      l.val[0] != r,
      l.val[1] != r,
      l.val[2] != r,
      l.val[3] != r,
      l.val[4] != r,
      l.val[5] != r,
      l.val[6] != r,
      l.val[7] != r
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator ! (const Vec8Simple<T> &l)
{
   return Vec8Simple<bool>(
      !l.val[0],
      !l.val[1],
      !l.val[2],
      !l.val[3],
      !l.val[4],
      !l.val[5],
      !l.val[6],
      !l.val[7]
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator > (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l.val[0] > r.val[0],
      l.val[1] > r.val[1],
      l.val[2] > r.val[2],
      l.val[3] > r.val[3],
      l.val[4] > r.val[4],
      l.val[5] > r.val[5],
      l.val[6] > r.val[6],
      l.val[7] > r.val[7]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator > (const Vec8Simple<T> &l, const S r)
{
   return Vec8Simple<bool>(
      l.val[0] > r,
      l.val[1] > r,
      l.val[2] > r,
      l.val[3] > r,
      l.val[4] > r,
      l.val[5] > r,
      l.val[6] > r,
      l.val[7] > r

   );
}



template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator > (const S l,const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l > r.val[0],
      l > r.val[1],
      l > r.val[2],
      l > r.val[3],
      l > r.val[4],
      l > r.val[5],
      l > r.val[6],
      l > r.val[7]

   );
}


template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator >= (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l.val[0] >= r.val[0],
      l.val[1] >= r.val[1],
      l.val[2] >= r.val[2],
      l.val[3] >= r.val[3],
      l.val[4] >= r.val[4],
      l.val[5] >= r.val[5],
      l.val[6] >= r.val[6],
      l.val[7] >= r.val[7]

   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator >= (const Vec8Simple<T> &l, const S r)
{
   return Vec8Simple<bool>(
      l.val[0] >= r,
      l.val[1] >= r,
      l.val[2] >= r,
      l.val[3] >= r,
      l.val[4] >= r,
      l.val[5] >= r,
      l.val[6] >= r,
      l.val[7] >= r

   );
}



template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator >= (const S l,const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l >= r.val[0],
      l >= r.val[1],
      l >= r.val[2],
      l >= r.val[3],
      l >= r.val[4],
      l >= r.val[5],
      l >= r.val[6],
      l >= r.val[7]

   );
}



template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator < (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l.val[0] < r.val[0],
      l.val[1] < r.val[1],
      l.val[2] < r.val[2],
      l.val[3] < r.val[3],
      l.val[4] < r.val[4],
      l.val[5] < r.val[5],
      l.val[6] < r.val[6],
      l.val[7] < r.val[7]

   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator < (const Vec8Simple<T> &l,const S &r)
{
   return Vec8Simple<bool>(
      l.val[0] < r,
      l.val[1] < r,
      l.val[2] < r,
      l.val[3] < r,
      l.val[4] < r,
      l.val[5] < r,
      l.val[6] < r,
      l.val[7] < r

   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator < (const S l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l < r.val[0],
      l < r.val[1],
      l < r.val[2],
      l < r.val[3],
      l < r.val[4],
      l < r.val[5],
      l < r.val[6],
      l < r.val[7]

   );
}



template <class T>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator <= (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l.val[0] <= r.val[0],
      l.val[1] <= r.val[1],
      l.val[2] <= r.val[2],
      l.val[3] <= r.val[3],
      l.val[4] <= r.val[4],
      l.val[5] <= r.val[5],
      l.val[6] <= r.val[6],
      l.val[7] <= r.val[7]

   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator <= (const Vec8Simple<T> &l,const S &r)
{
   return Vec8Simple<bool>(
      l.val[0] <= r,
      l.val[1] <= r,
      l.val[2] <= r,
      l.val[3] <= r,
      l.val[4] <= r,
      l.val[5] <= r,
      l.val[6] <= r,
      l.val[7] <= r

   );
}

template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<bool> operator <= (const S l, const Vec8Simple<T> &r)
{
   return Vec8Simple<bool>(
      l <= r.val[0],
      l <= r.val[1],
      l <= r.val[2],
      l <= r.val[3],
      l <= r.val[4],
      l <= r.val[5],
      l <= r.val[6],
      l <= r.val[7]

   );
}


template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> min(Vec8Simple<T> const & l, Vec8Simple<T> const & r){
   return Vec8Simple<T>(
      l.val[0]<r.val[0]?l.val[0]:r.val[0],
      l.val[1]<r.val[1]?l.val[1]:r.val[1],
      l.val[2]<r.val[2]?l.val[2]:r.val[2],
      l.val[3]<r.val[3]?l.val[3]:r.val[3],
      l.val[4]<r.val[4]?l.val[4]:r.val[4],
      l.val[5]<r.val[5]?l.val[5]:r.val[5],
      l.val[6]<r.val[6]?l.val[6]:r.val[6],
      l.val[7]<r.val[7]?l.val[7]:r.val[7]

   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> min(S const & l, Vec8Simple<T> const & r){
   return Vec8Simple<T>(
      l<r.val[0]?l:r.val[0],
      l<r.val[1]?l:r.val[1],
      l<r.val[2]?l:r.val[2],
      l<r.val[3]?l:r.val[3],
      l<r.val[4]?l:r.val[4],
      l<r.val[5]?l:r.val[5],
      l<r.val[6]?l:r.val[6],
      l<r.val[7]?l:r.val[7]
   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> max(Vec8Simple<T> const & l, Vec8Simple<T> const & r){
   return Vec8Simple<T>(
      l.val[0]>r.val[0]?l.val[0]:r.val[0],
      l.val[1]>r.val[1]?l.val[1]:r.val[1],
      l.val[2]>r.val[2]?l.val[2]:r.val[2],
      l.val[3]>r.val[3]?l.val[3]:r.val[3],
      l.val[4]>r.val[4]?l.val[4]:r.val[4],
      l.val[5]>r.val[5]?l.val[5]:r.val[5],
      l.val[6]>r.val[6]?l.val[6]:r.val[6],
      l.val[7]>r.val[7]?l.val[7]:r.val[7]

   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> max(Vec8Simple<T> const & l, S const & r){
   return Vec8Simple<T>(
      l.val[0] > r ? l.val[0] : r,
      l.val[1] > r ? l.val[1] : r,
      l.val[2] > r ? l.val[2] : r,
      l.val[3] > r ? l.val[3] : r,
      l.val[4] > r ? l.val[4] : r,
      l.val[5] > r ? l.val[5] : r,
      l.val[6] > r ? l.val[6] : r,
      l.val[7] > r ? l.val[7] : r

   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> max(S const & l, Vec8Simple<T> const & r){
  return Vec8Simple<T>(
     r.val[0] > l ? r.val[0] : l,
     r.val[1] > l ? r.val[1] : l,
     r.val[2] > l ? r.val[2] : l,
     r.val[3] > l ? r.val[3] : l,
     r.val[4] > l ? r.val[4] : l,
     r.val[5] > l ? r.val[5] : l,
     r.val[6] > l ? r.val[6] : l,
     r.val[7] > l ? r.val[7] : l

  );
}



template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> select(Vec8Simple<bool> const & a, Vec8Simple<T> const & b, Vec8Simple<T> const & c){
   return Vec8Simple<T>(
      a.val[0] ? b.val[0] : c.val[0],
      a.val[1] ? b.val[1] : c.val[1],
      a.val[2] ? b.val[2] : c.val[2],
      a.val[3] ? b.val[3] : c.val[3],
      a.val[4] ? b.val[4] : c.val[4],
      a.val[5] ? b.val[5] : c.val[5],
      a.val[6] ? b.val[6] : c.val[6],
      a.val[7] ? b.val[7] : c.val[7]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> select(Vec8Simple<bool> const & a, S const & b, Vec8Simple<T> const & c){
   return Vec8Simple<T>(
      a.val[0] ? b : c.val[0],
      a.val[1] ? b : c.val[1],
      a.val[2] ? b : c.val[2],
      a.val[3] ? b : c.val[3],
      a.val[4] ? b : c.val[4],
      a.val[5] ? b : c.val[5],
      a.val[6] ? b : c.val[6],
      a.val[7] ? b : c.val[7]
   );
}


template <class T, class S>
static CUDA_HOSTDEV inline Vec8Simple<T> select(Vec8Simple<bool> const & a, Vec8Simple<T> const & b, S const & c){
   return Vec8Simple<T>(
      a.val[0] ? b.val[0] : c,
      a.val[1] ? b.val[1] : c,
      a.val[2] ? b.val[2] : c,
      a.val[3] ? b.val[3] : c,
      a.val[4] ? b.val[4] : c,
      a.val[5] ? b.val[5] : c,
      a.val[6] ? b.val[6] : c,
      a.val[7] ? b.val[7] : c

   );
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<T> select(Vec8Simple<bool> const & a, T const & b, T const & c){
   return Vec8Simple<T>(
      a.val[0] ? b : c,
      a.val[1] ? b : c,
      a.val[2] ? b : c,
      a.val[3] ? b : c,
      a.val[4] ? b : c,
      a.val[5] ? b : c,
      a.val[6] ? b : c,
      a.val[7] ? b : c
      );
}

template <class T>
static CUDA_HOSTDEV inline bool horizontal_or(Vec8Simple<T> const & a){
  return a.val[0] || a.val[1] || a.val[2] || a.val[3] ||
     a.val[4] || a.val[5] || a.val[6] || a.val[7];
}


template <class T>
static CUDA_HOSTDEV inline bool horizontal_and(Vec8Simple<T> const & a){
   return a.val[0] && a.val[1] && a.val[2] && a.val[3] &&
      a.val[4] && a.val[5] && a.val[6] && a.val[7];
}



template <class T>
static CUDA_HOSTDEV inline Vec8Simple<int> truncate_to_int(Vec8Simple<T> const & a){
   return Vec8Simple<int>(a.val[0], a.val[1], a.val[2], a.val[3],
                          a.val[4], a.val[5], a.val[6], a.val[7]);

}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<double> to_double(Vec8Simple<T> const & a){
   return Vec8Simple<double>(a.val[0], a.val[1], a.val[2], a.val[3],
                             a.val[4], a.val[5], a.val[6], a.val[7]);
}

template <class T>
static CUDA_HOSTDEV inline Vec8Simple<float> to_float(Vec8Simple<T> const & a){
   return Vec8Simple<float>(a.val[0], a.val[1], a.val[2], a.val[3],
                            a.val[4], a.val[5], a.val[6], a.val[7]);
}
/*
template <typename T>
static CUDA_HOSTDEV Vec8Simple<T> truncate_to_int(Vec8Simple<T> const & a);
template <typename T>
static CUDA_HOSTDEV Vec8Simple<float> to_float(Vec8Simple<T> const & a);
*/
