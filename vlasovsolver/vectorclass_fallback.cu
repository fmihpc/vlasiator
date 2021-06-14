
/*

template <class T>
static CUDA_HOSTDEV inline Vec4Simple<T> abs(const Vec4Simple<T> &l){
   return Vec4Simple<T>(
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
*/
// Dummy functions that Agner vectorclass has.
// These are here to suppress compiler error messages only
//static CUDA_HOSTDEV void no_subnormals() {}

