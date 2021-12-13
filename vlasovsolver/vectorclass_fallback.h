/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef VECTORCLASS_PORTABLE_H
#define VECTORCLASS_PORTABLE_H
#include <math.h>
#include <iostream>
#include <initializer_list>

/*! \file vectorclass_fallback.h
  \brief Simple class for implementing a vector with 4 real values

*/

// Prefetching does nothing in the fallback vectorclass, if no system implementation
// is available
#ifndef _mm_prefetch
#define _mm_prefetch(...)
#endif



template <class T>
class VecSimple {
public:
   T val[VECL] __attribute__((aligned(32)));
   // donot initi v
   VecSimple() { }
   // Replicate scalar x across v.
   VecSimple(T x){
      for(unsigned int i=0;i<VECL;i++)
         val[i]=x;
   }

   // Pass vector values as an initializer list instead of a bunch of arguments.
   // usage: VecSimple<atype> newVec = {a, b, c, d} (for VECL values, this is checked against)
   // or   : Vec({a,b,c,d})
   VecSimple(std::initializer_list<T> list){
      // The rest of this is just sanity checking
      if(list.size() != VECL){
	 std::cerr <<  __FILE__ << ":" << __LINE__ <<
		 "Constructing a vector with a number of elements not equal to VECL = " << VECL <<
		 " (you had initializer_list size = "<<list.size()<<")";
	 abort();
      }
      else{
	 uint i = 0;
	 for(auto it = list.begin(); it != list.end(); ++it){
            val[i] = *it;
	    ++i;
	 }
	
      }
   }

   // The variadic-templated version of the above - this is a bit suspect with narrowing conversion warnings.
   //template<typename... Ts>
   //VecSimple(Ts... ts) : val{ts...}{}
   
   // Replicate VECL values across v.   
   VecSimple(T a,T b,T c,T d){
      if(VECL != 4) {
         std::cerr << __FILE__ << ":" << __LINE__ << ": Wrong function call for VECL=" << VECL << ", this function is for VECL=4!";
         abort();
      }
      val[0]=a;
      val[1]=b;
      val[2]=c;
      val[3]=d;
   }
   
   // Replicate VECL values across v.   
   VecSimple(T a,T b,T c,T d,T e,T f,T g,T h){
      if(VECL != 8) {
         std::cerr << __FILE__ << ":" << __LINE__ << ": Wrong function call for VECL=" << VECL << ", this function is for VECL=8!";
         abort();
      }
      val[0]=a;
      val[1]=b;
      val[2]=c;
      val[3]=d;
      val[4]=e;
      val[5]=f;
      val[6]=g;
      val[7]=h;
   }
   
   // Replicate VECL values across v.   
   VecSimple(T a,T b,T c,T d,T e,T f,T g,T h,T i,T j,T k,T l,T m,T n,T o,T p){
      if(VECL != 16) {
         std::cerr << __FILE__ << ":" << __LINE__ << ": Wrong function call for VECL=" << VECL << ", this function is for VECL=16!";
         abort();
      }
      val[0]=a;
      val[1]=b;
      val[2]=c;
      val[3]=d;
      val[4]=e;
      val[5]=f;
      val[6]=g;
      val[7]=h;
      val[8]=i;
      val[9]=j;
      val[10]=k;
      val[11]=l;
      val[12]=m;
      val[13]=n;
      val[14]=o;
      val[15]=p;
   }
   
   // Replicate VECL values across v.   
   VecSimple(T a,T b,T c,T d,T e,T f,T g,T h,T i,T j,T k,T l,T m,T n,T o,T p,T q,T r,T s,T t,T u,T v,T w,T x,T y,T z,T aa,T bb,T cc,T dd,T ee,T ff){
      if(VECL != 32) {
         std::cerr << __FILE__ << ":" << __LINE__ << ": Wrong function call for VECL=" << VECL << ", this function is for VECL=32!";
         abort();
      }
      val[0]=a;
      val[1]=b;
      val[2]=c;
      val[3]=d;
      val[4]=e;
      val[5]=f;
      val[6]=g;
      val[7]=h;
      val[8]=i;
      val[9]=j;
      val[10]=k;
      val[11]=l;
      val[12]=m;
      val[13]=n;
      val[14]=o;
      val[15]=p;
      val[16]=q;
      val[17]=r;
      val[18]=s;
      val[19]=t;
      val[20]=u;
      val[21]=v;
      val[22]=w;
      val[23]=x;
      val[24]=y;
      val[25]=z;
      val[26]=aa;
      val[27]=bb;
      val[28]=cc;
      val[29]=dd;
      val[30]=ee;
      val[31]=ff;
   }

   // Copy vector v.
   VecSimple(VecSimple const &x){
      for(unsigned int i=0;i<VECL;i++)
         val[i]=x.val[i];
   }
   
   // Member function to load from array (unaligned)
   VecSimple & load(T const * p)  {
      for(unsigned int i=0;i<VECL;i++)
         val[i]=p[i];
      return *this;
   }
   // Member function to load from array, aligned by 32
   VecSimple & load_a(T const * p){
      return this->load(p);
   }
   
   VecSimple & insert(int i,T const &x) {
      val[i]=x;
      return *this;
   }
   
   
   // Member function to store into array (unaligned)
   void store(T * p) const {
      for(unsigned int i=0;i<VECL;i++)
         p[i]=val[i];
   }
   // Member function to store into array, aligned by 32
   void store_a(T * p) const {
      this->store(p);
   }
   
   VecSimple & operator = (VecSimple const & r){
      for(unsigned int i=0;i<VECL;i++)
         val[i]=r.val[i];
      
      return *this;
   }
   
   T operator [](int i) const{
      return val[i];
   }
   
   VecSimple operator++ (int)
   {
      VecSimple<T> temp (*this);
      for(unsigned int i=0;i<VECL;i++)
         val[i]++;
      return temp;
   }
};

template <class T>
static inline VecSimple<T> abs(const VecSimple<T> &l){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, fabs(l.val[i]));
   return temp;
}

template <class T>
static inline VecSimple<T> sqrt(const VecSimple<T> &l){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, sqrt(l.val[i]));
   return temp;
}



template <class T>
static inline VecSimple<T> operator + (const VecSimple<T> &l, const VecSimple<T> &r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i]+r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> operator + (const S &l, const VecSimple<T> &r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l+r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> operator + (const VecSimple<T> &l, const S &r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i]+r);
   return temp;
}

template <class T>
static inline VecSimple<T> operator - (const VecSimple<T> &r)
{
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, -r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<T> operator - (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i]-r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> operator - (const S &l, const VecSimple<T> &r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l-r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> operator - (const VecSimple<T> &l, const S &r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i]-r);
   return temp;
}

template <class T>
static inline VecSimple<T> operator * (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i]*r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> operator * (const VecSimple<T> &l, const S &r)
{
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i]*r);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> operator * (const S &l,const VecSimple<T> &r)
{
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l*r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<T> operator / (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i]/r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> operator / (const VecSimple<T> &l, const S &r)
{
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i]/r);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> operator / (const S &l, const VecSimple<T> &r )
{
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l/r.val[i]);
   return temp;
}

template <class T>
static inline  VecSimple<T> & operator += (VecSimple<T> &l, const VecSimple<T> &r){
   l=l+r;
   return l;
}

template <class T, class S>
static inline  VecSimple<T> & operator += (VecSimple<T> &l, const S &r){
   l = l+r;
   return l;
}

template <class T>
static inline VecSimple<T> & operator -= (VecSimple<T> &l, const VecSimple<T> &r){
   l=l-r;
   return l;
}

template <class T, class S>   
static inline VecSimple<T> & operator -= (VecSimple<T> &l, const S &r){
   l = l - r;
   return l;
}

template <class T>
static inline VecSimple<bool> operator || (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] || r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<bool> operator && (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] && r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<bool> operator == (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] == r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator == (const VecSimple<T> &l, const S& r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] == r);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator != (const VecSimple<T> &l, const S& r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] != r);
   return temp;
}

template <class T>
static inline VecSimple<bool> operator ! (const VecSimple<T> &l)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, !l.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<bool> operator > (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] > r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator > (const VecSimple<T> &l, const S r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] > r);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator > (const S l,const VecSimple<T> &r) 
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l > r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<bool> operator >= (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] >= r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator >= (const VecSimple<T> &l, const S r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] >= r);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator >= (const S l,const VecSimple<T> &r) 
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l >= r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<bool> operator < (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] < r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator < (const VecSimple<T> &l,const S &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] < r);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator < (const S l, const VecSimple<T> &r) 
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l < r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<bool> operator <= (const VecSimple<T> &l, const VecSimple<T> &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] <= r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator <= (const VecSimple<T> &l,const S &r)
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] <= r);
   return temp;
}

template <class T, class S>
static inline VecSimple<bool> operator <= (const S l, const VecSimple<T> &r) 
{
   VecSimple<bool> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l <= r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<T> min(VecSimple<T> const & l, VecSimple<T> const & r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] < r.val[i] ? l.val[i] : r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> min(S const & l, VecSimple<T> const & r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l < r.val[i] ? l : r.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<T> max(VecSimple<T> const & l, VecSimple<T> const & r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] > r.val[i] ? l.val[i] : r.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> max(VecSimple<T> const & l, S const & r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, l.val[i] > r ? l.val[i] : r);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> max(S const & l, VecSimple<T> const & r){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, r.val[i] > l ? r.val[i] : l);
   return temp;
}

template <class T>
static inline VecSimple<T> select(VecSimple<bool> const & a, VecSimple<T> const & b, VecSimple<T> const & c){ 
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, a.val[i] ? b.val[i] : c.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> select(VecSimple<bool> const & a, S const & b, VecSimple<T> const & c){ 
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, a.val[i] ? b : c.val[i]);
   return temp;
}

template <class T, class S>
static inline VecSimple<T> select(VecSimple<bool> const & a, VecSimple<T> const & b, S const & c){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, a.val[i] ? b.val[i] : c);
   return temp;
}

template <class T>
static inline VecSimple<T> select(VecSimple<bool> const & a, T const & b, T const & c){
   VecSimple<T> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, a.val[i] ? b : c);
   return temp;
}

template <class T>
static inline bool horizontal_or(VecSimple<T> const & a){ 
   bool temp = a.val[0];
   for(unsigned int i=1;i<VECL;i++)
      temp = temp || a.val[i];
   return temp;
}

template <class T>
static inline bool horizontal_and(VecSimple<T> const & a){ 
   bool temp = a.val[0];
   for(unsigned int i=1;i<VECL;i++)
      temp = temp && a.val[i];
   return temp;
}

template <class T>
static inline VecSimple<int> truncate_to_int(VecSimple<T> const & a){ 
   VecSimple<int> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, (int)a.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<double> to_double(VecSimple<T> const & a){ 
   VecSimple<double> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, (double)a.val[i]);
   return temp;
}

template <class T>
static inline VecSimple<float> to_float(VecSimple<T> const & a){ 
   VecSimple<float> temp;
   for(unsigned int i=0;i<VECL;i++)
      temp.insert(i, (float)a.val[i]);
   return temp;
}














template <class T>
class Vec4Simple {
public:
   T val[4] __attribute__((aligned(32)));
   // donot initi v
   Vec4Simple() { }
   // Replicate scalar x across v.
   Vec4Simple(T x){
      for(unsigned int i=0;i<4;i++)
         val[i]=x;
   }
   
   // Replicate 4 values across v.   
   Vec4Simple(T a,T b,T c,T d){
      val[0]=a;
      val[1]=b;
      val[2]=c;
      val[3]=d;
      
   }
   // Copy vector v.
   Vec4Simple(Vec4Simple const &x){
      for(unsigned int i=0;i<4;i++)
         val[i]=x.val[i];
   }

   // Member function to load from array (unaligned)
   Vec4Simple & load(T const * p)  {
      for(unsigned int i=0;i<4;i++)
         val[i]=p[i];
      return *this;
   }
   // Member function to load from array, aligned by 32
   Vec4Simple & load_a(T const * p){
      return this->load(p);
   }
   
   Vec4Simple & insert(int i,T const &x) {
      val[i]=x;
      return *this;
   }


// Member function to store into array (unaligned)
   void store(T * p) const {
      for(unsigned int i=0;i<4;i++)
         p[i]=val[i];
   }
   // Member function to store into array, aligned by 32
   void store_a(T * p) const {
      this->store(p);
   }

   Vec4Simple & operator = (Vec4Simple const & r){
      for(unsigned int i=0;i<4;i++)
         val[i]=r.val[i];
      
      return *this;
   }

   T operator [](int i) const{
      return val[i];
   }

   Vec4Simple operator++ (int)
   {
      Vec4Simple<T> temp (*this);
      for(unsigned int i=0;i<4;i++)
         val[i]++;
      return temp;
   }
};

template <class T>
static inline Vec4Simple<T> abs(const Vec4Simple<T> &l){
   return Vec4Simple<T>(
      fabs(l.val[0]),
      fabs(l.val[1]),
      fabs(l.val[2]),
      fabs(l.val[3])
   );
}

template <class T>
static inline Vec4Simple<T> sqrt(const Vec4Simple<T> &l){
   return Vec4Simple<T>(
      sqrt(l.val[0]),
      sqrt(l.val[1]),
      sqrt(l.val[2]),
      sqrt(l.val[3])
   );
}



template <class T>
static inline Vec4Simple<T> operator + (const Vec4Simple<T> &l, const Vec4Simple<T> &r){
   return Vec4Simple<T>(
      l.val[0]+r.val[0],
      l.val[1]+r.val[1],
      l.val[2]+r.val[2],
      l.val[3]+r.val[3]
   );
}

template <class T, class S>
static inline Vec4Simple<T> operator + (const S &l, const Vec4Simple<T> &r){
   return Vec4Simple<T>(
      l+r.val[0],
      l+r.val[1],
      l+r.val[2],
      l+r.val[3]
   );
}

template <class T, class S>
static inline Vec4Simple<T> operator + (const Vec4Simple<T> &l, const S &r){
   return Vec4Simple<T>(
      l.val[0]+r,
      l.val[1]+r,
      l.val[2]+r,
      l.val[3]+r
   );
}
template <class T>
static inline Vec4Simple<T> operator - (const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      -r.val[0],
      -r.val[1],
      -r.val[2],
      -r.val[3]
   );
}




template <class T>
static inline Vec4Simple<T> operator - (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      l.val[0]-r.val[0],
      l.val[1]-r.val[1],
      l.val[2]-r.val[2],
      l.val[3]-r.val[3]
   );
}

template <class T, class S>
static inline Vec4Simple<T> operator - (const S &l, const Vec4Simple<T> &r){
   return Vec4Simple<T>(
      l-r.val[0],
      l-r.val[1],
      l-r.val[2],
      l-r.val[3]
   );
}

template <class T, class S>
static inline Vec4Simple<T> operator - (const Vec4Simple<T> &l, const S &r){
   return Vec4Simple<T>(
      l.val[0]-r,
      l.val[1]-r,
      l.val[2]-r,
      l.val[3]-r
   );
}

template <class T>
static inline Vec4Simple<T> operator * (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      l.val[0]*r.val[0],
      l.val[1]*r.val[1],
      l.val[2]*r.val[2],
      l.val[3]*r.val[3]
   );
}

template <class T, class S>
static inline Vec4Simple<T> operator * (const Vec4Simple<T> &l, const S &r)
{
   return Vec4Simple<T>(
      l.val[0]*r,
      l.val[1]*r,
      l.val[2]*r,
      l.val[3]*r
   );
}

template <class T, class S>
static inline Vec4Simple<T> operator * (const S &l,const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      l*r.val[0],
      l*r.val[1],
      l*r.val[2],
      l*r.val[3]
   );
}



template <class T>
static inline Vec4Simple<T> operator / (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<T>(
      l.val[0]/r.val[0],
      l.val[1]/r.val[1],
      l.val[2]/r.val[2],
      l.val[3]/r.val[3]
   );
}

template <class T, class S>
static inline Vec4Simple<T> operator / (const Vec4Simple<T> &l, const S &r)
{
   return Vec4Simple<T>(
      l.val[0]/r,
      l.val[1]/r,
      l.val[2]/r,
      l.val[3]/r
   );
}

template <class T, class S>
static inline Vec4Simple<T> operator / (const S &l, const Vec4Simple<T> &r )
{
   return Vec4Simple<T>(
      l/r.val[0],
      l/r.val[1],
      l/r.val[2],
      l/r.val[3]
   );
}

template <class T>
static inline  Vec4Simple<T> & operator += (Vec4Simple<T> &l, const Vec4Simple<T> &r){
   l=l+r;
   return l;
}

template <class T, class S>
static inline  Vec4Simple<T> & operator += (Vec4Simple<T> &l, const S &r){
   l = l+r;
   return l;
}

template <class T>
static inline Vec4Simple<T> & operator -= (Vec4Simple<T> &l, const Vec4Simple<T> &r){
   l=l-r;
   return l;
}

template <class T, class S>   
static inline Vec4Simple<T> & operator -= (Vec4Simple<T> &l, const S &r){
   l = l - r;
   return l;
}

template <class T>
static inline Vec4Simple<bool> operator || (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] || r.val[0],
      l.val[1] || r.val[1],
      l.val[2] || r.val[2],
      l.val[3] || r.val[3]
   );
}



template <class T>
static inline Vec4Simple<bool> operator && (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] && r.val[0],
      l.val[1] && r.val[1],
      l.val[2] && r.val[2],
      l.val[3] && r.val[3]
   );
}

template <class T>
static inline Vec4Simple<bool> operator == (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] == r.val[0],
      l.val[1] == r.val[1],
      l.val[2] == r.val[2],
      l.val[3] == r.val[3]
   );
}

template <class T, class S>
static inline Vec4Simple<bool> operator == (const Vec4Simple<T> &l, const S& r)
{
   return Vec4Simple<bool>(
      l.val[0] == r,
      l.val[1] == r,
      l.val[2] == r,
      l.val[3] == r
   );
}

template <class T, class S>
static inline Vec4Simple<bool> operator != (const Vec4Simple<T> &l, const S& r)
{
   return Vec4Simple<bool>(
      l.val[0] != r,
      l.val[1] != r,
      l.val[2] != r,
      l.val[3] != r
   );
}

template <class T>
static inline Vec4Simple<bool> operator ! (const Vec4Simple<T> &l)
{
   return Vec4Simple<bool>(
      !l.val[0],
      !l.val[1],
      !l.val[2],
      !l.val[3]
   );
}


template <class T>
static inline Vec4Simple<bool> operator > (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] > r.val[0],
      l.val[1] > r.val[1],
      l.val[2] > r.val[2],
      l.val[3] > r.val[3]
   );
}


template <class T, class S>
static inline Vec4Simple<bool> operator > (const Vec4Simple<T> &l, const S r)
{
   return Vec4Simple<bool>(
      l.val[0] > r,
      l.val[1] > r,
      l.val[2] > r,
      l.val[3] > r
   );
}



template <class T, class S>
  static inline Vec4Simple<bool> operator > (const S l,const Vec4Simple<T> &r) 
{
   return Vec4Simple<bool>(
      l > r.val[0],
      l > r.val[1],
      l > r.val[2],
      l > r.val[3]
   );
}


template <class T>
static inline Vec4Simple<bool> operator >= (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] >= r.val[0],
      l.val[1] >= r.val[1],
      l.val[2] >= r.val[2],
      l.val[3] >= r.val[3]
   );
}


template <class T, class S>
static inline Vec4Simple<bool> operator >= (const Vec4Simple<T> &l, const S r)
{
   return Vec4Simple<bool>(
      l.val[0] >= r,
      l.val[1] >= r,
      l.val[2] >= r,
      l.val[3] >= r
   );
}



template <class T, class S>
  static inline Vec4Simple<bool> operator >= (const S l,const Vec4Simple<T> &r) 
{
   return Vec4Simple<bool>(
      l >= r.val[0],
      l >= r.val[1],
      l >= r.val[2],
      l >= r.val[3]
   );
}



template <class T>
static inline Vec4Simple<bool> operator < (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] < r.val[0],
      l.val[1] < r.val[1],
      l.val[2] < r.val[2],
      l.val[3] < r.val[3]
   );
}


template <class T, class S>
static inline Vec4Simple<bool> operator < (const Vec4Simple<T> &l,const S &r)
{
   return Vec4Simple<bool>(
      l.val[0] < r,
      l.val[1] < r,
      l.val[2] < r,
      l.val[3] < r
   );
}

template <class T, class S>
  static inline Vec4Simple<bool> operator < (const S l, const Vec4Simple<T> &r) 
{
   return Vec4Simple<bool>(
      l < r.val[0],
      l < r.val[1],
      l < r.val[2],
      l < r.val[3]
   );
}



template <class T>
static inline Vec4Simple<bool> operator <= (const Vec4Simple<T> &l, const Vec4Simple<T> &r)
{
   return Vec4Simple<bool>(
      l.val[0] <= r.val[0],
      l.val[1] <= r.val[1],
      l.val[2] <= r.val[2],
      l.val[3] <= r.val[3]
   );
}


template <class T, class S>
static inline Vec4Simple<bool> operator <= (const Vec4Simple<T> &l,const S &r)
{
   return Vec4Simple<bool>(
      l.val[0] <= r,
      l.val[1] <= r,
      l.val[2] <= r,
      l.val[3] <= r
   );
}

template <class T, class S>
  static inline Vec4Simple<bool> operator <= (const S l, const Vec4Simple<T> &r) 
{
   return Vec4Simple<bool>(
      l <= r.val[0],
      l <= r.val[1],
      l <= r.val[2],
      l <= r.val[3]
   );
}




template <class T>
static inline Vec4Simple<T> min(Vec4Simple<T> const & l, Vec4Simple<T> const & r){
   return Vec4Simple<T>(
      l.val[0]<r.val[0]?l.val[0]:r.val[0],
      l.val[1]<r.val[1]?l.val[1]:r.val[1],
      l.val[2]<r.val[2]?l.val[2]:r.val[2],
      l.val[3]<r.val[3]?l.val[3]:r.val[3]
   );
}

template <class T, class S>
static inline Vec4Simple<T> min(S const & l, Vec4Simple<T> const & r){
   return Vec4Simple<T>(
      l <r.val[0] ? l:r.val[0],
      l <r.val[1] ? l:r.val[1],
      l <r.val[2] ? l:r.val[2],
      l <r.val[3] ? l:r.val[3]
   );
}

template <class T>
static inline Vec4Simple<T> max(Vec4Simple<T> const & l, Vec4Simple<T> const & r){
   return Vec4Simple<T>(
      l.val[0]>r.val[0]?l.val[0]:r.val[0],
      l.val[1]>r.val[1]?l.val[1]:r.val[1],
      l.val[2]>r.val[2]?l.val[2]:r.val[2],
      l.val[3]>r.val[3]?l.val[3]:r.val[3]
   );
}


template <class T, class S>
static inline Vec4Simple<T> max(Vec4Simple<T> const & l, S const & r){
   return Vec4Simple<T>(
      l.val[0] > r ? l.val[0] : r,
      l.val[1] > r ? l.val[1] : r,
      l.val[2] > r ? l.val[2] : r,
      l.val[3] > r ? l.val[3] : r
   );
}


template <class T, class S>
  static inline Vec4Simple<T> max(S const & l, Vec4Simple<T> const & r){
  return Vec4Simple<T>(
     r.val[0] > l ? r.val[0] : l,
     r.val[1] > l ? r.val[1] : l,
     r.val[2] > l ? r.val[2] : l,
     r.val[3] > l ? r.val[3] : l
  );
}



template <class T>
static inline Vec4Simple<T> select(Vec4Simple<bool> const & a, Vec4Simple<T> const & b, Vec4Simple<T> const & c){ 
   return Vec4Simple<T>(
      a.val[0] ? b.val[0] : c.val[0],
      a.val[1] ? b.val[1] : c.val[1],
      a.val[2] ? b.val[2] : c.val[2],
      a.val[3] ? b.val[3] : c.val[3]
   );
}


template <class T, class S>
static inline Vec4Simple<T> select(Vec4Simple<bool> const & a, S const & b, Vec4Simple<T> const & c){ 
   return Vec4Simple<T>(
      a.val[0] ? b : c.val[0],
      a.val[1] ? b : c.val[1],
      a.val[2] ? b : c.val[2],
      a.val[3] ? b : c.val[3]
   );
}


template <class T, class S>
  static inline Vec4Simple<T> select(Vec4Simple<bool> const & a, Vec4Simple<T> const & b, S const & c){
   return Vec4Simple<T>(
      a.val[0] ? b.val[0] : c,
      a.val[1] ? b.val[1] : c,
      a.val[2] ? b.val[2] : c,
      a.val[3] ? b.val[3] : c
   );
}


template <class T>
  static inline Vec4Simple<T> select(Vec4Simple<bool> const & a, T const & b, T const & c){
   return Vec4Simple<T>(
      a.val[0] ? b : c,
      a.val[1] ? b : c,
      a.val[2] ? b : c,
      a.val[3] ? b : c
   );
}

template <class T>
static inline bool horizontal_or(Vec4Simple<T> const & a){ 
  return a.val[0] || a.val[1] || a.val[2] || a.val[3];
}


template <class T>
static inline bool horizontal_and(Vec4Simple<T> const & a){ 
  return a.val[0] && a.val[1] && a.val[2] && a.val[3];
}



template <class T>
static inline Vec4Simple<int> truncate_to_int(Vec4Simple<T> const & a){ 
  return Vec4Simple<int>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

template <class T>
static inline Vec4Simple<double> to_double(Vec4Simple<T> const & a){ 
  return Vec4Simple<double>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

template <class T>
static inline Vec4Simple<float> to_float(Vec4Simple<T> const & a){ 
  return Vec4Simple<float>(a.val[0], a.val[1], a.val[2], a.val[3]);
}

// Dummy functions that Agner vectorclass has.
// These are here to suppress compiler error messages only
static void no_subnormals() 
{
}




template <class T>
class Vec8Simple {
public:
   T val[8] __attribute__((aligned(32)));
   // donot initi v
   Vec8Simple() { }
   // Replicate scalar x across v.
   Vec8Simple(T x){
      for(unsigned int i=0; i<8; i++)
         val[i]=x;
   }
   
   // Replicate 4 values across v.   
   Vec8Simple(T a,T b,T c,T d, T e,T f,T g,T h ){
      val[0]=a;
      val[1]=b;
      val[2]=c;
      val[3]=d;
      val[4]=e;
      val[5]=f;
      val[6]=g;
      val[7]=h;
   }
   // Copy vector v.
   Vec8Simple(Vec8Simple const &x){
      for(unsigned int i=0;i<8;i++)
         val[i]=x.val[i];
   }

   // Member function to load from array (unaligned)
   Vec8Simple & load(T const * p)  {
      for(unsigned int i=0;i<8;i++)
         val[i]=p[i];
      return *this;
   }
   // Member function to load from array, aligned by 32
   Vec8Simple & load_a(T const * p){
      return this->load(p);
   }
   
   Vec8Simple & insert(int i,T const &x) {
      val[i]=x;
      return *this;
   }


// Member function to store into array (unaligned)
   void store(T * p) const {
      for(unsigned int i=0;i<8;i++)
         p[i]=val[i];
   }
   // Member function to store into array, aligned by 32
   void store_a(T * p) const {
      this->store(p);
   }

   Vec8Simple & operator = (Vec8Simple const & r){
      for(unsigned int i=0;i<8;i++)
         val[i]=r.val[i];
      
      return *this;
   }

   T operator [](int i) const{
      return val[i];
   }

   Vec8Simple operator++ (int)
   {
      Vec8Simple<T> temp (*this);
      for(unsigned int i=0;i<8;i++)
         val[i]++;
      return temp;
   }
};

template <class T>
static inline Vec8Simple<T> abs(const Vec8Simple<T> &l){
   return Vec8Simple<T>(
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
static inline Vec8Simple<T> sqrt(const Vec8Simple<T> &l){
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
static inline Vec8Simple<T> operator + (const Vec8Simple<T> &l, const Vec8Simple<T> &r){
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
static inline Vec8Simple<T> operator + (const S &l, const Vec8Simple<T> &r){
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
static inline Vec8Simple<T> operator + (const Vec8Simple<T> &l, const S &r){
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
static inline Vec8Simple<T> operator - (const Vec8Simple<T> &r)
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
static inline Vec8Simple<T> operator - (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<T> operator - (const S &l, const Vec8Simple<T> &r){
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
static inline Vec8Simple<T> operator - (const Vec8Simple<T> &l, const S &r){
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
static inline Vec8Simple<T> operator * (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<T> operator * (const Vec8Simple<T> &l, const S &r)
{
   return Vec8Simple<T>(
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
static inline Vec8Simple<T> operator * (const S &l,const Vec8Simple<T> &r)
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
static inline Vec8Simple<T> operator / (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<T> operator / (const Vec8Simple<T> &l, const S &r)
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
static inline Vec8Simple<T> operator / (const S &l, const Vec8Simple<T> &r)
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
static inline  Vec8Simple<T> & operator += (Vec8Simple<T> &l, const Vec8Simple<T> &r){
   l=l+r;
   return l;
}

template <class T, class S>
static inline  Vec8Simple<T> & operator += (Vec8Simple<T> &l, const S &r){
   l = l+r;
   return l;
}

template <class T>
static inline Vec8Simple<T> & operator -= (Vec8Simple<T> &l, const Vec8Simple<T> &r){
   l=l-r;
   return l;
}

template <class T, class S>   
static inline Vec8Simple<T> & operator -= (Vec8Simple<T> &l, const S &r){
   l = l - r;
   return l;
}

template <class T>
static inline Vec8Simple<bool> operator || (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<bool> operator && (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<bool> operator == (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<bool> operator == (const Vec8Simple<T> &l, const S &r)
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
static inline Vec8Simple<bool> operator != (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<bool> operator != (const Vec8Simple<T> &l, const S &r)
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
static inline Vec8Simple<bool> operator ! (const Vec8Simple<T> &l)
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
static inline Vec8Simple<bool> operator > (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<bool> operator > (const Vec8Simple<T> &l, const S r)
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
  static inline Vec8Simple<bool> operator > (const S l,const Vec8Simple<T> &r) 
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
static inline Vec8Simple<bool> operator >= (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<bool> operator >= (const Vec8Simple<T> &l, const S r)
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
  static inline Vec8Simple<bool> operator >= (const S l,const Vec8Simple<T> &r) 
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
static inline Vec8Simple<bool> operator < (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<bool> operator < (const Vec8Simple<T> &l,const S &r)
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
  static inline Vec8Simple<bool> operator < (const S l, const Vec8Simple<T> &r) 
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
static inline Vec8Simple<bool> operator <= (const Vec8Simple<T> &l, const Vec8Simple<T> &r)
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
static inline Vec8Simple<bool> operator <= (const Vec8Simple<T> &l,const S &r)
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
  static inline Vec8Simple<bool> operator <= (const S l, const Vec8Simple<T> &r) 
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
static inline Vec8Simple<T> min(Vec8Simple<T> const & l, Vec8Simple<T> const & r){
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
static inline Vec8Simple<T> min(S const & l, Vec8Simple<T> const & r){
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
static inline Vec8Simple<T> max(Vec8Simple<T> const & l, Vec8Simple<T> const & r){
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
static inline Vec8Simple<T> max(Vec8Simple<T> const & l, S const & r){
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
  static inline Vec8Simple<T> max(S const & l, Vec8Simple<T> const & r){
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
static inline Vec8Simple<T> select(Vec8Simple<bool> const & a, Vec8Simple<T> const & b, Vec8Simple<T> const & c){ 
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
static inline Vec8Simple<T> select(Vec8Simple<bool> const & a, S const & b, Vec8Simple<T> const & c){ 
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
  static inline Vec8Simple<T> select(Vec8Simple<bool> const & a, Vec8Simple<T> const & b, S const & c){
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
  static inline Vec8Simple<T> select(Vec8Simple<bool> const & a, T const & b, T const & c){
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
static inline bool horizontal_or(Vec8Simple<T> const & a){ 
  return a.val[0] || a.val[1] || a.val[2] || a.val[3] ||
     a.val[4] || a.val[5] || a.val[6] || a.val[7];
}


template <class T>
static inline bool horizontal_and(Vec8Simple<T> const & a){ 
   return a.val[0] && a.val[1] && a.val[2] && a.val[3] && 
      a.val[4] && a.val[5] && a.val[6] && a.val[7];
}



template <class T>
static inline Vec8Simple<int> truncate_to_int(Vec8Simple<T> const & a){ 
   return Vec8Simple<int>(a.val[0], a.val[1], a.val[2], a.val[3],
                          a.val[4], a.val[5], a.val[6], a.val[7]);
   
}

template <class T>
static inline Vec8Simple<double> to_double(Vec8Simple<T> const & a){ 
   return Vec8Simple<double>(a.val[0], a.val[1], a.val[2], a.val[3],
                             a.val[4], a.val[5], a.val[6], a.val[7]);
}

template <class T>
static inline Vec8Simple<float> to_float(Vec8Simple<T> const & a){ 
   return Vec8Simple<float>(a.val[0], a.val[1], a.val[2], a.val[3],
                            a.val[4], a.val[5], a.val[6], a.val[7]);
}


#endif
