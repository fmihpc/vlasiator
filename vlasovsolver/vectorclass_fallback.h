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

#include "../arch/arch_device_api.h"
#include <math.h>
#include <stdio.h>
#include <initializer_list>
// Prefetching does nothing in the fallback vectorclass, if no system implementation
// is available
//#ifndef _mm_prefetch
//#define _mm_prefetch(...)
//#endif

template <typename T>
class VecSimple
{
  public:
    T val[VECL] __attribute__((aligned(32)));
    ARCH_HOSTDEV VecSimple();
    ARCH_HOSTDEV VecSimple(T x);
    /*
    VecSimple(T a,T b,T c,T d);
    VecSimple(T a,T b,T c,T d, T e,T f,T g,T h);
    VecSimple(T a,T b,T c,T d,T e,T f,T g,T h,T i,T j,T k,T l,T m,T n,T o,T p);
    VecSimple(T a,T b,T c,T d,T e,T f,T g,T h,T i,T j,T k,T l,T m,T n,T o,T p,T q,T r,T s,T t,T u,T v,T w,T x,T y,T z,T aa,T bb,T cc,T dd,T ee,T ff);
    */
    ARCH_HOSTDEV VecSimple(VecSimple const &x);
    ARCH_HOSTDEV VecSimple<T> & load(T const * p);
    ARCH_HOSTDEV VecSimple<T> & load_a(T const * p);
    ARCH_HOSTDEV VecSimple<T> & insert(int i,T const &x);
    ARCH_HOSTDEV void store(T * p) const;
    ARCH_HOSTDEV void store_a(T * p) const;
    ARCH_HOSTDEV VecSimple<T> & operator = (VecSimple<T> const & r);
    ARCH_HOSTDEV T operator [](int i) const;
    ARCH_HOSTDEV T & operator [](int i);
    ARCH_HOSTDEV VecSimple<T> operator++ (int);
    ARCH_HOSTDEV VecSimple<T> operator-- (int);
    // Pass vector values as an initializer list instead of a bunch of arguments.
    ARCH_HOSTDEV VecSimple(std::initializer_list<T> list)
    {
       if (list.size() == 1) {
          for(int i=0; i<VECL; ++i) {
             val[i] = *(list.begin());
          }
       } else if (list.size() == VECL) {
          unsigned int i = 0;
          for(auto it = list.begin(); it != list.end(); ++it) {
             val[i] = *it;
             ++i;
          }
       } else {
          printf("Constructing a vector with a number of elements not equal to 1 or VECL = %d \nInitializer_list size = %lu\n", VECL, list.size());
       }
    }

};
static ARCH_HOSTDEV void no_subnormals(){};

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> abs(const VecSimple<T> &l)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, fabs(l.val[i]));
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> sqrt(const VecSimple<T> &l)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, sqrt(l.val[i]));
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> operator + (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i]+r.val[i]);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> operator + (const S l, const VecSimple<T> &r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l+r.val[i]);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> operator + (const VecSimple<T> &l, const S r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i]+r);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> operator - (const VecSimple<T> &r)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, -r.val[i]);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> operator - (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i]-r.val[i]);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> operator - (const S l, const VecSimple<T> &r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l-r.val[i]);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> operator - (const VecSimple<T> &l, const S r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i]-r);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> operator * (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i]*r.val[i]);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> operator * (const VecSimple<T> &l, const S r)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i]*r);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> operator * (const S l,const VecSimple<T> &r)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l*r.val[i]);
  }
  return temp;
}



template <class T>
static ARCH_HOSTDEV inline VecSimple<T> operator / (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i]/r.val[i]);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> operator / (const VecSimple<T> &l, const S r)
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i]/r);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> operator / (const S l, const VecSimple<T> &r )
{
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l/r.val[i]);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline  VecSimple<T> & operator += (VecSimple<T> &l, const VecSimple<T> &r){
  l=l+r;
  return l;
}

template <class T, class S>
static ARCH_HOSTDEV inline  VecSimple<T> & operator += (VecSimple<T> &l, const S r){
  l = l+r;
  return l;
}

template <class T>
static ARCH_HOSTDEV inline  VecSimple<T> & operator *= (VecSimple<T> &l, const VecSimple<T> &r){
  l=l*r;
  return l;
}

template <class T, class S>
static ARCH_HOSTDEV inline  VecSimple<T> & operator *= (VecSimple<T> &l, const S r){
  l = l*r;
  return l;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> & operator -= (VecSimple<T> &l, const VecSimple<T> &r){
   l=l-r;
   return l;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> & operator -= (VecSimple<T> &l, const S r){
   l = l - r;
   return l;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<bool> operator || (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] || r.val[i]);
  }
  return temp;
}


template <class T>
static ARCH_HOSTDEV inline VecSimple<bool> operator && (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] && r.val[i]);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<bool> operator == (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] == r.val[i]);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator == (const VecSimple<T> &l, const S r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] == r);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator != (const VecSimple<T> &l, const S r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] != r);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<bool> operator ! (const VecSimple<T> &l)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, !l.val[i]);
  }
  return temp;
}


template <class T>
static ARCH_HOSTDEV inline VecSimple<bool> operator > (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] > r.val[i]);
  }
  return temp;
}


template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator > (const VecSimple<T> &l, const S r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] > r);
  }
  return temp;
}



template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator > (const S l,const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l > r.val[i]);
  }
  return temp;
}


template <class T>
static ARCH_HOSTDEV inline VecSimple<bool> operator >= (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] >= r.val[i]);
  }
  return temp;
}


template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator >= (const VecSimple<T> &l, const S r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] >= r);
  }
  return temp;
}



template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator >= (const S l,const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l >= r.val[i]);
  }
  return temp;
}



template <class T>
static ARCH_HOSTDEV inline VecSimple<bool> operator < (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] < r.val[i]);
  }
  return temp;
}


template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator < (const VecSimple<T> &l,const S r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] < r);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator < (const S l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l < r.val[i]);
  }
  return temp;
}



template <class T>
static ARCH_HOSTDEV inline VecSimple<bool> operator <= (const VecSimple<T> &l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] <= r.val[i]);
  }
  return temp;
}


template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator <= (const VecSimple<T> &l,const S r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] <= r);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<bool> operator <= (const S l, const VecSimple<T> &r)
{
  VecSimple<bool> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l <= r.val[i]);
  }
  return temp;
}




template <class T>
static ARCH_HOSTDEV inline VecSimple<T> min(VecSimple<T> const & l, VecSimple<T> const & r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] < r.val[i] ? l.val[i] : r.val[i]);
  }
  return temp;
}

template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> min(S const l, VecSimple<T> const & r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l < r.val[i] ? l : r.val[i]);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> max(VecSimple<T> const & l, VecSimple<T> const & r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] > r.val[i] ? l.val[i] : r.val[i]);
  }
  return temp;
}


template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> max(VecSimple<T> const & l, S const r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, l.val[i] > r ? l.val[i] : r);
  }
  return temp;
}


template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> max(S const l, VecSimple<T> const & r){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, r.val[i] > l ? r.val[i] : l);
  }
  return temp;
}



template <class T>
static ARCH_HOSTDEV inline VecSimple<T> select(VecSimple<bool> const & a, VecSimple<T> const & b, VecSimple<T> const & c){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, a.val[i] ? b.val[i] : c.val[i]);
  }
  return temp;
}


template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> select(VecSimple<bool> const & a, S const b, VecSimple<T> const & c){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, a.val[i] ? b : c.val[i]);
  }
  return temp;
}


template <class T, class S>
static ARCH_HOSTDEV inline VecSimple<T> select(VecSimple<bool> const & a, VecSimple<T> const & b, S const c){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, a.val[i] ? b.val[i] : c);
  }
  return temp;
}


template <class T>
static ARCH_HOSTDEV inline VecSimple<T> select(VecSimple<bool> const & a, T const b, T const c){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, a.val[i] ? b : c);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline bool horizontal_or(VecSimple<T> const & a){
  bool temp = a.val[0];
  for(unsigned int i=1;i<VECL;i++) {
     temp = temp || a.val[i];
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline bool horizontal_and(VecSimple<T> const & a){
  bool temp = a.val[0];
  for(unsigned int i=1;i<VECL;i++) {
     temp = temp && a.val[i];
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline T horizontal_add(VecSimple<T> const & a){
  T temp = a.val[0];
  for(unsigned int i=1;i<VECL;i++) {
     temp += a.val[i];
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<int> truncate_to_int(VecSimple<T> const & a){
  VecSimple<int> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, (int)a.val[i]);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<int> roundi(VecSimple<T> const & a){
   // function roundi: round to nearest integer or even.
   //std::fesetround(FE_TONEAREST);
   VecSimple<int> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, std::nearbyint(a.val[i]));
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<T> floor(VecSimple<T> const & a){
  VecSimple<T> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, floor(a.val[i]));
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<double> to_double(VecSimple<T> const & a){
  VecSimple<double> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, (double)a.val[i]);
  }
  return temp;
}

template <class T>
static ARCH_HOSTDEV inline VecSimple<float> to_float(VecSimple<T> const & a){
  VecSimple<float> temp;
  for(unsigned int i=0;i<VECL;i++) {
     temp.insert(i, (float)a.val[i]);
  }
  return temp;
}




template <typename T>
ARCH_HOSTDEV VecSimple<T>::VecSimple() { }

template <typename T>
ARCH_HOSTDEV inline VecSimple<T>::VecSimple(T x)
{
  for(unsigned int i=0;i<VECL;i++) {
    val[i]=x;
  }
}



template <typename T>
ARCH_HOSTDEV inline VecSimple<T>::VecSimple(VecSimple const &x)
{
  for(unsigned int i=0;i<VECL;i++) {
     val[i]=x.val[i];
  }
}
// Member function to load from array (unaligned)
template <typename T>
ARCH_HOSTDEV inline VecSimple<T> & VecSimple<T>::load(T const * p)
{
  for(unsigned int i=0;i<VECL;i++) {
     val[i]=p[i];
  }
  return *this;
}
// Member function to load from array, aligned by 32
template <typename T>
ARCH_HOSTDEV inline VecSimple<T> & VecSimple<T>::load_a(T const * p)
{
   return this->load(p);
}
template <typename T>
ARCH_HOSTDEV inline VecSimple<T> & VecSimple<T>::insert(int i,T const &x)
{
   val[i]=x;
   return *this;
}
// Member function to store into array (unaligned)
template <typename T>
ARCH_HOSTDEV inline void VecSimple<T>::store(T * p) const
{
  for(unsigned int i=0;i<VECL;i++) {
       p[i]=val[i];
  }
}
// Member function to store into array, aligned by 32
template <typename T>
ARCH_HOSTDEV inline void VecSimple<T>::store_a(T * p) const
{
   this->store(p);
}
template <typename T>
ARCH_HOSTDEV inline VecSimple<T> & VecSimple<T>::operator = (VecSimple<T> const & r)
{
   for(unsigned int i=0;i<VECL;i++) {
      val[i]=r.val[i];
   }
   return *this;
}
template <typename T>
ARCH_HOSTDEV inline T VecSimple<T>::operator [](int i) const
{
   return val[i];
}
template <typename T>
ARCH_HOSTDEV inline T & VecSimple<T>::operator [](int i)
{
   return val[i];
}
template <typename T>
ARCH_HOSTDEV inline VecSimple<T> VecSimple<T>::operator++ (int)
{
   for(unsigned int i=0;i<VECL;i++) {
      val[i]++;
   }
   return *this;
}
template <typename T>
ARCH_HOSTDEV inline VecSimple<T> VecSimple<T>::operator-- (int)
{
   for(unsigned int i=0;i<VECL;i++)
      val[i]--;
   return *this;
}

#endif
