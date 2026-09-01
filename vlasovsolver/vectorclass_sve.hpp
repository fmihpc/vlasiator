/*
 * This file is part of Vlasiator.
 * Copyright 2025 SiPearl
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef VECTORCLASS_SVE_FIXED_H
#define VECTORCLASS_SVE_FIXED_H

#include <arm_sve.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <initializer_list>
#include <new>
#include <type_traits>

#define SVE_INLINE __attribute__((always_inline))

#ifndef _mm_prefetch

// Ideally, we would want to use the __pld or __pldw intrinsic,
// but afaik they're available only in armclang...
// This weird switch is to keep compatibility with X86 flags
template <uint32_t hint>
void SVE_INLINE aarch64_do_prefetch(void* addr) {
   switch (hint) {
      case /*_MM_HINT_T0*/ 1:
         __asm__ __volatile__("prfm PLDL1KEEP, [%0]\n" : : "r"(addr) : "memory");
         break;
      case /*_MM_HINT_T1*/ 2:
         __asm__ __volatile__("prfm PLDL2KEEP, [%0]\n" : : "r"(addr) : "memory");
         break;
      case /*_MM_HINT_T2*/ 3:
         __asm__ __volatile__("prfm PLDL3KEEP, [%0]\n" : : "r"(addr) : "memory");
         break;
      case /*_MM_HINT_NTA*/ 0:
         __asm__ __volatile__("prfm PLDL1STRM, [%0]\n" : : "r"(addr) : "memory");
         break;
      case /*_MM_HINT_ENTA*/ 4:
         __asm__ __volatile__("prfm PSTL1KEEP, [%0]\n" : : "r"(addr) : "memory");
         break;
      case /*_MM_HINT_ET0*/ 5:
         __asm__ __volatile__("prfm PSTL1KEEP, [%0]\n" : : "r"(addr) : "memory");
         break;
      case /*_MM_HINT_ET1*/ 6:
         __asm__ __volatile__("prfm PSTL2KEEP, [%0]\n" : : "r"(addr) : "memory");
         break;
      case /*_MM_HINT_ET0*/ 7:
         __asm__ __volatile__("prfm PSTL3KEEP, [%0]\n" : : "r"(addr) : "memory");
         break;
    }
}

#define _mm_prefetch(addr, hint) aarch64_do_prefetch<hint>(addr)
#endif

#if !defined(CACHE_LINE_SIZE)
#define CACHE_LINE_SIZE (64)
#endif

#define SVE_DISPATCH_1ARG(prefix, arg1)                                                                                \
   if constexpr (std::is_same_v<T, __fp16>) {                                                                          \
      return prefix##_f16(arg1);                                                                                       \
   } else if constexpr (std::is_same_v<T, __bf16>) {                                                                   \
      return prefix##_bf16(arg1);                                                                                      \
   } else if constexpr (std::is_same_v<T, float>) {                                                                    \
      return prefix##_f32(arg1);                                                                                       \
   } else if constexpr (std::is_same_v<T, double>) {                                                                   \
      return prefix##_f64(arg1);                                                                                       \
   } else if constexpr (std::is_same_v<T, int32_t>) {                                                                  \
      return prefix##_s32(arg1);                                                                                       \
   } else if constexpr (std::is_same_v<T, int64_t>) {                                                                  \
      return prefix##_s64(arg1);                                                                                       \
   } else {                                                                                                            \
      static_assert(false, "Unsupported type for SVE dispatching");                                                    \
   }

#define SVE_DISPATCH_1ARG_PREDICATED_SFX(prefix, suffix, predicate, arg1)                                              \
   if constexpr (std::is_same_v<T, __fp16>) {                                                                          \
      return prefix##_f16##suffix(predicate, arg1);                                                                    \
   } else if constexpr (std::is_same_v<T, __bf16>) {                                                                   \
      return prefix##_bf16##suffix(predicate, arg1);                                                                   \
   } else if constexpr (std::is_same_v<T, float>) {                                                                    \
      return prefix##_f32##suffix(predicate, arg1);                                                                    \
   } else if constexpr (std::is_same_v<T, double>) {                                                                   \
      return prefix##_f64##suffix(predicate, arg1);                                                                    \
   } else if constexpr (std::is_same_v<T, int32_t>) {                                                                  \
      return prefix##_s32##suffix(predicate, arg1);                                                                    \
   } else if constexpr (std::is_same_v<T, int64_t>) {                                                                  \
      return prefix##_s64##suffix(predicate, arg1);                                                                    \
   } else {                                                                                                            \
      static_assert(false, "Unsupported type for SVE dispatching");                                                    \
   }

#define SVE_DISPATCH_2ARG_PREDICATED_SFX(prefix, suffix, predicate, arg1, arg2)                                        \
   if constexpr (std::is_same_v<T, __fp16>) {                                                                          \
      return prefix##_f16##suffix(predicate, arg1, arg2);                                                              \
   } else if constexpr (std::is_same_v<T, __bf16>) {                                                                   \
      return prefix##_bf16##suffix(predicate, arg1, arg2);                                                             \
   } else if constexpr (std::is_same_v<T, float>) {                                                                    \
      return prefix##_f32##suffix(predicate, arg1, arg2);                                                              \
   } else if constexpr (std::is_same_v<T, double>) {                                                                   \
      return prefix##_f64##suffix(predicate, arg1, arg2);                                                              \
   } else if constexpr (std::is_same_v<T, int32_t>) {                                                                  \
      return prefix##_s32##suffix(predicate, arg1, arg2);                                                              \
   } else if constexpr (std::is_same_v<T, int64_t>) {                                                                  \
      return prefix##_s64##suffix(predicate, arg1, arg2);                                                              \
   } else {                                                                                                            \
      static_assert(false, "Unsupported type for SVE dispatching");                                                    \
   }

#define SVE_DISPATCH_1ARG_PREDICATED_NOSFX(prefix, predicate, arg1)                                                    \
   if constexpr (std::is_same_v<T, __fp16>) {                                                                          \
      return prefix##_f16(predicate, arg1);                                                                            \
   } else if constexpr (std::is_same_v<T, __bf16>) {                                                                   \
      return prefix##_bf16(predicate, arg1);                                                                           \
   } else if constexpr (std::is_same_v<T, float>) {                                                                    \
      return prefix##_f32(predicate, arg1);                                                                            \
   } else if constexpr (std::is_same_v<T, double>) {                                                                   \
      return prefix##_f64(predicate, arg1);                                                                            \
   } else if constexpr (std::is_same_v<T, int32_t>) {                                                                  \
      return prefix##_s32(predicate, arg1);                                                                            \
   } else if constexpr (std::is_same_v<T, int64_t>) {                                                                  \
      return prefix##_s64(predicate, arg1);                                                                            \
   } else {                                                                                                            \
      static_assert(false, "Unsupported type for SVE dispatching");                                                    \
   }

#define SVE_DISPATCH_2ARG_PREDICATED_NOSFX(prefix, predicate, arg1, arg2)                                              \
   if constexpr (std::is_same_v<T, __fp16>) {                                                                          \
      return prefix##_f16(predicate, arg1, arg2);                                                                      \
   } else if constexpr (std::is_same_v<T, __bf16>) {                                                                   \
      return prefix##_bf16(predicate, arg1, arg2);                                                                     \
   } else if constexpr (std::is_same_v<T, float>) {                                                                    \
      return prefix##_f32(predicate, arg1, arg2);                                                                      \
   } else if constexpr (std::is_same_v<T, double>) {                                                                   \
      return prefix##_f64(predicate, arg1, arg2);                                                                      \
   } else if constexpr (std::is_same_v<T, int32_t>) {                                                                  \
      return prefix##_s32(predicate, arg1, arg2);                                                                      \
   } else if constexpr (std::is_same_v<T, int64_t>) {                                                                  \
      return prefix##_s64(predicate, arg1, arg2);                                                                      \
   } else {                                                                                                            \
      static_assert(false, "Unsupported type for SVE dispatching");                                                    \
   }

namespace sve {
template <typename T> struct type {
   using scalar = T;
};

template <> struct type<__fp16> {
   using scalar = __fp16;
   using integral = int16_t;
   using vec = svfloat16_t __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
};

template <> struct type<__bf16> {
   using scalar = __bf16;
   using vec = svbfloat16_t __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
};

template <> struct type<float> {
   using scalar = float;
   using integral = int32_t;
   using vec = svfloat32_t __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
};

template <> struct type<int32_t> {
   using scalar = int32_t;
   using integral = int32_t;
   using vec = svint32_t __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
};

template <> struct type<double> {
   using scalar = double;
   using integral = int64_t;
   using vec = svfloat64_t __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
};

template <> struct type<int64_t> {
   using scalar = int64_t;
   using integral = int64_t;
   using vec = svint64_t __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
};

namespace instructions {

// Yes, I'm aware that there are intrinsics that
// encapsulate this _dispatch_ across multiple types
// Also, there are a ton of intrinsics that do not, so for consistency,
// all the intrinsics that we use have this form.

// Memory Operations
template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE void svst(const svbool_t& predicate, T* location, const sve_type_t& value) {
   // SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svst1, predicate, location, value);
   ::svst1(predicate, location, value);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svld(const svbool_t& predicate, const T* ptr) {
   return ::svld1(predicate, ptr);
}

// Utility Operations
template <typename T, typename sve_type_t = typename sve::type<T>::vec> SVE_INLINE sve_type_t svdup(T value) {
   SVE_DISPATCH_1ARG(svdup_n, value);
}

template <typename T> SVE_INLINE svbool_t svtrue() {
   if constexpr (sizeof(T) == 2)
      return svptrue_b16();
   else if constexpr (sizeof(T) == 4)
      return svptrue_b32();
   else if constexpr (sizeof(T) == 8)
      return svptrue_b64();
   else
      static_assert(false, "Unsupported type for SVE dispatching");
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE T svextract(const svbool_t& predicate, const sve_type_t& vec) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svclastb_n, predicate, T(0), vec);
}

template <typename T> SVE_INLINE svbool_t pred_from_index(const std::size_t index) {
   if constexpr (sizeof(T) == 2) {
      return svwhilele_b16_s32(0, index);
   } else if constexpr (sizeof(T) == 4) {
      return svwhilele_b32_s32(0, index);
   } else if constexpr (sizeof(T) == 8) {
      return svwhilele_b64_s64(0, index);
   } else {
      static_assert(false, "Unsupported type for SVE dispatching");
   }
}

// Logic Operations
template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svorr(const svbool_t& predicate, const sve_type_t& v1, const sve_type_t& v2) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svorr, _x, predicate, v1, v2);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svsqrt(const svbool_t& predicate, const sve_type_t& v1) {
   SVE_DISPATCH_1ARG_PREDICATED_SFX(svsqrt, _x, predicate, v1);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svand(const svbool_t& predicate, const sve_type_t& v1, const sve_type_t& v2) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svand, _x, predicate, v1, v2);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t sveq(const svbool_t& mask, const sve_type_t& v1, const sve_type_t& v2) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpeq, mask, v1, v2);
}

// Math Operations
template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svadd(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svadd, _x, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svadd_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svadd_n, _m, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE T reduce(const svbool_t& predicate, const sve_type_t& vec) {
   SVE_DISPATCH_1ARG_PREDICATED_NOSFX(svaddv, predicate, vec);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svmul(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svmul, _x, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svmul_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svmul_n, _m, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svsub(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svsub, _x, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svsub_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svsub_n, _x, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svneg(const svbool_t& predicate, const sve_type_t& operand) {
   SVE_DISPATCH_1ARG_PREDICATED_SFX(svneg, _x, predicate, operand);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svdiv(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svdiv, _x, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svdiv_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svdiv_n, _x, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svabs(const svbool_t& predicate, const sve_type_t& operand) {
   SVE_DISPATCH_1ARG_PREDICATED_SFX(svabs, _x, predicate, operand);
}

// Comparison Operations
template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t eq(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpeq, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t neq(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpne, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t eq_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpeq_n, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t neq_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpne_n, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t lt(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmplt, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t lt_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmplt_n, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t leq(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmple, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t leq_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmple_n, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t gt(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpgt, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t gt_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpgt_n, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t geq(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpge, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE svbool_t geq_n(const svbool_t& predicate, const sve_type_t& lhs, const T rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svcmpge_n, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svsel(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_NOSFX(svsel, predicate, lhs, rhs);
}

// Min and Max
template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svmin(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svmin, _x, predicate, lhs, rhs);
}

template <typename T, typename sve_type_t = typename sve::type<T>::vec>
SVE_INLINE sve_type_t svmax(const svbool_t& predicate, const sve_type_t& lhs, const sve_type_t& rhs) {
   SVE_DISPATCH_2ARG_PREDICATED_SFX(svmax, _x, predicate, lhs, rhs);
}
} // namespace instructions

namespace functions {

template <typename sve_wrapper_vector_t, typename predicate_t, typename lhs_sve_vector_t, typename rhs_sve_vector_t,
          std::size_t len, typename return_sve_vector_t = typename sve_wrapper_vector_t::sve_vec_t>
sve_wrapper_vector_t transform(const predicate_t& predicate, const lhs_sve_vector_t* lhs, const rhs_sve_vector_t* rhs,
                               return_sve_vector_t (*binary_operation)(const svbool_t&, const lhs_sve_vector_t&,
                                                                       const rhs_sve_vector_t&)) {
   return_sve_vector_t ret[len];

   for (std::size_t i = 0; i < len; ++i) {
      ret[i] = binary_operation(predicate.at(i), lhs[i], rhs[i]);
   }
   return sve_wrapper_vector_t(ret);
}

template <typename sve_wrapper_vector_t, typename lhs_sve_vector_t, typename rhs_sve_vector_t, std::size_t len,
          typename return_sve_vector_t = typename sve_wrapper_vector_t::sve_vec_t>
sve_wrapper_vector_t transform(const lhs_sve_vector_t* lhs, const rhs_sve_vector_t* rhs,
                               return_sve_vector_t (*binary_operation)(const svbool_t&, const lhs_sve_vector_t&,
                                                                       const rhs_sve_vector_t&)) {
   return_sve_vector_t ret[len];
   auto all_lanes = sve::instructions::svtrue<typename sve_wrapper_vector_t::sve_scalar_t>();

   for (std::size_t i = 0; i < len; ++i) {
      ret[i] = binary_operation(all_lanes, lhs[i], rhs[i]);
   }
   return sve_wrapper_vector_t(ret);
}

template <typename sve_wrapper_vector_t, typename lhs_sve_vector_t, std::size_t len,
          typename return_sve_vector_t = typename sve_wrapper_vector_t::sve_vec_t,
          typename scalar_t = typename sve_wrapper_vector_t::sve_scalar_t>
sve_wrapper_vector_t transform(const lhs_sve_vector_t* lhs, const scalar_t rhs,
                               return_sve_vector_t (*binary_operation_with_scalar)(const svbool_t&,
                                                                                   const lhs_sve_vector_t&,
                                                                                   const scalar_t)) {
   return_sve_vector_t ret[len];
   auto all_lanes = sve::instructions::svtrue<typename sve_wrapper_vector_t::sve_scalar_t>();

   for (std::size_t i = 0; i < len; ++i) {
      ret[i] = binary_operation_with_scalar(all_lanes, lhs[i], rhs);
   }
   return sve_wrapper_vector_t(ret);
}

template <typename sve_wrapper_vector_t, typename operand_vector_t, std::size_t len,
          typename return_sve_vector_t = typename sve_wrapper_vector_t::sve_vec_t>
sve_wrapper_vector_t transform(const operand_vector_t* operand,
                               return_sve_vector_t (*unary_operation)(const svbool_t&, const operand_vector_t&)) {
   return_sve_vector_t ret[len];
   auto all_lanes = sve::instructions::svtrue<typename sve_wrapper_vector_t::sve_scalar_t>();

   for (std::size_t i = 0; i < len; ++i) {
      ret[i] = unary_operation(all_lanes, operand[i]);
   }
   return sve_wrapper_vector_t(ret);
}

template <typename sve_wrapper_vector_t, typename predicate_t, typename lhs_sve_vector_t, std::size_t len,
          typename return_sve_vector_t = typename sve_wrapper_vector_t::sve_vec_t>
sve_wrapper_vector_t transform(const lhs_sve_vector_t* operand, const predicate_t* predicate,
                               return_sve_vector_t (*unary_operation)(const svbool_t&, const lhs_sve_vector_t&)) {
   return_sve_vector_t ret[len];

   for (std::size_t i = 0; i < len; ++i) {
      ret[i] = unary_operation(predicate[i], operand[i]);
   }
   return sve_wrapper_vector_t(ret);
}

// Version for boolean ops
template <typename sve_wrapper_vector_t, typename predicate_t, typename lhs_sve_vector_t, typename rhs_sve_vector_t,
          std::size_t len, typename return_sve_vector_t = typename sve_wrapper_vector_t::sve_vec_t,
          typename return_scalar_t = typename sve_wrapper_vector_t::sve_scalar_t>
sve_wrapper_vector_t transform(const lhs_sve_vector_t* lhs, const rhs_sve_vector_t* rhs, const predicate_t* predicate,
                               svbool_t (*binary_operation)(const svbool_t&, const lhs_sve_vector_t&,
                                                            const rhs_sve_vector_t&)) {
   return_sve_vector_t ret[len];
   auto ones = sve::instructions::svdup(return_scalar_t(1));
   auto zeroes = sve::instructions::svdup(return_scalar_t(0));

   for (std::size_t i = 0; i < len; ++i) {
      auto result_mask = binary_operation(predicate[i], lhs[i], rhs[i]);
      ret[i] = sve::instructions::svsel<return_scalar_t>(result_mask, ones, zeroes);
   }
   return sve_wrapper_vector_t(ret);
}

// Unpredicated
template <typename sve_wrapper_vector_t, typename lhs_sve_vector_t, typename rhs_sve_vector_t, std::size_t len,
          typename return_sve_vector_t = typename sve_wrapper_vector_t::sve_vec_t,
          typename return_scalar_t = typename sve_wrapper_vector_t::sve_scalar_t>
sve_wrapper_vector_t transform(const lhs_sve_vector_t* lhs, const rhs_sve_vector_t* rhs,
                               svbool_t (*binary_operation)(const svbool_t&, const lhs_sve_vector_t&,
                                                            const rhs_sve_vector_t&)) {
   return_sve_vector_t ret[len];
   auto ones = sve::instructions::svdup(return_scalar_t(1));
   auto zeroes = sve::instructions::svdup(return_scalar_t(0));
   auto all_lanes = sve::instructions::svtrue<return_scalar_t>();

   for (std::size_t i = 0; i < len; ++i) {
      auto result_mask = binary_operation(all_lanes, lhs[i], rhs[i]);
      ret[i] = sve::instructions::svsel<return_scalar_t>(result_mask, ones, zeroes);
   }
   return sve_wrapper_vector_t(ret);
}

template <typename sve_wrapper_vector_t, typename lhs_sve_vector_t, std::size_t len,
          typename scalar_t = typename sve_wrapper_vector_t::sve_scalar_t,
          typename return_sve_vector_t = typename sve_wrapper_vector_t::sve_vec_t,
          typename return_scalar_t = typename sve_wrapper_vector_t::sve_scalar_t>
sve_wrapper_vector_t transform(const lhs_sve_vector_t* lhs, const scalar_t rhs,
                               svbool_t (*binary_operation_with_scalar)(const svbool_t&, const lhs_sve_vector_t&,
                                                                        const scalar_t)) {
   return_sve_vector_t ret[len];
   auto ones = sve::instructions::svdup(return_scalar_t(1));
   auto zeroes = sve::instructions::svdup(return_scalar_t(0));
   auto all_lanes = sve::instructions::svtrue<return_scalar_t>();

   for (std::size_t i = 0; i < len; ++i) {
      auto result_mask = binary_operation_with_scalar(all_lanes, lhs[i], rhs);
      ret[i] = sve::instructions::svsel<return_scalar_t>(result_mask, ones, zeroes);
   }
   return sve_wrapper_vector_t(ret);
}

template <typename T, std::size_t len, std::size_t stride, typename vec_t = typename sve::type<T>::vec>
void load_into(vec_t* dst, const T* src) {
   auto all_lanes = sve::instructions::svtrue<T>();
   for (std::size_t i = 0; i < len; ++i) {
      dst[i] = sve::instructions::svld(all_lanes, src + i * stride);
   }
}

template <typename T, std::size_t len, std::size_t stride, typename vec_t = typename sve::type<T>::vec>
void store_into(T* dst, const vec_t* src) {
   auto all_lanes = sve::instructions::svtrue<T>();
   for (std::size_t i = 0; i < len; ++i) {
      sve::instructions::svst(all_lanes, dst + i * stride, src[i]);
   }
}
} // namespace functions
} // namespace sve

static void no_subnormals() {};

template <typename T, std::size_t N> class VecSVEBase {

public:
   using sve_vec_t = typename sve::type<T>::vec;
   using sve_integral_t = typename sve::type<T>::integral;
   using sve_scalar_t = typename sve::type<T>::scalar;
   static constexpr std::size_t sve_size_bytes = __ARM_FEATURE_SVE_BITS >> 3;
   static constexpr std::size_t len = N * sizeof(T) / sve_size_bytes;
   static constexpr std::size_t elems_per_vec = sve_size_bytes / sizeof(T);

protected:
   sve_vec_t m_data[len] __attribute__((aligned(CACHE_LINE_SIZE)));

public:
   T operator[](const std::size_t i) const {
      auto idx_in_data = i / elems_per_vec;
      auto idx_in_vec = i % elems_per_vec;
      auto pred = sve::instructions::pred_from_index<T>(idx_in_vec);
      auto ret = sve::instructions::svextract<T>(pred, this->m_data[idx_in_data]);
      return ret;
   }

   SVE_INLINE sve_scalar_t reduce() const {
      auto all_lanes = sve::instructions::svtrue<sve_scalar_t>();
      sve_scalar_t sum = 0;
      for (std::size_t i = 0; i < len; ++i) {
         sum += sve::instructions::reduce<sve_scalar_t>(all_lanes, this->m_data[i]);
      }
      return sum;
   }

   sve_scalar_t* leak_ptr() const {
      auto* ptr = reinterpret_cast<sve_scalar_t*>(std::aligned_alloc(CACHE_LINE_SIZE, sizeof(sve_scalar_t) * N));
      sve::functions::store_into<T, len, elems_per_vec>(ptr, m_data);
      return ptr;
   }
};

/*
 * Mask class, to be used for selecting values from VecSVE.
 * Ideally, we would want to use a lower precision datatype
 * However, _mixed_ operations in SVE are very tricky
 * (they are sparrigly available) and thus it requires conversion.
 * We prefer to use a bit more memory, for better performance
 */
template <typename base_t, std::size_t N>
class VecSVEMask : public VecSVEBase<typename sve::type<base_t>::integral, N> {

   using T = typename sve::type<base_t>::integral;

public:
   using sve_vec_t = typename VecSVEBase<T, N>::sve_vec_t;
   using sve_scalar_t = typename VecSVEBase<T, N>::sve_scalar_t;
   static constexpr std::size_t len = VecSVEBase<T, N>::len;
   static constexpr std::size_t elems_per_vec = VecSVEBase<T, N>::elems_per_vec;

   VecSVEMask() = default;

   SVE_INLINE VecSVEMask(sve_vec_t* data) {
      for (std::size_t i = 0; i < len; ++i) {
         this->m_data[i] = data[i];
      }
   }

   SVE_INLINE VecSVEMask(const T* data) { sve::functions::load_into<T, len, elems_per_vec>(this->m_data, data); }

   SVE_INLINE svbool_t at(std::size_t index) const {
      auto ones = sve::instructions::svdup(sve_scalar_t(1));
      auto all_lanes = sve::instructions::svtrue<sve_scalar_t>();
      auto sve_values = this->m_data[index];
      return sve::instructions::sveq<sve_scalar_t>(all_lanes, sve_values, ones);
   }

   SVE_INLINE bool any() const { return this->reduce() != 0; }

   SVE_INLINE bool all() const { return this->reduce() == N; }

   SVE_INLINE VecSVEMask operator||(const VecSVEMask& rhs) const {
      return sve::functions::transform<VecSVEMask, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                              sve::instructions::svorr<sve_scalar_t>);
   }

   SVE_INLINE VecSVEMask operator&&(const VecSVEMask& rhs) const {
      return sve::functions::transform<VecSVEMask, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                              sve::instructions::svand<sve_scalar_t>);
   }

   SVE_INLINE VecSVEMask operator!() const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      auto all_ones = sve::instructions::svdup(sve_scalar_t(1));
      auto all_zeros = sve::instructions::svdup(sve_scalar_t(0));
      auto all_lanes = sve::instructions::svtrue<sve_scalar_t>();
      sve_vec_t storage[len];

      for (std::size_t i = 0; i < len; ++i)
         storage[i] = sve::instructions::svsel<sve_scalar_t>(at(i), all_zeros, all_ones);

      return VecSVEMask(storage);
   }
};

/*
 * SVE-backed fixed-size vector
 * Exactly what it sounds like.
 */
template <typename T, std::size_t N> class VecSVE : public VecSVEBase<T, N> {

public:
   using sve_vec_t = typename VecSVEBase<T, N>::sve_vec_t;
   using sve_scalar_t = typename VecSVEBase<T, N>::sve_scalar_t;
   static constexpr std::size_t sve_size_bytes = VecSVEBase<T, N>::sve_size_bytes;
   static constexpr std::size_t len = VecSVEBase<T, N>::len;
   static constexpr std::size_t elems_per_vec = VecSVEBase<T, N>::elems_per_vec;

   SVE_INLINE VecSVE() = default;

   SVE_INLINE VecSVE(const T* data) { sve::functions::load_into<T, len, elems_per_vec>(this->m_data, data); }

   SVE_INLINE VecSVE(sve_vec_t* data) {
      for (std::size_t i = 0; i < len; ++i) {
         this->m_data[i] = data[i];
      }
   }

   SVE_INLINE VecSVE(const VecSVE& other) {
      for (std::size_t i = 0; i < len; ++i) {
         this->m_data[i] = other.m_data[i];
      }
   }

   SVE_INLINE VecSVE(std::initializer_list<T> values) {
      // This is a quirk of VectorclassFallback/Agner's Vectorclass
      // to allow initialization from an initializer list
      // with a single value.
      // It is kept here merely as compatibility.
      // This behaviour is very much _nonstandard_, if you want this,
      // use the VecSVE(value) constructor instead.
      // Also, here, a [[likely]] would be great,
      // but I'd rather keep C++17 compatibility
      if (values.size() != 1) {
         const auto* values_ptr = std::data(values);
         sve::functions::load_into<T, len, elems_per_vec>(this->m_data, values_ptr);
      } else {
         auto sve_value = sve::instructions::svdup(*(values.begin()));
         for (std::size_t i = 0; i < len; ++i) {
            this->m_data[i] = sve_value;
         }
      }
   }

   SVE_INLINE VecSVE(const T value) {
      auto sve_value = sve::instructions::svdup(value);
      for (std::size_t i = 0; i < len; ++i) {
         this->m_data[i] = sve_value;
      }
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   explicit VecSVE(U value) : VecSVE(T(value)) {}

   SVE_INLINE VecSVE& operator=(const VecSVE<T, N>& other) {
      for (std::size_t i = 0; i < len; ++i) {
         this->m_data[i] = other.m_data[i];
      }
      return *this;
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVE& operator=(const U value) {
      auto data = sve::instructions::svdup(T(value));
      for (std::size_t i = 0; i < len; ++i) {
         this->m_data[i] = data;
      }
      return *this;
   }

   SVE_INLINE VecSVE& load(const T* ptr) {
      sve::functions::load_into<T, len, elems_per_vec>(this->m_data, ptr);
      return *this;
   }

   SVE_INLINE VecSVE& load_a(const T* ptr) {
      const auto* ptr_aligned = static_cast<const T*>(__builtin_assume_aligned(ptr, CACHE_LINE_SIZE));
      return load(ptr_aligned);
   }

   SVE_INLINE void store(T* destination) const {
      sve::functions::store_into<T, len, elems_per_vec>(destination, this->m_data);
   }

   SVE_INLINE void store_a(T* destination) const {
      auto* destination_aligned = static_cast<T*>(__builtin_assume_aligned(destination, CACHE_LINE_SIZE));
      store(destination_aligned);
   }

   SVE_INLINE VecSVE sqrt() const {
      return sve::functions::transform<VecSVE, sve_vec_t, len>(this->m_data, sve::instructions::svsqrt<T>);
   }

   template <typename conv_t> VecSVE<conv_t, N> as() const {

      static_assert(sizeof(T) == sizeof(conv_t), "Size should match");

      using sve_conv_vec_t = typename VecSVE<conv_t, N>::sve_vec_t;

      sve_conv_vec_t storage[len];
      auto all_lanes = sve::instructions::svtrue<conv_t>();

      for (std::size_t i = 0; i < len; ++i) {
         if constexpr (std::is_same_v<conv_t, float>)
            storage[i] = svcvt_f32_s32_x(all_lanes, this->m_data[i]);
         else if constexpr (std::is_same_v<conv_t, double>)
            storage[i] = svcvt_f64_s64_x(all_lanes, this->m_data[i]);
         else if constexpr (std::is_same_v<conv_t, int16_t>)
            storage[i] = svcvt_s16_f16_x(all_lanes, this->m_data[i]);
         else if constexpr (std::is_same_v<conv_t, int32_t>)
            storage[i] = svcvt_s32_f32_x(all_lanes, this->m_data[i]);
         else if constexpr (std::is_same_v<conv_t, int64_t>)
            storage[i] = svcvt_s64_f64_x(all_lanes, this->m_data[i]);
         else if constexpr (std::is_same_v<T, __bf16>)
            static_assert(false, "Cannot cast from bfloat16 vector to s16 vector.");
         else
            static_assert(false, "Unsupported conversion!");
      }
      return VecSVE<conv_t, N>(storage);
   }

   SVE_INLINE VecSVE operator+(const VecSVE& rhs) const {
      return sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                          sve::instructions::svadd<T>);
   }

   SVE_INLINE VecSVE& operator+=(const VecSVE& rhs) {
      *this = sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                           sve::instructions::svadd<T>);
      return *this;
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVE operator+(const U value) const {
      auto converted_value = T(value);
      return sve::functions::transform<VecSVE, sve_vec_t, len>(this->m_data, converted_value,
                                                               sve::instructions::svadd_n<T>);
   }

   SVE_INLINE VecSVE operator*(const VecSVE& rhs) const {
      return sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                          sve::instructions::svmul<T>);
   }

   SVE_INLINE VecSVE& operator*=(const VecSVE& rhs) {
      *this = sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                           sve::instructions::svmul<T>);
      return *this;
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVE operator*(const U value) const {
      auto converted_value = T(value);
      return sve::functions::transform<VecSVE, sve_vec_t, len>(this->m_data, converted_value,
                                                               sve::instructions::svmul_n<T>);
   }

   SVE_INLINE VecSVE operator-(const VecSVE& rhs) const {
      return sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                          sve::instructions::svsub<T>);
   }

   SVE_INLINE VecSVE& operator-=(const VecSVE& rhs) const {
      *this = sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                           sve::instructions::svsub<T>);
      return *this;
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVE operator-(const U value) const {
      auto converted_value = T(value);
      return sve::functions::transform<VecSVE, sve_vec_t, len>(this->m_data, converted_value,
                                                               sve::instructions::svsub_n<T>);
   }

   SVE_INLINE VecSVE operator-() const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support negating a bf16!");
      return sve::functions::transform<VecSVE, sve_vec_t, len>(this->m_data, sve::instructions::svneg<T>);
   }

   SVE_INLINE VecSVE operator/(const VecSVE& rhs) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support division between bf16!");
      return sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                          sve::instructions::svdiv<T>);
   }

   SVE_INLINE VecSVE& operator/=(const VecSVE& rhs) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support division between bf16!");
      *this = sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                           sve::instructions::svdiv<T>);
      return *this;
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVE operator/(const U value) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support division between bf16!");
      auto converted_value = T(value);
      return sve::functions::transform<VecSVE, sve_vec_t, len>(this->m_data, converted_value,
                                                               sve::instructions::svdiv_n<T>);
   }

   SVE_INLINE VecSVE abs() const {
      return sve::functions::transform<VecSVE, sve_vec_t, len>(this->m_data, sve::instructions::svabs<T>);
   }

   SVE_INLINE VecSVE min(const VecSVE& rhs) const {
      return sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                          sve::instructions::svmin<T>);
   }

   SVE_INLINE VecSVE max(const VecSVE& rhs) const {
      return sve::functions::transform<VecSVE, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                          sve::instructions::svmax<T>);
   }

   SVE_INLINE VecSVE select(const VecSVEMask<T, N>& mask, const VecSVE& rhs) const {
      return sve::functions::transform<VecSVE, VecSVEMask<T, N>, sve_vec_t, sve_vec_t, len>(
          mask, this->m_data, rhs.m_data, sve::instructions::svsel<T>);
   }

   SVE_INLINE VecSVEMask<T, N> operator==(const VecSVE& rhs) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                                    sve::instructions::eq<T>);
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVEMask<T, N> operator==(const U& value) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      auto rhs = T(value);
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, len>(this->m_data, rhs, sve::instructions::eq_n<T>);
   }

   SVE_INLINE VecSVEMask<T, N> operator!=(const VecSVE& rhs) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                                    sve::instructions::neq<T>);
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVEMask<T, N> operator!=(const U& value) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      auto rhs = T(value);
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, len>(this->m_data, rhs,
                                                                         sve::instructions::neq_n<T>);
   }

   SVE_INLINE VecSVEMask<T, N> operator>(const VecSVE& rhs) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                                    sve::instructions::gt<T>);
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVEMask<T, N> operator>(const U& value) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      auto rhs = T(value);
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, len>(this->m_data, rhs, sve::instructions::gt_n<T>);
   }

   SVE_INLINE VecSVEMask<T, N> operator<(const VecSVE& rhs) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                                    sve::instructions::lt<T>);
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVEMask<T, N> operator<(const U& value) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      auto rhs = T(value);
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, len>(this->m_data, rhs, sve::instructions::lt_n<T>);
   }

   SVE_INLINE VecSVEMask<T, N> operator>=(const VecSVE& rhs) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                                    sve::instructions::geq<T>);
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVEMask<T, N> operator>=(const U& value) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      auto rhs = T(value);
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, len>(this->m_data, rhs,
                                                                         sve::instructions::geq_n<T>);
   }

   SVE_INLINE VecSVEMask<T, N> operator<=(const VecSVE& rhs) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, sve_vec_t, len>(this->m_data, rhs.m_data,
                                                                                    sve::instructions::leq<T>);
   }

   template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
   SVE_INLINE VecSVEMask<T, N> operator<=(const U& value) const {
      static_assert(!std::is_same_v<T, __bf16>, "Did you know? SVE doesn't support comparison between bf16!");
      auto rhs = T(value);
      return sve::functions::transform<VecSVEMask<T, N>, sve_vec_t, len>(this->m_data, rhs,
                                                                         sve::instructions::leq_n<T>);
   }
};

// Support functions
template <typename T, std::size_t N> SVE_INLINE VecSVE<T, N> abs(const VecSVE<T, N>& v) { return v.abs(); }

template <typename T, std::size_t N> SVE_INLINE VecSVE<T, N> sqrt(const VecSVE<T, N>& v) { return v.sqrt(); }

template <typename T, std::size_t N>
SVE_INLINE VecSVE<T, N> select(const VecSVEMask<T, N>& mask, const VecSVE<T, N>& v1, const VecSVE<T, N>& v2) {
   return v1.select(mask, v2);
}

template <typename T, std::size_t N, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
SVE_INLINE VecSVE<T, N> select(const VecSVEMask<T, N>& mask, const VecSVE<T, N>& v1, const U value) {
   // svsel doesn't have a _n option to select against a value
   // so we gotta create a new vector instead...
   auto v2 = VecSVE<T, N>(value);
   return v1.select(mask, v2);
}

template <typename T, std::size_t N, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
SVE_INLINE VecSVE<T, N> select(const VecSVEMask<T, N>& mask, const U value, const VecSVE<T, N>& v1) {
   // svsel doesn't have a _n option to select against a value
   // so we gotta create a new vector instead...
   auto v2 = VecSVE<T, N>(value);
   return v1.select(mask, v2);
}

template <typename T, std::size_t N> SVE_INLINE VecSVE<T, N> min(const VecSVE<T, N>& lhs, const VecSVE<T, N>& rhs) {
   return lhs.min(rhs);
}

template <typename T, std::size_t N> SVE_INLINE VecSVE<T, N> max(const VecSVE<T, N>& lhs, const VecSVE<T, N>& rhs) {
   return lhs.max(rhs);
}

template <typename T, std::size_t N> SVE_INLINE bool horizontal_or(const VecSVEMask<T, N>& vec) { return vec.any(); }

template <typename T, std::size_t N> SVE_INLINE bool horizontal_and(const VecSVEMask<T, N>& vec) { return vec.all(); }

template <typename other_t, std::size_t N> VecSVE<float, N> to_float(const VecSVE<other_t, N>& vec) {
   return vec.template as<float>();
}

template <typename other_t, std::size_t N> VecSVE<double, N> to_double(const VecSVE<other_t, N>& vec) {
   return vec.template as<double>();
}

template <typename real_t, std::size_t N, typename integral_t = typename sve::type<real_t>::integral>
VecSVE<integral_t, N> truncate_to_int(const VecSVE<real_t, N>& arg) {
   return arg.template as<integral_t>();
}

template <typename real_t, std::size_t N, typename integral_t = typename sve::type<real_t>::integral>
VecSVE<integral_t, N> roundi(const VecSVE<real_t, N>& arg) {
   auto res = arg + real_t(0.5);
   return arg.template as<integral_t>();
}

template <typename T, std::size_t N> VecSVE<T, N> floor(const VecSVE<T, N>& arg) { return arg; }

template <typename T, std::size_t N> T horizontal_add(const VecSVE<T, N>& vec) { return vec.reduce(); }

// Symmetric operators, because thanks C++ /s
template <typename T, std::size_t N, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
VecSVE<T, N> operator+(const U lhs, const VecSVE<T, N>& rhs) {
   return rhs + lhs;
}

template <typename T, std::size_t N, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
VecSVE<T, N> operator*(const U lhs, const VecSVE<T, N>& rhs) {
   return rhs * lhs;
}

template <typename T, std::size_t N, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
VecSVE<T, N> operator-(U lhs, const VecSVE<T, N>& rhs) {
   return rhs - lhs;
}

template <typename T, std::size_t N, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
VecSVE<T, N> operator/(U lhs, const VecSVE<T, N>& rhs) {
   return rhs / lhs;
}

#endif // VECTORCLASS_SVE_FIXED_H
