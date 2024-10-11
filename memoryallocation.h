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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef MEMORYALLOCATION_H
#define MEMORYALLOCATION_H

#include <cstdlib>
#include <cstdint>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string.h>

#ifdef USE_JEMALLOC
#include "jemalloc/jemalloc.h"
#endif

#ifndef NDEBUG
#ifndef INITIALIZE_ALIGNED_MALLOC_WITH_NAN
#define INITIALIZE_ALIGNED_MALLOC_WITH_NAN
#endif
#endif

#ifdef DEBUG_SPATIAL_CELL
#ifndef INITIALIZE_ALIGNED_MALLOC_WITH_NAN
#define INITIALIZE_ALIGNED_MALLOC_WITH_NAN
#endif
#endif

/*! Return the amount of free memory on the node in bytes*/
uint64_t get_node_free_memory();

/*! Purge jemalloc arenas to actually return memory to the system. If not using jemalloc, does nothing.*/
void memory_purge();

/*! Initialize memory allocator configuration.*/
void memory_configurator();

/*! Measures memory consumption and writes it into logfile. Collective
 *  operation on MPI_COMM_WORLD
 *  extra_bytes is used for additional buffer for the high water mark,
 *  for example when estimating refinement memory usage
 */
void report_process_memory_consumption(double extra_bytes = 0.0);

/*! Alligned malloc, could be done using aligned_alloc*/
inline void * aligned_malloc(size_t size,std::size_t align) {
   /* Allocate necessary memory area
    * client request - size parameter -
    * plus area to store the address
    * of the memory returned by standard
    * malloc(), or je_malloc().
    */
   void *ptr;
#ifdef USE_JEMALLOC
   void *p = je_malloc(size + align - 1 + sizeof(void*));
#else
   void *p = malloc(size + align - 1 + sizeof(void*));
#endif
#ifdef INITIALIZE_ALIGNED_MALLOC_WITH_NAN
   memset(p, ~0u, size + align - 1 + sizeof(void*));
#endif

   if (p != NULL) {
      /* Address of the aligned memory according to the align parameter*/
      ptr = (void*) (((unsigned long)p + sizeof(void*) + align -1) & ~(align-1));
      /* store the address of the malloc() above
       * at the beginning of our total memory area.
       * You can also use *((void **)ptr-1) = p
       * instead of the one below.
       */
      *((void**)((unsigned long)ptr - sizeof(void*))) = p;
      /* Return the address of aligned memory */
      return ptr;
   }
   return NULL;
}

inline void aligned_free(void *p) {
   /* Get the address of the memory, stored at the
    * start of our total memory area. Alternatively,
    * you can use void *ptr = *((void **)p-1) instead
    * of the one below.
    */
   void *ptr = *((void**)((unsigned long)p - sizeof(void*)));
#ifdef USE_JEMALLOC
   je_free(ptr);
#else
   free(ptr);
#endif
   return;
}





/**
 * Allocator for aligned data.
 *
 * Copied from https://gist.github.com/1471329
 * Modified from the Mallocator from Stephan T. Lavavej.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 */

template <typename T, std::size_t Alignment>
class aligned_allocator
{
public:

   // The following will be the same for virtually all allocators.
   typedef T * pointer;
   typedef const T * const_pointer;
   typedef T& reference;
   typedef const T& const_reference;
   typedef T value_type;
   typedef std::size_t size_type;
   typedef ptrdiff_t difference_type;

   T * address(T& r) const
      {
         return &r;
      }

   const T * address(const T& s) const
      {
         return &s;
      }

   std::size_t max_size() const
      {
         // The following has been carefully written to be independent of
         // the definition of size_t and to avoid signed/unsigned warnings.
         return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
      }


   // The following must be the same for all allocators.
   template <typename U>
   struct rebind
   {
      typedef aligned_allocator<U, Alignment> other;
   } ;

   bool operator!=(const aligned_allocator& other) const
      {
         return !(*this == other);
      }

   void construct(T * const p, const T& t) const
      {
         void * const pv = static_cast<void *>(p);

         new (pv) T(t);
      }

   void destroy(T * const p) const
      {
         p->~T();
      }

   // Returns true if and only if storage allocated from *this
   // can be deallocated from other, and vice versa.
   // Always returns true for stateless allocators.
   bool operator==(const aligned_allocator& other) const
      {
         return true;
      }


   // Default constructor, copy constructor, rebinding constructor, and destructor.
   // Empty for stateless allocators.
   aligned_allocator() { }

   aligned_allocator(const aligned_allocator&) { }

   template <typename U> aligned_allocator(const aligned_allocator<U, Alignment>&) { }

   ~aligned_allocator() { }


   // The following will be different for each allocator.
   T * allocate(const std::size_t n) const
      {
         // The return value of allocate(0) is unspecified.
         // Mallocator returns NULL in order to avoid depending
         // on malloc(0)'s implementation-defined behavior
         // (the implementation can define malloc(0) to return NULL,
         // in which case the bad_alloc check below would fire).
         // All allocators can return NULL in this case.
         if (n == 0) {
            return NULL;
         }

         // All allocators should contain an integer overflow check.
         // The Standardization Committee recommends that std::length_error
         // be thrown in the case of integer overflow.
         if (n > max_size())
         {
            throw std::length_error("aligned_allocator<T>::allocate() - Integer overflow.");
         }

         // Mallocator wraps malloc().
         void * const pv = aligned_malloc(n * sizeof(T), Alignment);

         // Allocators should throw std::bad_alloc in the case of memory allocation failure.
         if (pv == NULL)
         {
            throw std::bad_alloc();
         }

         return static_cast<T *>(pv);
      }

   void deallocate(T * const p, const std::size_t ) const
      {
         aligned_free(p);
      }


   // The following will be the same for all allocators that ignore hints.
   template <typename U>
   T * allocate(const std::size_t n, const U * /* const hint */) const
      {
         return allocate(n);
      }


   // Allocators are not required to be assignable, so
   // all allocators should have a private unimplemented
   // assignment operator. Note that this will trigger the
   // off-by-default (enabled under /Wall) warning C4626
   // "assignment operator could not be generated because a
   // base class assignment operator is inaccessible" within
   // the STL headers, but that warning is useless.
private:
   aligned_allocator& operator=(const aligned_allocator&);
};


#endif
