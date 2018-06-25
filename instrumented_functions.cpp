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

/*-------------------------------------------------------------------------
  Function call tracing hooks,
  based on etrace by N. Devillard and V. Chudnovsky.

  To activate this feature, compile with the -finstrument-functions flag.
  --------------------------------------------------------------------------*/

#if (__GNUC__>2) || ((__GNUC__ == 2) && (__GNUC_MINOR__ > 95))

/*---------------------------------------------------------------------------
                           Includes
 ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <execinfo.h>

/*---------------------------------------------------------------------------
                           Defines
 ---------------------------------------------------------------------------*/

#define FUNCTION_ENTRY   "enter"
#define FUNCTION_EXIT    "exit"

/*---------------------------------------------------------------------------
                     Function codes
 ---------------------------------------------------------------------------*/

static bool doOutput=false; 

void instrumentation_init(bool doOutputOnThisRank) {
   doOutput = doOutputOnThisRank;
}

extern "C" {
   /** Function called by every function event */
   void __attribute__((__no_instrument_function__)) gnu_ptrace(char * what, void * p) {
      char** fun_name = backtrace_symbols(&p,1);
      fprintf(stderr, "TRACE %ld> %s %s\n", (long)getpid(), what, fun_name[0]);
      free(fun_name);
      return ;
   }

   /** According to gcc documentation: called upon function entry */
   void __attribute__((__no_instrument_function__)) __cyg_profile_func_enter(void *this_fn, void *call_site) {
      if(doOutput) {
         gnu_ptrace(FUNCTION_ENTRY, this_fn);
         (void)call_site;
      }
   }

   /** According to gcc documentation: called upon function exit */
   void __attribute__((__no_instrument_function__)) __cyg_profile_func_exit(void *this_fn, void *call_site) {
      if(doOutput) {
         gnu_ptrace(FUNCTION_EXIT, this_fn);
         (void)call_site;
      }
   }
}

#else
#warning No instrumented function support!
#endif
