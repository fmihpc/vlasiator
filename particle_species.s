	.text
	.file	"particle_species.ll"
	.file	1 "/home/talgat/vlasiator/particle_species.cpp"
	.globl	_ZN7species7SpeciesC1Ev # -- Begin function _ZN7species7SpeciesC1Ev
	.p2align	4, 0x90
	.type	_ZN7species7SpeciesC1Ev,@function
_ZN7species7SpeciesC1Ev:                # @_ZN7species7SpeciesC1Ev
.Lfunc_begin0:
	.loc	1 29 0                  # particle_species.cpp:29:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: Species: <- $rdi
	#DEBUG_VALUE: Species: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp0:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	1 434 1 prologue_end    # particle_species.cpp:434:1
	leaq	16(%rdi), %rax
	movq	%rax, (%rdi)
	movq	$0, 8(%rdi)
	movb	$0, 16(%rdi)
	.loc	1 29 1                  # particle_species.cpp:29:1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp1:
.Lfunc_end0:
	.size	_ZN7species7SpeciesC1Ev, .Lfunc_end0-_ZN7species7SpeciesC1Ev
	.cfi_endproc
                                        # -- End function
	.globl	_ZN7species7SpeciesC2Ev # -- Begin function _ZN7species7SpeciesC2Ev
	.p2align	4, 0x90
	.type	_ZN7species7SpeciesC2Ev,@function
_ZN7species7SpeciesC2Ev:                # @_ZN7species7SpeciesC2Ev
.Lfunc_begin1:
	.loc	1 0 0                   # particle_species.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp2:
	.loc	1 29 1 prologue_end     # particle_species.cpp:29:1
	leaq	16(%rdi), %rax
	movq	%rax, (%rdi)
	movq	$0, 8(%rdi)
	movb	$0, 16(%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp3:
.Lfunc_end1:
	.size	_ZN7species7SpeciesC2Ev, .Lfunc_end1-_ZN7species7SpeciesC2Ev
	.cfi_endproc
                                        # -- End function
	.globl	_ZN7species7SpeciesC1ERKS0_ # -- Begin function _ZN7species7SpeciesC1ERKS0_
	.p2align	4, 0x90
	.type	_ZN7species7SpeciesC1ERKS0_,@function
_ZN7species7SpeciesC1ERKS0_:            # @_ZN7species7SpeciesC1ERKS0_
.Lfunc_begin2:
	.loc	1 31 0                  # particle_species.cpp:31:0
	.cfi_startproc
	.cfi_personality 3, __gxx_personality_v0
	.cfi_lsda 3, .Lexception0
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: Species: <- $rdi
	#DEBUG_VALUE: Species:other <- $rsi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%rbx
	pushq	%rax
	.cfi_offset %rbx, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
.Ltmp7:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: other <- $rsi
	movq	%rsi, %r14
	movq	%rdi, %rbx
.Ltmp8:
	#DEBUG_VALUE: other <- $r14
	#DEBUG_VALUE: Species:other <- $r14
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: Species: <- $rbx
	.loc	1 434 1 prologue_end    # particle_species.cpp:434:1
	leaq	16(%rdi), %r15
	movq	%r15, (%rdi)
	movq	$0, 8(%rdi)
	movb	$0, 16(%rdi)
.Ltmp4:
.Ltmp9:
	#DEBUG_VALUE: __str <- $r14
	#DEBUG_VALUE: operator=:__str <- $r14
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: operator=: <- $rbx
	.file	2 "/usr/include/c++/9/bits/basic_string.h"
	.loc	2 1366 1                # /usr/include/c++/9/bits/basic_string.h:1366:1
	callq	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_assignERKS4_
.Ltmp5:
.Ltmp10:
# %bb.1:                                # %L.B0220
	#DEBUG_VALUE: Species: <- $rbx
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: Species:other <- $r14
	#DEBUG_VALUE: other <- $r14
	.loc	1 33 1                  # particle_species.cpp:33:1
	movq	32(%r14), %rax
	movq	%rax, 32(%rbx)
	.loc	1 34 1                  # particle_species.cpp:34:1
	movq	40(%r14), %rax
	movq	%rax, 40(%rbx)
	.loc	1 35 1                  # particle_species.cpp:35:1
	movq	48(%r14), %rax
	movq	%rax, 48(%rbx)
	.loc	1 36 1                  # particle_species.cpp:36:1
	movq	56(%r14), %rax
	movq	%rax, 56(%rbx)
	.loc	1 37 1                  # particle_species.cpp:37:1
	addq	$8, %rsp
	popq	%rbx
.Ltmp11:
	popq	%r14
.Ltmp12:
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB2_2:                                # %L_T19107768_7575
	.cfi_def_cfa %rbp, 16
.Ltmp13:
	#DEBUG_VALUE: Species: <- $rbx
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: Species:other <- $r14
	#DEBUG_VALUE: other <- $r14
.Ltmp6:
	.loc	1 0 1 is_stmt 0         # particle_species.cpp:0:1
	movq	%rax, %r14
.Ltmp14:
	.loc	1 658 1 is_stmt 1       # particle_species.cpp:658:1
	movq	(%rbx), %rdi
	cmpq	%r15, %rdi
	je	.LBB2_4
.Ltmp15:
# %bb.3:                                # %L.B0221
	#DEBUG_VALUE: Species: <- $rbx
	#DEBUG_VALUE:  <- $rbx
	callq	_ZdlPv
.Ltmp16:
.LBB2_4:                                # %L..inline.10994
	#DEBUG_VALUE: Species: <- $rbx
	#DEBUG_VALUE:  <- $rbx
	.loc	1 0 1 is_stmt 0         # particle_species.cpp:0:1
	movq	%r14, %rdi
	callq	_Unwind_Resume
.Lfunc_end2:
	.size	_ZN7species7SpeciesC1ERKS0_, .Lfunc_end2-_ZN7species7SpeciesC1ERKS0_
	.cfi_endproc
	.section	.gcc_except_table,"a",@progbits
	.p2align	2
GCC_except_table2:
.Lexception0:
	.byte	255                     # @LPStart Encoding = omit
	.byte	255                     # @TType Encoding = omit
	.byte	1                       # Call site Encoding = uleb128
	.uleb128 .Lcst_end0-.Lcst_begin0
.Lcst_begin0:
	.uleb128 .Ltmp4-.Lfunc_begin2   # >> Call Site 1 <<
	.uleb128 .Ltmp5-.Ltmp4          #   Call between .Ltmp4 and .Ltmp5
	.uleb128 .Ltmp6-.Lfunc_begin2   #     jumps to .Ltmp6
	.byte	0                       #   On action: cleanup
	.uleb128 .Ltmp5-.Lfunc_begin2   # >> Call Site 2 <<
	.uleb128 .Lfunc_end2-.Ltmp5     #   Call between .Ltmp5 and .Lfunc_end2
	.byte	0                       #     has no landing pad
	.byte	0                       #   On action: cleanup
.Lcst_end0:
	.p2align	2
                                        # -- End function
	.text
	.globl	_ZN7species7SpeciesC2ERKS0_ # -- Begin function _ZN7species7SpeciesC2ERKS0_
	.p2align	4, 0x90
	.type	_ZN7species7SpeciesC2ERKS0_,@function
_ZN7species7SpeciesC2ERKS0_:            # @_ZN7species7SpeciesC2ERKS0_
.Lfunc_begin3:
	.loc	1 0 0 is_stmt 1         # particle_species.cpp:0:0
	.cfi_startproc
	.cfi_personality 3, __gxx_personality_v0
	.cfi_lsda 3, .Lexception1
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: Species: <- $rdi
	#DEBUG_VALUE: Species:other <- $rsi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%rbx
	pushq	%rax
	.cfi_offset %rbx, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
.Ltmp20:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: other <- $rsi
	movq	%rsi, %r14
	movq	%rdi, %rbx
.Ltmp21:
	#DEBUG_VALUE: other <- $r14
	#DEBUG_VALUE: Species:other <- $r14
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: Species: <- $rbx
	.loc	1 434 1 prologue_end    # particle_species.cpp:434:1
	leaq	16(%rdi), %r15
	movq	%r15, (%rdi)
	movq	$0, 8(%rdi)
	movb	$0, 16(%rdi)
.Ltmp17:
.Ltmp22:
	#DEBUG_VALUE: __str <- $r14
	#DEBUG_VALUE: operator=:__str <- $r14
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: operator=: <- $rbx
	.loc	2 1366 1                # /usr/include/c++/9/bits/basic_string.h:1366:1
	callq	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_assignERKS4_
.Ltmp18:
.Ltmp23:
# %bb.1:                                # %_ZN7species7SpeciesC1ERKS0_.exit
	#DEBUG_VALUE: Species: <- $rbx
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: Species:other <- $r14
	#DEBUG_VALUE: other <- $r14
	.loc	1 33 1                  # particle_species.cpp:33:1
	movq	32(%r14), %rax
	movq	%rax, 32(%rbx)
	.loc	1 34 1                  # particle_species.cpp:34:1
	movq	40(%r14), %rax
	movq	%rax, 40(%rbx)
	.loc	1 35 1                  # particle_species.cpp:35:1
	movq	48(%r14), %rax
	movq	%rax, 48(%rbx)
	.loc	1 36 1                  # particle_species.cpp:36:1
	movq	56(%r14), %rax
	movq	%rax, 56(%rbx)
.Ltmp24:
	.loc	1 36 1 is_stmt 0        # particle_species.cpp:36:1
	addq	$8, %rsp
	popq	%rbx
.Ltmp25:
	popq	%r14
.Ltmp26:
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB3_2:                                # %L_T19107768_7575.i
	.cfi_def_cfa %rbp, 16
.Ltmp27:
	#DEBUG_VALUE: Species: <- $rbx
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: Species:other <- $r14
	#DEBUG_VALUE: other <- $r14
.Ltmp19:
	.loc	1 0 1                   # particle_species.cpp:0:1
	movq	%rax, %r14
.Ltmp28:
	.loc	1 658 1 is_stmt 1       # particle_species.cpp:658:1
	movq	(%rbx), %rdi
	cmpq	%r15, %rdi
	je	.LBB3_4
.Ltmp29:
# %bb.3:                                # %L.B0221.i
	#DEBUG_VALUE: Species: <- $rbx
	#DEBUG_VALUE:  <- $rbx
	callq	_ZdlPv
.Ltmp30:
.LBB3_4:                                # %L..inline.10994.i
	#DEBUG_VALUE: Species: <- $rbx
	#DEBUG_VALUE:  <- $rbx
	.loc	1 0 1 is_stmt 0         # particle_species.cpp:0:1
	movq	%r14, %rdi
	callq	_Unwind_Resume
.Lfunc_end3:
	.size	_ZN7species7SpeciesC2ERKS0_, .Lfunc_end3-_ZN7species7SpeciesC2ERKS0_
	.cfi_endproc
	.section	.gcc_except_table,"a",@progbits
	.p2align	2
GCC_except_table3:
.Lexception1:
	.byte	255                     # @LPStart Encoding = omit
	.byte	255                     # @TType Encoding = omit
	.byte	1                       # Call site Encoding = uleb128
	.uleb128 .Lcst_end1-.Lcst_begin1
.Lcst_begin1:
	.uleb128 .Ltmp17-.Lfunc_begin3  # >> Call Site 1 <<
	.uleb128 .Ltmp18-.Ltmp17        #   Call between .Ltmp17 and .Ltmp18
	.uleb128 .Ltmp19-.Lfunc_begin3  #     jumps to .Ltmp19
	.byte	0                       #   On action: cleanup
	.uleb128 .Ltmp18-.Lfunc_begin3  # >> Call Site 2 <<
	.uleb128 .Lfunc_end3-.Ltmp18    #   Call between .Ltmp18 and .Lfunc_end3
	.byte	0                       #     has no landing pad
	.byte	0                       #   On action: cleanup
.Lcst_end1:
	.p2align	2
                                        # -- End function
	.text
	.globl	_ZN7species7SpeciesD1Ev # -- Begin function _ZN7species7SpeciesD1Ev
	.p2align	4, 0x90
	.type	_ZN7species7SpeciesD1Ev,@function
_ZN7species7SpeciesD1Ev:                # @_ZN7species7SpeciesD1Ev
.Lfunc_begin4:
	.loc	1 39 0 is_stmt 1        # particle_species.cpp:39:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~Species: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp31:
	#DEBUG_VALUE:  <- $rdi
	movq	%rdi, %rax
.Ltmp32:
	#DEBUG_VALUE:  <- $rax
	#DEBUG_VALUE: ~Species: <- $rax
	.loc	1 658 1 prologue_end    # particle_species.cpp:658:1
	movq	(%rdi), %rdi
	addq	$16, %rax
.Ltmp33:
	cmpq	%rax, %rdi
	je	.LBB4_1
# %bb.2:                                # %L.B0224
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPv                  # TAILCALL
.LBB4_1:                                # %L..inline.11140
	.cfi_def_cfa %rbp, 16
	.loc	1 39 1                  # particle_species.cpp:39:1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp34:
.Lfunc_end4:
	.size	_ZN7species7SpeciesD1Ev, .Lfunc_end4-_ZN7species7SpeciesD1Ev
	.cfi_endproc
                                        # -- End function
	.globl	_ZN7species7SpeciesD2Ev # -- Begin function _ZN7species7SpeciesD2Ev
	.p2align	4, 0x90
	.type	_ZN7species7SpeciesD2Ev,@function
_ZN7species7SpeciesD2Ev:                # @_ZN7species7SpeciesD2Ev
.Lfunc_begin5:
	.loc	1 0 0                   # particle_species.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	%rdi, %rax
.Ltmp35:
	.loc	1 39 1 prologue_end     # particle_species.cpp:39:1
	movq	(%rdi), %rdi
	addq	$16, %rax
	cmpq	%rax, %rdi
	je	.LBB5_1
# %bb.2:                                # %L.B0227
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPv                  # TAILCALL
.LBB5_1:                                # %L..inline.11252
	.cfi_def_cfa %rbp, 16
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp36:
.Lfunc_end5:
	.size	_ZN7species7SpeciesD2Ev, .Lfunc_end5-_ZN7species7SpeciesD2Ev
	.cfi_endproc
                                        # -- End function
	.globl	__sti___20_particle_species_cpp_29edf090 # -- Begin function __sti___20_particle_species_cpp_29edf090
	.p2align	4, 0x90
	.type	__sti___20_particle_species_cpp_29edf090,@function
__sti___20_particle_species_cpp_29edf090: # @__sti___20_particle_species_cpp_29edf090
.Lfunc_begin6:
	.loc	1 0 0                   # particle_species.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp37:
	.loc	1 74 1 prologue_end     # particle_species.cpp:74:1
	cmpl	$1, __I___20_particle_species_cpp_29edf090(%rip)
	jne	.LBB6_2
# %bb.1:                                # %L.B0012
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB6_2:                                # %L.B0234
	.cfi_def_cfa %rbp, 16
	movl	$1, __I___20_particle_species_cpp_29edf090(%rip)
	movl	$_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE, %edi
	callq	_ZNSt8ios_base4InitC1Ev
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	movl	$_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE, %esi
	movl	$__dso_handle, %edx
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	__cxa_atexit            # TAILCALL
.Ltmp38:
.Lfunc_end6:
	.size	__sti___20_particle_species_cpp_29edf090, .Lfunc_end6-__sti___20_particle_species_cpp_29edf090
	.cfi_endproc
                                        # -- End function
	.type	__I___20_particle_species_cpp_29edf090,@object # @__I___20_particle_species_cpp_29edf090
	.bss
	.globl	__I___20_particle_species_cpp_29edf090
	.p2align	2
__I___20_particle_species_cpp_29edf090:
	.long	0                       # 0x0
	.size	__I___20_particle_species_cpp_29edf090, 4

	.type	_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE,@object # @_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE
	.local	_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE
	.comm	_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE,1,1
	.section	.init_array,"aw",@init_array
	.p2align	3
	.quad	__sti___20_particle_species_cpp_29edf090
	.section	.debug_str,"MS",@progbits,1
.Linfo_string0:
	.asciz	" NVC++ 21.2-0"         # string offset=0
.Linfo_string1:
	.asciz	"particle_species.cpp"  # string offset=14
.Linfo_string2:
	.asciz	"/home/talgat/vlasiator" # string offset=35
.Linfo_string3:
	.asciz	"__I___20_particle_species_cpp_29edf090" # string offset=58
.Linfo_string4:
	.asciz	"int"                   # string offset=97
.Linfo_string5:
	.asciz	"_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE" # string offset=101
.Linfo_string6:
	.asciz	"signed char"           # string offset=161
.Linfo_string7:
	.asciz	"__ARRAY_SIZE_TYPE__"   # string offset=173
.Linfo_string8:
	.asciz	"_ZNSt8ios_base4InitE"  # string offset=193
.Linfo_string9:
	.asciz	"__dso_handle"          # string offset=214
.Linfo_string10:
	.asciz	"void"                  # string offset=227
.Linfo_string11:
	.asciz	"_ZN42_INTERNAL_20_particle_species_cpp_29edf0905vmesh15INVALID_LOCALIDE" # string offset=232
.Linfo_string12:
	.asciz	"unsigned"              # string offset=304
.Linfo_string13:
	.asciz	"species"               # string offset=313
.Linfo_string14:
	.asciz	"name"                  # string offset=321
.Linfo_string15:
	.asciz	"std"                   # string offset=326
.Linfo_string16:
	.asciz	"__cxx11"               # string offset=330
.Linfo_string17:
	.asciz	"_M_dataplus"           # string offset=338
.Linfo_string18:
	.asciz	"allocator"             # string offset=350
.Linfo_string19:
	.asciz	"_ZSaIcE"               # string offset=360
.Linfo_string20:
	.asciz	"_M_p"                  # string offset=368
.Linfo_string21:
	.asciz	"_Alloc_hider"          # string offset=373
.Linfo_string22:
	.asciz	"_M_string_length"      # string offset=386
.Linfo_string23:
	.asciz	"unsigned long"         # string offset=403
.Linfo_string24:
	.asciz	"_M_local_buf"          # string offset=417
.Linfo_string25:
	.asciz	"_M_allocated_capacity" # string offset=430
.Linfo_string26:
	.asciz	"basic_string"          # string offset=452
.Linfo_string27:
	.asciz	"charge"                # string offset=465
.Linfo_string28:
	.asciz	"double"                # string offset=472
.Linfo_string29:
	.asciz	"mass"                  # string offset=479
.Linfo_string30:
	.asciz	"sparseMinValue"        # string offset=484
.Linfo_string31:
	.asciz	"velocityMesh"          # string offset=499
.Linfo_string32:
	.asciz	"sparseBlockAddWidthV"  # string offset=512
.Linfo_string33:
	.asciz	"sparse_conserve_mass"  # string offset=533
.Linfo_string34:
	.asciz	"sparseDynamicAlgorithm" # string offset=554
.Linfo_string35:
	.asciz	"sparseDynamicBulkValue1" # string offset=577
.Linfo_string36:
	.asciz	"sparseDynamicBulkValue2" # string offset=601
.Linfo_string37:
	.asciz	"sparseDynamicMinValue1" # string offset=625
.Linfo_string38:
	.asciz	"sparseDynamicMinValue2" # string offset=648
.Linfo_string39:
	.asciz	"thermalRadius"         # string offset=671
.Linfo_string40:
	.asciz	"thermalV"              # string offset=685
.Linfo_string41:
	.asciz	"_M_elems"              # string offset=694
.Linfo_string42:
	.asciz	"_ZSt5arrayIdLm3EE"     # string offset=703
.Linfo_string43:
	.asciz	"EnergyDensityLimit1"   # string offset=721
.Linfo_string44:
	.asciz	"EnergyDensityLimit2"   # string offset=741
.Linfo_string45:
	.asciz	"SolarWindEnergy"       # string offset=761
.Linfo_string46:
	.asciz	"SolarWindSpeed"        # string offset=777
.Linfo_string47:
	.asciz	"precipitationNChannels" # string offset=792
.Linfo_string48:
	.asciz	"precipitationEmin"     # string offset=815
.Linfo_string49:
	.asciz	"precipitationEmax"     # string offset=833
.Linfo_string50:
	.asciz	"precipitationLossConeAngle" # string offset=851
.Linfo_string51:
	.asciz	"Species"               # string offset=878
.Linfo_string52:
	.asciz	"_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE" # string offset=886
.Linfo_string53:
	.asciz	"operator="             # string offset=941
.Linfo_string54:
	.asciz	"__str"                 # string offset=951
.Linfo_string55:
	.asciz	"_ZN7species7SpeciesE"  # string offset=957
.Linfo_string56:
	.asciz	"other"                 # string offset=978
.Linfo_string57:
	.asciz	"_ZN7species7SpeciesC2Ev" # string offset=984
.Linfo_string58:
	.asciz	"_ZN7species7SpeciesC2ERKS0_" # string offset=1008
.Linfo_string59:
	.asciz	"~Species"              # string offset=1036
.Linfo_string60:
	.asciz	"_ZN7species7SpeciesD2Ev" # string offset=1045
.Linfo_string61:
	.asciz	"__sti___20_particle_species_cpp_29edf090" # string offset=1069
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp8-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	.Lfunc_end2-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc1:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp8-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # DW_OP_reg4
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc2:
	.quad	.Ltmp7-.Lfunc_begin0
	.quad	.Ltmp8-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	.Lfunc_end2-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc3:
	.quad	.Ltmp7-.Lfunc_begin0
	.quad	.Ltmp8-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # DW_OP_reg4
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc4:
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc5:
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc6:
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc7:
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc8:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp21-.Lfunc_begin0
	.quad	.Ltmp25-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	.Ltmp27-.Lfunc_begin0
	.quad	.Lfunc_end3-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc9:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # DW_OP_reg4
	.quad	.Ltmp21-.Lfunc_begin0
	.quad	.Ltmp26-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	.Ltmp27-.Lfunc_begin0
	.quad	.Ltmp28-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc10:
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp21-.Lfunc_begin0
	.quad	.Ltmp25-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	.Ltmp27-.Lfunc_begin0
	.quad	.Lfunc_end3-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc11:
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # DW_OP_reg4
	.quad	.Ltmp21-.Lfunc_begin0
	.quad	.Ltmp26-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	.Ltmp27-.Lfunc_begin0
	.quad	.Ltmp28-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc12:
	.quad	.Ltmp22-.Lfunc_begin0
	.quad	.Ltmp23-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc13:
	.quad	.Ltmp22-.Lfunc_begin0
	.quad	.Ltmp23-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc14:
	.quad	.Ltmp22-.Lfunc_begin0
	.quad	.Ltmp23-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc15:
	.quad	.Ltmp22-.Lfunc_begin0
	.quad	.Ltmp23-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc16:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp32-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp33-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	80                      # DW_OP_reg0
	.quad	0
	.quad	0
.Ldebug_loc17:
	.quad	.Ltmp31-.Lfunc_begin0
	.quad	.Ltmp32-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp33-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	80                      # DW_OP_reg0
	.quad	0
	.quad	0
	.section	.debug_abbrev,"",@progbits
	.byte	1                       # Abbreviation Code
	.byte	17                      # DW_TAG_compile_unit
	.byte	1                       # DW_CHILDREN_yes
	.byte	37                      # DW_AT_producer
	.byte	14                      # DW_FORM_strp
	.byte	19                      # DW_AT_language
	.byte	5                       # DW_FORM_data2
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	16                      # DW_AT_stmt_list
	.byte	6                       # DW_FORM_data4
	.byte	27                      # DW_AT_comp_dir
	.byte	14                      # DW_FORM_strp
	.ascii	"\264B"                 # DW_AT_GNU_pubnames
	.byte	12                      # DW_FORM_flag
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	2                       # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	63                      # DW_AT_external
	.byte	12                      # DW_FORM_flag
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	5                       # DW_FORM_data2
	.byte	2                       # DW_AT_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	3                       # Abbreviation Code
	.byte	36                      # DW_TAG_base_type
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	62                      # DW_AT_encoding
	.byte	11                      # DW_FORM_data1
	.byte	11                      # DW_AT_byte_size
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	4                       # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	2                       # DW_AT_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	5                       # Abbreviation Code
	.byte	19                      # DW_TAG_structure_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	11                      # DW_AT_byte_size
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	6                       # Abbreviation Code
	.byte	13                      # DW_TAG_member
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	56                      # DW_AT_data_member_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	7                       # Abbreviation Code
	.byte	1                       # DW_TAG_array_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	8                       # Abbreviation Code
	.byte	33                      # DW_TAG_subrange_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	55                      # DW_AT_count
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	9                       # Abbreviation Code
	.byte	36                      # DW_TAG_base_type
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	11                      # DW_AT_byte_size
	.byte	11                      # DW_FORM_data1
	.byte	62                      # DW_AT_encoding
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	10                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	63                      # DW_AT_external
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	11                      # Abbreviation Code
	.byte	15                      # DW_TAG_pointer_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	12                      # Abbreviation Code
	.byte	59                      # DW_TAG_unspecified_type
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	13                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	14                      # Abbreviation Code
	.byte	57                      # DW_TAG_namespace
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	15                      # Abbreviation Code
	.byte	19                      # DW_TAG_structure_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	11                      # DW_AT_byte_size
	.byte	11                      # DW_FORM_data1
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.ascii	"\346\177"              # DW_AT_APPLE_runtime_class
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	16                      # Abbreviation Code
	.byte	13                      # DW_TAG_member
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	56                      # DW_AT_data_member_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	17                      # Abbreviation Code
	.byte	46                      # DW_TAG_subprogram
	.byte	1                       # DW_CHILDREN_yes
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
	.byte	64                      # DW_AT_frame_base
	.byte	10                      # DW_FORM_block1
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.byte	63                      # DW_AT_external
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	18                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	10                      # DW_FORM_block1
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	52                      # DW_AT_artificial
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	19                      # Abbreviation Code
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	20                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	10                      # DW_FORM_block1
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	52                      # DW_AT_artificial
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	21                      # Abbreviation Code
	.byte	2                       # DW_TAG_class_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	11                      # DW_AT_byte_size
	.byte	11                      # DW_FORM_data1
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.ascii	"\346\177"              # DW_AT_APPLE_runtime_class
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	22                      # Abbreviation Code
	.byte	28                      # DW_TAG_inheritance
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	56                      # DW_AT_data_member_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	23                      # Abbreviation Code
	.byte	13                      # DW_TAG_member
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	56                      # DW_AT_data_member_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	24                      # Abbreviation Code
	.byte	23                      # DW_TAG_union_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	11                      # DW_AT_byte_size
	.byte	11                      # DW_FORM_data1
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.ascii	"\346\177"              # DW_AT_APPLE_runtime_class
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	25                      # Abbreviation Code
	.byte	13                      # DW_TAG_member
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	56                      # DW_AT_data_member_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	26                      # Abbreviation Code
	.byte	46                      # DW_TAG_subprogram
	.byte	0                       # DW_CHILDREN_no
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
	.byte	64                      # DW_AT_frame_base
	.byte	10                      # DW_FORM_block1
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	63                      # DW_AT_external
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	27                      # Abbreviation Code
	.byte	46                      # DW_TAG_subprogram
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	5                       # DW_FORM_data2
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	63                      # DW_AT_external
	.byte	12                      # DW_FORM_flag
	.byte	32                      # DW_AT_inline
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	28                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	52                      # DW_AT_artificial
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	29                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	30                      # Abbreviation Code
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	31                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	52                      # DW_AT_artificial
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	32                      # Abbreviation Code
	.byte	16                      # DW_TAG_reference_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	33                      # Abbreviation Code
	.byte	46                      # DW_TAG_subprogram
	.byte	1                       # DW_CHILDREN_yes
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
	.byte	64                      # DW_AT_frame_base
	.byte	10                      # DW_FORM_block1
	.byte	49                      # DW_AT_abstract_origin
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	34                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	6                       # DW_FORM_data4
	.byte	49                      # DW_AT_abstract_origin
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	35                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	6                       # DW_FORM_data4
	.byte	49                      # DW_AT_abstract_origin
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	36                      # Abbreviation Code
	.byte	29                      # DW_TAG_inlined_subroutine
	.byte	1                       # DW_CHILDREN_yes
	.byte	49                      # DW_AT_abstract_origin
	.byte	19                      # DW_FORM_ref4
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
	.byte	88                      # DW_AT_call_file
	.byte	11                      # DW_FORM_data1
	.byte	89                      # DW_AT_call_line
	.byte	11                      # DW_FORM_data1
	.byte	87                      # DW_AT_call_column
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	37                      # Abbreviation Code
	.byte	46                      # DW_TAG_subprogram
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.byte	63                      # DW_AT_external
	.byte	12                      # DW_FORM_flag
	.byte	32                      # DW_AT_inline
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	38                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	6                       # DW_FORM_data4
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	52                      # DW_AT_artificial
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	39                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	6                       # DW_FORM_data4
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	52                      # DW_AT_artificial
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	40                      # Abbreviation Code
	.byte	46                      # DW_TAG_subprogram
	.byte	1                       # DW_CHILDREN_yes
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
	.byte	64                      # DW_AT_frame_base
	.byte	10                      # DW_FORM_block1
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	63                      # DW_AT_external
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	41                      # Abbreviation Code
	.byte	29                      # DW_TAG_inlined_subroutine
	.byte	1                       # DW_CHILDREN_yes
	.byte	49                      # DW_AT_abstract_origin
	.byte	19                      # DW_FORM_ref4
	.byte	85                      # DW_AT_ranges
	.byte	6                       # DW_FORM_data4
	.byte	88                      # DW_AT_call_file
	.byte	11                      # DW_FORM_data1
	.byte	89                      # DW_AT_call_line
	.byte	11                      # DW_FORM_data1
	.byte	87                      # DW_AT_call_column
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	42                      # Abbreviation Code
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	85                      # DW_AT_ranges
	.byte	6                       # DW_FORM_data4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	0                       # EOM(3)
	.section	.debug_info,"",@progbits
.Lcu_begin0:
	.long	.Ldebug_info_end0-.Ldebug_info_start0 # Length of Unit
.Ldebug_info_start0:
	.short	2                       # DWARF version number
	.long	.debug_abbrev           # Offset Into Abbrev. Section
	.byte	8                       # Address Size (in bytes)
	.byte	1                       # Abbrev [1] 0xb:0x6af DW_TAG_compile_unit
	.long	.Linfo_string0          # DW_AT_producer
	.short	4                       # DW_AT_language
	.long	.Linfo_string1          # DW_AT_name
	.long	.Lline_table_start0     # DW_AT_stmt_list
	.long	.Linfo_string2          # DW_AT_comp_dir
	.byte	1                       # DW_AT_GNU_pubnames
	.quad	.Lfunc_begin0           # DW_AT_low_pc
	.quad	.Lfunc_end6             # DW_AT_high_pc
	.byte	2                       # Abbrev [2] 0x2f:0x17 DW_TAG_variable
	.long	.Linfo_string3          # DW_AT_name
	.long	70                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	1                       # DW_AT_decl_file
	.short	7591                    # DW_AT_decl_line
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	__I___20_particle_species_cpp_29edf090
	.byte	3                       # Abbrev [3] 0x46:0x7 DW_TAG_base_type
	.long	.Linfo_string4          # DW_AT_name
	.byte	5                       # DW_AT_encoding
	.byte	4                       # DW_AT_byte_size
	.byte	4                       # Abbrev [4] 0x4d:0x13 DW_TAG_variable
	.long	.Linfo_string5          # DW_AT_name
	.long	96                      # DW_AT_type
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE
	.byte	5                       # Abbrev [5] 0x60:0x11 DW_TAG_structure_type
	.long	.Linfo_string8          # DW_AT_name
	.byte	1                       # DW_AT_byte_size
	.byte	1                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x67:0x9 DW_TAG_member
	.long	113                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x71:0xc DW_TAG_array_type
	.long	125                     # DW_AT_type
	.byte	8                       # Abbrev [8] 0x76:0x6 DW_TAG_subrange_type
	.long	132                     # DW_AT_type
	.byte	0                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x7d:0x7 DW_TAG_base_type
	.long	.Linfo_string6          # DW_AT_name
	.byte	6                       # DW_AT_encoding
	.byte	1                       # DW_AT_byte_size
	.byte	9                       # Abbrev [9] 0x84:0x7 DW_TAG_base_type
	.long	.Linfo_string7          # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	7                       # DW_AT_encoding
	.byte	10                      # Abbrev [10] 0x8b:0xa DW_TAG_variable
	.long	.Linfo_string9          # DW_AT_name
	.long	149                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	11                      # Abbrev [11] 0x95:0x5 DW_TAG_pointer_type
	.long	154                     # DW_AT_type
	.byte	12                      # Abbrev [12] 0x9a:0x5 DW_TAG_unspecified_type
	.long	.Linfo_string10         # DW_AT_name
	.byte	13                      # Abbrev [13] 0x9f:0x9 DW_TAG_variable
	.long	.Linfo_string11         # DW_AT_name
	.long	168                     # DW_AT_type
	.byte	3                       # Abbrev [3] 0xa8:0x7 DW_TAG_base_type
	.long	.Linfo_string12         # DW_AT_name
	.byte	7                       # DW_AT_encoding
	.byte	4                       # DW_AT_byte_size
	.byte	14                      # Abbrev [14] 0xaf:0x1a0 DW_TAG_namespace
	.long	.Linfo_string13         # DW_AT_name
	.byte	15                      # Abbrev [15] 0xb4:0x19a DW_TAG_structure_type
	.long	.Linfo_string51         # DW_AT_name
	.byte	208                     # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	16                      # Abbrev [16] 0xbe:0xf DW_TAG_member
	.long	.Linfo_string14         # DW_AT_name
	.long	601                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	16                      # Abbrev [16] 0xcd:0xf DW_TAG_member
	.long	.Linfo_string27         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	16                      # Abbrev [16] 0xdc:0xf DW_TAG_member
	.long	.Linfo_string29         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	16                      # Abbrev [16] 0xeb:0xf DW_TAG_member
	.long	.Linfo_string30         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	16                      # Abbrev [16] 0xfa:0xf DW_TAG_member
	.long	.Linfo_string31         # DW_AT_name
	.long	755                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	56
	.byte	16                      # Abbrev [16] 0x109:0xf DW_TAG_member
	.long	.Linfo_string32         # DW_AT_name
	.long	70                      # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	64
	.byte	16                      # Abbrev [16] 0x118:0xf DW_TAG_member
	.long	.Linfo_string33         # DW_AT_name
	.long	125                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	68
	.byte	16                      # Abbrev [16] 0x127:0xf DW_TAG_member
	.long	.Linfo_string34         # DW_AT_name
	.long	70                      # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	72
	.byte	16                      # Abbrev [16] 0x136:0xf DW_TAG_member
	.long	.Linfo_string35         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	80
	.byte	16                      # Abbrev [16] 0x145:0xf DW_TAG_member
	.long	.Linfo_string36         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	88
	.byte	16                      # Abbrev [16] 0x154:0xf DW_TAG_member
	.long	.Linfo_string37         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	96
	.byte	16                      # Abbrev [16] 0x163:0xf DW_TAG_member
	.long	.Linfo_string38         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	104
	.byte	16                      # Abbrev [16] 0x172:0xf DW_TAG_member
	.long	.Linfo_string39         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	112
	.byte	16                      # Abbrev [16] 0x181:0xf DW_TAG_member
	.long	.Linfo_string40         # DW_AT_name
	.long	781                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	120
	.byte	16                      # Abbrev [16] 0x190:0x10 DW_TAG_member
	.long	.Linfo_string43         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\220\001"
	.byte	16                      # Abbrev [16] 0x1a0:0x10 DW_TAG_member
	.long	.Linfo_string44         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\230\001"
	.byte	16                      # Abbrev [16] 0x1b0:0x10 DW_TAG_member
	.long	.Linfo_string45         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\240\001"
	.byte	16                      # Abbrev [16] 0x1c0:0x10 DW_TAG_member
	.long	.Linfo_string46         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\250\001"
	.byte	16                      # Abbrev [16] 0x1d0:0x10 DW_TAG_member
	.long	.Linfo_string47         # DW_AT_name
	.long	70                      # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\260\001"
	.byte	16                      # Abbrev [16] 0x1e0:0x10 DW_TAG_member
	.long	.Linfo_string48         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\270\001"
	.byte	16                      # Abbrev [16] 0x1f0:0x10 DW_TAG_member
	.long	.Linfo_string49         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\300\001"
	.byte	16                      # Abbrev [16] 0x200:0x10 DW_TAG_member
	.long	.Linfo_string50         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\310\001"
	.byte	17                      # Abbrev [17] 0x210:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin0           # DW_AT_low_pc
	.quad	.Lfunc_end0             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string51         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	29                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	18                      # Abbrev [18] 0x22a:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1504                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	19                      # Abbrev [19] 0x232:0x1a DW_TAG_lexical_block
	.quad	.Ltmp0                  # DW_AT_low_pc
	.quad	.Ltmp1                  # DW_AT_high_pc
	.byte	20                      # Abbrev [20] 0x243:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1504                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	14                      # Abbrev [14] 0x24f:0x8e DW_TAG_namespace
	.long	.Linfo_string15         # DW_AT_name
	.byte	14                      # Abbrev [14] 0x254:0x88 DW_TAG_namespace
	.long	.Linfo_string16         # DW_AT_name
	.byte	21                      # Abbrev [21] 0x259:0x82 DW_TAG_class_type
	.long	.Linfo_string26         # DW_AT_name
	.byte	32                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	16                      # Abbrev [16] 0x263:0xf DW_TAG_member
	.long	.Linfo_string17         # DW_AT_name
	.long	626                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	15                      # Abbrev [15] 0x272:0x29 DW_TAG_structure_type
	.long	.Linfo_string21         # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	22                      # Abbrev [22] 0x27c:0xf DW_TAG_inheritance
	.long	.Linfo_string18         # DW_AT_name
	.long	733                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	16                      # Abbrev [16] 0x28b:0xf DW_TAG_member
	.long	.Linfo_string20         # DW_AT_name
	.long	750                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	16                      # Abbrev [16] 0x29b:0xf DW_TAG_member
	.long	.Linfo_string22         # DW_AT_name
	.long	755                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	23                      # Abbrev [23] 0x2aa:0xb DW_TAG_member
	.long	693                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	24                      # Abbrev [24] 0x2b5:0x25 DW_TAG_union_type
	.byte	16                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	16                      # Abbrev [16] 0x2bb:0xf DW_TAG_member
	.long	.Linfo_string24         # DW_AT_name
	.long	762                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	16                      # Abbrev [16] 0x2ca:0xf DW_TAG_member
	.long	.Linfo_string25         # DW_AT_name
	.long	755                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x2dd:0x11 DW_TAG_structure_type
	.long	.Linfo_string19         # DW_AT_name
	.byte	1                       # DW_AT_byte_size
	.byte	1                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x2e4:0x9 DW_TAG_member
	.long	113                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	11                      # Abbrev [11] 0x2ee:0x5 DW_TAG_pointer_type
	.long	125                     # DW_AT_type
	.byte	3                       # Abbrev [3] 0x2f3:0x7 DW_TAG_base_type
	.long	.Linfo_string23         # DW_AT_name
	.byte	7                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	7                       # Abbrev [7] 0x2fa:0xc DW_TAG_array_type
	.long	125                     # DW_AT_type
	.byte	8                       # Abbrev [8] 0x2ff:0x6 DW_TAG_subrange_type
	.long	132                     # DW_AT_type
	.byte	16                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x306:0x7 DW_TAG_base_type
	.long	.Linfo_string28         # DW_AT_name
	.byte	4                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	5                       # Abbrev [5] 0x30d:0x15 DW_TAG_structure_type
	.long	.Linfo_string42         # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	25                      # Abbrev [25] 0x314:0xd DW_TAG_member
	.long	.Linfo_string41         # DW_AT_name
	.long	802                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x322:0xc DW_TAG_array_type
	.long	774                     # DW_AT_type
	.byte	8                       # Abbrev [8] 0x327:0x6 DW_TAG_subrange_type
	.long	132                     # DW_AT_type
	.byte	3                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	26                      # Abbrev [26] 0x32e:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin1           # DW_AT_low_pc
	.quad	.Lfunc_end1             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string57         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	5                       # Abbrev [5] 0x346:0x5a DW_TAG_structure_type
	.long	.Linfo_string52         # DW_AT_name
	.byte	32                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	25                      # Abbrev [25] 0x34d:0xd DW_TAG_member
	.long	.Linfo_string17         # DW_AT_name
	.long	626                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	25                      # Abbrev [25] 0x35a:0xd DW_TAG_member
	.long	.Linfo_string22         # DW_AT_name
	.long	755                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	6                       # Abbrev [6] 0x367:0x9 DW_TAG_member
	.long	693                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	27                      # Abbrev [27] 0x370:0x2f DW_TAG_subprogram
	.long	.Linfo_string53         # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.short	666                     # DW_AT_decl_line
	.long	928                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	1                       # DW_AT_inline
	.byte	28                      # Abbrev [28] 0x37e:0x6 DW_TAG_formal_parameter
	.long	928                     # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0x384:0x9 DW_TAG_formal_parameter
	.long	.Linfo_string54         # DW_AT_name
	.long	933                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x38d:0x11 DW_TAG_lexical_block
	.byte	13                      # Abbrev [13] 0x38e:0x9 DW_TAG_variable
	.long	.Linfo_string54         # DW_AT_name
	.long	933                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x397:0x6 DW_TAG_variable
	.long	928                     # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	11                      # Abbrev [11] 0x3a0:0x5 DW_TAG_pointer_type
	.long	601                     # DW_AT_type
	.byte	32                      # Abbrev [32] 0x3a5:0x5 DW_TAG_reference_type
	.long	838                     # DW_AT_type
	.byte	5                       # Abbrev [5] 0x3aa:0x236 DW_TAG_structure_type
	.long	.Linfo_string55         # DW_AT_name
	.byte	208                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	25                      # Abbrev [25] 0x3b1:0xd DW_TAG_member
	.long	.Linfo_string14         # DW_AT_name
	.long	838                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	25                      # Abbrev [25] 0x3be:0xd DW_TAG_member
	.long	.Linfo_string27         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	25                      # Abbrev [25] 0x3cb:0xd DW_TAG_member
	.long	.Linfo_string29         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	25                      # Abbrev [25] 0x3d8:0xd DW_TAG_member
	.long	.Linfo_string30         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	25                      # Abbrev [25] 0x3e5:0xd DW_TAG_member
	.long	.Linfo_string31         # DW_AT_name
	.long	755                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	56
	.byte	25                      # Abbrev [25] 0x3f2:0xd DW_TAG_member
	.long	.Linfo_string32         # DW_AT_name
	.long	70                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	64
	.byte	25                      # Abbrev [25] 0x3ff:0xd DW_TAG_member
	.long	.Linfo_string33         # DW_AT_name
	.long	125                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	68
	.byte	25                      # Abbrev [25] 0x40c:0xd DW_TAG_member
	.long	.Linfo_string34         # DW_AT_name
	.long	70                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	72
	.byte	25                      # Abbrev [25] 0x419:0xd DW_TAG_member
	.long	.Linfo_string35         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	80
	.byte	25                      # Abbrev [25] 0x426:0xd DW_TAG_member
	.long	.Linfo_string36         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	88
	.byte	25                      # Abbrev [25] 0x433:0xd DW_TAG_member
	.long	.Linfo_string37         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	96
	.byte	25                      # Abbrev [25] 0x440:0xd DW_TAG_member
	.long	.Linfo_string38         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	104
	.byte	25                      # Abbrev [25] 0x44d:0xd DW_TAG_member
	.long	.Linfo_string39         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	112
	.byte	25                      # Abbrev [25] 0x45a:0xd DW_TAG_member
	.long	.Linfo_string40         # DW_AT_name
	.long	781                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	120
	.byte	25                      # Abbrev [25] 0x467:0xe DW_TAG_member
	.long	.Linfo_string43         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\220\001"
	.byte	25                      # Abbrev [25] 0x475:0xe DW_TAG_member
	.long	.Linfo_string44         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\230\001"
	.byte	25                      # Abbrev [25] 0x483:0xe DW_TAG_member
	.long	.Linfo_string45         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\240\001"
	.byte	25                      # Abbrev [25] 0x491:0xe DW_TAG_member
	.long	.Linfo_string46         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\250\001"
	.byte	25                      # Abbrev [25] 0x49f:0xe DW_TAG_member
	.long	.Linfo_string47         # DW_AT_name
	.long	70                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\260\001"
	.byte	25                      # Abbrev [25] 0x4ad:0xe DW_TAG_member
	.long	.Linfo_string48         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\270\001"
	.byte	25                      # Abbrev [25] 0x4bb:0xe DW_TAG_member
	.long	.Linfo_string49         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\300\001"
	.byte	25                      # Abbrev [25] 0x4c9:0xe DW_TAG_member
	.long	.Linfo_string50         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\310\001"
	.byte	33                      # Abbrev [33] 0x4d7:0x9d DW_TAG_subprogram
	.quad	.Lfunc_begin2           # DW_AT_low_pc
	.quad	.Lfunc_end2             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	1396                    # DW_AT_abstract_origin
	.byte	34                      # Abbrev [34] 0x4ee:0x9 DW_TAG_formal_parameter
	.long	.Ldebug_loc0            # DW_AT_location
	.long	1405                    # DW_AT_abstract_origin
	.byte	34                      # Abbrev [34] 0x4f7:0x9 DW_TAG_formal_parameter
	.long	.Ldebug_loc1            # DW_AT_location
	.long	1411                    # DW_AT_abstract_origin
	.byte	19                      # Abbrev [19] 0x500:0x73 DW_TAG_lexical_block
	.quad	.Ltmp8                  # DW_AT_low_pc
	.quad	.Ltmp16                 # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0x511:0x9 DW_TAG_variable
	.long	.Ldebug_loc2            # DW_AT_location
	.long	1421                    # DW_AT_abstract_origin
	.byte	35                      # Abbrev [35] 0x51a:0x9 DW_TAG_variable
	.long	.Ldebug_loc3            # DW_AT_location
	.long	1427                    # DW_AT_abstract_origin
	.byte	36                      # Abbrev [36] 0x523:0x4f DW_TAG_inlined_subroutine
	.long	880                     # DW_AT_abstract_origin
	.quad	.Ltmp9                  # DW_AT_low_pc
	.quad	.Ltmp10                 # DW_AT_high_pc
	.byte	1                       # DW_AT_call_file
	.byte	32                      # DW_AT_call_line
	.byte	1                       # DW_AT_call_column
	.byte	34                      # Abbrev [34] 0x53b:0x9 DW_TAG_formal_parameter
	.long	.Ldebug_loc7            # DW_AT_location
	.long	894                     # DW_AT_abstract_origin
	.byte	34                      # Abbrev [34] 0x544:0x9 DW_TAG_formal_parameter
	.long	.Ldebug_loc5            # DW_AT_location
	.long	900                     # DW_AT_abstract_origin
	.byte	19                      # Abbrev [19] 0x54d:0x24 DW_TAG_lexical_block
	.quad	.Ltmp9                  # DW_AT_low_pc
	.quad	.Ltmp10                 # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0x55e:0x9 DW_TAG_variable
	.long	.Ldebug_loc4            # DW_AT_location
	.long	910                     # DW_AT_abstract_origin
	.byte	35                      # Abbrev [35] 0x567:0x9 DW_TAG_variable
	.long	.Ldebug_loc6            # DW_AT_location
	.long	919                     # DW_AT_abstract_origin
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	37                      # Abbrev [37] 0x574:0x2a DW_TAG_subprogram
	.long	.Linfo_string51         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	31                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	1                       # DW_AT_inline
	.byte	28                      # Abbrev [28] 0x57d:0x6 DW_TAG_formal_parameter
	.long	1504                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0x583:0x9 DW_TAG_formal_parameter
	.long	.Linfo_string56         # DW_AT_name
	.long	1509                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x58c:0x11 DW_TAG_lexical_block
	.byte	31                      # Abbrev [31] 0x58d:0x6 DW_TAG_variable
	.long	1504                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	13                      # Abbrev [13] 0x593:0x9 DW_TAG_variable
	.long	.Linfo_string56         # DW_AT_name
	.long	1509                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	17                      # Abbrev [17] 0x59e:0x41 DW_TAG_subprogram
	.quad	.Lfunc_begin4           # DW_AT_low_pc
	.quad	.Lfunc_end4             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string59         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	39                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	38                      # Abbrev [38] 0x5b8:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc16           # DW_AT_location
	.long	1504                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	19                      # Abbrev [19] 0x5c2:0x1c DW_TAG_lexical_block
	.quad	.Ltmp32                 # DW_AT_low_pc
	.quad	.Ltmp34                 # DW_AT_high_pc
	.byte	39                      # Abbrev [39] 0x5d3:0xa DW_TAG_variable
	.long	.Ldebug_loc17           # DW_AT_location
	.long	1504                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	11                      # Abbrev [11] 0x5e0:0x5 DW_TAG_pointer_type
	.long	180                     # DW_AT_type
	.byte	32                      # Abbrev [32] 0x5e5:0x5 DW_TAG_reference_type
	.long	938                     # DW_AT_type
	.byte	40                      # Abbrev [40] 0x5ea:0x9f DW_TAG_subprogram
	.quad	.Lfunc_begin3           # DW_AT_low_pc
	.quad	.Lfunc_end3             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string58         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	41                      # Abbrev [41] 0x602:0x86 DW_TAG_inlined_subroutine
	.long	1396                    # DW_AT_abstract_origin
	.long	.Ldebug_ranges0         # DW_AT_ranges
	.byte	1                       # DW_AT_call_file
	.byte	36                      # DW_AT_call_line
	.byte	1                       # DW_AT_call_column
	.byte	34                      # Abbrev [34] 0x60e:0x9 DW_TAG_formal_parameter
	.long	.Ldebug_loc8            # DW_AT_location
	.long	1405                    # DW_AT_abstract_origin
	.byte	34                      # Abbrev [34] 0x617:0x9 DW_TAG_formal_parameter
	.long	.Ldebug_loc9            # DW_AT_location
	.long	1411                    # DW_AT_abstract_origin
	.byte	42                      # Abbrev [42] 0x620:0x67 DW_TAG_lexical_block
	.long	.Ldebug_ranges1         # DW_AT_ranges
	.byte	35                      # Abbrev [35] 0x625:0x9 DW_TAG_variable
	.long	.Ldebug_loc10           # DW_AT_location
	.long	1421                    # DW_AT_abstract_origin
	.byte	35                      # Abbrev [35] 0x62e:0x9 DW_TAG_variable
	.long	.Ldebug_loc11           # DW_AT_location
	.long	1427                    # DW_AT_abstract_origin
	.byte	36                      # Abbrev [36] 0x637:0x4f DW_TAG_inlined_subroutine
	.long	880                     # DW_AT_abstract_origin
	.quad	.Ltmp22                 # DW_AT_low_pc
	.quad	.Ltmp23                 # DW_AT_high_pc
	.byte	1                       # DW_AT_call_file
	.byte	32                      # DW_AT_call_line
	.byte	1                       # DW_AT_call_column
	.byte	34                      # Abbrev [34] 0x64f:0x9 DW_TAG_formal_parameter
	.long	.Ldebug_loc15           # DW_AT_location
	.long	894                     # DW_AT_abstract_origin
	.byte	34                      # Abbrev [34] 0x658:0x9 DW_TAG_formal_parameter
	.long	.Ldebug_loc13           # DW_AT_location
	.long	900                     # DW_AT_abstract_origin
	.byte	19                      # Abbrev [19] 0x661:0x24 DW_TAG_lexical_block
	.quad	.Ltmp22                 # DW_AT_low_pc
	.quad	.Ltmp23                 # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0x672:0x9 DW_TAG_variable
	.long	.Ldebug_loc12           # DW_AT_location
	.long	910                     # DW_AT_abstract_origin
	.byte	35                      # Abbrev [35] 0x67b:0x9 DW_TAG_variable
	.long	.Ldebug_loc14           # DW_AT_location
	.long	919                     # DW_AT_abstract_origin
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	26                      # Abbrev [26] 0x689:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin5           # DW_AT_low_pc
	.quad	.Lfunc_end5             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string60         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	26                      # Abbrev [26] 0x6a1:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin6           # DW_AT_low_pc
	.quad	.Lfunc_end6             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string61         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	0                       # End Of Children Mark
.Ldebug_info_end0:
	.section	.debug_ranges,"",@progbits
.Ldebug_ranges0:
	.quad	.Ltmp21-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.quad	.Ltmp28-.Lfunc_begin0
	.quad	.Ltmp30-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges1:
	.quad	.Ltmp21-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.quad	.Ltmp28-.Lfunc_begin0
	.quad	.Ltmp30-.Lfunc_begin0
	.quad	0
	.quad	0
	.section	.debug_pubnames,"",@progbits
	.long	.LpubNames_end0-.LpubNames_begin0 # Length of Public Names Info
.LpubNames_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	1722                    # Compilation Unit Length
	.long	1697                    # DIE offset
	.asciz	"__sti___20_particle_species_cpp_29edf090" # External Name
	.long	1514                    # DIE offset
	.asciz	"_ZN7species7SpeciesC2ERKS0_" # External Name
	.long	77                      # DIE offset
	.asciz	"_ZN42_INTERNAL_20_particle_species_cpp_29edf090St8__ioinitE" # External Name
	.long	591                     # DIE offset
	.asciz	"std"                   # External Name
	.long	528                     # DIE offset
	.asciz	"species::Species::Species" # External Name
	.long	175                     # DIE offset
	.asciz	"species"               # External Name
	.long	139                     # DIE offset
	.asciz	"__dso_handle"          # External Name
	.long	1438                    # DIE offset
	.asciz	"_ZN7species7SpeciesE::~Species" # External Name
	.long	596                     # DIE offset
	.asciz	"std::__cxx11"          # External Name
	.long	1396                    # DIE offset
	.asciz	"_ZN7species7SpeciesE::Species" # External Name
	.long	47                      # DIE offset
	.asciz	"__I___20_particle_species_cpp_29edf090" # External Name
	.long	880                     # DIE offset
	.asciz	"_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE::operator=" # External Name
	.long	159                     # DIE offset
	.asciz	"_ZN42_INTERNAL_20_particle_species_cpp_29edf0905vmesh15INVALID_LOCALIDE" # External Name
	.long	814                     # DIE offset
	.asciz	"_ZN7species7SpeciesC2Ev" # External Name
	.long	1673                    # DIE offset
	.asciz	"_ZN7species7SpeciesD2Ev" # External Name
	.long	0                       # End Mark
.LpubNames_end0:
	.section	.debug_pubtypes,"",@progbits
	.long	.LpubTypes_end0-.LpubTypes_begin0 # Length of Public Types Info
.LpubTypes_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	1722                    # Compilation Unit Length
	.long	70                      # DIE offset
	.asciz	"int"                   # External Name
	.long	180                     # DIE offset
	.asciz	"species::Species"      # External Name
	.long	755                     # DIE offset
	.asciz	"unsigned long"         # External Name
	.long	154                     # DIE offset
	.asciz	"void"                  # External Name
	.long	938                     # DIE offset
	.asciz	"_ZN7species7SpeciesE"  # External Name
	.long	781                     # DIE offset
	.asciz	"_ZSt5arrayIdLm3EE"     # External Name
	.long	96                      # DIE offset
	.asciz	"_ZNSt8ios_base4InitE"  # External Name
	.long	125                     # DIE offset
	.asciz	"signed char"           # External Name
	.long	774                     # DIE offset
	.asciz	"double"                # External Name
	.long	838                     # DIE offset
	.asciz	"_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE" # External Name
	.long	168                     # DIE offset
	.asciz	"unsigned"              # External Name
	.long	601                     # DIE offset
	.asciz	"std::__cxx11::basic_string" # External Name
	.long	733                     # DIE offset
	.asciz	"_ZSaIcE"               # External Name
	.long	0                       # End Mark
.LpubTypes_end0:
	.section	".note.GNU-stack","",@progbits
	.section	.debug_line,"",@progbits
.Lline_table_start0:
