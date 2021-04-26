	.text
	.file	"linedipole.ll"
	.file	1 "/home/talgat/vlasiator/backgroundfield/linedipole.cpp"
	.weak	_ZN11T3DFunctionD1Ev    # -- Begin function _ZN11T3DFunctionD1Ev
	.p2align	4, 0x90
	.type	_ZN11T3DFunctionD1Ev,@function
_ZN11T3DFunctionD1Ev:                   # @_ZN11T3DFunctionD1Ev
.Lfunc_begin0:
	.file	2 "/home/talgat/vlasiator/backgroundfield/functions.hpp"
	.loc	2 31 0                  # backgroundfield/functions.hpp:31:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T3DFunction: <- $rdi
	#DEBUG_VALUE: ~T3DFunction: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp0:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 31 1 prologue_end     # backgroundfield/functions.hpp:31:1
	movq	$_ZTV11T3DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp1:
.Lfunc_end0:
	.size	_ZN11T3DFunctionD1Ev, .Lfunc_end0-_ZN11T3DFunctionD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T3DFunctionD0Ev    # -- Begin function _ZN11T3DFunctionD0Ev
	.p2align	4, 0x90
	.type	_ZN11T3DFunctionD0Ev,@function
_ZN11T3DFunctionD0Ev:                   # @_ZN11T3DFunctionD0Ev
.Lfunc_begin1:
	.loc	1 0 0                   # backgroundfield/linedipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp2:
	.loc	1 31 1 prologue_end     # backgroundfield/linedipole.cpp:31:1
	movl	$8, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp3:
.Lfunc_end1:
	.size	_ZN11T3DFunctionD0Ev, .Lfunc_end1-_ZN11T3DFunctionD0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN13FieldFunctionD1Ev  # -- Begin function _ZN13FieldFunctionD1Ev
	.p2align	4, 0x90
	.type	_ZN13FieldFunctionD1Ev,@function
_ZN13FieldFunctionD1Ev:                 # @_ZN13FieldFunctionD1Ev
.Lfunc_begin2:
	.file	3 "/home/talgat/vlasiator/backgroundfield/fieldfunction.hpp"
	.loc	3 0 0                   # backgroundfield/fieldfunction.hpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~FieldFunction: <- $rdi
	#DEBUG_VALUE: ~FieldFunction: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp4:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	3 31 1 prologue_end     # backgroundfield/fieldfunction.hpp:31:1
	movq	$_ZTV11T3DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp5:
.Lfunc_end2:
	.size	_ZN13FieldFunctionD1Ev, .Lfunc_end2-_ZN13FieldFunctionD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN13FieldFunctionD0Ev  # -- Begin function _ZN13FieldFunctionD0Ev
	.p2align	4, 0x90
	.type	_ZN13FieldFunctionD0Ev,@function
_ZN13FieldFunctionD0Ev:                 # @_ZN13FieldFunctionD0Ev
.Lfunc_begin3:
	.loc	1 0 0                   # backgroundfield/linedipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp6:
	.loc	1 31 1 prologue_end     # backgroundfield/linedipole.cpp:31:1
	movl	$24, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp7:
.Lfunc_end3:
	.size	_ZN13FieldFunctionD0Ev, .Lfunc_end3-_ZN13FieldFunctionD0Ev
	.cfi_endproc
                                        # -- End function
	.globl	_ZN10LineDipole10initializeEdddd # -- Begin function _ZN10LineDipole10initializeEdddd
	.p2align	4, 0x90
	.type	_ZN10LineDipole10initializeEdddd,@function
_ZN10LineDipole10initializeEdddd:       # @_ZN10LineDipole10initializeEdddd
.Lfunc_begin4:
	.loc	1 32 0                  # backgroundfield/linedipole.cpp:32:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: initialize: <- $rdi
	#DEBUG_VALUE: initialize: <- $rdi
	#DEBUG_VALUE: initialize:moment <- $xmm0
	#DEBUG_VALUE: initialize:moment <- $xmm0
	#DEBUG_VALUE: initialize:center_x <- $xmm1
	#DEBUG_VALUE: initialize:center_x <- $xmm1
	#DEBUG_VALUE: initialize:center_y <- $xmm2
	#DEBUG_VALUE: initialize:center_y <- $xmm2
	#DEBUG_VALUE: initialize:center_z <- $xmm3
	#DEBUG_VALUE: initialize:center_z <- $xmm3
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp8:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: moment <- $xmm0
	#DEBUG_VALUE: moment <- $xmm0
	#DEBUG_VALUE: center_x <- $xmm1
	#DEBUG_VALUE: center_x <- $xmm1
	#DEBUG_VALUE: center_y <- $xmm2
	#DEBUG_VALUE: center_y <- $xmm2
	#DEBUG_VALUE: center_z <- $xmm3
	#DEBUG_VALUE: center_z <- $xmm3
	.loc	1 33 1 prologue_end     # backgroundfield/linedipole.cpp:33:1
	movb	$1, 20(%rdi)
	.loc	1 35 1                  # backgroundfield/linedipole.cpp:35:1
	vxorps	%xmm4, %xmm4, %xmm4
	vmovups	%xmm4, 24(%rdi)
	.loc	1 36 1                  # backgroundfield/linedipole.cpp:36:1
	vmovsd	%xmm0, 40(%rdi)
	.loc	1 37 1                  # backgroundfield/linedipole.cpp:37:1
	vmovsd	%xmm1, 48(%rdi)
	.loc	1 38 1                  # backgroundfield/linedipole.cpp:38:1
	vmovsd	%xmm2, 56(%rdi)
	.loc	1 39 1                  # backgroundfield/linedipole.cpp:39:1
	vmovsd	%xmm3, 64(%rdi)
	.loc	1 40 1                  # backgroundfield/linedipole.cpp:40:1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp9:
.Lfunc_end4:
	.size	_ZN10LineDipole10initializeEdddd, .Lfunc_end4-_ZN10LineDipole10initializeEdddd
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _ZNK10LineDipole4callEddd
.LCPI5_0:
	.quad	4720717001298091704     # double 40592189.439999998
.LCPI5_2:
	.quad	4613937818241073152     # double 3
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4
.LCPI5_1:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
	.text
	.globl	_ZNK10LineDipole4callEddd
	.p2align	4, 0x90
	.type	_ZNK10LineDipole4callEddd,@function
_ZNK10LineDipole4callEddd:              # @_ZNK10LineDipole4callEddd
.Lfunc_begin5:
	.loc	1 45 0                  # backgroundfield/linedipole.cpp:45:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:x <- $xmm0
	#DEBUG_VALUE: call:z <- $xmm2
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp10:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: z <- $xmm2
	vmovapd	%xmm0, %xmm1
.Ltmp11:
	#DEBUG_VALUE: y <- undef
	#DEBUG_VALUE: call:y <- undef
	vxorpd	%xmm0, %xmm0, %xmm0
.Ltmp12:
	#DEBUG_VALUE: z <- $xmm2
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: x <- $xmm1
	#DEBUG_VALUE: call:x <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 47 1 prologue_end     # backgroundfield/linedipole.cpp:47:1
	cmpb	$0, 20(%rdi)
	je	.LBB5_13
.Ltmp13:
# %bb.1:                                # %L.B0000
	#DEBUG_VALUE: call:x <- $xmm1
	#DEBUG_VALUE: x <- $xmm1
	#DEBUG_VALUE: z <- $xmm2
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 51 1                  # backgroundfield/linedipole.cpp:51:1
	vsubsd	48(%rdi), %xmm1, %xmm3
.Ltmp14:
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 64 64] undef
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	.loc	1 53 1                  # backgroundfield/linedipole.cpp:53:1
	vsubsd	64(%rdi), %xmm2, %xmm1
.Ltmp15:
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	.loc	1 55 1                  # backgroundfield/linedipole.cpp:55:1
	vmulsd	%xmm3, %xmm3, %xmm5
	vmovapd	%xmm1, %xmm6
	vfmadd213sd	%xmm5, %xmm1, %xmm6 # xmm6 = (xmm1 * xmm6) + xmm5
.Ltmp16:
	#DEBUG_VALUE: r2 <- $xmm6
	.loc	1 0 1 is_stmt 0         # backgroundfield/linedipole.cpp:0:1
	vmovsd	.LCPI5_0(%rip), %xmm2   # xmm2 = mem[0],zero
.Ltmp17:
	.loc	1 57 1 is_stmt 1        # backgroundfield/linedipole.cpp:57:1
	vucomisd	%xmm6, %xmm2
	ja	.LBB5_13
.Ltmp18:
# %bb.2:                                # %L.B0001
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 61 1                  # backgroundfield/linedipole.cpp:61:1
	vmulsd	%xmm6, %xmm6, %xmm4
.Ltmp19:
	#DEBUG_VALUE: r6 <- undef
	.loc	1 63 1                  # backgroundfield/linedipole.cpp:63:1
	vmovsd	40(%rdi), %xmm7         # xmm7 = mem[0],zero
	vxorpd	.LCPI5_1(%rip), %xmm7, %xmm2
.Ltmp20:
	#DEBUG_VALUE: DerivativeDiffComponent <- undef
	#DEBUG_VALUE: DerivativeSameComponent <- undef
	#DEBUG_VALUE: D <- $xmm2
	.loc	1 70 1                  # backgroundfield/linedipole.cpp:70:1
	movl	16(%rdi), %eax
	cmpl	$1, %eax
	je	.LBB5_9
.Ltmp21:
# %bb.3:                                # %L.B0001
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: D <- $xmm2
	testl	%eax, %eax
	jne	.LBB5_13
.Ltmp22:
# %bb.4:                                # %L.B0047
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: D <- $xmm2
	.loc	1 71 1                  # backgroundfield/linedipole.cpp:71:1
	movl	8(%rdi), %eax
	cmpl	$2, %eax
	je	.LBB5_7
.Ltmp23:
# %bb.5:                                # %L.B0047
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	testl	%eax, %eax
	jne	.LBB5_13
.Ltmp24:
# %bb.6:                                # %L.B0048
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 72 1                  # backgroundfield/linedipole.cpp:72:1
	vsubsd	%xmm7, %xmm2, %xmm0
	vmulsd	%xmm0, %xmm3, %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	jmp	.LBB5_8
.Ltmp25:
.LBB5_9:                                # %L.B0054
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: r6 <- undef
	#DEBUG_VALUE: DerivativeSameComponent <- undef
	#DEBUG_VALUE: DerivativeDiffComponent <- undef
	.loc	1 80 1                  # backgroundfield/linedipole.cpp:80:1
	movl	12(%rdi), %eax
	cmpl	$1, %eax
	je	.LBB5_13
.Ltmp26:
# %bb.10:                               # %L.B0055
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: r6 <- undef
	#DEBUG_VALUE: DerivativeSameComponent <- undef
	#DEBUG_VALUE: DerivativeDiffComponent <- undef
	movl	8(%rdi), %ecx
	cmpl	$1, %ecx
	je	.LBB5_13
.Ltmp27:
# %bb.11:                               # %L.B0008
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/linedipole.cpp:0:1
	vmulsd	%xmm4, %xmm6, %xmm4
.Ltmp28:
	#DEBUG_VALUE: r6 <- $xmm4
	vmovsd	.LCPI5_2(%rip), %xmm5   # xmm5 = mem[0],zero
.Ltmp29:
	#DEBUG_VALUE: DerivativeDiffComponent <- undef
	#DEBUG_VALUE: DerivativeSameComponent <- undef
	.loc	1 83 1 is_stmt 1        # backgroundfield/linedipole.cpp:83:1
	cmpl	%ecx, %eax
	jne	.LBB5_12
.Ltmp30:
# %bb.14:                               # %L.B0057
	#DEBUG_VALUE: r6 <- $xmm4
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 0 is_stmt 0         # backgroundfield/linedipole.cpp:0:0
	vmulsd	%xmm5, %xmm3, %xmm5
	vmulsd	%xmm5, %xmm3, %xmm3
.Ltmp31:
	vfmsub231sd	%xmm1, %xmm1, %xmm3 # xmm3 = (xmm1 * xmm1) - xmm3
	vaddsd	%xmm1, %xmm1, %xmm1
.Ltmp32:
	vmulsd	%xmm3, %xmm1, %xmm1
	vmulsd	%xmm2, %xmm1, %xmm1
	vdivsd	%xmm4, %xmm1, %xmm1
.Ltmp33:
	#DEBUG_VALUE: DerivativeSameComponent <- $xmm1
	.loc	1 84 1 is_stmt 1        # backgroundfield/linedipole.cpp:84:1
	testl	%eax, %eax
	je	.LBB5_15
.Ltmp34:
# %bb.16:                               # %L.B0057
	#DEBUG_VALUE: DerivativeSameComponent <- $xmm1
	#DEBUG_VALUE: r6 <- $xmm4
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	cmpl	$2, %eax
	jne	.LBB5_13
.Ltmp35:
# %bb.17:                               # %L.B0060
	#DEBUG_VALUE: DerivativeSameComponent <- $xmm1
	#DEBUG_VALUE: r6 <- $xmm4
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 88 1                  # backgroundfield/linedipole.cpp:88:1
	vxorpd	.LCPI5_1(%rip), %xmm1, %xmm0
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp36:
.LBB5_7:                                # %L.B0050
	.cfi_def_cfa %rbp, 16
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 74 1                  # backgroundfield/linedipole.cpp:74:1
	vfmsub213sd	%xmm5, %xmm1, %xmm1 # xmm1 = (xmm1 * xmm1) - xmm5
.Ltmp37:
	vmulsd	%xmm2, %xmm1, %xmm0
.Ltmp38:
.LBB5_8:                                # %L.B0050
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 0 is_stmt 0         # backgroundfield/linedipole.cpp:0:0
	vdivsd	%xmm4, %xmm0, %xmm0
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp39:
.LBB5_12:
	.cfi_def_cfa %rbp, 16
	#DEBUG_VALUE: r6 <- $xmm4
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 128 64] $xmm1
	#DEBUG_VALUE: r <- [DW_OP_LLVM_fragment 0 64] $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	vmulsd	%xmm5, %xmm1, %xmm0
	vmulsd	%xmm0, %xmm1, %xmm0
	vfmsub231sd	%xmm3, %xmm3, %xmm0 # xmm0 = (xmm3 * xmm3) - xmm0
	vaddsd	%xmm3, %xmm3, %xmm1
.Ltmp40:
	vmulsd	%xmm0, %xmm1, %xmm0
	vmulsd	%xmm2, %xmm0, %xmm0
	vdivsd	%xmm4, %xmm0, %xmm0
.Ltmp41:
	#DEBUG_VALUE: DerivativeDiffComponent <- $xmm0
.LBB5_13:                               # %L.B0043
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 48 1 is_stmt 1        # backgroundfield/linedipole.cpp:48:1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp42:
.LBB5_15:
	.cfi_def_cfa %rbp, 16
	#DEBUG_VALUE: DerivativeSameComponent <- $xmm1
	#DEBUG_VALUE: r6 <- $xmm4
	#DEBUG_VALUE: D <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm6
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/linedipole.cpp:0:1
	vmovapd	%xmm1, %xmm0
	.loc	1 48 1                  # backgroundfield/linedipole.cpp:48:1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp43:
.Lfunc_end5:
	.size	_ZNK10LineDipole4callEddd, .Lfunc_end5-_ZNK10LineDipole4callEddd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN10LineDipoleD1Ev     # -- Begin function _ZN10LineDipoleD1Ev
	.p2align	4, 0x90
	.type	_ZN10LineDipoleD1Ev,@function
_ZN10LineDipoleD1Ev:                    # @_ZN10LineDipoleD1Ev
.Lfunc_begin6:
	.file	4 "/home/talgat/vlasiator/backgroundfield/linedipole.hpp"
	.loc	4 47 0 is_stmt 1        # backgroundfield/linedipole.hpp:47:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~LineDipole: <- $rdi
	#DEBUG_VALUE: ~LineDipole: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp44:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	4 47 1 prologue_end     # backgroundfield/linedipole.hpp:47:1
	movq	$_ZTV11T3DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp45:
.Lfunc_end6:
	.size	_ZN10LineDipoleD1Ev, .Lfunc_end6-_ZN10LineDipoleD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN10LineDipoleD0Ev     # -- Begin function _ZN10LineDipoleD0Ev
	.p2align	4, 0x90
	.type	_ZN10LineDipoleD0Ev,@function
_ZN10LineDipoleD0Ev:                    # @_ZN10LineDipoleD0Ev
.Lfunc_begin7:
	.loc	1 0 0                   # backgroundfield/linedipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp46:
	.loc	1 47 1 prologue_end     # backgroundfield/linedipole.cpp:47:1
	movl	$72, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp47:
.Lfunc_end7:
	.size	_ZN10LineDipoleD0Ev, .Lfunc_end7-_ZN10LineDipoleD0Ev
	.cfi_endproc
                                        # -- End function
	.globl	__sti___30_backgroundfield_linedipole_cpp_4a71163a # -- Begin function __sti___30_backgroundfield_linedipole_cpp_4a71163a
	.p2align	4, 0x90
	.type	__sti___30_backgroundfield_linedipole_cpp_4a71163a,@function
__sti___30_backgroundfield_linedipole_cpp_4a71163a: # @__sti___30_backgroundfield_linedipole_cpp_4a71163a
.Lfunc_begin8:
	.loc	1 0 0                   # backgroundfield/linedipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp48:
	.loc	1 74 1 prologue_end     # backgroundfield/linedipole.cpp:74:1
	cmpl	$1, __I___30_backgroundfield_linedipole_cpp_4a71163a(%rip)
	jne	.LBB8_2
# %bb.1:                                # %L.B0016
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB8_2:                                # %L.B0069
	.cfi_def_cfa %rbp, 16
	movl	$1, __I___30_backgroundfield_linedipole_cpp_4a71163a(%rip)
	movl	$_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE, %edi
	callq	_ZNSt8ios_base4InitC1Ev
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	movl	$_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE, %esi
	movl	$__dso_handle, %edx
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	__cxa_atexit            # TAILCALL
.Ltmp49:
.Lfunc_end8:
	.size	__sti___30_backgroundfield_linedipole_cpp_4a71163a, .Lfunc_end8-__sti___30_backgroundfield_linedipole_cpp_4a71163a
	.cfi_endproc
                                        # -- End function
	.type	_ZTV11T3DFunction,@object # @_ZTV11T3DFunction
	.data
	.weak	_ZTV11T3DFunction
	.p2align	4
_ZTV11T3DFunction:
	.quad	0
	.quad	_ZTI11T3DFunction
	.quad	__cxa_pure_virtual
	.quad	_ZN11T3DFunctionD1Ev
	.quad	_ZN11T3DFunctionD0Ev
	.size	_ZTV11T3DFunction, 40

	.type	_ZTV13FieldFunction,@object # @_ZTV13FieldFunction
	.weak	_ZTV13FieldFunction
	.p2align	4
_ZTV13FieldFunction:
	.quad	0
	.quad	_ZTI13FieldFunction
	.quad	__cxa_pure_virtual
	.quad	_ZN13FieldFunctionD1Ev
	.quad	_ZN13FieldFunctionD0Ev
	.size	_ZTV13FieldFunction, 40

	.type	_ZTV10LineDipole,@object # @_ZTV10LineDipole
	.weak	_ZTV10LineDipole
	.p2align	4
_ZTV10LineDipole:
	.quad	0
	.quad	_ZTI10LineDipole
	.quad	_ZNK10LineDipole4callEddd
	.quad	_ZN10LineDipoleD1Ev
	.quad	_ZN10LineDipoleD0Ev
	.size	_ZTV10LineDipole, 40

	.type	_ZTI11T3DFunction,@object # @_ZTI11T3DFunction
	.weak	_ZTI11T3DFunction
	.p2align	4
_ZTI11T3DFunction:
	.quad	_ZTVN10__cxxabiv117__class_type_infoE+16
	.quad	_ZTS11T3DFunction
	.size	_ZTI11T3DFunction, 16

	.type	_ZTI13FieldFunction,@object # @_ZTI13FieldFunction
	.weak	_ZTI13FieldFunction
	.p2align	4
_ZTI13FieldFunction:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS13FieldFunction
	.quad	_ZTI11T3DFunction
	.size	_ZTI13FieldFunction, 24

	.type	_ZTI10LineDipole,@object # @_ZTI10LineDipole
	.weak	_ZTI10LineDipole
	.p2align	4
_ZTI10LineDipole:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS10LineDipole
	.quad	_ZTI13FieldFunction
	.size	_ZTI10LineDipole, 24

	.type	_ZTS11T3DFunction,@object # @_ZTS11T3DFunction
	.weak	_ZTS11T3DFunction
	.p2align	3
_ZTS11T3DFunction:
	.asciz	"11T3DFunction"
	.size	_ZTS11T3DFunction, 14

	.type	_ZTS13FieldFunction,@object # @_ZTS13FieldFunction
	.weak	_ZTS13FieldFunction
	.p2align	4
_ZTS13FieldFunction:
	.asciz	"13FieldFunction"
	.size	_ZTS13FieldFunction, 16

	.type	_ZTS10LineDipole,@object # @_ZTS10LineDipole
	.weak	_ZTS10LineDipole
	.p2align	3
_ZTS10LineDipole:
	.asciz	"10LineDipole"
	.size	_ZTS10LineDipole, 13

	.type	__I___30_backgroundfield_linedipole_cpp_4a71163a,@object # @__I___30_backgroundfield_linedipole_cpp_4a71163a
	.bss
	.globl	__I___30_backgroundfield_linedipole_cpp_4a71163a
	.p2align	2
__I___30_backgroundfield_linedipole_cpp_4a71163a:
	.long	0                       # 0x0
	.size	__I___30_backgroundfield_linedipole_cpp_4a71163a, 4

	.type	_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE,@object # @_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE
	.local	_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE
	.comm	_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE,1,1
	.section	.init_array,"aw",@init_array
	.p2align	3
	.quad	__sti___30_backgroundfield_linedipole_cpp_4a71163a
	.section	.debug_str,"MS",@progbits,1
.Linfo_string0:
	.asciz	" NVC++ 21.2-0"         # string offset=0
.Linfo_string1:
	.asciz	"backgroundfield/linedipole.cpp" # string offset=14
.Linfo_string2:
	.asciz	"/home/talgat/vlasiator" # string offset=45
.Linfo_string3:
	.asciz	"_ZTV11T3DFunction"     # string offset=68
.Linfo_string4:
	.asciz	"int"                   # string offset=86
.Linfo_string5:
	.asciz	"__ARRAY_SIZE_TYPE__"   # string offset=90
.Linfo_string6:
	.asciz	"_ZTV13FieldFunction"   # string offset=110
.Linfo_string7:
	.asciz	"_ZTV10LineDipole"      # string offset=130
.Linfo_string8:
	.asciz	"__I___30_backgroundfield_linedipole_cpp_4a71163a" # string offset=147
.Linfo_string9:
	.asciz	"_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE" # string offset=196
.Linfo_string10:
	.asciz	"signed char"           # string offset=266
.Linfo_string11:
	.asciz	"_ZNSt8ios_base4InitE"  # string offset=278
.Linfo_string12:
	.asciz	"__dso_handle"          # string offset=299
.Linfo_string13:
	.asciz	"void"                  # string offset=312
.Linfo_string14:
	.asciz	"_ZTI11T3DFunction"     # string offset=317
.Linfo_string15:
	.asciz	"__class_type_info"     # string offset=335
.Linfo_string16:
	.asciz	"_ZTI13FieldFunction"   # string offset=353
.Linfo_string17:
	.asciz	"__si_class_type_info"  # string offset=373
.Linfo_string18:
	.asciz	"WID"                   # string offset=394
.Linfo_string19:
	.asciz	"_ZTI10LineDipole"      # string offset=398
.Linfo_string20:
	.asciz	"WID2"                  # string offset=415
.Linfo_string21:
	.asciz	"WID3"                  # string offset=420
.Linfo_string22:
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # string offset=425
.Linfo_string23:
	.asciz	"_ZTS11T3DFunction"     # string offset=463
.Linfo_string24:
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # string offset=481
.Linfo_string25:
	.asciz	"_ZTS13FieldFunction"   # string offset=522
.Linfo_string26:
	.asciz	"_ZTS10LineDipole"      # string offset=542
.Linfo_string27:
	.asciz	"_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163a5vmesh15INVALID_LOCALIDE" # string offset=559
.Linfo_string28:
	.asciz	"unsigned"              # string offset=641
.Linfo_string29:
	.asciz	"_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163a17physicalconstants3R_EE" # string offset=650
.Linfo_string30:
	.asciz	"double"                # string offset=732
.Linfo_string31:
	.asciz	"T3DFunction"           # string offset=739
.Linfo_string32:
	.asciz	"FieldFunction"         # string offset=751
.Linfo_string33:
	.asciz	"_fComponent"           # string offset=765
.Linfo_string34:
	.asciz	"X"                     # string offset=777
.Linfo_string35:
	.asciz	"Y"                     # string offset=779
.Linfo_string36:
	.asciz	"Z"                     # string offset=781
.Linfo_string37:
	.asciz	"coordinate"            # string offset=783
.Linfo_string38:
	.asciz	"_dComponent"           # string offset=794
.Linfo_string39:
	.asciz	"_derivative"           # string offset=806
.Linfo_string40:
	.asciz	"initialized"           # string offset=818
.Linfo_string41:
	.asciz	"q"                     # string offset=830
.Linfo_string42:
	.asciz	"center"                # string offset=832
.Linfo_string43:
	.asciz	"LineDipole"            # string offset=839
.Linfo_string44:
	.asciz	"~T3DFunction"          # string offset=850
.Linfo_string45:
	.asciz	"_ZN11T3DFunctionD0Ev"  # string offset=863
.Linfo_string46:
	.asciz	"~FieldFunction"        # string offset=884
.Linfo_string47:
	.asciz	"_ZN13FieldFunctionD0Ev" # string offset=899
.Linfo_string48:
	.asciz	"initialize"            # string offset=922
.Linfo_string49:
	.asciz	"call"                  # string offset=933
.Linfo_string50:
	.asciz	"~LineDipole"           # string offset=938
.Linfo_string51:
	.asciz	"_ZN10LineDipoleD0Ev"   # string offset=950
.Linfo_string52:
	.asciz	"__sti___30_backgroundfield_linedipole_cpp_4a71163a" # string offset=970
.Linfo_string53:
	.asciz	"moment"                # string offset=1021
.Linfo_string54:
	.asciz	"center_x"              # string offset=1028
.Linfo_string55:
	.asciz	"center_y"              # string offset=1037
.Linfo_string56:
	.asciz	"center_z"              # string offset=1046
.Linfo_string57:
	.asciz	"x"                     # string offset=1055
.Linfo_string58:
	.asciz	"z"                     # string offset=1057
.Linfo_string59:
	.asciz	"y"                     # string offset=1059
.Linfo_string60:
	.asciz	"r"                     # string offset=1061
.Linfo_string61:
	.asciz	"r2"                    # string offset=1063
.Linfo_string62:
	.asciz	"r6"                    # string offset=1066
.Linfo_string63:
	.asciz	"DerivativeDiffComponent" # string offset=1069
.Linfo_string64:
	.asciz	"DerivativeSameComponent" # string offset=1093
.Linfo_string65:
	.asciz	"D"                     # string offset=1117
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
	.quad	.Lfunc_begin5-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp12-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc1:
	.quad	.Lfunc_begin5-.Lfunc_begin0
	.quad	.Ltmp17-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc2:
	.quad	.Ltmp10-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp12-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc3:
	.quad	.Ltmp10-.Lfunc_begin0
	.quad	.Ltmp17-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc4:
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Ltmp31-.Lfunc_begin0
	.short	8                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	98                      # DW_OP_reg18
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp31-.Lfunc_begin0
	.quad	.Ltmp32-.Lfunc_begin0
	.short	5                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	16                      # 16
	.byte	98                      # DW_OP_reg18
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp36-.Lfunc_begin0
	.quad	.Ltmp37-.Lfunc_begin0
	.short	8                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	98                      # DW_OP_reg18
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp37-.Lfunc_begin0
	.quad	.Ltmp39-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp39-.Lfunc_begin0
	.quad	.Ltmp40-.Lfunc_begin0
	.short	8                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	98                      # DW_OP_reg18
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp40-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	0
	.quad	0
.Ldebug_loc5:
	.quad	.Ltmp16-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp42-.Lfunc_begin0
	.quad	.Lfunc_end5-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	0
	.quad	0
.Ldebug_loc6:
	.quad	.Ltmp28-.Lfunc_begin0
	.quad	.Ltmp36-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp39-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp42-.Lfunc_begin0
	.quad	.Lfunc_end5-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	0
	.quad	0
.Ldebug_loc7:
	.quad	.Ltmp33-.Lfunc_begin0
	.quad	.Ltmp36-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp42-.Lfunc_begin0
	.quad	.Lfunc_end5-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc8:
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp42-.Lfunc_begin0
	.quad	.Lfunc_end5-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
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
	.byte	2                       # DW_AT_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	3                       # Abbreviation Code
	.byte	1                       # DW_TAG_array_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	4                       # Abbreviation Code
	.byte	33                      # DW_TAG_subrange_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	55                      # DW_AT_count
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	5                       # Abbreviation Code
	.byte	15                      # DW_TAG_pointer_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	6                       # Abbreviation Code
	.byte	21                      # DW_TAG_subroutine_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	7                       # Abbreviation Code
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
	.byte	8                       # Abbreviation Code
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
	.byte	9                       # Abbreviation Code
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
	.byte	10                      # Abbreviation Code
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
	.byte	11                      # Abbreviation Code
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
	.byte	12                      # Abbreviation Code
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
	.byte	13                      # Abbreviation Code
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
	.byte	14                      # Abbreviation Code
	.byte	59                      # DW_TAG_unspecified_type
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	15                      # Abbreviation Code
	.byte	19                      # DW_TAG_structure_type
	.byte	0                       # DW_CHILDREN_no
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
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	17                      # Abbreviation Code
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
	.byte	18                      # Abbreviation Code
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
	.byte	19                      # Abbreviation Code
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
	.byte	20                      # Abbreviation Code
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	21                      # Abbreviation Code
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
	.byte	22                      # Abbreviation Code
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
	.byte	23                      # Abbreviation Code
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
	.byte	24                      # Abbreviation Code
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
	.byte	25                      # Abbreviation Code
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
	.byte	26                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	10                      # DW_FORM_block1
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	27                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	10                      # DW_FORM_block1
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	28                      # Abbreviation Code
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
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	63                      # DW_AT_external
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	29                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	6                       # DW_FORM_data4
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	30                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	31                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	2                       # DW_AT_location
	.byte	6                       # DW_FORM_data4
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	32                      # Abbreviation Code
	.byte	4                       # DW_TAG_enumeration_type
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
	.byte	33                      # Abbreviation Code
	.byte	40                      # DW_TAG_enumerator
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	28                      # DW_AT_const_value
	.byte	13                      # DW_FORM_sdata
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
	.byte	1                       # Abbrev [1] 0xb:0x548 DW_TAG_compile_unit
	.long	.Linfo_string0          # DW_AT_producer
	.short	4                       # DW_AT_language
	.long	.Linfo_string1          # DW_AT_name
	.long	.Lline_table_start0     # DW_AT_stmt_list
	.long	.Linfo_string2          # DW_AT_comp_dir
	.byte	1                       # DW_AT_GNU_pubnames
	.quad	.Lfunc_begin0           # DW_AT_low_pc
	.quad	.Lfunc_end8             # DW_AT_high_pc
	.byte	2                       # Abbrev [2] 0x2f:0x14 DW_TAG_variable
	.long	.Linfo_string3          # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV11T3DFunction
	.byte	3                       # Abbrev [3] 0x43:0xc DW_TAG_array_type
	.long	79                      # DW_AT_type
	.byte	4                       # Abbrev [4] 0x48:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	5                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x4f:0x5 DW_TAG_pointer_type
	.long	84                      # DW_AT_type
	.byte	6                       # Abbrev [6] 0x54:0x5 DW_TAG_subroutine_type
	.long	89                      # DW_AT_type
	.byte	5                       # Abbrev [5] 0x59:0x5 DW_TAG_pointer_type
	.long	94                      # DW_AT_type
	.byte	6                       # Abbrev [6] 0x5e:0x5 DW_TAG_subroutine_type
	.long	99                      # DW_AT_type
	.byte	7                       # Abbrev [7] 0x63:0x7 DW_TAG_base_type
	.long	.Linfo_string4          # DW_AT_name
	.byte	5                       # DW_AT_encoding
	.byte	4                       # DW_AT_byte_size
	.byte	8                       # Abbrev [8] 0x6a:0x7 DW_TAG_base_type
	.long	.Linfo_string5          # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	7                       # DW_AT_encoding
	.byte	2                       # Abbrev [2] 0x71:0x14 DW_TAG_variable
	.long	.Linfo_string6          # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV13FieldFunction
	.byte	2                       # Abbrev [2] 0x85:0x14 DW_TAG_variable
	.long	.Linfo_string7          # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV10LineDipole
	.byte	9                       # Abbrev [9] 0x99:0x17 DW_TAG_variable
	.long	.Linfo_string8          # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	1                       # DW_AT_decl_file
	.short	8110                    # DW_AT_decl_line
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	__I___30_backgroundfield_linedipole_cpp_4a71163a
	.byte	10                      # Abbrev [10] 0xb0:0x13 DW_TAG_variable
	.long	.Linfo_string9          # DW_AT_name
	.long	195                     # DW_AT_type
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE
	.byte	11                      # Abbrev [11] 0xc3:0x11 DW_TAG_structure_type
	.long	.Linfo_string11         # DW_AT_name
	.byte	1                       # DW_AT_byte_size
	.byte	1                       # DW_AT_alignment
	.byte	12                      # Abbrev [12] 0xca:0x9 DW_TAG_member
	.long	212                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0xd4:0xc DW_TAG_array_type
	.long	224                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0xd9:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	0                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0xe0:0x7 DW_TAG_base_type
	.long	.Linfo_string10         # DW_AT_name
	.byte	6                       # DW_AT_encoding
	.byte	1                       # DW_AT_byte_size
	.byte	13                      # Abbrev [13] 0xe7:0xa DW_TAG_variable
	.long	.Linfo_string12         # DW_AT_name
	.long	241                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	5                       # Abbrev [5] 0xf1:0x5 DW_TAG_pointer_type
	.long	246                     # DW_AT_type
	.byte	14                      # Abbrev [14] 0xf6:0x5 DW_TAG_unspecified_type
	.long	.Linfo_string13         # DW_AT_name
	.byte	2                       # Abbrev [2] 0xfb:0x14 DW_TAG_variable
	.long	.Linfo_string14         # DW_AT_name
	.long	271                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI11T3DFunction
	.byte	15                      # Abbrev [15] 0x10f:0xa DW_TAG_structure_type
	.long	.Linfo_string15         # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	2                       # Abbrev [2] 0x119:0x14 DW_TAG_variable
	.long	.Linfo_string16         # DW_AT_name
	.long	301                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI13FieldFunction
	.byte	15                      # Abbrev [15] 0x12d:0xa DW_TAG_structure_type
	.long	.Linfo_string17         # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	16                      # Abbrev [16] 0x137:0x9 DW_TAG_variable
	.long	.Linfo_string18         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	2                       # Abbrev [2] 0x140:0x14 DW_TAG_variable
	.long	.Linfo_string19         # DW_AT_name
	.long	301                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI10LineDipole
	.byte	16                      # Abbrev [16] 0x154:0x9 DW_TAG_variable
	.long	.Linfo_string20         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	16                      # Abbrev [16] 0x15d:0x9 DW_TAG_variable
	.long	.Linfo_string21         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	13                      # Abbrev [13] 0x166:0xa DW_TAG_variable
	.long	.Linfo_string22         # DW_AT_name
	.long	368                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	3                       # Abbrev [3] 0x170:0xc DW_TAG_array_type
	.long	79                      # DW_AT_type
	.byte	4                       # Abbrev [4] 0x175:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	0                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	2                       # Abbrev [2] 0x17c:0x14 DW_TAG_variable
	.long	.Linfo_string23         # DW_AT_name
	.long	400                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS11T3DFunction
	.byte	3                       # Abbrev [3] 0x190:0xc DW_TAG_array_type
	.long	224                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x195:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	14                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	13                      # Abbrev [13] 0x19c:0xa DW_TAG_variable
	.long	.Linfo_string24         # DW_AT_name
	.long	368                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	2                       # Abbrev [2] 0x1a6:0x14 DW_TAG_variable
	.long	.Linfo_string25         # DW_AT_name
	.long	442                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS13FieldFunction
	.byte	3                       # Abbrev [3] 0x1ba:0xc DW_TAG_array_type
	.long	224                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x1bf:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	16                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	2                       # Abbrev [2] 0x1c6:0x14 DW_TAG_variable
	.long	.Linfo_string26         # DW_AT_name
	.long	474                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS10LineDipole
	.byte	3                       # Abbrev [3] 0x1da:0xc DW_TAG_array_type
	.long	224                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x1df:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	13                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	16                      # Abbrev [16] 0x1e6:0x9 DW_TAG_variable
	.long	.Linfo_string27         # DW_AT_name
	.long	495                     # DW_AT_type
	.byte	7                       # Abbrev [7] 0x1ef:0x7 DW_TAG_base_type
	.long	.Linfo_string28         # DW_AT_name
	.byte	7                       # DW_AT_encoding
	.byte	4                       # DW_AT_byte_size
	.byte	16                      # Abbrev [16] 0x1f6:0x9 DW_TAG_variable
	.long	.Linfo_string29         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	7                       # Abbrev [7] 0x1ff:0x7 DW_TAG_base_type
	.long	.Linfo_string30         # DW_AT_name
	.byte	4                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	17                      # Abbrev [17] 0x206:0x48 DW_TAG_class_type
	.long	.Linfo_string31         # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	18                      # Abbrev [18] 0x210:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin0           # DW_AT_low_pc
	.quad	.Lfunc_end0             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string44         # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	31                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x22a:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1347                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x232:0x1a DW_TAG_lexical_block
	.quad	.Ltmp0                  # DW_AT_low_pc
	.quad	.Ltmp1                  # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x243:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1347                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	22                      # Abbrev [22] 0x24e:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin1           # DW_AT_low_pc
	.quad	.Lfunc_end1             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string45         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	23                      # Abbrev [23] 0x266:0x3b DW_TAG_subprogram
	.quad	.Lfunc_begin2           # DW_AT_low_pc
	.quad	.Lfunc_end2             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string46         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x27e:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1352                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x286:0x1a DW_TAG_lexical_block
	.quad	.Ltmp4                  # DW_AT_low_pc
	.quad	.Ltmp5                  # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x297:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1352                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	22                      # Abbrev [22] 0x2a1:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin3           # DW_AT_low_pc
	.quad	.Lfunc_end3             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string47         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	17                      # Abbrev [17] 0x2b9:0x1ea DW_TAG_class_type
	.long	.Linfo_string43         # DW_AT_name
	.byte	72                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	24                      # Abbrev [24] 0x2c3:0xf DW_TAG_inheritance
	.long	.Linfo_string32         # DW_AT_name
	.long	1187                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	25                      # Abbrev [25] 0x2d2:0xf DW_TAG_member
	.long	.Linfo_string40         # DW_AT_name
	.long	224                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	20
	.byte	25                      # Abbrev [25] 0x2e1:0xf DW_TAG_member
	.long	.Linfo_string41         # DW_AT_name
	.long	1287                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	25                      # Abbrev [25] 0x2f0:0xf DW_TAG_member
	.long	.Linfo_string42         # DW_AT_name
	.long	1287                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	18                      # Abbrev [18] 0x2ff:0x95 DW_TAG_subprogram
	.quad	.Lfunc_begin4           # DW_AT_low_pc
	.quad	.Lfunc_end4             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string48         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	32                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x319:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1357                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	26                      # Abbrev [26] 0x321:0xb DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	97
	.long	.Linfo_string53         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	26                      # Abbrev [26] 0x32c:0xb DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	98
	.long	.Linfo_string54         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	26                      # Abbrev [26] 0x337:0xb DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	99
	.long	.Linfo_string55         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	26                      # Abbrev [26] 0x342:0xb DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	100
	.long	.Linfo_string56         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	20                      # Abbrev [20] 0x34d:0x46 DW_TAG_lexical_block
	.quad	.Ltmp8                  # DW_AT_low_pc
	.quad	.Ltmp9                  # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x35e:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1357                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0x366:0xb DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	97
	.long	.Linfo_string53         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x371:0xb DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	98
	.long	.Linfo_string54         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x37c:0xb DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	99
	.long	.Linfo_string55         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x387:0xb DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	100
	.long	.Linfo_string56         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	28                      # Abbrev [28] 0x394:0xd1 DW_TAG_subprogram
	.quad	.Lfunc_begin5           # DW_AT_low_pc
	.quad	.Lfunc_end5             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string49         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	45                      # DW_AT_decl_line
	.long	511                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x3b2:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1357                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0x3ba:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc0            # DW_AT_location
	.long	.Linfo_string57         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x3c7:0x9 DW_TAG_formal_parameter
	.long	.Linfo_string59         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	29                      # Abbrev [29] 0x3d0:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc1            # DW_AT_location
	.long	.Linfo_string58         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	20                      # Abbrev [20] 0x3dd:0x87 DW_TAG_lexical_block
	.quad	.Ltmp12                 # DW_AT_low_pc
	.quad	.Ltmp43                 # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x3ee:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1357                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	31                      # Abbrev [31] 0x3f6:0xd DW_TAG_variable
	.long	.Ldebug_loc2            # DW_AT_location
	.long	.Linfo_string57         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x403:0xd DW_TAG_variable
	.long	.Ldebug_loc3            # DW_AT_location
	.long	.Linfo_string58         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	16                      # Abbrev [16] 0x410:0x9 DW_TAG_variable
	.long	.Linfo_string59         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x419:0xd DW_TAG_variable
	.long	.Ldebug_loc4            # DW_AT_location
	.long	.Linfo_string60         # DW_AT_name
	.long	1287                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x426:0xd DW_TAG_variable
	.long	.Ldebug_loc5            # DW_AT_location
	.long	.Linfo_string61         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x433:0xd DW_TAG_variable
	.long	.Ldebug_loc6            # DW_AT_location
	.long	.Linfo_string62         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	16                      # Abbrev [16] 0x440:0x9 DW_TAG_variable
	.long	.Linfo_string63         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x449:0xd DW_TAG_variable
	.long	.Ldebug_loc7            # DW_AT_location
	.long	.Linfo_string64         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x456:0xd DW_TAG_variable
	.long	.Ldebug_loc8            # DW_AT_location
	.long	.Linfo_string65         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	18                      # Abbrev [18] 0x465:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin6           # DW_AT_low_pc
	.quad	.Lfunc_end6             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string50         # DW_AT_name
	.byte	4                       # DW_AT_decl_file
	.byte	47                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x47f:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1357                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x487:0x1a DW_TAG_lexical_block
	.quad	.Ltmp44                 # DW_AT_low_pc
	.quad	.Ltmp45                 # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x498:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1357                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	17                      # Abbrev [17] 0x4a3:0x47 DW_TAG_class_type
	.long	.Linfo_string32         # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	24                      # Abbrev [24] 0x4ad:0xf DW_TAG_inheritance
	.long	.Linfo_string31         # DW_AT_name
	.long	518                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	25                      # Abbrev [25] 0x4bc:0xf DW_TAG_member
	.long	.Linfo_string33         # DW_AT_name
	.long	1258                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	25                      # Abbrev [25] 0x4cb:0xf DW_TAG_member
	.long	.Linfo_string38         # DW_AT_name
	.long	1258                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	12
	.byte	25                      # Abbrev [25] 0x4da:0xf DW_TAG_member
	.long	.Linfo_string39         # DW_AT_name
	.long	495                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	0                       # End Of Children Mark
	.byte	32                      # Abbrev [32] 0x4ea:0x1d DW_TAG_enumeration_type
	.long	.Linfo_string37         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	4                       # DW_AT_alignment
	.byte	33                      # Abbrev [33] 0x4f4:0x6 DW_TAG_enumerator
	.long	.Linfo_string34         # DW_AT_name
	.byte	0                       # DW_AT_const_value
	.byte	33                      # Abbrev [33] 0x4fa:0x6 DW_TAG_enumerator
	.long	.Linfo_string35         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	33                      # Abbrev [33] 0x500:0x6 DW_TAG_enumerator
	.long	.Linfo_string36         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x507:0xc DW_TAG_array_type
	.long	511                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x50c:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	3                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	22                      # Abbrev [22] 0x513:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin7           # DW_AT_low_pc
	.quad	.Lfunc_end7             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string51         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	22                      # Abbrev [22] 0x52b:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin8           # DW_AT_low_pc
	.quad	.Lfunc_end8             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string52         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	5                       # Abbrev [5] 0x543:0x5 DW_TAG_pointer_type
	.long	518                     # DW_AT_type
	.byte	5                       # Abbrev [5] 0x548:0x5 DW_TAG_pointer_type
	.long	1187                    # DW_AT_type
	.byte	5                       # Abbrev [5] 0x54d:0x5 DW_TAG_pointer_type
	.long	697                     # DW_AT_type
	.byte	0                       # End Of Children Mark
.Ldebug_info_end0:
	.section	.debug_pubnames,"",@progbits
	.long	.LpubNames_end0-.LpubNames_begin0 # Length of Public Names Info
.LpubNames_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	1363                    # Compilation Unit Length
	.long	311                     # DIE offset
	.asciz	"WID"                   # External Name
	.long	251                     # DIE offset
	.asciz	"_ZTI11T3DFunction"     # External Name
	.long	590                     # DIE offset
	.asciz	"_ZN11T3DFunctionD0Ev"  # External Name
	.long	1323                    # DIE offset
	.asciz	"__sti___30_backgroundfield_linedipole_cpp_4a71163a" # External Name
	.long	231                     # DIE offset
	.asciz	"__dso_handle"          # External Name
	.long	113                     # DIE offset
	.asciz	"_ZTV13FieldFunction"   # External Name
	.long	380                     # DIE offset
	.asciz	"_ZTS11T3DFunction"     # External Name
	.long	528                     # DIE offset
	.asciz	"T3DFunction::~T3DFunction" # External Name
	.long	916                     # DIE offset
	.asciz	"LineDipole::call"      # External Name
	.long	412                     # DIE offset
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # External Name
	.long	1268                    # DIE offset
	.asciz	"X"                     # External Name
	.long	1280                    # DIE offset
	.asciz	"Z"                     # External Name
	.long	320                     # DIE offset
	.asciz	"_ZTI10LineDipole"      # External Name
	.long	1125                    # DIE offset
	.asciz	"LineDipole::~LineDipole" # External Name
	.long	1274                    # DIE offset
	.asciz	"Y"                     # External Name
	.long	358                     # DIE offset
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # External Name
	.long	281                     # DIE offset
	.asciz	"_ZTI13FieldFunction"   # External Name
	.long	673                     # DIE offset
	.asciz	"_ZN13FieldFunctionD0Ev" # External Name
	.long	454                     # DIE offset
	.asciz	"_ZTS10LineDipole"      # External Name
	.long	153                     # DIE offset
	.asciz	"__I___30_backgroundfield_linedipole_cpp_4a71163a" # External Name
	.long	614                     # DIE offset
	.asciz	"~FieldFunction"        # External Name
	.long	133                     # DIE offset
	.asciz	"_ZTV10LineDipole"      # External Name
	.long	767                     # DIE offset
	.asciz	"LineDipole::initialize" # External Name
	.long	486                     # DIE offset
	.asciz	"_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163a5vmesh15INVALID_LOCALIDE" # External Name
	.long	422                     # DIE offset
	.asciz	"_ZTS13FieldFunction"   # External Name
	.long	340                     # DIE offset
	.asciz	"WID2"                  # External Name
	.long	47                      # DIE offset
	.asciz	"_ZTV11T3DFunction"     # External Name
	.long	349                     # DIE offset
	.asciz	"WID3"                  # External Name
	.long	176                     # DIE offset
	.asciz	"_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163aSt8__ioinitE" # External Name
	.long	502                     # DIE offset
	.asciz	"_ZN52_INTERNAL_30_backgroundfield_linedipole_cpp_4a71163a17physicalconstants3R_EE" # External Name
	.long	1299                    # DIE offset
	.asciz	"_ZN10LineDipoleD0Ev"   # External Name
	.long	0                       # End Mark
.LpubNames_end0:
	.section	.debug_pubtypes,"",@progbits
	.long	.LpubTypes_end0-.LpubTypes_begin0 # Length of Public Types Info
.LpubTypes_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	1363                    # Compilation Unit Length
	.long	271                     # DIE offset
	.asciz	"__class_type_info"     # External Name
	.long	518                     # DIE offset
	.asciz	"T3DFunction"           # External Name
	.long	246                     # DIE offset
	.asciz	"void"                  # External Name
	.long	697                     # DIE offset
	.asciz	"LineDipole"            # External Name
	.long	1258                    # DIE offset
	.asciz	"coordinate"            # External Name
	.long	195                     # DIE offset
	.asciz	"_ZNSt8ios_base4InitE"  # External Name
	.long	224                     # DIE offset
	.asciz	"signed char"           # External Name
	.long	1187                    # DIE offset
	.asciz	"FieldFunction"         # External Name
	.long	99                      # DIE offset
	.asciz	"int"                   # External Name
	.long	301                     # DIE offset
	.asciz	"__si_class_type_info"  # External Name
	.long	495                     # DIE offset
	.asciz	"unsigned"              # External Name
	.long	511                     # DIE offset
	.asciz	"double"                # External Name
	.long	0                       # End Mark
.LpubTypes_end0:
	.section	".note.GNU-stack","",@progbits
	.section	.debug_line,"",@progbits
.Lline_table_start0:
