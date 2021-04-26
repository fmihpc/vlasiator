	.text
	.file	"dipole.ll"
	.file	1 "/home/talgat/vlasiator/backgroundfield/dipole.cpp"
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
	.loc	1 0 0                   # backgroundfield/dipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp2:
	.loc	1 31 1 prologue_end     # backgroundfield/dipole.cpp:31:1
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
	.loc	1 0 0                   # backgroundfield/dipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp6:
	.loc	1 31 1 prologue_end     # backgroundfield/dipole.cpp:31:1
	movl	$24, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp7:
.Lfunc_end3:
	.size	_ZN13FieldFunctionD0Ev, .Lfunc_end3-_ZN13FieldFunctionD0Ev
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4               # -- Begin function _ZN6Dipole10initializeEddddd
.LCPI4_0:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
	.text
	.globl	_ZN6Dipole10initializeEddddd
	.p2align	4, 0x90
	.type	_ZN6Dipole10initializeEddddd,@function
_ZN6Dipole10initializeEddddd:           # @_ZN6Dipole10initializeEddddd
.Lfunc_begin4:
	.loc	1 32 0                  # backgroundfield/dipole.cpp:32:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: initialize: <- $rdi
	#DEBUG_VALUE: initialize:moment <- $xmm0
	#DEBUG_VALUE: initialize:center_x <- $xmm1
	#DEBUG_VALUE: initialize:center_y <- $xmm2
	#DEBUG_VALUE: initialize:center_z <- $xmm3
	#DEBUG_VALUE: initialize:tilt_angle <- $xmm4
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%rbx
	subq	$40, %rsp
	.cfi_offset %rbx, -24
.Ltmp8:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: moment <- $xmm0
	#DEBUG_VALUE: center_x <- $xmm1
	#DEBUG_VALUE: center_y <- $xmm2
	#DEBUG_VALUE: center_z <- $xmm3
	#DEBUG_VALUE: tilt_angle <- $xmm4
	vmovsd	%xmm3, -40(%rbp)        # 8-byte Spill
	vmovsd	%xmm2, -32(%rbp)        # 8-byte Spill
	vmovsd	%xmm1, -24(%rbp)        # 8-byte Spill
	vmovsd	%xmm0, -16(%rbp)        # 8-byte Spill
	movq	%rdi, %rbx
	#DEBUG_VALUE: tilt_angle <- $xmm4
	#DEBUG_VALUE: initialize:tilt_angle <- $xmm4
.Ltmp9:
	#DEBUG_VALUE: center_z <- [DW_OP_constu 40, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:center_z <- [DW_OP_constu 40, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: center_y <- [DW_OP_constu 32, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:center_y <- [DW_OP_constu 32, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: center_x <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:center_x <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: moment <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:moment <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: initialize: <- $rbx
	.loc	1 33 1 prologue_end     # backgroundfield/dipole.cpp:33:1
	movb	$1, 20(%rdi)
	.loc	1 34 1                  # backgroundfield/dipole.cpp:34:1
	vmovapd	%xmm4, %xmm0
.Ltmp10:
	#DEBUG_VALUE: tilt_angle <- $xmm0
	#DEBUG_VALUE: initialize:tilt_angle <- $xmm0
	callq	__fd_sincos_1
.Ltmp11:
	.loc	1 0 1 is_stmt 0         # backgroundfield/dipole.cpp:0:1
	vmovsd	-16(%rbp), %xmm3        # 8-byte Reload
                                        # xmm3 = mem[0],zero
.Ltmp12:
	#DEBUG_VALUE: moment <- $xmm3
	#DEBUG_VALUE: initialize:moment <- $xmm3
	.loc	1 34 1                  # backgroundfield/dipole.cpp:34:1
	vmulsd	%xmm3, %xmm0, %xmm0
	vmovapd	.LCPI4_0(%rip), %xmm2   # xmm2 = [-0.0E+0,-0.0E+0]
	vxorpd	%xmm2, %xmm0, %xmm0
	vmovlpd	%xmm0, 24(%rbx)
	.loc	1 35 1 is_stmt 1        # backgroundfield/dipole.cpp:35:1
	movq	$0, 32(%rbx)
	.loc	1 36 1                  # backgroundfield/dipole.cpp:36:1
	vmulsd	%xmm3, %xmm1, %xmm0
	vxorpd	%xmm2, %xmm0, %xmm0
	vmovlpd	%xmm0, 40(%rbx)
	vmovsd	-24(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp13:
	#DEBUG_VALUE: center_x <- $xmm0
	#DEBUG_VALUE: initialize:center_x <- $xmm0
	.loc	1 37 1                  # backgroundfield/dipole.cpp:37:1
	vmovsd	%xmm0, 48(%rbx)
	vmovsd	-32(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp14:
	#DEBUG_VALUE: center_y <- $xmm0
	#DEBUG_VALUE: initialize:center_y <- $xmm0
	.loc	1 38 1                  # backgroundfield/dipole.cpp:38:1
	vmovsd	%xmm0, 56(%rbx)
	vmovsd	-40(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp15:
	#DEBUG_VALUE: center_z <- $xmm0
	#DEBUG_VALUE: initialize:center_z <- $xmm0
	.loc	1 39 1                  # backgroundfield/dipole.cpp:39:1
	vmovsd	%xmm0, 64(%rbx)
	.loc	1 40 1                  # backgroundfield/dipole.cpp:40:1
	addq	$40, %rsp
	popq	%rbx
.Ltmp16:
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp17:
.Lfunc_end4:
	.size	_ZN6Dipole10initializeEddddd, .Lfunc_end4-_ZN6Dipole10initializeEddddd
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _ZNK6Dipole4callEddd
.LCPI5_0:
	.quad	4720717001298091704     # double 40592189.439999998
.LCPI5_1:
	.quad	4613937818241073152     # double 3
.LCPI5_2:
	.quad	4617315517961601024     # double 5
.LCPI5_3:
	.quad	4607182418800017408     # double 1
	.text
	.globl	_ZNK6Dipole4callEddd
	.p2align	4, 0x90
	.type	_ZNK6Dipole4callEddd,@function
_ZNK6Dipole4callEddd:                   # @_ZNK6Dipole4callEddd
.Lfunc_begin5:
	.loc	1 45 0                  # backgroundfield/dipole.cpp:45:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:x <- $xmm0
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE: call:z <- $xmm2
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	vmovapd	%xmm0, %xmm3
.Ltmp18:
	#DEBUG_VALUE: z <- undef
	#DEBUG_VALUE: y <- undef
	#DEBUG_VALUE: x <- undef
	vxorpd	%xmm0, %xmm0, %xmm0
.Ltmp19:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE: call:x <- $xmm3
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 47 1 prologue_end     # backgroundfield/dipole.cpp:47:1
	cmpb	$0, 20(%rdi)
	je	.LBB5_4
.Ltmp20:
# %bb.1:                                # %L.B0000
	#DEBUG_VALUE: call:x <- $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 51 1                  # backgroundfield/dipole.cpp:51:1
	vunpcklpd	%xmm1, %xmm3, %xmm1 # xmm1 = xmm3[0],xmm1[0]
.Ltmp21:
	vsubpd	48(%rdi), %xmm1, %xmm3
.Ltmp22:
	vmovapd	%xmm3, -32(%rbp)
	.loc	1 53 1                  # backgroundfield/dipole.cpp:53:1
	vsubsd	64(%rdi), %xmm2, %xmm4
	vmovsd	%xmm4, -16(%rbp)
	.loc	1 55 1                  # backgroundfield/dipole.cpp:55:1
	vmulsd	%xmm3, %xmm3, %xmm1
	vpermilpd	$1, %xmm3, %xmm5 # xmm5 = xmm3[1,0]
	vfmadd231sd	%xmm5, %xmm5, %xmm1 # xmm1 = (xmm5 * xmm5) + xmm1
	vfmadd231sd	%xmm4, %xmm4, %xmm1 # xmm1 = (xmm4 * xmm4) + xmm1
.Ltmp23:
	#DEBUG_VALUE: r2 <- $xmm1
	.loc	1 0 1 is_stmt 0         # backgroundfield/dipole.cpp:0:1
	vmovsd	.LCPI5_0(%rip), %xmm2   # xmm2 = mem[0],zero
.Ltmp24:
	.loc	1 57 1 is_stmt 1        # backgroundfield/dipole.cpp:57:1
	vucomisd	%xmm1, %xmm2
	ja	.LBB5_4
.Ltmp25:
# %bb.2:                                # %L.B0001
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 61 1                  # backgroundfield/dipole.cpp:61:1
	vsqrtsd	%xmm1, %xmm1, %xmm2
	vmulsd	%xmm1, %xmm1, %xmm6
	vmulsd	%xmm6, %xmm2, %xmm2
.Ltmp26:
	#DEBUG_VALUE: r5 <- $xmm2
	.loc	1 62 1                  # backgroundfield/dipole.cpp:62:1
	vmulsd	24(%rdi), %xmm3, %xmm3
	vfmadd231sd	32(%rdi), %xmm5, %xmm3 # xmm3 = (xmm5 * mem) + xmm3
	vfmadd231sd	40(%rdi), %xmm4, %xmm3 # xmm3 = (xmm4 * mem) + xmm3
.Ltmp27:
	#DEBUG_VALUE: rdotq <- $xmm3
	.loc	1 64 1                  # backgroundfield/dipole.cpp:64:1
	movl	8(%rdi), %eax
	vmovsd	-32(%rbp,%rax,8), %xmm4 # xmm4 = mem[0],zero
	vmulsd	.LCPI5_1(%rip), %xmm4, %xmm6
	.loc	1 66 1                  # backgroundfield/dipole.cpp:66:1
	movl	16(%rdi), %ecx
	.loc	1 64 1                  # backgroundfield/dipole.cpp:64:1
	vmovsd	24(%rdi,%rax,8), %xmm5  # xmm5 = mem[0],zero
	vmulsd	%xmm5, %xmm1, %xmm7
	vfmsub231sd	%xmm6, %xmm3, %xmm7 # xmm7 = (xmm3 * xmm6) - xmm7
	vdivsd	%xmm2, %xmm7, %xmm6
.Ltmp28:
	#DEBUG_VALUE: B <- $xmm6
	.loc	1 66 1                  # backgroundfield/dipole.cpp:66:1
	testl	%ecx, %ecx
	je	.LBB5_3
.Ltmp29:
# %bb.5:                                # %L.B0001
	#DEBUG_VALUE: B <- $xmm6
	#DEBUG_VALUE: rdotq <- $xmm3
	#DEBUG_VALUE: r5 <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	cmpl	$1, %ecx
	jne	.LBB5_4
.Ltmp30:
# %bb.6:                                # %L.B0044
	#DEBUG_VALUE: B <- $xmm6
	#DEBUG_VALUE: rdotq <- $xmm3
	#DEBUG_VALUE: r5 <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 78 1                  # backgroundfield/dipole.cpp:78:1
	vmulsd	.LCPI5_2(%rip), %xmm6, %xmm0
	movl	12(%rdi), %ecx
	vmovsd	-32(%rbp,%rcx,8), %xmm6 # xmm6 = mem[0],zero
.Ltmp31:
	.loc	1 0 1 is_stmt 0         # backgroundfield/dipole.cpp:0:1
	vmovsd	.LCPI5_1(%rip), %xmm8   # xmm8 = mem[0],zero
	.loc	1 78 1                  # backgroundfield/dipole.cpp:78:1
	vmulsd	24(%rdi,%rcx,8), %xmm8, %xmm7
	vmulsd	%xmm6, %xmm0, %xmm0
	vaddsd	%xmm5, %xmm5, %xmm5
	vmulsd	%xmm6, %xmm5, %xmm5
	vfmsub231sd	%xmm7, %xmm4, %xmm5 # xmm5 = (xmm4 * xmm7) - xmm5
	vmulsd	%xmm8, %xmm3, %xmm3
.Ltmp32:
	cmpl	%eax, %ecx
	je	.LBB5_7
.Ltmp33:
# %bb.8:                                # %L.B0044
	#DEBUG_VALUE: r5 <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1                   # backgroundfield/dipole.cpp:0:1
	vxorpd	%xmm4, %xmm4, %xmm4
	jmp	.LBB5_9
.Ltmp34:
.LBB5_3:
	#DEBUG_VALUE: B <- $xmm6
	#DEBUG_VALUE: rdotq <- $xmm3
	#DEBUG_VALUE: r5 <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	vmovapd	%xmm6, %xmm0
.Ltmp35:
.LBB5_4:                                # %L.B0038
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 48 1 is_stmt 1        # backgroundfield/dipole.cpp:48:1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp36:
.LBB5_7:
	.cfi_def_cfa %rbp, 16
	#DEBUG_VALUE: r5 <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/dipole.cpp:0:1
	vmovsd	.LCPI5_3(%rip), %xmm4   # xmm4 = mem[0],zero
.Ltmp37:
.LBB5_9:                                # %L.B0044
	#DEBUG_VALUE: r5 <- $xmm2
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 78 1 is_stmt 1        # backgroundfield/dipole.cpp:78:1
	vfmadd213sd	%xmm5, %xmm3, %xmm4 # xmm4 = (xmm3 * xmm4) + xmm5
	vunpcklpd	%xmm0, %xmm4, %xmm0 # xmm0 = xmm4[0],xmm0[0]
	vunpcklpd	%xmm1, %xmm2, %xmm1 # xmm1 = xmm2[0],xmm1[0]
.Ltmp38:
	vdivpd	%xmm1, %xmm0, %xmm0
	vpermilpd	$1, %xmm0, %xmm1 # xmm1 = xmm0[1,0]
	vsubsd	%xmm1, %xmm0, %xmm0
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp39:
.Lfunc_end5:
	.size	_ZNK6Dipole4callEddd, .Lfunc_end5-_ZNK6Dipole4callEddd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN6DipoleD1Ev          # -- Begin function _ZN6DipoleD1Ev
	.p2align	4, 0x90
	.type	_ZN6DipoleD1Ev,@function
_ZN6DipoleD1Ev:                         # @_ZN6DipoleD1Ev
.Lfunc_begin6:
	.file	4 "/home/talgat/vlasiator/backgroundfield/dipole.hpp"
	.loc	4 44 0                  # backgroundfield/dipole.hpp:44:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~Dipole: <- $rdi
	#DEBUG_VALUE: ~Dipole: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp40:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	4 44 1 prologue_end     # backgroundfield/dipole.hpp:44:1
	movq	$_ZTV11T3DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp41:
.Lfunc_end6:
	.size	_ZN6DipoleD1Ev, .Lfunc_end6-_ZN6DipoleD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN6DipoleD0Ev          # -- Begin function _ZN6DipoleD0Ev
	.p2align	4, 0x90
	.type	_ZN6DipoleD0Ev,@function
_ZN6DipoleD0Ev:                         # @_ZN6DipoleD0Ev
.Lfunc_begin7:
	.loc	1 0 0                   # backgroundfield/dipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp42:
	.loc	1 44 1 prologue_end     # backgroundfield/dipole.cpp:44:1
	movl	$72, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp43:
.Lfunc_end7:
	.size	_ZN6DipoleD0Ev, .Lfunc_end7-_ZN6DipoleD0Ev
	.cfi_endproc
                                        # -- End function
	.globl	__sti___26_backgroundfield_dipole_cpp_41c6b7aa # -- Begin function __sti___26_backgroundfield_dipole_cpp_41c6b7aa
	.p2align	4, 0x90
	.type	__sti___26_backgroundfield_dipole_cpp_41c6b7aa,@function
__sti___26_backgroundfield_dipole_cpp_41c6b7aa: # @__sti___26_backgroundfield_dipole_cpp_41c6b7aa
.Lfunc_begin8:
	.loc	1 0 0                   # backgroundfield/dipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp44:
	.loc	1 74 1 prologue_end     # backgroundfield/dipole.cpp:74:1
	cmpl	$1, __I___26_backgroundfield_dipole_cpp_41c6b7aa(%rip)
	jne	.LBB8_2
# %bb.1:                                # %L.B0011
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB8_2:                                # %L.B0052
	.cfi_def_cfa %rbp, 16
	movl	$1, __I___26_backgroundfield_dipole_cpp_41c6b7aa(%rip)
	movl	$_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE, %edi
	callq	_ZNSt8ios_base4InitC1Ev
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	movl	$_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE, %esi
	movl	$__dso_handle, %edx
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	__cxa_atexit            # TAILCALL
.Ltmp45:
.Lfunc_end8:
	.size	__sti___26_backgroundfield_dipole_cpp_41c6b7aa, .Lfunc_end8-__sti___26_backgroundfield_dipole_cpp_41c6b7aa
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

	.type	_ZTV6Dipole,@object     # @_ZTV6Dipole
	.weak	_ZTV6Dipole
	.p2align	4
_ZTV6Dipole:
	.quad	0
	.quad	_ZTI6Dipole
	.quad	_ZNK6Dipole4callEddd
	.quad	_ZN6DipoleD1Ev
	.quad	_ZN6DipoleD0Ev
	.size	_ZTV6Dipole, 40

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

	.type	_ZTI6Dipole,@object     # @_ZTI6Dipole
	.weak	_ZTI6Dipole
	.p2align	4
_ZTI6Dipole:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS6Dipole
	.quad	_ZTI13FieldFunction
	.size	_ZTI6Dipole, 24

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

	.type	_ZTS6Dipole,@object     # @_ZTS6Dipole
	.weak	_ZTS6Dipole
	.p2align	3
_ZTS6Dipole:
	.asciz	"6Dipole"
	.size	_ZTS6Dipole, 8

	.type	__I___26_backgroundfield_dipole_cpp_41c6b7aa,@object # @__I___26_backgroundfield_dipole_cpp_41c6b7aa
	.bss
	.globl	__I___26_backgroundfield_dipole_cpp_41c6b7aa
	.p2align	2
__I___26_backgroundfield_dipole_cpp_41c6b7aa:
	.long	0                       # 0x0
	.size	__I___26_backgroundfield_dipole_cpp_41c6b7aa, 4

	.type	_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE,@object # @_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE
	.local	_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE
	.comm	_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE,1,1
	.section	.init_array,"aw",@init_array
	.p2align	3
	.quad	__sti___26_backgroundfield_dipole_cpp_41c6b7aa
	.section	.debug_str,"MS",@progbits,1
.Linfo_string0:
	.asciz	" NVC++ 21.2-0"         # string offset=0
.Linfo_string1:
	.asciz	"backgroundfield/dipole.cpp" # string offset=14
.Linfo_string2:
	.asciz	"/home/talgat/vlasiator" # string offset=41
.Linfo_string3:
	.asciz	"_ZTV11T3DFunction"     # string offset=64
.Linfo_string4:
	.asciz	"int"                   # string offset=82
.Linfo_string5:
	.asciz	"__ARRAY_SIZE_TYPE__"   # string offset=86
.Linfo_string6:
	.asciz	"_ZTV13FieldFunction"   # string offset=106
.Linfo_string7:
	.asciz	"_ZTV6Dipole"           # string offset=126
.Linfo_string8:
	.asciz	"__I___26_backgroundfield_dipole_cpp_41c6b7aa" # string offset=138
.Linfo_string9:
	.asciz	"_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE" # string offset=183
.Linfo_string10:
	.asciz	"signed char"           # string offset=249
.Linfo_string11:
	.asciz	"_ZNSt8ios_base4InitE"  # string offset=261
.Linfo_string12:
	.asciz	"__dso_handle"          # string offset=282
.Linfo_string13:
	.asciz	"void"                  # string offset=295
.Linfo_string14:
	.asciz	"_ZTI11T3DFunction"     # string offset=300
.Linfo_string15:
	.asciz	"__class_type_info"     # string offset=318
.Linfo_string16:
	.asciz	"_ZTI13FieldFunction"   # string offset=336
.Linfo_string17:
	.asciz	"__si_class_type_info"  # string offset=356
.Linfo_string18:
	.asciz	"WID"                   # string offset=377
.Linfo_string19:
	.asciz	"_ZTI6Dipole"           # string offset=381
.Linfo_string20:
	.asciz	"WID2"                  # string offset=393
.Linfo_string21:
	.asciz	"WID3"                  # string offset=398
.Linfo_string22:
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # string offset=403
.Linfo_string23:
	.asciz	"_ZTS11T3DFunction"     # string offset=441
.Linfo_string24:
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # string offset=459
.Linfo_string25:
	.asciz	"_ZTS13FieldFunction"   # string offset=500
.Linfo_string26:
	.asciz	"_ZTS6Dipole"           # string offset=520
.Linfo_string27:
	.asciz	"_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aa5vmesh15INVALID_LOCALIDE" # string offset=532
.Linfo_string28:
	.asciz	"unsigned"              # string offset=610
.Linfo_string29:
	.asciz	"_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aa17physicalconstants3R_EE" # string offset=619
.Linfo_string30:
	.asciz	"double"                # string offset=697
.Linfo_string31:
	.asciz	"T3DFunction"           # string offset=704
.Linfo_string32:
	.asciz	"FieldFunction"         # string offset=716
.Linfo_string33:
	.asciz	"_fComponent"           # string offset=730
.Linfo_string34:
	.asciz	"X"                     # string offset=742
.Linfo_string35:
	.asciz	"Y"                     # string offset=744
.Linfo_string36:
	.asciz	"Z"                     # string offset=746
.Linfo_string37:
	.asciz	"coordinate"            # string offset=748
.Linfo_string38:
	.asciz	"_dComponent"           # string offset=759
.Linfo_string39:
	.asciz	"_derivative"           # string offset=771
.Linfo_string40:
	.asciz	"initialized"           # string offset=783
.Linfo_string41:
	.asciz	"q"                     # string offset=795
.Linfo_string42:
	.asciz	"center"                # string offset=797
.Linfo_string43:
	.asciz	"Dipole"                # string offset=804
.Linfo_string44:
	.asciz	"~T3DFunction"          # string offset=811
.Linfo_string45:
	.asciz	"_ZN11T3DFunctionD0Ev"  # string offset=824
.Linfo_string46:
	.asciz	"~FieldFunction"        # string offset=845
.Linfo_string47:
	.asciz	"_ZN13FieldFunctionD0Ev" # string offset=860
.Linfo_string48:
	.asciz	"initialize"            # string offset=883
.Linfo_string49:
	.asciz	"call"                  # string offset=894
.Linfo_string50:
	.asciz	"~Dipole"               # string offset=899
.Linfo_string51:
	.asciz	"_ZN6DipoleD0Ev"        # string offset=907
.Linfo_string52:
	.asciz	"__sti___26_backgroundfield_dipole_cpp_41c6b7aa" # string offset=922
.Linfo_string53:
	.asciz	"moment"                # string offset=969
.Linfo_string54:
	.asciz	"center_x"              # string offset=976
.Linfo_string55:
	.asciz	"center_y"              # string offset=985
.Linfo_string56:
	.asciz	"center_z"              # string offset=994
.Linfo_string57:
	.asciz	"tilt_angle"            # string offset=1003
.Linfo_string58:
	.asciz	"r"                     # string offset=1014
.Linfo_string59:
	.asciz	"x"                     # string offset=1016
.Linfo_string60:
	.asciz	"y"                     # string offset=1018
.Linfo_string61:
	.asciz	"z"                     # string offset=1020
.Linfo_string62:
	.asciz	"r2"                    # string offset=1022
.Linfo_string63:
	.asciz	"r5"                    # string offset=1025
.Linfo_string64:
	.asciz	"rdotq"                 # string offset=1028
.Linfo_string65:
	.asciz	"B"                     # string offset=1034
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp16-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc1:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	112                     # -16
	.quad	.Ltmp12-.Lfunc_begin0
	.quad	.Lfunc_end4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc2:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp13-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	104                     # -24
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc3:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	96                      # -32
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc4:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	88                      # -40
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Lfunc_end4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc5:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp10-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc6:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp16-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc7:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	112                     # -16
	.quad	.Ltmp12-.Lfunc_begin0
	.quad	.Lfunc_end4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc8:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp13-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	104                     # -24
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc9:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	96                      # -32
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc10:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	88                      # -40
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Lfunc_end4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc11:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp10-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc12:
	.quad	.Lfunc_begin5-.Lfunc_begin0
	.quad	.Ltmp19-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp19-.Lfunc_begin0
	.quad	.Ltmp22-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc13:
	.quad	.Lfunc_begin5-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc14:
	.quad	.Lfunc_begin5-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc15:
	.quad	.Ltmp23-.Lfunc_begin0
	.quad	.Ltmp35-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp36-.Lfunc_begin0
	.quad	.Ltmp38-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc16:
	.quad	.Ltmp26-.Lfunc_begin0
	.quad	.Ltmp35-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp36-.Lfunc_begin0
	.quad	.Lfunc_end5-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc17:
	.quad	.Ltmp27-.Lfunc_begin0
	.quad	.Ltmp32-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp34-.Lfunc_begin0
	.quad	.Ltmp35-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc18:
	.quad	.Ltmp28-.Lfunc_begin0
	.quad	.Ltmp31-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp34-.Lfunc_begin0
	.quad	.Ltmp35-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
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
	.byte	6                       # DW_FORM_data4
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	52                      # DW_AT_artificial
	.byte	12                      # DW_FORM_flag
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	27                      # Abbreviation Code
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
	.byte	28                      # Abbreviation Code
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
	.byte	29                      # Abbreviation Code
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
	.byte	30                      # Abbreviation Code
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
	.byte	31                      # Abbreviation Code
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
	.byte	1                       # Abbrev [1] 0xb:0x568 DW_TAG_compile_unit
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
	.quad	_ZTV6Dipole
	.byte	9                       # Abbrev [9] 0x99:0x17 DW_TAG_variable
	.long	.Linfo_string8          # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	1                       # DW_AT_decl_file
	.short	8115                    # DW_AT_decl_line
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	__I___26_backgroundfield_dipole_cpp_41c6b7aa
	.byte	10                      # Abbrev [10] 0xb0:0x13 DW_TAG_variable
	.long	.Linfo_string9          # DW_AT_name
	.long	195                     # DW_AT_type
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE
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
	.quad	_ZTI6Dipole
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
	.quad	_ZTS6Dipole
	.byte	3                       # Abbrev [3] 0x1da:0xc DW_TAG_array_type
	.long	224                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x1df:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	8                       # DW_AT_count
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
	.long	1379                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x232:0x1a DW_TAG_lexical_block
	.quad	.Ltmp0                  # DW_AT_low_pc
	.quad	.Ltmp1                  # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x243:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1379                    # DW_AT_type
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
	.long	1384                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x286:0x1a DW_TAG_lexical_block
	.quad	.Ltmp4                  # DW_AT_low_pc
	.quad	.Ltmp5                  # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x297:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1384                    # DW_AT_type
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
	.byte	17                      # Abbrev [17] 0x2b9:0x20a DW_TAG_class_type
	.long	.Linfo_string43         # DW_AT_name
	.byte	72                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	24                      # Abbrev [24] 0x2c3:0xf DW_TAG_inheritance
	.long	.Linfo_string32         # DW_AT_name
	.long	1219                    # DW_AT_type
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
	.long	1319                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	25                      # Abbrev [25] 0x2f0:0xf DW_TAG_member
	.long	.Linfo_string42         # DW_AT_name
	.long	1319                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	18                      # Abbrev [18] 0x2ff:0xc3 DW_TAG_subprogram
	.quad	.Lfunc_begin4           # DW_AT_low_pc
	.quad	.Lfunc_end4             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string48         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	32                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	26                      # Abbrev [26] 0x319:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc0            # DW_AT_location
	.long	1389                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0x323:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc1            # DW_AT_location
	.long	.Linfo_string53         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x330:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc2            # DW_AT_location
	.long	.Linfo_string54         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x33d:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc3            # DW_AT_location
	.long	.Linfo_string55         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x34a:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc4            # DW_AT_location
	.long	.Linfo_string56         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x357:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc5            # DW_AT_location
	.long	.Linfo_string57         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	20                      # Abbrev [20] 0x364:0x5d DW_TAG_lexical_block
	.quad	.Ltmp9                  # DW_AT_low_pc
	.quad	.Ltmp17                 # DW_AT_high_pc
	.byte	28                      # Abbrev [28] 0x375:0xa DW_TAG_variable
	.long	.Ldebug_loc6            # DW_AT_location
	.long	1389                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0x37f:0xd DW_TAG_variable
	.long	.Ldebug_loc7            # DW_AT_location
	.long	.Linfo_string53         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	29                      # Abbrev [29] 0x38c:0xd DW_TAG_variable
	.long	.Ldebug_loc8            # DW_AT_location
	.long	.Linfo_string54         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	29                      # Abbrev [29] 0x399:0xd DW_TAG_variable
	.long	.Ldebug_loc9            # DW_AT_location
	.long	.Linfo_string55         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	29                      # Abbrev [29] 0x3a6:0xd DW_TAG_variable
	.long	.Ldebug_loc10           # DW_AT_location
	.long	.Linfo_string56         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	29                      # Abbrev [29] 0x3b3:0xd DW_TAG_variable
	.long	.Ldebug_loc11           # DW_AT_location
	.long	.Linfo_string57         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	30                      # Abbrev [30] 0x3c2:0xc3 DW_TAG_subprogram
	.quad	.Lfunc_begin5           # DW_AT_low_pc
	.quad	.Lfunc_end5             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string49         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	45                      # DW_AT_decl_line
	.long	511                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x3e0:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1389                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0x3e8:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc12           # DW_AT_location
	.long	.Linfo_string59         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x3f5:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc13           # DW_AT_location
	.long	.Linfo_string60         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x402:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc14           # DW_AT_location
	.long	.Linfo_string61         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x40f:0xc DW_TAG_variable
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	96
	.long	.Linfo_string58         # DW_AT_name
	.long	1319                    # DW_AT_type
	.byte	20                      # Abbrev [20] 0x41b:0x69 DW_TAG_lexical_block
	.quad	.Ltmp19                 # DW_AT_low_pc
	.quad	.Ltmp39                 # DW_AT_high_pc
	.byte	16                      # Abbrev [16] 0x42c:0x9 DW_TAG_variable
	.long	.Linfo_string61         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	16                      # Abbrev [16] 0x435:0x9 DW_TAG_variable
	.long	.Linfo_string60         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	16                      # Abbrev [16] 0x43e:0x9 DW_TAG_variable
	.long	.Linfo_string59         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	21                      # Abbrev [21] 0x447:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1389                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0x44f:0xd DW_TAG_variable
	.long	.Ldebug_loc15           # DW_AT_location
	.long	.Linfo_string62         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	29                      # Abbrev [29] 0x45c:0xd DW_TAG_variable
	.long	.Ldebug_loc16           # DW_AT_location
	.long	.Linfo_string63         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	29                      # Abbrev [29] 0x469:0xd DW_TAG_variable
	.long	.Ldebug_loc17           # DW_AT_location
	.long	.Linfo_string64         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	29                      # Abbrev [29] 0x476:0xd DW_TAG_variable
	.long	.Ldebug_loc18           # DW_AT_location
	.long	.Linfo_string65         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	18                      # Abbrev [18] 0x485:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin6           # DW_AT_low_pc
	.quad	.Lfunc_end6             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string50         # DW_AT_name
	.byte	4                       # DW_AT_decl_file
	.byte	44                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x49f:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1389                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x4a7:0x1a DW_TAG_lexical_block
	.quad	.Ltmp40                 # DW_AT_low_pc
	.quad	.Ltmp41                 # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x4b8:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1389                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	17                      # Abbrev [17] 0x4c3:0x47 DW_TAG_class_type
	.long	.Linfo_string32         # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	24                      # Abbrev [24] 0x4cd:0xf DW_TAG_inheritance
	.long	.Linfo_string31         # DW_AT_name
	.long	518                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	25                      # Abbrev [25] 0x4dc:0xf DW_TAG_member
	.long	.Linfo_string33         # DW_AT_name
	.long	1290                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	25                      # Abbrev [25] 0x4eb:0xf DW_TAG_member
	.long	.Linfo_string38         # DW_AT_name
	.long	1290                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	12
	.byte	25                      # Abbrev [25] 0x4fa:0xf DW_TAG_member
	.long	.Linfo_string39         # DW_AT_name
	.long	495                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	0                       # End Of Children Mark
	.byte	32                      # Abbrev [32] 0x50a:0x1d DW_TAG_enumeration_type
	.long	.Linfo_string37         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	4                       # DW_AT_alignment
	.byte	33                      # Abbrev [33] 0x514:0x6 DW_TAG_enumerator
	.long	.Linfo_string34         # DW_AT_name
	.byte	0                       # DW_AT_const_value
	.byte	33                      # Abbrev [33] 0x51a:0x6 DW_TAG_enumerator
	.long	.Linfo_string35         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	33                      # Abbrev [33] 0x520:0x6 DW_TAG_enumerator
	.long	.Linfo_string36         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x527:0xc DW_TAG_array_type
	.long	511                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x52c:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	3                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	22                      # Abbrev [22] 0x533:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin7           # DW_AT_low_pc
	.quad	.Lfunc_end7             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string51         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	22                      # Abbrev [22] 0x54b:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin8           # DW_AT_low_pc
	.quad	.Lfunc_end8             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string52         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	5                       # Abbrev [5] 0x563:0x5 DW_TAG_pointer_type
	.long	518                     # DW_AT_type
	.byte	5                       # Abbrev [5] 0x568:0x5 DW_TAG_pointer_type
	.long	1219                    # DW_AT_type
	.byte	5                       # Abbrev [5] 0x56d:0x5 DW_TAG_pointer_type
	.long	697                     # DW_AT_type
	.byte	0                       # End Of Children Mark
.Ldebug_info_end0:
	.section	.debug_pubnames,"",@progbits
	.long	.LpubNames_end0-.LpubNames_begin0 # Length of Public Names Info
.LpubNames_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	1395                    # Compilation Unit Length
	.long	311                     # DIE offset
	.asciz	"WID"                   # External Name
	.long	251                     # DIE offset
	.asciz	"_ZTI11T3DFunction"     # External Name
	.long	320                     # DIE offset
	.asciz	"_ZTI6Dipole"           # External Name
	.long	590                     # DIE offset
	.asciz	"_ZN11T3DFunctionD0Ev"  # External Name
	.long	962                     # DIE offset
	.asciz	"Dipole::call"          # External Name
	.long	231                     # DIE offset
	.asciz	"__dso_handle"          # External Name
	.long	486                     # DIE offset
	.asciz	"_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aa5vmesh15INVALID_LOCALIDE" # External Name
	.long	113                     # DIE offset
	.asciz	"_ZTV13FieldFunction"   # External Name
	.long	1157                    # DIE offset
	.asciz	"Dipole::~Dipole"       # External Name
	.long	380                     # DIE offset
	.asciz	"_ZTS11T3DFunction"     # External Name
	.long	454                     # DIE offset
	.asciz	"_ZTS6Dipole"           # External Name
	.long	528                     # DIE offset
	.asciz	"T3DFunction::~T3DFunction" # External Name
	.long	1355                    # DIE offset
	.asciz	"__sti___26_backgroundfield_dipole_cpp_41c6b7aa" # External Name
	.long	412                     # DIE offset
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # External Name
	.long	1306                    # DIE offset
	.asciz	"Y"                     # External Name
	.long	1312                    # DIE offset
	.asciz	"Z"                     # External Name
	.long	1300                    # DIE offset
	.asciz	"X"                     # External Name
	.long	176                     # DIE offset
	.asciz	"_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aaSt8__ioinitE" # External Name
	.long	1331                    # DIE offset
	.asciz	"_ZN6DipoleD0Ev"        # External Name
	.long	502                     # DIE offset
	.asciz	"_ZN48_INTERNAL_26_backgroundfield_dipole_cpp_41c6b7aa17physicalconstants3R_EE" # External Name
	.long	358                     # DIE offset
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # External Name
	.long	281                     # DIE offset
	.asciz	"_ZTI13FieldFunction"   # External Name
	.long	153                     # DIE offset
	.asciz	"__I___26_backgroundfield_dipole_cpp_41c6b7aa" # External Name
	.long	767                     # DIE offset
	.asciz	"Dipole::initialize"    # External Name
	.long	673                     # DIE offset
	.asciz	"_ZN13FieldFunctionD0Ev" # External Name
	.long	614                     # DIE offset
	.asciz	"~FieldFunction"        # External Name
	.long	422                     # DIE offset
	.asciz	"_ZTS13FieldFunction"   # External Name
	.long	340                     # DIE offset
	.asciz	"WID2"                  # External Name
	.long	47                      # DIE offset
	.asciz	"_ZTV11T3DFunction"     # External Name
	.long	349                     # DIE offset
	.asciz	"WID3"                  # External Name
	.long	133                     # DIE offset
	.asciz	"_ZTV6Dipole"           # External Name
	.long	0                       # End Mark
.LpubNames_end0:
	.section	.debug_pubtypes,"",@progbits
	.long	.LpubTypes_end0-.LpubTypes_begin0 # Length of Public Types Info
.LpubTypes_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	1395                    # Compilation Unit Length
	.long	271                     # DIE offset
	.asciz	"__class_type_info"     # External Name
	.long	518                     # DIE offset
	.asciz	"T3DFunction"           # External Name
	.long	246                     # DIE offset
	.asciz	"void"                  # External Name
	.long	697                     # DIE offset
	.asciz	"Dipole"                # External Name
	.long	1290                    # DIE offset
	.asciz	"coordinate"            # External Name
	.long	195                     # DIE offset
	.asciz	"_ZNSt8ios_base4InitE"  # External Name
	.long	224                     # DIE offset
	.asciz	"signed char"           # External Name
	.long	1219                    # DIE offset
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
