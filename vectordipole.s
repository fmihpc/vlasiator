	.text
	.file	"vectordipole.ll"
	.file	1 "/home/talgat/vlasiator/backgroundfield/vectordipole.cpp"
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
	.loc	1 0 0                   # backgroundfield/vectordipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp2:
	.loc	1 31 1 prologue_end     # backgroundfield/vectordipole.cpp:31:1
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
	.loc	1 0 0                   # backgroundfield/vectordipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp6:
	.loc	1 31 1 prologue_end     # backgroundfield/vectordipole.cpp:31:1
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
	.p2align	4               # -- Begin function _ZN12VectorDipole10initializeEddddddddddd
.LCPI4_0:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
	.text
	.globl	_ZN12VectorDipole10initializeEddddddddddd
	.p2align	4, 0x90
	.type	_ZN12VectorDipole10initializeEddddddddddd,@function
_ZN12VectorDipole10initializeEddddddddddd: # @_ZN12VectorDipole10initializeEddddddddddd
.Lfunc_begin4:
	.loc	1 34 0                  # backgroundfield/vectordipole.cpp:34:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: initialize: <- $rdi
	#DEBUG_VALUE: initialize:moment <- $xmm0
	#DEBUG_VALUE: initialize:center_x <- $xmm1
	#DEBUG_VALUE: initialize:center_y <- $xmm2
	#DEBUG_VALUE: initialize:center_z <- $xmm3
	#DEBUG_VALUE: initialize:tilt_angle_phi <- $xmm4
	#DEBUG_VALUE: initialize:tilt_angle_theta <- $xmm5
	#DEBUG_VALUE: initialize:xlimit_f <- $xmm6
	#DEBUG_VALUE: initialize:xlimit_z <- $xmm7
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%rbx
	subq	$104, %rsp
	.cfi_offset %rbx, -24
.Ltmp8:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: moment <- $xmm0
	#DEBUG_VALUE: center_x <- $xmm1
	#DEBUG_VALUE: center_y <- $xmm2
	#DEBUG_VALUE: center_z <- $xmm3
	#DEBUG_VALUE: tilt_angle_phi <- $xmm4
	#DEBUG_VALUE: tilt_angle_theta <- $xmm5
	#DEBUG_VALUE: xlimit_f <- $xmm6
	#DEBUG_VALUE: xlimit_z <- $xmm7
	vmovsd	%xmm7, -56(%rbp)        # 8-byte Spill
	vmovsd	%xmm6, -48(%rbp)        # 8-byte Spill
	vmovsd	%xmm4, -16(%rbp)        # 8-byte Spill
	vmovsd	%xmm3, -40(%rbp)        # 8-byte Spill
	vmovsd	%xmm2, -32(%rbp)        # 8-byte Spill
	vmovsd	%xmm1, -24(%rbp)        # 8-byte Spill
	vmovaps	%xmm0, -112(%rbp)       # 16-byte Spill
	movq	%rdi, %rbx
.Ltmp9:
	#DEBUG_VALUE: IMF_Bz <- [DW_OP_plus_uconst 32] [$rbp+0]
	#DEBUG_VALUE: initialize:IMF_Bz <- [DW_OP_plus_uconst 32] [$rbp+0]
	#DEBUG_VALUE: IMF_By <- [DW_OP_plus_uconst 24] [$rbp+0]
	#DEBUG_VALUE: initialize:IMF_By <- [DW_OP_plus_uconst 24] [$rbp+0]
	#DEBUG_VALUE: IMF_Bx <- [DW_OP_plus_uconst 16] [$rbp+0]
	#DEBUG_VALUE: initialize:IMF_Bx <- [DW_OP_plus_uconst 16] [$rbp+0]
	#DEBUG_VALUE: xlimit_z <- [DW_OP_constu 56, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:xlimit_z <- [DW_OP_constu 56, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: xlimit_f <- [DW_OP_constu 48, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:xlimit_f <- [DW_OP_constu 48, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: tilt_angle_theta <- $xmm5
	#DEBUG_VALUE: initialize:tilt_angle_theta <- $xmm5
	#DEBUG_VALUE: tilt_angle_phi <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:tilt_angle_phi <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: center_z <- [DW_OP_constu 40, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:center_z <- [DW_OP_constu 40, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: center_y <- [DW_OP_constu 32, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:center_y <- [DW_OP_constu 32, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: center_x <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:center_x <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: moment <- [DW_OP_constu 112, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: initialize:moment <- [DW_OP_constu 112, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE:  <- $rbx
	#DEBUG_VALUE: initialize: <- $rbx
	.loc	1 35 1 prologue_end     # backgroundfield/vectordipole.cpp:35:1
	movb	$1, 20(%rdi)
	.loc	1 37 1                  # backgroundfield/vectordipole.cpp:37:1
	vmovaps	%xmm5, %xmm0
.Ltmp10:
	#DEBUG_VALUE: tilt_angle_theta <- $xmm0
	#DEBUG_VALUE: initialize:tilt_angle_theta <- $xmm0
	callq	__fd_sincos_1
.Ltmp11:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovaps	%xmm0, -96(%rbp)        # 16-byte Spill
	vmovapd	%xmm1, -80(%rbp)        # 16-byte Spill
	vmovsd	-16(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp12:
	#DEBUG_VALUE: tilt_angle_phi <- $xmm0
	#DEBUG_VALUE: initialize:tilt_angle_phi <- $xmm0
	.loc	1 37 1                  # backgroundfield/vectordipole.cpp:37:1
	callq	__fd_sincos_1
.Ltmp13:
	.loc	1 0 1                   # backgroundfield/vectordipole.cpp:0:1
	vmovapd	-80(%rbp), %xmm2        # 16-byte Reload
	.loc	1 37 1                  # backgroundfield/vectordipole.cpp:37:1
	vunpcklpd	-96(%rbp), %xmm2, %xmm2 # 16-byte Folded Reload
                                        # xmm2 = xmm2[0],mem[0]
	vmovddup	%xmm0, %xmm0    # xmm0 = xmm0[0,0]
	vmulpd	%xmm0, %xmm2, %xmm0
	vmovapd	-112(%rbp), %xmm3       # 16-byte Reload
.Ltmp14:
	#DEBUG_VALUE: moment <- $xmm3
	#DEBUG_VALUE: initialize:moment <- $xmm3
	vmovddup	%xmm3, %xmm2    # xmm2 = xmm3[0,0]
	vmulpd	%xmm2, %xmm0, %xmm0
	vmovapd	.LCPI4_0(%rip), %xmm2   # xmm2 = [-0.0E+0,-0.0E+0]
	vxorpd	%xmm2, %xmm0, %xmm0
	vmovupd	%xmm0, 24(%rbx)
	.loc	1 39 1 is_stmt 1        # backgroundfield/vectordipole.cpp:39:1
	vmulsd	%xmm3, %xmm1, %xmm0
	vxorpd	%xmm2, %xmm0, %xmm0
	vmovlpd	%xmm0, 40(%rbx)
	vmovsd	-24(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp15:
	#DEBUG_VALUE: center_x <- $xmm0
	#DEBUG_VALUE: initialize:center_x <- $xmm0
	.loc	1 41 1                  # backgroundfield/vectordipole.cpp:41:1
	vmovsd	%xmm0, 48(%rbx)
	vmovsd	-32(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp16:
	#DEBUG_VALUE: center_y <- $xmm0
	#DEBUG_VALUE: initialize:center_y <- $xmm0
	.loc	1 42 1                  # backgroundfield/vectordipole.cpp:42:1
	vmovsd	%xmm0, 56(%rbx)
	vmovsd	-40(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp17:
	#DEBUG_VALUE: center_z <- $xmm0
	#DEBUG_VALUE: initialize:center_z <- $xmm0
	.loc	1 43 1                  # backgroundfield/vectordipole.cpp:43:1
	vmovsd	%xmm0, 64(%rbx)
	vmovsd	-48(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp18:
	#DEBUG_VALUE: xlimit_f <- $xmm0
	#DEBUG_VALUE: initialize:xlimit_f <- $xmm0
	.loc	1 46 1                  # backgroundfield/vectordipole.cpp:46:1
	vmovsd	%xmm0, 72(%rbx)
	vmovsd	-56(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp19:
	#DEBUG_VALUE: xlimit_z <- $xmm0
	#DEBUG_VALUE: initialize:xlimit_z <- $xmm0
	.loc	1 47 1                  # backgroundfield/vectordipole.cpp:47:1
	vmovsd	%xmm0, 80(%rbx)
	.loc	1 50 1                  # backgroundfield/vectordipole.cpp:50:1
	vmovaps	16(%rbp), %xmm0
.Ltmp20:
	vmovups	%xmm0, 88(%rbx)
	.loc	1 52 1                  # backgroundfield/vectordipole.cpp:52:1
	vmovsd	32(%rbp), %xmm0         # xmm0 = mem[0],zero
	vmovsd	%xmm0, 104(%rbx)
	.loc	1 55 1                  # backgroundfield/vectordipole.cpp:55:1
	addq	$104, %rsp
	popq	%rbx
.Ltmp21:
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp22:
.Lfunc_end4:
	.size	_ZN12VectorDipole10initializeEddddddddddd, .Lfunc_end4-_ZN12VectorDipole10initializeEddddddddddd
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _ZNK12VectorDipole4callEddd
.LCPI5_0:
	.quad	4720717001298091704     # double 40592189.439999998
.LCPI5_1:
	.quad	4613937818241073152     # double 3
.LCPI5_2:
	.quad	4617315517961601024     # double 5
.LCPI5_3:
	.quad	4607182418800017408     # double 1
.LCPI5_5:
	.quad	4621819117588971520     # double 10
.LCPI5_6:
	.quad	4618441417868443648     # double 6
.LCPI5_7:
	.quad	4624633867356078080     # double 15
.LCPI5_8:
	.quad	4629137466983448576     # double 30
.LCPI5_9:
	.quad	4633641066610819072     # double 60
.LCPI5_11:
	.quad	-4609434218613702656    # double -3
.LCPI5_12:
	.quad	4602678819172646912     # double 0.5
.LCPI5_13:
	.quad	-4620693217682128896    # double -0.5
.LCPI5_15:
	.quad	-9223372036854775808    # double -0
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4
.LCPI5_4:
	.quad	4602678819172646912     # double 0.5
	.quad	4602678819172646912     # double 0.5
.LCPI5_10:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
.LCPI5_14:
	.quad	-9223372036854775808    # double -0
	.zero	8
	.text
	.globl	_ZNK12VectorDipole4callEddd
	.p2align	4, 0x90
	.type	_ZNK12VectorDipole4callEddd,@function
_ZNK12VectorDipole4callEddd:            # @_ZNK12VectorDipole4callEddd
.Lfunc_begin5:
	.loc	1 60 0                  # backgroundfield/vectordipole.cpp:60:0
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
	subq	$448, %rsp              # imm = 0x1C0
	vmovapd	%xmm0, %xmm3
.Ltmp23:
	#DEBUG_VALUE: z <- undef
	#DEBUG_VALUE: y <- undef
	#DEBUG_VALUE: x <- undef
	vxorpd	%xmm0, %xmm0, %xmm0
.Ltmp24:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE: call:x <- $xmm3
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 62 1 prologue_end     # backgroundfield/vectordipole.cpp:62:1
	cmpb	$0, 20(%rdi)
	je	.LBB5_23
.Ltmp25:
# %bb.1:                                # %L.B0000
	#DEBUG_VALUE: call:x <- $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 66 1                  # backgroundfield/vectordipole.cpp:66:1
	vpermilpd	$1, 48(%rdi), %xmm4 # xmm4 = mem[1,0]
	.loc	1 67 1                  # backgroundfield/vectordipole.cpp:67:1
	vunpcklpd	%xmm3, %xmm1, %xmm1 # xmm1 = xmm1[0],xmm3[0]
.Ltmp26:
	vsubpd	%xmm4, %xmm1, %xmm6
	.loc	1 66 1                  # backgroundfield/vectordipole.cpp:66:1
	vpermilpd	$1, %xmm6, %xmm4 # xmm4 = xmm6[1,0]
	vmovapd	%xmm4, -112(%rbp)
	.loc	1 68 1                  # backgroundfield/vectordipole.cpp:68:1
	vsubsd	64(%rdi), %xmm2, %xmm7
	vmovsd	%xmm7, -96(%rbp)
	.loc	1 70 1                  # backgroundfield/vectordipole.cpp:70:1
	vmulsd	%xmm4, %xmm4, %xmm1
	vfmadd231sd	%xmm6, %xmm6, %xmm1 # xmm1 = (xmm6 * xmm6) + xmm1
	vfmadd231sd	%xmm7, %xmm7, %xmm1 # xmm1 = (xmm7 * xmm7) + xmm1
.Ltmp27:
	#DEBUG_VALUE: r2 <- $xmm1
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	.LCPI5_0(%rip), %xmm2   # xmm2 = mem[0],zero
.Ltmp28:
	.loc	1 72 1 is_stmt 1        # backgroundfield/vectordipole.cpp:72:1
	vucomisd	%xmm1, %xmm2
	ja	.LBB5_23
.Ltmp29:
# %bb.2:                                # %L.B0001
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE: call:x <- $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 76 1                  # backgroundfield/vectordipole.cpp:76:1
	vmovsd	80(%rdi), %xmm13        # xmm13 = mem[0],zero
	vucomisd	%xmm13, %xmm4
	jae	.LBB5_3
.Ltmp30:
# %bb.5:                                # %L.B0002
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE: call:x <- $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 87 1                  # backgroundfield/vectordipole.cpp:87:1
	vsqrtsd	%xmm1, %xmm1, %xmm2
.Ltmp31:
	#DEBUG_VALUE: r1 <- $xmm2
	.loc	1 88 1                  # backgroundfield/vectordipole.cpp:88:1
	vmulsd	%xmm1, %xmm1, %xmm3
.Ltmp32:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	%xmm2, -64(%rbp)        # 8-byte Spill
.Ltmp33:
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	.loc	1 88 1                  # backgroundfield/vectordipole.cpp:88:1
	vmulsd	%xmm3, %xmm2, %xmm9
.Ltmp34:
	#DEBUG_VALUE: r5 <- $xmm9
	.loc	1 89 1 is_stmt 1        # backgroundfield/vectordipole.cpp:89:1
	vmovsd	40(%rdi), %xmm10        # xmm10 = mem[0],zero
	vmovsd	32(%rdi), %xmm11        # xmm11 = mem[0],zero
	vmovsd	24(%rdi), %xmm2         # xmm2 = mem[0],zero
	vmulsd	%xmm2, %xmm4, %xmm15
	vfmadd231sd	%xmm11, %xmm6, %xmm15 # xmm15 = (xmm6 * xmm11) + xmm15
	vfmadd231sd	%xmm10, %xmm7, %xmm15 # xmm15 = (xmm7 * xmm10) + xmm15
.Ltmp35:
	#DEBUG_VALUE: rdotq <- $xmm15
	.loc	1 90 1                  # backgroundfield/vectordipole.cpp:90:1
	movl	8(%rdi), %eax
	vmovsd	-112(%rbp,%rax,8), %xmm14 # xmm14 = mem[0],zero
	vmulsd	.LCPI5_1(%rip), %xmm14, %xmm3
	.loc	1 92 1                  # backgroundfield/vectordipole.cpp:92:1
	movl	16(%rdi), %ecx
	.loc	1 90 1                  # backgroundfield/vectordipole.cpp:90:1
	vmovsd	24(%rdi,%rax,8), %xmm12 # xmm12 = mem[0],zero
	vmulsd	%xmm12, %xmm1, %xmm5
	vfmsub231sd	%xmm3, %xmm15, %xmm5 # xmm5 = (xmm15 * xmm3) - xmm5
	vdivsd	%xmm9, %xmm5, %xmm3
.Ltmp36:
	#DEBUG_VALUE: B <- $xmm3
	.loc	1 92 1                  # backgroundfield/vectordipole.cpp:92:1
	testl	%ecx, %ecx
	je	.LBB5_6
.Ltmp37:
# %bb.8:                                # %L.B0005
	#DEBUG_VALUE: B <- $xmm3
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 0 is_stmt 0         # backgroundfield/vectordipole.cpp:0:0
	vmovsd	72(%rdi), %xmm8         # xmm8 = mem[0],zero
	movb	$1, %dl
.Ltmp38:
	.loc	1 96 1 is_stmt 1        # backgroundfield/vectordipole.cpp:96:1
	cmpl	$1, %ecx
	jne	.LBB5_11
.Ltmp39:
# %bb.9:                                # %L.B0060
	#DEBUG_VALUE: B <- $xmm3
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	vucomisd	%xmm4, %xmm8
	jae	.LBB5_14
.Ltmp40:
# %bb.10:
	#DEBUG_VALUE: B <- $xmm3
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	xorl	%edx, %edx
	jmp	.LBB5_11
.Ltmp41:
.LBB5_3:                                # %L.B0053
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE: call:x <- $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 78 1 is_stmt 1        # backgroundfield/vectordipole.cpp:78:1
	cmpl	$0, 16(%rdi)
	jne	.LBB5_23
.Ltmp42:
# %bb.4:                                # %L.B0054
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE: call:x <- $xmm3
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 79 1                  # backgroundfield/vectordipole.cpp:79:1
	movl	8(%rdi), %eax
	vmovsd	88(%rdi,%rax,8), %xmm0  # xmm0 = mem[0],zero
	jmp	.LBB5_23
.Ltmp43:
.LBB5_6:                                # %L.B0057
	#DEBUG_VALUE: B <- $xmm3
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 92 1                  # backgroundfield/vectordipole.cpp:92:1
	vmovsd	72(%rdi), %xmm8         # xmm8 = mem[0],zero
	movb	$1, %dl
	vucomisd	%xmm4, %xmm8
	jae	.LBB5_7
.Ltmp44:
.LBB5_11:                               # %L.B0006
	#DEBUG_VALUE: B <- $xmm3
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	%xmm12, -88(%rbp)       # 8-byte Spill
	vmovsd	%xmm14, -48(%rbp)       # 8-byte Spill
	vmovsd	%xmm3, -8(%rbp)         # 8-byte Spill
.Ltmp45:
	#DEBUG_VALUE: B <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	vmovapd	%xmm9, -144(%rbp)       # 16-byte Spill
.Ltmp46:
	#DEBUG_VALUE: r5 <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	.loc	1 120 1 is_stmt 1       # backgroundfield/vectordipole.cpp:120:1
	vunpcklpd	%xmm7, %xmm4, %xmm9 # xmm9 = xmm4[0],xmm7[0]
	vmovapd	%xmm11, -336(%rbp)      # 16-byte Spill
	vunpcklpd	%xmm2, %xmm11, %xmm5 # xmm5 = xmm11[0],xmm2[0]
	vmulpd	%xmm5, %xmm9, %xmm5
	vmovapd	%xmm10, -352(%rbp)      # 16-byte Spill
	vmovapd	%xmm2, -128(%rbp)       # 16-byte Spill
	vunpcklpd	%xmm10, %xmm2, %xmm2 # xmm2 = xmm2[0],xmm10[0]
	vfmsub213pd	%xmm5, %xmm6, %xmm2 # xmm2 = (xmm6 * xmm2) - xmm5
	vmovapd	%xmm2, -192(%rbp)       # 16-byte Spill
	.loc	1 124 1                 # backgroundfield/vectordipole.cpp:124:1
	vmovsd	88(%rdi), %xmm2         # xmm2 = mem[0],zero
	.loc	1 123 1                 # backgroundfield/vectordipole.cpp:123:1
	vmovsd	96(%rdi), %xmm5         # xmm5 = mem[0],zero
	vmovsd	104(%rdi), %xmm3        # xmm3 = mem[0],zero
	vmovapd	%xmm5, -304(%rbp)       # 16-byte Spill
	.loc	1 125 1                 # backgroundfield/vectordipole.cpp:125:1
	vunpcklpd	%xmm2, %xmm5, %xmm5 # xmm5 = xmm5[0],xmm2[0]
	vmulpd	%xmm5, %xmm9, %xmm10
	vmovapd	%xmm2, -320(%rbp)       # 16-byte Spill
	vmovapd	%xmm3, -288(%rbp)       # 16-byte Spill
	vunpcklpd	%xmm3, %xmm2, %xmm14 # xmm14 = xmm2[0],xmm3[0]
	.loc	1 129 1                 # backgroundfield/vectordipole.cpp:129:1
	vsubsd	%xmm4, %xmm13, %xmm5
	vsubsd	%xmm8, %xmm13, %xmm9
	vdivsd	%xmm9, %xmm5, %xmm2
.Ltmp47:
	#DEBUG_VALUE: s <- $xmm2
	.loc	1 130 1                 # backgroundfield/vectordipole.cpp:130:1
	vmulsd	%xmm2, %xmm2, %xmm11
.Ltmp48:
	#DEBUG_VALUE: ss <- $xmm11
	.loc	1 132 1                 # backgroundfield/vectordipole.cpp:132:1
	vmulsd	.LCPI5_5(%rip), %xmm11, %xmm12
	vmulsd	.LCPI5_6(%rip), %xmm11, %xmm5
	.loc	1 125 1                 # backgroundfield/vectordipole.cpp:125:1
	vfmsub213pd	%xmm10, %xmm6, %xmm14 # xmm14 = (xmm6 * xmm14) - xmm10
	.loc	1 132 1                 # backgroundfield/vectordipole.cpp:132:1
	vmulsd	%xmm5, %xmm11, %xmm5
	vmulsd	.LCPI5_7(%rip), %xmm11, %xmm3
	vmulsd	%xmm3, %xmm11, %xmm3
	vfmsub231sd	%xmm5, %xmm2, %xmm3 # xmm3 = (xmm2 * xmm5) - xmm3
	vfmadd231sd	%xmm12, %xmm2, %xmm3 # xmm3 = (xmm2 * xmm12) + xmm3
.Ltmp49:
	#DEBUG_VALUE: S2 <- $xmm3
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	%xmm3, -24(%rbp)        # 8-byte Spill
.Ltmp50:
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	vmovsd	.LCPI5_8(%rip), %xmm10  # xmm10 = mem[0],zero
	.loc	1 133 1 is_stmt 1       # backgroundfield/vectordipole.cpp:133:1
	vmulsd	%xmm10, %xmm11, %xmm12
	vmovsd	.LCPI5_9(%rip), %xmm3   # xmm3 = mem[0],zero
	vmulsd	%xmm3, %xmm11, %xmm5
	vmovsd	%xmm2, -72(%rbp)        # 8-byte Spill
.Ltmp51:
	#DEBUG_VALUE: s <- [DW_OP_constu 72, DW_OP_minus] [$rbp+0]
	vmulsd	%xmm5, %xmm2, %xmm5
	vfmsub231sd	%xmm12, %xmm11, %xmm5 # xmm5 = (xmm11 * xmm12) - xmm5
	vfmadd231sd	%xmm10, %xmm11, %xmm5 # xmm5 = (xmm11 * xmm10) + xmm5
	vdivsd	%xmm9, %xmm5, %xmm2
	vmovapd	%xmm2, -160(%rbp)       # 16-byte Spill
	vxorpd	.LCPI5_10(%rip), %xmm2, %xmm12
.Ltmp52:
	#DEBUG_VALUE: dS2dx <- $xmm12
	.loc	1 136 1                 # backgroundfield/vectordipole.cpp:136:1
	vsubsd	%xmm8, %xmm4, %xmm2
	vdivsd	%xmm9, %xmm2, %xmm10
.Ltmp53:
	#DEBUG_VALUE: IMFs <- $xmm10
	.loc	1 137 1                 # backgroundfield/vectordipole.cpp:137:1
	vmulsd	%xmm10, %xmm10, %xmm13
.Ltmp54:
	#DEBUG_VALUE: IMFS2 <- undef
	#DEBUG_VALUE: IMFss <- $xmm13
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	.LCPI5_8(%rip), %xmm2   # xmm2 = mem[0],zero
	.loc	1 140 1 is_stmt 1       # backgroundfield/vectordipole.cpp:140:1
	vmulsd	%xmm2, %xmm13, %xmm5
	vmulsd	%xmm3, %xmm13, %xmm3
	vmulsd	%xmm3, %xmm10, %xmm3
	vfmsub231sd	%xmm5, %xmm13, %xmm3 # xmm3 = (xmm13 * xmm5) - xmm3
	.loc	1 126 1                 # backgroundfield/vectordipole.cpp:126:1
	vmovsd	88(%rdi,%rax,8), %xmm5  # xmm5 = mem[0],zero
.Ltmp55:
	#DEBUG_VALUE: IMFB <- $xmm5
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	%xmm5, -16(%rbp)        # 8-byte Spill
.Ltmp56:
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	.loc	1 140 1 is_stmt 1       # backgroundfield/vectordipole.cpp:140:1
	vfmadd231sd	%xmm2, %xmm13, %xmm3 # xmm3 = (xmm13 * xmm2) + xmm3
	vmovsd	%xmm9, -80(%rbp)        # 8-byte Spill
	vdivsd	%xmm9, %xmm3, %xmm9
.Ltmp57:
	#DEBUG_VALUE: IMFdS2dx <- $xmm9
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	%xmm12, -176(%rbp)      # 16-byte Spill
.Ltmp58:
	#DEBUG_VALUE: dS2dx <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	.loc	1 144 1 is_stmt 1       # backgroundfield/vectordipole.cpp:144:1
	vmovlpd	%xmm12, -424(%rbp)
	.loc	1 146 1                 # backgroundfield/vectordipole.cpp:146:1
	vxorpd	%xmm3, %xmm3, %xmm3
	vmovupd	%xmm3, -416(%rbp)
	.loc	1 150 1                 # backgroundfield/vectordipole.cpp:150:1
	vmovsd	%xmm9, -400(%rbp)
	.loc	1 152 1                 # backgroundfield/vectordipole.cpp:152:1
	vmovupd	%xmm3, -392(%rbp)
	vmovsd	-64(%rbp), %xmm5        # 8-byte Reload
                                        # xmm5 = mem[0],zero
	.loc	1 118 1                 # backgroundfield/vectordipole.cpp:118:1
	vmulsd	%xmm5, %xmm1, %xmm12
.Ltmp59:
	#DEBUG_VALUE: A <- [DW_OP_LLVM_fragment 0 64] undef
	.loc	1 120 1                 # backgroundfield/vectordipole.cpp:120:1
	vmovddup	%xmm12, %xmm3   # xmm3 = xmm12[0,0]
	vmovapd	-192(%rbp), %xmm2       # 16-byte Reload
	vdivpd	%xmm3, %xmm2, %xmm2
.Ltmp60:
	#DEBUG_VALUE: IMFA <- [DW_OP_LLVM_fragment 0 64] undef
	.loc	1 125 1                 # backgroundfield/vectordipole.cpp:125:1
	vmulpd	.LCPI5_4(%rip), %xmm14, %xmm14
	vmovapd	%xmm8, %xmm3
.Ltmp61:
	.loc	1 154 1                 # backgroundfield/vectordipole.cpp:154:1
	vucomisd	%xmm8, %xmm5
	jbe	.LBB5_18
.Ltmp62:
# %bb.12:                               # %L.B0006
	#DEBUG_VALUE: dS2dx <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFdS2dx <- $xmm9
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFss <- $xmm13
	#DEBUG_VALUE: IMFs <- $xmm10
	#DEBUG_VALUE: s <- [DW_OP_constu 72, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: ss <- $xmm11
	#DEBUG_VALUE: r5 <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: B <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	testl	%ecx, %ecx
	jne	.LBB5_18
.Ltmp63:
# %bb.13:                               # %L.B0064
	#DEBUG_VALUE: dS2dx <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFdS2dx <- $xmm9
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFss <- $xmm13
	#DEBUG_VALUE: IMFs <- $xmm10
	#DEBUG_VALUE: s <- [DW_OP_constu 72, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: ss <- $xmm11
	#DEBUG_VALUE: r5 <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: B <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 139 1                 # backgroundfield/vectordipole.cpp:139:1
	vmulsd	.LCPI5_5(%rip), %xmm13, %xmm0
	vmulsd	.LCPI5_6(%rip), %xmm13, %xmm1
.Ltmp64:
	vmulsd	.LCPI5_7(%rip), %xmm13, %xmm3
	vmulsd	%xmm1, %xmm13, %xmm1
	vmulsd	%xmm3, %xmm13, %xmm3
	vfmsub231sd	%xmm1, %xmm10, %xmm3 # xmm3 = (xmm10 * xmm1) - xmm3
.Ltmp65:
	.loc	1 192 1                 # backgroundfield/vectordipole.cpp:192:1
	movq	$0, -376(%rbp)
	vmovapd	-160(%rbp), %xmm1       # 16-byte Reload
	.loc	1 193 1                 # backgroundfield/vectordipole.cpp:193:1
	vunpcklpd	-176(%rbp), %xmm1, %xmm1 # 16-byte Folded Reload
                                        # xmm1 = xmm1[0],mem[0]
	vmulpd	%xmm1, %xmm2, %xmm1
	vmovupd	%xmm1, -368(%rbp)
	.loc	1 198 1                 # backgroundfield/vectordipole.cpp:198:1
	movq	$0, -216(%rbp)
	.loc	1 199 1                 # backgroundfield/vectordipole.cpp:199:1
	vmulsd	%xmm9, %xmm14, %xmm1
	vxorpd	.LCPI5_10(%rip), %xmm1, %xmm1
	vmovlpd	%xmm1, -208(%rbp)
	.loc	1 200 1                 # backgroundfield/vectordipole.cpp:200:1
	vpermilpd	$1, %xmm14, %xmm1 # xmm1 = xmm14[1,0]
	vmulsd	%xmm9, %xmm1, %xmm1
	vmovsd	%xmm1, -200(%rbp)
.Ltmp66:
	.loc	1 139 1                 # backgroundfield/vectordipole.cpp:139:1
	vfmadd231sd	%xmm0, %xmm10, %xmm3 # xmm3 = (xmm10 * xmm0) + xmm3
.Ltmp67:
	#DEBUG_VALUE: IMFS2 <- $xmm3
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	-8(%rbp), %xmm0         # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp68:
	#DEBUG_VALUE: B <- $xmm0
	vmovsd	-24(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
.Ltmp69:
	.loc	1 203 1 is_stmt 1       # backgroundfield/vectordipole.cpp:203:1
	vfmadd213sd	-376(%rbp,%rax,8), %xmm0, %xmm1 # xmm1 = (xmm0 * xmm1) + mem
	vfmadd231sd	-16(%rbp), %xmm3, %xmm1 # 8-byte Folded Reload
                                        # xmm1 = (xmm3 * mem) + xmm1
	vaddsd	-216(%rbp,%rax,8), %xmm1, %xmm0
.Ltmp70:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	jmp	.LBB5_23
.Ltmp71:
.LBB5_18:                               # %L.B0011
	#DEBUG_VALUE: dS2dx <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFdS2dx <- $xmm9
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFss <- $xmm13
	#DEBUG_VALUE: IMFs <- $xmm10
	#DEBUG_VALUE: s <- [DW_OP_constu 72, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: ss <- $xmm11
	#DEBUG_VALUE: r5 <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: B <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 154 1 is_stmt 1       # backgroundfield/vectordipole.cpp:154:1
	vucomisd	%xmm3, %xmm5
	setbe	%cl
.Ltmp72:
	.loc	1 206 1                 # backgroundfield/vectordipole.cpp:206:1
	orb	%cl, %dl
	jne	.LBB5_23
.Ltmp73:
# %bb.19:                               # %L.B0067
	#DEBUG_VALUE: dS2dx <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFdS2dx <- $xmm9
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFss <- $xmm13
	#DEBUG_VALUE: IMFs <- $xmm10
	#DEBUG_VALUE: s <- [DW_OP_constu 72, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: ss <- $xmm11
	#DEBUG_VALUE: r5 <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: B <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	%xmm2, -256(%rbp)       # 16-byte Spill
	vmovapd	%xmm9, -64(%rbp)        # 16-byte Spill
.Ltmp74:
	#DEBUG_VALUE: IMFdS2dx <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	vmovsd	-8(%rbp), %xmm0         # 8-byte Reload
                                        # xmm0 = mem[0],zero
.Ltmp75:
	#DEBUG_VALUE: B <- $xmm0
	.loc	1 238 1 is_stmt 1       # backgroundfield/vectordipole.cpp:238:1
	vmulsd	.LCPI5_2(%rip), %xmm0, %xmm0
.Ltmp76:
	movl	12(%rdi), %ecx
	vmovsd	-112(%rbp,%rcx,8), %xmm2 # xmm2 = mem[0],zero
	vmulsd	%xmm2, %xmm0, %xmm0
	vmovapd	%xmm0, -240(%rbp)       # 16-byte Spill
	vmovsd	.LCPI5_1(%rip), %xmm5   # xmm5 = mem[0],zero
	vmulsd	24(%rdi,%rcx,8), %xmm5, %xmm3
	vmovsd	-88(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vaddsd	%xmm0, %xmm0, %xmm0
	vmulsd	%xmm2, %xmm0, %xmm0
	vfmsub231sd	-48(%rbp), %xmm3, %xmm0 # 8-byte Folded Reload
                                        # xmm0 = (xmm3 * mem) - xmm0
	vmulsd	%xmm5, %xmm15, %xmm2
	vxorpd	%xmm15, %xmm15, %xmm15
.Ltmp77:
	cmpl	%eax, %ecx
	vmovapd	%xmm14, -272(%rbp)      # 16-byte Spill
	je	.LBB5_20
.Ltmp78:
# %bb.21:                               # %L.B0067
	#DEBUG_VALUE: IMFdS2dx <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dS2dx <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFss <- $xmm13
	#DEBUG_VALUE: IMFs <- $xmm10
	#DEBUG_VALUE: s <- [DW_OP_constu 72, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: ss <- $xmm11
	#DEBUG_VALUE: r5 <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	vxorpd	%xmm3, %xmm3, %xmm3
	jmp	.LBB5_22
.Ltmp79:
.LBB5_7:
	#DEBUG_VALUE: B <- $xmm3
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	%xmm3, %xmm0
	jmp	.LBB5_23
.Ltmp80:
.LBB5_14:                               # %L.B0061
	#DEBUG_VALUE: B <- $xmm3
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 106 1 is_stmt 1       # backgroundfield/vectordipole.cpp:106:1
	vmulsd	.LCPI5_2(%rip), %xmm3, %xmm0
	movl	12(%rdi), %ecx
	vmovsd	-112(%rbp,%rcx,8), %xmm2 # xmm2 = mem[0],zero
	vmovsd	.LCPI5_1(%rip), %xmm3   # xmm3 = mem[0],zero
.Ltmp81:
	vmulsd	24(%rdi,%rcx,8), %xmm3, %xmm4
	vmulsd	%xmm2, %xmm0, %xmm0
	vaddsd	%xmm12, %xmm12, %xmm5
	vmulsd	%xmm2, %xmm5, %xmm2
	vfmsub231sd	%xmm4, %xmm14, %xmm2 # xmm2 = (xmm14 * xmm4) - xmm2
	vmulsd	%xmm3, %xmm15, %xmm3
	cmpl	%eax, %ecx
	je	.LBB5_15
.Ltmp82:
# %bb.16:                               # %L.B0061
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vxorpd	%xmm4, %xmm4, %xmm4
	jmp	.LBB5_17
.Ltmp83:
.LBB5_20:
	#DEBUG_VALUE: IMFdS2dx <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dS2dx <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFss <- $xmm13
	#DEBUG_VALUE: IMFs <- $xmm10
	#DEBUG_VALUE: s <- [DW_OP_constu 72, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: ss <- $xmm11
	#DEBUG_VALUE: r5 <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	vmovsd	.LCPI5_3(%rip), %xmm3   # xmm3 = mem[0],zero
.Ltmp84:
.LBB5_22:                               # %L.B0067
	#DEBUG_VALUE: IMFdS2dx <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dS2dx <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFB <- [DW_OP_constu 16, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: IMFss <- $xmm13
	#DEBUG_VALUE: IMFs <- $xmm10
	#DEBUG_VALUE: s <- [DW_OP_constu 72, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: S2 <- [DW_OP_constu 24, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: ss <- $xmm11
	#DEBUG_VALUE: r5 <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	vmovapd	-128(%rbp), %xmm13      # 16-byte Reload
.Ltmp85:
	.loc	1 238 1 is_stmt 1       # backgroundfield/vectordipole.cpp:238:1
	vfmadd213sd	%xmm0, %xmm2, %xmm3 # xmm3 = (xmm2 * xmm3) + xmm0
	vmovapd	%xmm3, -48(%rbp)        # 16-byte Spill
	vmovsd	.LCPI5_11(%rip), %xmm9  # xmm9 = mem[0],zero
	.loc	1 249 1                 # backgroundfield/vectordipole.cpp:249:1
	vdivsd	-144(%rbp), %xmm9, %xmm3 # 16-byte Folded Reload
	vmovapd	-192(%rbp), %xmm2       # 16-byte Reload
	vpermilpd	$1, %xmm2, %xmm0 # xmm0 = xmm2[1,0]
	vmulsd	%xmm0, %xmm3, %xmm8
	.loc	1 252 1                 # backgroundfield/vectordipole.cpp:252:1
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovapd	-352(%rbp), %xmm3       # 16-byte Reload
	.loc	1 249 1                 # backgroundfield/vectordipole.cpp:249:1
	vdivsd	%xmm12, %xmm3, %xmm3
	.loc	1 251 1                 # backgroundfield/vectordipole.cpp:251:1
	vdivsd	%xmm12, %xmm13, %xmm14
	vmovapd	-336(%rbp), %xmm0       # 16-byte Reload
	.loc	1 252 1                 # backgroundfield/vectordipole.cpp:252:1
	vdivsd	%xmm12, %xmm0, %xmm5
	.loc	1 249 1                 # backgroundfield/vectordipole.cpp:249:1
	vfmadd231sd	%xmm8, %xmm4, %xmm3 # xmm3 = (xmm4 * xmm8) + xmm3
.Ltmp86:
	#DEBUG_VALUE: delAy <- [DW_OP_LLVM_fragment 0 64] $xmm3
	.loc	1 252 1                 # backgroundfield/vectordipole.cpp:252:1
	vfmsub213sd	%xmm5, %xmm2, %xmm4 # xmm4 = (xmm2 * xmm4) - xmm5
.Ltmp87:
	#DEBUG_VALUE: delAz <- [DW_OP_LLVM_fragment 0 64] $xmm4
	.loc	1 250 1                 # backgroundfield/vectordipole.cpp:250:1
	vmulsd	%xmm8, %xmm6, %xmm13
.Ltmp88:
	#DEBUG_VALUE: delAy <- [DW_OP_LLVM_fragment 64 64] $xmm13
	.loc	1 251 1                 # backgroundfield/vectordipole.cpp:251:1
	vfmsub213sd	%xmm14, %xmm7, %xmm8 # xmm8 = (xmm7 * xmm8) - xmm14
.Ltmp89:
	#DEBUG_VALUE: delAy <- [DW_OP_LLVM_fragment 128 64] $xmm8
	.loc	1 253 1                 # backgroundfield/vectordipole.cpp:253:1
	vfmadd231sd	%xmm2, %xmm6, %xmm14 # xmm14 = (xmm6 * xmm2) + xmm14
.Ltmp90:
	#DEBUG_VALUE: delAz <- [DW_OP_LLVM_fragment 64 64] $xmm14
	.loc	1 254 1                 # backgroundfield/vectordipole.cpp:254:1
	vmulsd	%xmm2, %xmm7, %xmm2
.Ltmp91:
	#DEBUG_VALUE: delAz <- [DW_OP_LLVM_fragment 128 64] $xmm2
	.loc	1 351 1                 # backgroundfield/vectordipole.cpp:351:1
	vmulsd	%xmm9, %xmm11, %xmm5
.Ltmp92:
	.loc	1 277 1                 # backgroundfield/vectordipole.cpp:277:1
	vaddsd	%xmm11, %xmm11, %xmm6
	#DEBUG_VALUE: delAz <- [DW_OP_LLVM_fragment 128 64] $xmm2
.Ltmp93:
	#DEBUG_VALUE: delAz <- [DW_OP_LLVM_fragment 64 64] $xmm14
	#DEBUG_VALUE: delAy <- [DW_OP_LLVM_fragment 128 64] $xmm8
	#DEBUG_VALUE: delAy <- [DW_OP_LLVM_fragment 64 64] $xmm13
	#DEBUG_VALUE: delAz <- [DW_OP_LLVM_fragment 0 64] $xmm4
	#DEBUG_VALUE: delAy <- [DW_OP_LLVM_fragment 0 64] $xmm3
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	-72(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	.loc	1 277 1                 # backgroundfield/vectordipole.cpp:277:1
	vfmadd231sd	%xmm6, %xmm0, %xmm5 # xmm5 = (xmm0 * xmm6) + xmm5
	vaddsd	%xmm5, %xmm0, %xmm5
	vmulsd	.LCPI5_9(%rip), %xmm5, %xmm5
	.loc	1 323 1 is_stmt 1       # backgroundfield/vectordipole.cpp:323:1
	vunpcklpd	%xmm3, %xmm2, %xmm9 # xmm9 = xmm2[0],xmm3[0]
	vmovsd	-80(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	.loc	1 277 1                 # backgroundfield/vectordipole.cpp:277:1
	vmulsd	%xmm0, %xmm0, %xmm2
.Ltmp94:
	vdivsd	%xmm2, %xmm5, %xmm5
.Ltmp95:
	#DEBUG_VALUE: deldS2dx <- [DW_OP_LLVM_fragment 128 64] 0.000000e+00
	#DEBUG_VALUE: deldS2dx <- [DW_OP_LLVM_fragment 64 64] 0.000000e+00
	#DEBUG_VALUE: deldS2dx <- [DW_OP_LLVM_fragment 0 64] $xmm5
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	-256(%rbp), %xmm3       # 16-byte Reload
.Ltmp96:
	.loc	1 321 1 is_stmt 1       # backgroundfield/vectordipole.cpp:321:1
	vmulsd	%xmm5, %xmm3, %xmm10
.Ltmp97:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	-160(%rbp), %xmm6       # 16-byte Reload
	.loc	1 321 1                 # backgroundfield/vectordipole.cpp:321:1
	vfmsub231sd	%xmm4, %xmm6, %xmm10 # xmm10 = (xmm6 * xmm4) - xmm10
	.loc	1 322 1 is_stmt 1       # backgroundfield/vectordipole.cpp:322:1
	vmovapd	.LCPI5_14(%rip), %xmm2  # xmm2 = <-0.0E+0,u>
	vunpcklpd	%xmm5, %xmm2, %xmm7 # xmm7 = xmm2[0],xmm5[0]
	vmulpd	%xmm7, %xmm3, %xmm4
.Ltmp98:
	vfmadd213sd	%xmm4, %xmm6, %xmm14 # xmm14 = (xmm6 * xmm14) + xmm4
.Ltmp99:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	-176(%rbp), %xmm0       # 16-byte Reload
	.loc	1 323 1 is_stmt 1       # backgroundfield/vectordipole.cpp:323:1
	vunpcklpd	%xmm0, %xmm6, %xmm2 # xmm2 = xmm6[0],xmm0[0]
	vfmadd213pd	%xmm4, %xmm9, %xmm2 # xmm2 = (xmm9 * xmm2) + xmm4
	.loc	1 326 1                 # backgroundfield/vectordipole.cpp:326:1
	vpermilpd	$1, %xmm3, %xmm3 # xmm3 = xmm3[1,0]
	vmulsd	%xmm15, %xmm3, %xmm3
	vfmadd213sd	%xmm3, %xmm0, %xmm13 # xmm13 = (xmm0 * xmm13) + xmm3
.Ltmp100:
	.loc	1 327 1                 # backgroundfield/vectordipole.cpp:327:1
	vfmadd213sd	%xmm3, %xmm0, %xmm8 # xmm8 = (xmm0 * xmm8) + xmm3
.Ltmp101:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	.LCPI5_12(%rip), %xmm3  # xmm3 = mem[0],zero
.Ltmp102:
	#DEBUG_VALUE: IMFdelAx <- [DW_OP_LLVM_fragment 128 64] undef
	#DEBUG_VALUE: IMFdelAx <- [DW_OP_LLVM_fragment 64 64] undef
	#DEBUG_VALUE: IMFdelAx <- [DW_OP_LLVM_fragment 0 64] 0.000000e+00
	.loc	1 268 1 is_stmt 1       # backgroundfield/vectordipole.cpp:268:1
	vmulsd	-288(%rbp), %xmm3, %xmm9 # 16-byte Folded Reload
.Ltmp103:
	#DEBUG_VALUE: IMFdelAy <- [DW_OP_LLVM_fragment 0 64] $xmm9
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	.LCPI5_13(%rip), %xmm4  # xmm4 = mem[0],zero
.Ltmp104:
	#DEBUG_VALUE: IMFdelAy <- [DW_OP_LLVM_fragment 64 64] 0.000000e+00
	.loc	1 271 1 is_stmt 1       # backgroundfield/vectordipole.cpp:271:1
	vmulsd	-304(%rbp), %xmm4, %xmm6 # 16-byte Folded Reload
.Ltmp105:
	#DEBUG_VALUE: IMFdelAz <- [DW_OP_LLVM_fragment 0 64] $xmm6
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	-320(%rbp), %xmm0       # 16-byte Reload
	.loc	1 270 1 is_stmt 1       # backgroundfield/vectordipole.cpp:270:1
	vmulsd	%xmm4, %xmm0, %xmm4
.Ltmp106:
	#DEBUG_VALUE: IMFdelAy <- [DW_OP_LLVM_fragment 128 64] $xmm4
	.loc	1 272 1                 # backgroundfield/vectordipole.cpp:272:1
	vmulsd	%xmm3, %xmm0, %xmm3
.Ltmp107:
	#DEBUG_VALUE: IMFdelAz <- [DW_OP_LLVM_fragment 128 64] 0.000000e+00
	#DEBUG_VALUE: IMFdelAz <- [DW_OP_LLVM_fragment 64 64] $xmm3
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	-272(%rbp), %xmm12      # 16-byte Reload
	.loc	1 337 1 is_stmt 1       # backgroundfield/vectordipole.cpp:337:1
	vmulsd	%xmm5, %xmm12, %xmm5
.Ltmp108:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	-64(%rbp), %xmm0        # 16-byte Reload
.Ltmp109:
	#DEBUG_VALUE: IMFdS2dx <- $xmm0
	.loc	1 337 1                 # backgroundfield/vectordipole.cpp:337:1
	vfnmsub231sd	%xmm6, %xmm0, %xmm5 # xmm5 = -(xmm0 * xmm6) - xmm5
.Ltmp110:
	.loc	1 351 1 is_stmt 1       # backgroundfield/vectordipole.cpp:351:1
	vmovddup	.LCPI5_15(%rip), %xmm6 # xmm6 = [-0.0E+0,-0.0E+0]
                                        # xmm6 = mem[0,0]
.Ltmp111:
	vxorpd	%xmm6, %xmm0, %xmm6
.Ltmp112:
	.loc	1 338 1                 # backgroundfield/vectordipole.cpp:338:1
	vmulpd	%xmm7, %xmm12, %xmm7
.Ltmp113:
	#DEBUG_VALUE: IMFdelAz <- [DW_OP_LLVM_fragment 128 64] 0.000000e+00
	#DEBUG_VALUE: IMFdelAz <- [DW_OP_LLVM_fragment 64 64] $xmm3
	#DEBUG_VALUE: IMFdelAy <- [DW_OP_LLVM_fragment 128 64] $xmm4
	#DEBUG_VALUE: IMFdelAy <- [DW_OP_LLVM_fragment 64 64] 0.000000e+00
	#DEBUG_VALUE: IMFdelAy <- [DW_OP_LLVM_fragment 0 64] $xmm9
	#DEBUG_VALUE: IMFdelAx <- [DW_OP_LLVM_fragment 0 64] 0.000000e+00
	#DEBUG_VALUE: deldS2dx <- [DW_OP_LLVM_fragment 128 64] 0.000000e+00
	#DEBUG_VALUE: deldS2dx <- [DW_OP_LLVM_fragment 64 64] 0.000000e+00
	.loc	1 339 1                 # backgroundfield/vectordipole.cpp:339:1
	vunpcklpd	%xmm0, %xmm6, %xmm6 # xmm6 = xmm6[0],xmm0[0]
	.loc	1 318 1                 # backgroundfield/vectordipole.cpp:318:1
	vxorpd	%xmm11, %xmm11, %xmm11
.Ltmp114:
	.loc	1 339 1                 # backgroundfield/vectordipole.cpp:339:1
	vunpcklpd	%xmm9, %xmm11, %xmm9 # xmm9 = xmm11[0],xmm9[0]
.Ltmp115:
	vfmadd213pd	%xmm7, %xmm6, %xmm9 # xmm9 = (xmm6 * xmm9) + xmm7
	.loc	1 338 1                 # backgroundfield/vectordipole.cpp:338:1
	vfnmadd213sd	%xmm7, %xmm0, %xmm3 # xmm3 = -(xmm0 * xmm3) + xmm7
.Ltmp116:
	.loc	1 342 1                 # backgroundfield/vectordipole.cpp:342:1
	vpermilpd	$1, %xmm12, %xmm6 # xmm6 = xmm12[1,0]
	vmulsd	%xmm15, %xmm6, %xmm6
	vfmadd213sd	%xmm6, %xmm0, %xmm15 # xmm15 = (xmm0 * xmm15) + xmm6
	.loc	1 343 1                 # backgroundfield/vectordipole.cpp:343:1
	vfmadd213sd	%xmm6, %xmm0, %xmm4 # xmm4 = (xmm0 * xmm4) + xmm6
.Ltmp117:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovapd	-48(%rbp), %xmm0        # 16-byte Reload
.Ltmp118:
	.loc	1 238 1 is_stmt 1       # backgroundfield/vectordipole.cpp:238:1
	vunpcklpd	-240(%rbp), %xmm0, %xmm6 # 16-byte Folded Reload
                                        # xmm6 = xmm0[0],mem[0]
	vmovapd	-144(%rbp), %xmm0       # 16-byte Reload
.Ltmp119:
	#DEBUG_VALUE: r5 <- $xmm0
	vunpcklpd	%xmm1, %xmm0, %xmm1 # xmm1 = xmm0[0],xmm1[0]
.Ltmp120:
	.loc	1 318 1                 # backgroundfield/vectordipole.cpp:318:1
	vmovapd	%xmm11, -576(%rbp)
	movq	$0, -560(%rbp)
	.loc	1 321 1                 # backgroundfield/vectordipole.cpp:321:1
	vmovsd	%xmm10, -552(%rbp)
	.loc	1 322 1                 # backgroundfield/vectordipole.cpp:322:1
	vmovsd	%xmm14, -544(%rbp)
	.loc	1 323 1                 # backgroundfield/vectordipole.cpp:323:1
	vmovupd	%xmm2, -536(%rbp)
	.loc	1 326 1                 # backgroundfield/vectordipole.cpp:326:1
	vmovsd	%xmm13, -520(%rbp)
	.loc	1 238 1                 # backgroundfield/vectordipole.cpp:238:1
	vdivpd	%xmm1, %xmm6, %xmm1
	.loc	1 327 1                 # backgroundfield/vectordipole.cpp:327:1
	vmovsd	%xmm8, -512(%rbp)
	.loc	1 334 1                 # backgroundfield/vectordipole.cpp:334:1
	vmovapd	%xmm11, -496(%rbp)
	.loc	1 238 1                 # backgroundfield/vectordipole.cpp:238:1
	vpermilpd	$1, %xmm1, %xmm0 # xmm0 = xmm1[1,0]
.Ltmp121:
	vsubsd	%xmm0, %xmm1, %xmm0
.Ltmp122:
	#DEBUG_VALUE: delB <- $xmm0
	.loc	1 334 1                 # backgroundfield/vectordipole.cpp:334:1
	movq	$0, -480(%rbp)
	.loc	1 337 1                 # backgroundfield/vectordipole.cpp:337:1
	vmovsd	%xmm5, -472(%rbp)
	.loc	1 338 1                 # backgroundfield/vectordipole.cpp:338:1
	vmovsd	%xmm3, -464(%rbp)
	.loc	1 339 1                 # backgroundfield/vectordipole.cpp:339:1
	vmovupd	%xmm9, -456(%rbp)
	.loc	1 342 1                 # backgroundfield/vectordipole.cpp:342:1
	vmovsd	%xmm15, -440(%rbp)
	vmovsd	-8(%rbp), %xmm1         # 8-byte Reload
                                        # xmm1 = mem[0],zero
	.loc	1 346 1                 # backgroundfield/vectordipole.cpp:346:1
	vmulsd	-424(%rbp,%rcx,8), %xmm1, %xmm1
	vfmadd231sd	-24(%rbp), %xmm0, %xmm1 # 8-byte Folded Reload
                                        # xmm1 = (xmm0 * mem) + xmm1
	.loc	1 343 1                 # backgroundfield/vectordipole.cpp:343:1
	vmovsd	%xmm4, -432(%rbp)
	.loc	1 346 1                 # backgroundfield/vectordipole.cpp:346:1
	leaq	(%rax,%rax,2), %rax
	shlq	$3, %rax
	leaq	(%rax,%rcx,8), %rax
	vaddsd	-576(%rbp,%rax), %xmm1, %xmm0
.Ltmp123:
	.loc	1 0 1 is_stmt 0         # backgroundfield/vectordipole.cpp:0:1
	vmovsd	-16(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	.loc	1 346 1                 # backgroundfield/vectordipole.cpp:346:1
	vfmadd231sd	-400(%rbp,%rcx,8), %xmm1, %xmm0 # xmm0 = (xmm1 * mem) + xmm0
	vaddsd	-496(%rbp,%rax), %xmm0, %xmm0
	jmp	.LBB5_23
.Ltmp124:
.LBB5_15:
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 1                   # backgroundfield/vectordipole.cpp:0:1
	vmovsd	.LCPI5_3(%rip), %xmm4   # xmm4 = mem[0],zero
.Ltmp125:
.LBB5_17:                               # %L.B0061
	#DEBUG_VALUE: rdotq <- $xmm15
	#DEBUG_VALUE: r5 <- $xmm9
	#DEBUG_VALUE: r1 <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: r2 <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 106 1 is_stmt 1       # backgroundfield/vectordipole.cpp:106:1
	vfmadd213sd	%xmm2, %xmm3, %xmm4 # xmm4 = (xmm3 * xmm4) + xmm2
	vunpcklpd	%xmm0, %xmm4, %xmm0 # xmm0 = xmm4[0],xmm0[0]
	vunpcklpd	%xmm1, %xmm9, %xmm1 # xmm1 = xmm9[0],xmm1[0]
.Ltmp126:
	vdivpd	%xmm1, %xmm0, %xmm0
	vpermilpd	$1, %xmm0, %xmm1 # xmm1 = xmm0[1,0]
	vsubsd	%xmm1, %xmm0, %xmm0
.Ltmp127:
.LBB5_23:                               # %L.B0049
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	1 0 0 is_stmt 0         # backgroundfield/vectordipole.cpp:0:0
	addq	$448, %rsp              # imm = 0x1C0
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp128:
.Lfunc_end5:
	.size	_ZNK12VectorDipole4callEddd, .Lfunc_end5-_ZNK12VectorDipole4callEddd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN12VectorDipoleD1Ev   # -- Begin function _ZN12VectorDipoleD1Ev
	.p2align	4, 0x90
	.type	_ZN12VectorDipoleD1Ev,@function
_ZN12VectorDipoleD1Ev:                  # @_ZN12VectorDipoleD1Ev
.Lfunc_begin6:
	.file	4 "/home/talgat/vlasiator/backgroundfield/vectordipole.hpp"
	.loc	4 47 0 is_stmt 1        # backgroundfield/vectordipole.hpp:47:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~VectorDipole: <- $rdi
	#DEBUG_VALUE: ~VectorDipole: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp129:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	4 47 1 prologue_end     # backgroundfield/vectordipole.hpp:47:1
	movq	$_ZTV11T3DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp130:
.Lfunc_end6:
	.size	_ZN12VectorDipoleD1Ev, .Lfunc_end6-_ZN12VectorDipoleD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN12VectorDipoleD0Ev   # -- Begin function _ZN12VectorDipoleD0Ev
	.p2align	4, 0x90
	.type	_ZN12VectorDipoleD0Ev,@function
_ZN12VectorDipoleD0Ev:                  # @_ZN12VectorDipoleD0Ev
.Lfunc_begin7:
	.loc	1 0 0                   # backgroundfield/vectordipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp131:
	.loc	1 47 1 prologue_end     # backgroundfield/vectordipole.cpp:47:1
	movl	$112, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp132:
.Lfunc_end7:
	.size	_ZN12VectorDipoleD0Ev, .Lfunc_end7-_ZN12VectorDipoleD0Ev
	.cfi_endproc
                                        # -- End function
	.globl	__sti___32_backgroundfield_vectordipole_cpp_83390dc3 # -- Begin function __sti___32_backgroundfield_vectordipole_cpp_83390dc3
	.p2align	4, 0x90
	.type	__sti___32_backgroundfield_vectordipole_cpp_83390dc3,@function
__sti___32_backgroundfield_vectordipole_cpp_83390dc3: # @__sti___32_backgroundfield_vectordipole_cpp_83390dc3
.Lfunc_begin8:
	.loc	1 0 0                   # backgroundfield/vectordipole.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp133:
	.loc	1 74 1 prologue_end     # backgroundfield/vectordipole.cpp:74:1
	cmpl	$1, __I___32_backgroundfield_vectordipole_cpp_83390dc3(%rip)
	jne	.LBB8_2
# %bb.1:                                # %L.B0022
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB8_2:                                # %L.B0075
	.cfi_def_cfa %rbp, 16
	movl	$1, __I___32_backgroundfield_vectordipole_cpp_83390dc3(%rip)
	movl	$_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE, %edi
	callq	_ZNSt8ios_base4InitC1Ev
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	movl	$_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE, %esi
	movl	$__dso_handle, %edx
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	__cxa_atexit            # TAILCALL
.Ltmp134:
.Lfunc_end8:
	.size	__sti___32_backgroundfield_vectordipole_cpp_83390dc3, .Lfunc_end8-__sti___32_backgroundfield_vectordipole_cpp_83390dc3
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

	.type	_ZTV12VectorDipole,@object # @_ZTV12VectorDipole
	.weak	_ZTV12VectorDipole
	.p2align	4
_ZTV12VectorDipole:
	.quad	0
	.quad	_ZTI12VectorDipole
	.quad	_ZNK12VectorDipole4callEddd
	.quad	_ZN12VectorDipoleD1Ev
	.quad	_ZN12VectorDipoleD0Ev
	.size	_ZTV12VectorDipole, 40

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

	.type	_ZTI12VectorDipole,@object # @_ZTI12VectorDipole
	.weak	_ZTI12VectorDipole
	.p2align	4
_ZTI12VectorDipole:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS12VectorDipole
	.quad	_ZTI13FieldFunction
	.size	_ZTI12VectorDipole, 24

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

	.type	_ZTS12VectorDipole,@object # @_ZTS12VectorDipole
	.weak	_ZTS12VectorDipole
	.p2align	3
_ZTS12VectorDipole:
	.asciz	"12VectorDipole"
	.size	_ZTS12VectorDipole, 15

	.type	__I___32_backgroundfield_vectordipole_cpp_83390dc3,@object # @__I___32_backgroundfield_vectordipole_cpp_83390dc3
	.bss
	.globl	__I___32_backgroundfield_vectordipole_cpp_83390dc3
	.p2align	2
__I___32_backgroundfield_vectordipole_cpp_83390dc3:
	.long	0                       # 0x0
	.size	__I___32_backgroundfield_vectordipole_cpp_83390dc3, 4

	.type	_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE,@object # @_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE
	.local	_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE
	.comm	_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE,1,1
	.section	.init_array,"aw",@init_array
	.p2align	3
	.quad	__sti___32_backgroundfield_vectordipole_cpp_83390dc3
	.section	.debug_str,"MS",@progbits,1
.Linfo_string0:
	.asciz	" NVC++ 21.2-0"         # string offset=0
.Linfo_string1:
	.asciz	"backgroundfield/vectordipole.cpp" # string offset=14
.Linfo_string2:
	.asciz	"/home/talgat/vlasiator" # string offset=47
.Linfo_string3:
	.asciz	"_ZTV11T3DFunction"     # string offset=70
.Linfo_string4:
	.asciz	"int"                   # string offset=88
.Linfo_string5:
	.asciz	"__ARRAY_SIZE_TYPE__"   # string offset=92
.Linfo_string6:
	.asciz	"_ZTV13FieldFunction"   # string offset=112
.Linfo_string7:
	.asciz	"_ZTV12VectorDipole"    # string offset=132
.Linfo_string8:
	.asciz	"__I___32_backgroundfield_vectordipole_cpp_83390dc3" # string offset=151
.Linfo_string9:
	.asciz	"_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE" # string offset=202
.Linfo_string10:
	.asciz	"signed char"           # string offset=274
.Linfo_string11:
	.asciz	"_ZNSt8ios_base4InitE"  # string offset=286
.Linfo_string12:
	.asciz	"__dso_handle"          # string offset=307
.Linfo_string13:
	.asciz	"void"                  # string offset=320
.Linfo_string14:
	.asciz	"_ZTI11T3DFunction"     # string offset=325
.Linfo_string15:
	.asciz	"__class_type_info"     # string offset=343
.Linfo_string16:
	.asciz	"_ZTI13FieldFunction"   # string offset=361
.Linfo_string17:
	.asciz	"__si_class_type_info"  # string offset=381
.Linfo_string18:
	.asciz	"WID"                   # string offset=402
.Linfo_string19:
	.asciz	"_ZTI12VectorDipole"    # string offset=406
.Linfo_string20:
	.asciz	"WID2"                  # string offset=425
.Linfo_string21:
	.asciz	"WID3"                  # string offset=430
.Linfo_string22:
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # string offset=435
.Linfo_string23:
	.asciz	"_ZTS11T3DFunction"     # string offset=473
.Linfo_string24:
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # string offset=491
.Linfo_string25:
	.asciz	"_ZTS13FieldFunction"   # string offset=532
.Linfo_string26:
	.asciz	"_ZTS12VectorDipole"    # string offset=552
.Linfo_string27:
	.asciz	"_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc35vmesh15INVALID_LOCALIDE" # string offset=571
.Linfo_string28:
	.asciz	"unsigned"              # string offset=655
.Linfo_string29:
	.asciz	"_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc317physicalconstants3R_EE" # string offset=664
.Linfo_string30:
	.asciz	"double"                # string offset=748
.Linfo_string31:
	.asciz	"T3DFunction"           # string offset=755
.Linfo_string32:
	.asciz	"FieldFunction"         # string offset=767
.Linfo_string33:
	.asciz	"_fComponent"           # string offset=781
.Linfo_string34:
	.asciz	"X"                     # string offset=793
.Linfo_string35:
	.asciz	"Y"                     # string offset=795
.Linfo_string36:
	.asciz	"Z"                     # string offset=797
.Linfo_string37:
	.asciz	"coordinate"            # string offset=799
.Linfo_string38:
	.asciz	"_dComponent"           # string offset=810
.Linfo_string39:
	.asciz	"_derivative"           # string offset=822
.Linfo_string40:
	.asciz	"initialized"           # string offset=834
.Linfo_string41:
	.asciz	"q"                     # string offset=846
.Linfo_string42:
	.asciz	"center"                # string offset=848
.Linfo_string43:
	.asciz	"xlimit"                # string offset=855
.Linfo_string44:
	.asciz	"IMF"                   # string offset=862
.Linfo_string45:
	.asciz	"VectorDipole"          # string offset=866
.Linfo_string46:
	.asciz	"~T3DFunction"          # string offset=879
.Linfo_string47:
	.asciz	"_ZN11T3DFunctionD0Ev"  # string offset=892
.Linfo_string48:
	.asciz	"~FieldFunction"        # string offset=913
.Linfo_string49:
	.asciz	"_ZN13FieldFunctionD0Ev" # string offset=928
.Linfo_string50:
	.asciz	"initialize"            # string offset=951
.Linfo_string51:
	.asciz	"call"                  # string offset=962
.Linfo_string52:
	.asciz	"~VectorDipole"         # string offset=967
.Linfo_string53:
	.asciz	"_ZN12VectorDipoleD0Ev" # string offset=981
.Linfo_string54:
	.asciz	"__sti___32_backgroundfield_vectordipole_cpp_83390dc3" # string offset=1003
.Linfo_string55:
	.asciz	"moment"                # string offset=1056
.Linfo_string56:
	.asciz	"center_x"              # string offset=1063
.Linfo_string57:
	.asciz	"center_y"              # string offset=1072
.Linfo_string58:
	.asciz	"center_z"              # string offset=1081
.Linfo_string59:
	.asciz	"tilt_angle_phi"        # string offset=1090
.Linfo_string60:
	.asciz	"tilt_angle_theta"      # string offset=1105
.Linfo_string61:
	.asciz	"xlimit_f"              # string offset=1122
.Linfo_string62:
	.asciz	"xlimit_z"              # string offset=1131
.Linfo_string63:
	.asciz	"IMF_Bz"                # string offset=1140
.Linfo_string64:
	.asciz	"IMF_By"                # string offset=1147
.Linfo_string65:
	.asciz	"IMF_Bx"                # string offset=1154
.Linfo_string66:
	.asciz	"r"                     # string offset=1161
.Linfo_string67:
	.asciz	"dS2cart"               # string offset=1163
.Linfo_string68:
	.asciz	"IMFdS2cart"            # string offset=1171
.Linfo_string69:
	.asciz	"delS2crossA"           # string offset=1182
.Linfo_string70:
	.asciz	"IMFdelS2crossA"        # string offset=1194
.Linfo_string71:
	.asciz	"ddS2crossA"            # string offset=1209
.Linfo_string72:
	.asciz	"IMFddS2crossA"         # string offset=1220
.Linfo_string73:
	.asciz	"x"                     # string offset=1234
.Linfo_string74:
	.asciz	"y"                     # string offset=1236
.Linfo_string75:
	.asciz	"z"                     # string offset=1238
.Linfo_string76:
	.asciz	"r2"                    # string offset=1240
.Linfo_string77:
	.asciz	"r1"                    # string offset=1243
.Linfo_string78:
	.asciz	"r5"                    # string offset=1246
.Linfo_string79:
	.asciz	"rdotq"                 # string offset=1249
.Linfo_string80:
	.asciz	"B"                     # string offset=1255
.Linfo_string81:
	.asciz	"s"                     # string offset=1257
.Linfo_string82:
	.asciz	"ss"                    # string offset=1259
.Linfo_string83:
	.asciz	"S2"                    # string offset=1262
.Linfo_string84:
	.asciz	"dS2dx"                 # string offset=1265
.Linfo_string85:
	.asciz	"IMFs"                  # string offset=1271
.Linfo_string86:
	.asciz	"IMFS2"                 # string offset=1276
.Linfo_string87:
	.asciz	"IMFss"                 # string offset=1282
.Linfo_string88:
	.asciz	"IMFB"                  # string offset=1288
.Linfo_string89:
	.asciz	"IMFdS2dx"              # string offset=1293
.Linfo_string90:
	.asciz	"A"                     # string offset=1302
.Linfo_string91:
	.asciz	"IMFA"                  # string offset=1304
.Linfo_string92:
	.asciz	"delAy"                 # string offset=1309
.Linfo_string93:
	.asciz	"delAz"                 # string offset=1315
.Linfo_string94:
	.asciz	"deldS2dx"              # string offset=1321
.Linfo_string95:
	.asciz	"IMFdelAx"              # string offset=1330
.Linfo_string96:
	.asciz	"IMFdelAy"              # string offset=1339
.Linfo_string97:
	.asciz	"IMFdelAz"              # string offset=1348
.Linfo_string98:
	.asciz	"delB"                  # string offset=1357
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
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
	.quad	.Ltmp14-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	144                     # -112
	.byte	127                     # 
	.quad	.Ltmp14-.Lfunc_begin0
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
	.quad	.Ltmp15-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	104                     # -24
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Ltmp16-.Lfunc_begin0
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
	.quad	.Ltmp16-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	96                      # -32
	.quad	.Ltmp16-.Lfunc_begin0
	.quad	.Ltmp17-.Lfunc_begin0
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
	.quad	.Ltmp17-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	88                      # -40
	.quad	.Ltmp17-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc5:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	112                     # -16
	.quad	.Ltmp12-.Lfunc_begin0
	.quad	.Ltmp13-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc6:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	102                     # DW_OP_reg22
	.quad	.Ltmp10-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc7:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	80                      # -48
	.quad	.Ltmp18-.Lfunc_begin0
	.quad	.Ltmp19-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc8:
	.quad	.Lfunc_begin4-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	104                     # DW_OP_reg24
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp19-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	72                      # -56
	.quad	.Ltmp19-.Lfunc_begin0
	.quad	.Ltmp20-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc9:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	0
	.quad	0
.Ldebug_loc10:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	144                     # -112
	.byte	127                     # 
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Lfunc_end4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc11:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	104                     # -24
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Ltmp16-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc12:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp16-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	96                      # -32
	.quad	.Ltmp16-.Lfunc_begin0
	.quad	.Ltmp17-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc13:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp17-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	88                      # -40
	.quad	.Ltmp17-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc14:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	112                     # -16
	.quad	.Ltmp12-.Lfunc_begin0
	.quad	.Ltmp13-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc15:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	102                     # DW_OP_reg22
	.quad	.Ltmp10-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc16:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	80                      # -48
	.quad	.Ltmp18-.Lfunc_begin0
	.quad	.Ltmp19-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc17:
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	104                     # DW_OP_reg24
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp19-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	72                      # -56
	.quad	.Ltmp19-.Lfunc_begin0
	.quad	.Ltmp20-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc18:
	.quad	.Lfunc_begin5-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp24-.Lfunc_begin0
	.quad	.Ltmp32-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp41-.Lfunc_begin0
	.quad	.Ltmp43-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc19:
	.quad	.Lfunc_begin5-.Lfunc_begin0
	.quad	.Ltmp26-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc20:
	.quad	.Lfunc_begin5-.Lfunc_begin0
	.quad	.Ltmp28-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc21:
	.quad	.Ltmp27-.Lfunc_begin0
	.quad	.Ltmp64-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp71-.Lfunc_begin0
	.quad	.Ltmp120-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp124-.Lfunc_begin0
	.quad	.Ltmp126-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc22:
	.quad	.Ltmp31-.Lfunc_begin0
	.quad	.Ltmp33-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp33-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	64                      # -64
	.quad	.Ltmp43-.Lfunc_begin0
	.quad	.Ltmp127-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	64                      # -64
	.quad	0
	.quad	0
.Ldebug_loc23:
	.quad	.Ltmp34-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp43-.Lfunc_begin0
	.quad	.Ltmp46-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp46-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	240                     # -144
	.byte	126                     # 
	.quad	.Ltmp79-.Lfunc_begin0
	.quad	.Ltmp83-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp119-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	240                     # -144
	.byte	126                     # 
	.quad	.Ltmp119-.Lfunc_begin0
	.quad	.Ltmp121-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp124-.Lfunc_begin0
	.quad	.Ltmp127-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	0
	.quad	0
.Ldebug_loc24:
	.quad	.Ltmp35-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	144                     # DW_OP_regx
	.byte	32                      # 32
	.quad	.Ltmp43-.Lfunc_begin0
	.quad	.Ltmp77-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	144                     # DW_OP_regx
	.byte	32                      # 32
	.quad	.Ltmp79-.Lfunc_begin0
	.quad	.Ltmp83-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	144                     # DW_OP_regx
	.byte	32                      # 32
	.quad	.Ltmp124-.Lfunc_begin0
	.quad	.Ltmp127-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	144                     # DW_OP_regx
	.byte	32                      # 32
	.quad	0
	.quad	0
.Ldebug_loc25:
	.quad	.Ltmp36-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp43-.Lfunc_begin0
	.quad	.Ltmp45-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp45-.Lfunc_begin0
	.quad	.Ltmp68-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	120                     # -8
	.quad	.Ltmp68-.Lfunc_begin0
	.quad	.Ltmp70-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp71-.Lfunc_begin0
	.quad	.Ltmp75-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	120                     # -8
	.quad	.Ltmp75-.Lfunc_begin0
	.quad	.Ltmp76-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp79-.Lfunc_begin0
	.quad	.Ltmp81-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc26:
	.quad	.Ltmp47-.Lfunc_begin0
	.quad	.Ltmp51-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	184                     # -72
	.byte	127                     # 
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	184                     # -72
	.byte	127                     # 
	.quad	0
	.quad	0
.Ldebug_loc27:
	.quad	.Ltmp48-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	108                     # DW_OP_reg28
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp114-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	108                     # DW_OP_reg28
	.quad	0
	.quad	0
.Ldebug_loc28:
	.quad	.Ltmp49-.Lfunc_begin0
	.quad	.Ltmp50-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp50-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	104                     # -24
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	104                     # -24
	.quad	0
	.quad	0
.Ldebug_loc29:
	.quad	.Ltmp52-.Lfunc_begin0
	.quad	.Ltmp58-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	109                     # DW_OP_reg29
	.quad	.Ltmp58-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	208                     # -176
	.byte	126                     # 
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	208                     # -176
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc30:
	.quad	.Ltmp53-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	107                     # DW_OP_reg27
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp97-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	107                     # DW_OP_reg27
	.quad	0
	.quad	0
.Ldebug_loc31:
	.quad	.Ltmp67-.Lfunc_begin0
	.quad	.Ltmp71-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc32:
	.quad	.Ltmp54-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	110                     # DW_OP_reg30
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp85-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	110                     # DW_OP_reg30
	.quad	0
	.quad	0
.Ldebug_loc33:
	.quad	.Ltmp55-.Lfunc_begin0
	.quad	.Ltmp56-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	102                     # DW_OP_reg22
	.quad	.Ltmp56-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	112                     # -16
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	112                     # -16
	.quad	0
	.quad	0
.Ldebug_loc34:
	.quad	.Ltmp57-.Lfunc_begin0
	.quad	.Ltmp74-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp74-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	64                      # -64
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp109-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	64                      # -64
	.quad	.Ltmp109-.Lfunc_begin0
	.quad	.Ltmp118-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc35:
	.quad	.Ltmp86-.Lfunc_begin0
	.quad	.Ltmp88-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp88-.Lfunc_begin0
	.quad	.Ltmp89-.Lfunc_begin0
	.short	6                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	110                     # DW_OP_reg30
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp89-.Lfunc_begin0
	.quad	.Ltmp96-.Lfunc_begin0
	.short	9                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	110                     # DW_OP_reg30
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	105                     # DW_OP_reg25
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp96-.Lfunc_begin0
	.quad	.Ltmp100-.Lfunc_begin0
	.short	8                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	110                     # DW_OP_reg30
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	105                     # DW_OP_reg25
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp100-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	5                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	16                      # 16
	.byte	105                     # DW_OP_reg25
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	0
	.quad	0
.Ldebug_loc36:
	.quad	.Ltmp87-.Lfunc_begin0
	.quad	.Ltmp90-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp90-.Lfunc_begin0
	.quad	.Ltmp91-.Lfunc_begin0
	.short	6                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	111                     # DW_OP_reg31
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp91-.Lfunc_begin0
	.quad	.Ltmp94-.Lfunc_begin0
	.short	9                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	111                     # DW_OP_reg31
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	99                      # DW_OP_reg19
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp94-.Lfunc_begin0
	.quad	.Ltmp98-.Lfunc_begin0
	.short	6                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	111                     # DW_OP_reg31
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp98-.Lfunc_begin0
	.quad	.Ltmp99-.Lfunc_begin0
	.short	5                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	111                     # DW_OP_reg31
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	0
	.quad	0
.Ldebug_loc37:
	.quad	.Ltmp95-.Lfunc_begin0
	.quad	.Ltmp108-.Lfunc_begin0
	.short	9                       # Loc expr size
	.byte	102                     # DW_OP_reg22
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp108-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	8                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	0
	.quad	0
.Ldebug_loc38:
	.quad	.Ltmp102-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	0
	.quad	0
.Ldebug_loc39:
	.quad	.Ltmp103-.Lfunc_begin0
	.quad	.Ltmp104-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp104-.Lfunc_begin0
	.quad	.Ltmp106-.Lfunc_begin0
	.short	6                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp106-.Lfunc_begin0
	.quad	.Ltmp115-.Lfunc_begin0
	.short	9                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	101                     # DW_OP_reg21
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp115-.Lfunc_begin0
	.quad	.Ltmp117-.Lfunc_begin0
	.short	8                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	101                     # DW_OP_reg21
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp117-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	5                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	0
	.quad	0
.Ldebug_loc40:
	.quad	.Ltmp105-.Lfunc_begin0
	.quad	.Ltmp107-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp107-.Lfunc_begin0
	.quad	.Ltmp111-.Lfunc_begin0
	.short	9                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp111-.Lfunc_begin0
	.quad	.Ltmp116-.Lfunc_begin0
	.short	8                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	100                     # DW_OP_reg20
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	.Ltmp116-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	5                       # Loc expr size
	.byte	147                     # DW_OP_piece
	.byte	16                      # 16
	.byte	48                      # DW_OP_lit0
	.byte	147                     # DW_OP_piece
	.byte	8                       # 8
	.quad	0
	.quad	0
.Ldebug_loc41:
	.quad	.Ltmp122-.Lfunc_begin0
	.quad	.Ltmp123-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
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
	.byte	29                      # Abbreviation Code
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
	.byte	30                      # Abbreviation Code
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
	.byte	33                      # Abbreviation Code
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	85                      # DW_AT_ranges
	.byte	6                       # DW_FORM_data4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	34                      # Abbreviation Code
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
	.byte	35                      # Abbreviation Code
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
	.byte	1                       # Abbrev [1] 0xb:0x784 DW_TAG_compile_unit
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
	.quad	_ZTV12VectorDipole
	.byte	9                       # Abbrev [9] 0x99:0x17 DW_TAG_variable
	.long	.Linfo_string8          # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	1                       # DW_AT_decl_file
	.short	8117                    # DW_AT_decl_line
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	__I___32_backgroundfield_vectordipole_cpp_83390dc3
	.byte	10                      # Abbrev [10] 0xb0:0x13 DW_TAG_variable
	.long	.Linfo_string9          # DW_AT_name
	.long	195                     # DW_AT_type
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE
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
	.quad	_ZTI12VectorDipole
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
	.quad	_ZTS12VectorDipole
	.byte	3                       # Abbrev [3] 0x1da:0xc DW_TAG_array_type
	.long	224                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x1df:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	15                      # DW_AT_count
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
	.long	.Linfo_string46         # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	31                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x22a:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1901                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x232:0x1a DW_TAG_lexical_block
	.quad	.Ltmp0                  # DW_AT_low_pc
	.quad	.Ltmp1                  # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x243:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1901                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	22                      # Abbrev [22] 0x24e:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin1           # DW_AT_low_pc
	.quad	.Lfunc_end1             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string47         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	23                      # Abbrev [23] 0x266:0x3b DW_TAG_subprogram
	.quad	.Lfunc_begin2           # DW_AT_low_pc
	.quad	.Lfunc_end2             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string48         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x27e:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1906                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x286:0x1a DW_TAG_lexical_block
	.quad	.Ltmp4                  # DW_AT_low_pc
	.quad	.Ltmp5                  # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x297:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1906                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	22                      # Abbrev [22] 0x2a1:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin3           # DW_AT_low_pc
	.quad	.Lfunc_end3             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string49         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	17                      # Abbrev [17] 0x2b9:0x408 DW_TAG_class_type
	.long	.Linfo_string45         # DW_AT_name
	.byte	112                     # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	24                      # Abbrev [24] 0x2c3:0xf DW_TAG_inheritance
	.long	.Linfo_string32         # DW_AT_name
	.long	1729                    # DW_AT_type
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
	.long	1829                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	25                      # Abbrev [25] 0x2f0:0xf DW_TAG_member
	.long	.Linfo_string42         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	25                      # Abbrev [25] 0x2ff:0xf DW_TAG_member
	.long	.Linfo_string43         # DW_AT_name
	.long	1841                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	72
	.byte	25                      # Abbrev [25] 0x30e:0xf DW_TAG_member
	.long	.Linfo_string44         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	88
	.byte	18                      # Abbrev [18] 0x31d:0x159 DW_TAG_subprogram
	.quad	.Lfunc_begin4           # DW_AT_low_pc
	.quad	.Lfunc_end4             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string50         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	34                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	26                      # Abbrev [26] 0x337:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc0            # DW_AT_location
	.long	1911                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0x341:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc1            # DW_AT_location
	.long	.Linfo_string55         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x34e:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc2            # DW_AT_location
	.long	.Linfo_string56         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x35b:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc3            # DW_AT_location
	.long	.Linfo_string57         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x368:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc4            # DW_AT_location
	.long	.Linfo_string58         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x375:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc5            # DW_AT_location
	.long	.Linfo_string59         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x382:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc6            # DW_AT_location
	.long	.Linfo_string60         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x38f:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc7            # DW_AT_location
	.long	.Linfo_string61         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x39c:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc8            # DW_AT_location
	.long	.Linfo_string62         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	28                      # Abbrev [28] 0x3a9:0xc DW_TAG_formal_parameter
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	16
	.long	.Linfo_string65         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	28                      # Abbrev [28] 0x3b5:0xc DW_TAG_formal_parameter
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	24
	.long	.Linfo_string64         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	28                      # Abbrev [28] 0x3c1:0xc DW_TAG_formal_parameter
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	32
	.long	.Linfo_string63         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	20                      # Abbrev [20] 0x3cd:0xa8 DW_TAG_lexical_block
	.quad	.Ltmp9                  # DW_AT_low_pc
	.quad	.Ltmp22                 # DW_AT_high_pc
	.byte	29                      # Abbrev [29] 0x3de:0xa DW_TAG_variable
	.long	.Ldebug_loc9            # DW_AT_location
	.long	1911                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	30                      # Abbrev [30] 0x3e8:0xd DW_TAG_variable
	.long	.Ldebug_loc10           # DW_AT_location
	.long	.Linfo_string55         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x3f5:0xd DW_TAG_variable
	.long	.Ldebug_loc11           # DW_AT_location
	.long	.Linfo_string56         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x402:0xd DW_TAG_variable
	.long	.Ldebug_loc12           # DW_AT_location
	.long	.Linfo_string57         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x40f:0xd DW_TAG_variable
	.long	.Ldebug_loc13           # DW_AT_location
	.long	.Linfo_string58         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x41c:0xd DW_TAG_variable
	.long	.Ldebug_loc14           # DW_AT_location
	.long	.Linfo_string59         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x429:0xd DW_TAG_variable
	.long	.Ldebug_loc15           # DW_AT_location
	.long	.Linfo_string60         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x436:0xd DW_TAG_variable
	.long	.Ldebug_loc16           # DW_AT_location
	.long	.Linfo_string61         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x443:0xd DW_TAG_variable
	.long	.Ldebug_loc17           # DW_AT_location
	.long	.Linfo_string62         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x450:0xc DW_TAG_variable
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	32
	.long	.Linfo_string63         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x45c:0xc DW_TAG_variable
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	24
	.long	.Linfo_string64         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x468:0xc DW_TAG_variable
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	16
	.long	.Linfo_string65         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	32                      # Abbrev [32] 0x476:0x20d DW_TAG_subprogram
	.quad	.Lfunc_begin5           # DW_AT_low_pc
	.quad	.Lfunc_end5             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string51         # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	60                      # DW_AT_decl_line
	.long	511                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x494:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1911                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0x49c:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc18           # DW_AT_location
	.long	.Linfo_string73         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x4a9:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc19           # DW_AT_location
	.long	.Linfo_string74         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	27                      # Abbrev [27] 0x4b6:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc20           # DW_AT_location
	.long	.Linfo_string75         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x4c3:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\220\177"
	.long	.Linfo_string66         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x4d0:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\330|"
	.long	.Linfo_string67         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x4dd:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\360|"
	.long	.Linfo_string68         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	20                      # Abbrev [20] 0x4ea:0x198 DW_TAG_lexical_block
	.quad	.Ltmp24                 # DW_AT_low_pc
	.quad	.Ltmp128                # DW_AT_high_pc
	.byte	16                      # Abbrev [16] 0x4fb:0x9 DW_TAG_variable
	.long	.Linfo_string75         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	16                      # Abbrev [16] 0x504:0x9 DW_TAG_variable
	.long	.Linfo_string74         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	16                      # Abbrev [16] 0x50d:0x9 DW_TAG_variable
	.long	.Linfo_string73         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	21                      # Abbrev [21] 0x516:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1911                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	30                      # Abbrev [30] 0x51e:0xd DW_TAG_variable
	.long	.Ldebug_loc21           # DW_AT_location
	.long	.Linfo_string76         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x52b:0xd DW_TAG_variable
	.long	.Ldebug_loc22           # DW_AT_location
	.long	.Linfo_string77         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x538:0xd DW_TAG_variable
	.long	.Ldebug_loc23           # DW_AT_location
	.long	.Linfo_string78         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x545:0xd DW_TAG_variable
	.long	.Ldebug_loc24           # DW_AT_location
	.long	.Linfo_string79         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x552:0xd DW_TAG_variable
	.long	.Ldebug_loc25           # DW_AT_location
	.long	.Linfo_string80         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x55f:0xd DW_TAG_variable
	.long	.Ldebug_loc26           # DW_AT_location
	.long	.Linfo_string81         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x56c:0xd DW_TAG_variable
	.long	.Ldebug_loc27           # DW_AT_location
	.long	.Linfo_string82         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x579:0xd DW_TAG_variable
	.long	.Ldebug_loc28           # DW_AT_location
	.long	.Linfo_string83         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x586:0xd DW_TAG_variable
	.long	.Ldebug_loc29           # DW_AT_location
	.long	.Linfo_string84         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x593:0xd DW_TAG_variable
	.long	.Ldebug_loc30           # DW_AT_location
	.long	.Linfo_string85         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x5a0:0xd DW_TAG_variable
	.long	.Ldebug_loc31           # DW_AT_location
	.long	.Linfo_string86         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x5ad:0xd DW_TAG_variable
	.long	.Ldebug_loc32           # DW_AT_location
	.long	.Linfo_string87         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x5ba:0xd DW_TAG_variable
	.long	.Ldebug_loc33           # DW_AT_location
	.long	.Linfo_string88         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	30                      # Abbrev [30] 0x5c7:0xd DW_TAG_variable
	.long	.Ldebug_loc34           # DW_AT_location
	.long	.Linfo_string89         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	16                      # Abbrev [16] 0x5d4:0x9 DW_TAG_variable
	.long	.Linfo_string90         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	16                      # Abbrev [16] 0x5dd:0x9 DW_TAG_variable
	.long	.Linfo_string91         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	33                      # Abbrev [33] 0x5e6:0x20 DW_TAG_lexical_block
	.long	.Ldebug_ranges0         # DW_AT_ranges
	.byte	31                      # Abbrev [31] 0x5eb:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\210}"
	.long	.Linfo_string69         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x5f8:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\250~"
	.long	.Linfo_string70         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	33                      # Abbrev [33] 0x606:0x7b DW_TAG_lexical_block
	.long	.Ldebug_ranges1         # DW_AT_ranges
	.byte	31                      # Abbrev [31] 0x60b:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\300{"
	.long	.Linfo_string71         # DW_AT_name
	.long	1916                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x618:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\220|"
	.long	.Linfo_string72         # DW_AT_name
	.long	1916                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x625:0xd DW_TAG_variable
	.long	.Ldebug_loc35           # DW_AT_location
	.long	.Linfo_string92         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x632:0xd DW_TAG_variable
	.long	.Ldebug_loc36           # DW_AT_location
	.long	.Linfo_string93         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x63f:0xd DW_TAG_variable
	.long	.Ldebug_loc37           # DW_AT_location
	.long	.Linfo_string94         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x64c:0xd DW_TAG_variable
	.long	.Ldebug_loc38           # DW_AT_location
	.long	.Linfo_string95         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x659:0xd DW_TAG_variable
	.long	.Ldebug_loc39           # DW_AT_location
	.long	.Linfo_string96         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x666:0xd DW_TAG_variable
	.long	.Ldebug_loc40           # DW_AT_location
	.long	.Linfo_string97         # DW_AT_name
	.long	1829                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x673:0xd DW_TAG_variable
	.long	.Ldebug_loc41           # DW_AT_location
	.long	.Linfo_string98         # DW_AT_name
	.long	511                     # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	18                      # Abbrev [18] 0x683:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin6           # DW_AT_low_pc
	.quad	.Lfunc_end6             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string52         # DW_AT_name
	.byte	4                       # DW_AT_decl_file
	.byte	47                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	19                      # Abbrev [19] 0x69d:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1911                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	20                      # Abbrev [20] 0x6a5:0x1a DW_TAG_lexical_block
	.quad	.Ltmp129                # DW_AT_low_pc
	.quad	.Ltmp130                # DW_AT_high_pc
	.byte	21                      # Abbrev [21] 0x6b6:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	1911                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	17                      # Abbrev [17] 0x6c1:0x47 DW_TAG_class_type
	.long	.Linfo_string32         # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	24                      # Abbrev [24] 0x6cb:0xf DW_TAG_inheritance
	.long	.Linfo_string31         # DW_AT_name
	.long	518                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	25                      # Abbrev [25] 0x6da:0xf DW_TAG_member
	.long	.Linfo_string33         # DW_AT_name
	.long	1800                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	25                      # Abbrev [25] 0x6e9:0xf DW_TAG_member
	.long	.Linfo_string38         # DW_AT_name
	.long	1800                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	12
	.byte	25                      # Abbrev [25] 0x6f8:0xf DW_TAG_member
	.long	.Linfo_string39         # DW_AT_name
	.long	495                     # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	0                       # End Of Children Mark
	.byte	34                      # Abbrev [34] 0x708:0x1d DW_TAG_enumeration_type
	.long	.Linfo_string37         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	4                       # DW_AT_alignment
	.byte	35                      # Abbrev [35] 0x712:0x6 DW_TAG_enumerator
	.long	.Linfo_string34         # DW_AT_name
	.byte	0                       # DW_AT_const_value
	.byte	35                      # Abbrev [35] 0x718:0x6 DW_TAG_enumerator
	.long	.Linfo_string35         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	35                      # Abbrev [35] 0x71e:0x6 DW_TAG_enumerator
	.long	.Linfo_string36         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x725:0xc DW_TAG_array_type
	.long	511                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x72a:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	3                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x731:0xc DW_TAG_array_type
	.long	511                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x736:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	2                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	22                      # Abbrev [22] 0x73d:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin7           # DW_AT_low_pc
	.quad	.Lfunc_end7             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string53         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	22                      # Abbrev [22] 0x755:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin8           # DW_AT_low_pc
	.quad	.Lfunc_end8             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string54         # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	5                       # Abbrev [5] 0x76d:0x5 DW_TAG_pointer_type
	.long	518                     # DW_AT_type
	.byte	5                       # Abbrev [5] 0x772:0x5 DW_TAG_pointer_type
	.long	1729                    # DW_AT_type
	.byte	5                       # Abbrev [5] 0x777:0x5 DW_TAG_pointer_type
	.long	697                     # DW_AT_type
	.byte	3                       # Abbrev [3] 0x77c:0x12 DW_TAG_array_type
	.long	511                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x781:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	3                       # DW_AT_count
	.byte	4                       # Abbrev [4] 0x787:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	3                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
.Ldebug_info_end0:
	.section	.debug_ranges,"",@progbits
.Ldebug_ranges0:
	.quad	.Ltmp61-.Lfunc_begin0
	.quad	.Ltmp63-.Lfunc_begin0
	.quad	.Ltmp65-.Lfunc_begin0
	.quad	.Ltmp66-.Lfunc_begin0
	.quad	.Ltmp69-.Lfunc_begin0
	.quad	.Ltmp72-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges1:
	.quad	.Ltmp72-.Lfunc_begin0
	.quad	.Ltmp79-.Lfunc_begin0
	.quad	.Ltmp85-.Lfunc_begin0
	.quad	.Ltmp91-.Lfunc_begin0
	.quad	.Ltmp92-.Lfunc_begin0
	.quad	.Ltmp110-.Lfunc_begin0
	.quad	.Ltmp112-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.quad	0
	.quad	0
	.section	.debug_pubnames,"",@progbits
	.long	.LpubNames_end0-.LpubNames_begin0 # Length of Public Names Info
.LpubNames_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	1935                    # Compilation Unit Length
	.long	1142                    # DIE offset
	.asciz	"VectorDipole::call"    # External Name
	.long	311                     # DIE offset
	.asciz	"WID"                   # External Name
	.long	1877                    # DIE offset
	.asciz	"__sti___32_backgroundfield_vectordipole_cpp_83390dc3" # External Name
	.long	251                     # DIE offset
	.asciz	"_ZTI11T3DFunction"     # External Name
	.long	320                     # DIE offset
	.asciz	"_ZTI12VectorDipole"    # External Name
	.long	590                     # DIE offset
	.asciz	"_ZN11T3DFunctionD0Ev"  # External Name
	.long	231                     # DIE offset
	.asciz	"__dso_handle"          # External Name
	.long	113                     # DIE offset
	.asciz	"_ZTV13FieldFunction"   # External Name
	.long	1667                    # DIE offset
	.asciz	"VectorDipole::~VectorDipole" # External Name
	.long	380                     # DIE offset
	.asciz	"_ZTS11T3DFunction"     # External Name
	.long	528                     # DIE offset
	.asciz	"T3DFunction::~T3DFunction" # External Name
	.long	454                     # DIE offset
	.asciz	"_ZTS12VectorDipole"    # External Name
	.long	133                     # DIE offset
	.asciz	"_ZTV12VectorDipole"    # External Name
	.long	412                     # DIE offset
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # External Name
	.long	1816                    # DIE offset
	.asciz	"Y"                     # External Name
	.long	1822                    # DIE offset
	.asciz	"Z"                     # External Name
	.long	1810                    # DIE offset
	.asciz	"X"                     # External Name
	.long	797                     # DIE offset
	.asciz	"VectorDipole::initialize" # External Name
	.long	153                     # DIE offset
	.asciz	"__I___32_backgroundfield_vectordipole_cpp_83390dc3" # External Name
	.long	358                     # DIE offset
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # External Name
	.long	281                     # DIE offset
	.asciz	"_ZTI13FieldFunction"   # External Name
	.long	486                     # DIE offset
	.asciz	"_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc35vmesh15INVALID_LOCALIDE" # External Name
	.long	673                     # DIE offset
	.asciz	"_ZN13FieldFunctionD0Ev" # External Name
	.long	614                     # DIE offset
	.asciz	"~FieldFunction"        # External Name
	.long	1853                    # DIE offset
	.asciz	"_ZN12VectorDipoleD0Ev" # External Name
	.long	422                     # DIE offset
	.asciz	"_ZTS13FieldFunction"   # External Name
	.long	176                     # DIE offset
	.asciz	"_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc3St8__ioinitE" # External Name
	.long	502                     # DIE offset
	.asciz	"_ZN54_INTERNAL_32_backgroundfield_vectordipole_cpp_83390dc317physicalconstants3R_EE" # External Name
	.long	47                      # DIE offset
	.asciz	"_ZTV11T3DFunction"     # External Name
	.long	349                     # DIE offset
	.asciz	"WID3"                  # External Name
	.long	340                     # DIE offset
	.asciz	"WID2"                  # External Name
	.long	0                       # End Mark
.LpubNames_end0:
	.section	.debug_pubtypes,"",@progbits
	.long	.LpubTypes_end0-.LpubTypes_begin0 # Length of Public Types Info
.LpubTypes_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	1935                    # Compilation Unit Length
	.long	271                     # DIE offset
	.asciz	"__class_type_info"     # External Name
	.long	518                     # DIE offset
	.asciz	"T3DFunction"           # External Name
	.long	246                     # DIE offset
	.asciz	"void"                  # External Name
	.long	697                     # DIE offset
	.asciz	"VectorDipole"          # External Name
	.long	1800                    # DIE offset
	.asciz	"coordinate"            # External Name
	.long	195                     # DIE offset
	.asciz	"_ZNSt8ios_base4InitE"  # External Name
	.long	224                     # DIE offset
	.asciz	"signed char"           # External Name
	.long	1729                    # DIE offset
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
