	.text
	.file	"integratefunction.ll"
	.file	1 "/home/talgat/vlasiator/backgroundfield/integratefunction.cpp"
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _Z11lineAverageRK11T3DFunction10coordinatedPKdd
.LCPI0_0:
	.quad	4607182418800017408     # double 1
	.text
	.globl	_Z11lineAverageRK11T3DFunction10coordinatedPKdd
	.p2align	4, 0x90
	.type	_Z11lineAverageRK11T3DFunction10coordinatedPKdd,@function
_Z11lineAverageRK11T3DFunction10coordinatedPKdd: # @_Z11lineAverageRK11T3DFunction10coordinatedPKdd
.Lfunc_begin0:
	.loc	1 43 0                  # backgroundfield/integratefunction.cpp:43:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: lineAverage:f1 <- $rdi
	#DEBUG_VALUE: lineAverage:line <- $esi
	#DEBUG_VALUE: lineAverage:accuracy <- $xmm0
	#DEBUG_VALUE: lineAverage:r1 <- $rdx
	#DEBUG_VALUE: lineAverage:L <- $xmm1
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$112, %rsp
.Ltmp0:
	#DEBUG_VALUE: f1 <- undef
	vmovsd	.LCPI0_0(%rip), %xmm2   # xmm2 = mem[0],zero
.Ltmp1:
	#DEBUG_VALUE: L <- $xmm1
	#DEBUG_VALUE: r1 <- $rdx
	#DEBUG_VALUE: accuracy <- $xmm0
	#DEBUG_VALUE: line <- $esi
	#DEBUG_VALUE: lineAverage:L <- $xmm1
	#DEBUG_VALUE: lineAverage:r1 <- $rdx
	#DEBUG_VALUE: lineAverage:accuracy <- $xmm0
	#DEBUG_VALUE: lineAverage:line <- $esi
	#DEBUG_VALUE: lineAverage:f1 <- $rdi
	.loc	1 48 1 prologue_end     # backgroundfield/integratefunction.cpp:48:1
	vdivsd	%xmm1, %xmm2, %xmm3
.Ltmp2:
	#DEBUG_VALUE: norm <- $xmm3
	.loc	1 49 1                  # backgroundfield/integratefunction.cpp:49:1
	vmulsd	%xmm1, %xmm0, %xmm2
.Ltmp3:
	#DEBUG_VALUE: acc <- $xmm2
	.loc	1 50 1                  # backgroundfield/integratefunction.cpp:50:1
	movl	%esi, %eax
	vmovsd	(%rdx,%rax,8), %xmm0    # xmm0 = mem[0],zero
.Ltmp4:
	#DEBUG_VALUE: a <- $xmm0
	.loc	1 51 1                  # backgroundfield/integratefunction.cpp:51:1
	vaddsd	%xmm1, %xmm0, %xmm1
.Ltmp5:
	#DEBUG_VALUE: b <- $xmm1
	.loc	1 53 1                  # backgroundfield/integratefunction.cpp:53:1
	testl	%esi, %esi
	je	.LBB0_4
.Ltmp6:
# %bb.1:                                # %L.entry
	#DEBUG_VALUE: b <- $xmm1
	#DEBUG_VALUE: a <- $xmm0
	#DEBUG_VALUE: acc <- $xmm2
	#DEBUG_VALUE: norm <- $xmm3
	#DEBUG_VALUE: line <- $esi
	#DEBUG_VALUE: r1 <- $rdx
	#DEBUG_VALUE: lineAverage:r1 <- $rdx
	#DEBUG_VALUE: lineAverage:line <- $esi
	#DEBUG_VALUE: lineAverage:f1 <- $rdi
	cmpl	$1, %esi
	je	.LBB0_5
.Ltmp7:
# %bb.2:                                # %L.entry
	#DEBUG_VALUE: b <- $xmm1
	#DEBUG_VALUE: a <- $xmm0
	#DEBUG_VALUE: acc <- $xmm2
	#DEBUG_VALUE: norm <- $xmm3
	#DEBUG_VALUE: line <- $esi
	#DEBUG_VALUE: r1 <- $rdx
	#DEBUG_VALUE: lineAverage:r1 <- $rdx
	#DEBUG_VALUE: lineAverage:line <- $esi
	#DEBUG_VALUE: lineAverage:f1 <- $rdi
	cmpl	$2, %esi
	jne	.LBB0_6
.Ltmp8:
# %bb.3:                                # %L.B0008
	#DEBUG_VALUE: b <- $xmm1
	#DEBUG_VALUE: a <- $xmm0
	#DEBUG_VALUE: acc <- $xmm2
	#DEBUG_VALUE: norm <- $xmm3
	#DEBUG_VALUE: line <- $esi
	#DEBUG_VALUE: r1 <- $rdx
	#DEBUG_VALUE: lineAverage:r1 <- $rdx
	#DEBUG_VALUE: lineAverage:line <- $esi
	#DEBUG_VALUE: lineAverage:f1 <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/integratefunction.cpp:0:1
	vmovsd	%xmm3, -8(%rbp)         # 8-byte Spill
.Ltmp9:
	#DEBUG_VALUE: norm <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	.loc	1 68 1 is_stmt 1        # backgroundfield/integratefunction.cpp:68:1
	vmovupd	(%rdx), %xmm3
.Ltmp10:
	.loc	1 94 1                  # backgroundfield/integratefunction.cpp:94:1
	movq	$_ZTV9T3D_fix12+16, -104(%rbp)
	movq	%rdi, -96(%rbp)
	vmovupd	%xmm3, -88(%rbp)
	leaq	-104(%rbp), %rdi
.Ltmp11:
	.loc	1 69 1                  # backgroundfield/integratefunction.cpp:69:1
	callq	_Z7RombergRK11T1DFunctionddd
.Ltmp12:
	vmulsd	-8(%rbp), %xmm0, %xmm0  # 8-byte Folded Reload
.Ltmp13:
	#DEBUG_VALUE: value <- $xmm0
	.loc	1 96 1                  # backgroundfield/integratefunction.cpp:96:1
	movq	$_ZTV11T1DFunction+16, -104(%rbp)
	jmp	.LBB0_7
.Ltmp14:
.LBB0_4:                                # %L.B0002
	#DEBUG_VALUE: b <- $xmm1
	#DEBUG_VALUE: a <- $xmm0
	#DEBUG_VALUE: acc <- $xmm2
	#DEBUG_VALUE: norm <- $xmm3
	#DEBUG_VALUE: line <- $esi
	#DEBUG_VALUE: r1 <- $rdx
	#DEBUG_VALUE: lineAverage:r1 <- $rdx
	#DEBUG_VALUE: lineAverage:line <- $esi
	#DEBUG_VALUE: lineAverage:f1 <- $rdi
	.loc	1 0 1 is_stmt 0         # backgroundfield/integratefunction.cpp:0:1
	vmovsd	%xmm3, -8(%rbp)         # 8-byte Spill
.Ltmp15:
	#DEBUG_VALUE: norm <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	.loc	1 56 1 is_stmt 1        # backgroundfield/integratefunction.cpp:56:1
	vmovupd	8(%rdx), %xmm3
.Ltmp16:
	.loc	1 114 1                 # backgroundfield/integratefunction.cpp:114:1
	movq	$_ZTV9T3D_fix23+16, -72(%rbp)
	movq	%rdi, -64(%rbp)
	vmovupd	%xmm3, -56(%rbp)
	leaq	-72(%rbp), %rdi
.Ltmp17:
	.loc	1 57 1                  # backgroundfield/integratefunction.cpp:57:1
	callq	_Z7RombergRK11T1DFunctionddd
.Ltmp18:
	vmulsd	-8(%rbp), %xmm0, %xmm0  # 8-byte Folded Reload
.Ltmp19:
	#DEBUG_VALUE: value <- $xmm0
	.loc	1 116 1                 # backgroundfield/integratefunction.cpp:116:1
	movq	$_ZTV11T1DFunction+16, -72(%rbp)
	jmp	.LBB0_7
.Ltmp20:
.LBB0_5:                                # %L.B0005
	#DEBUG_VALUE: b <- $xmm1
	#DEBUG_VALUE: a <- $xmm0
	#DEBUG_VALUE: acc <- $xmm2
	#DEBUG_VALUE: norm <- $xmm3
	#DEBUG_VALUE: line <- $esi
	#DEBUG_VALUE: r1 <- $rdx
	#DEBUG_VALUE: lineAverage:r1 <- $rdx
	#DEBUG_VALUE: lineAverage:line <- $esi
	#DEBUG_VALUE: lineAverage:f1 <- $rdi
	.loc	1 62 1                  # backgroundfield/integratefunction.cpp:62:1
	movq	(%rdx), %rax
	movq	16(%rdx), %rcx
.Ltmp21:
	.loc	1 104 1                 # backgroundfield/integratefunction.cpp:104:1
	movq	$_ZTV9T3D_fix13+16, -40(%rbp)
	movq	%rdi, -32(%rbp)
	movq	%rax, -24(%rbp)
	movq	%rcx, -16(%rbp)
	leaq	-40(%rbp), %rdi
.Ltmp22:
	.loc	1 0 1 is_stmt 0         # backgroundfield/integratefunction.cpp:0:1
	vmovsd	%xmm3, -8(%rbp)         # 8-byte Spill
.Ltmp23:
	#DEBUG_VALUE: norm <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	.loc	1 63 1 is_stmt 1        # backgroundfield/integratefunction.cpp:63:1
	callq	_Z7RombergRK11T1DFunctionddd
.Ltmp24:
	vmulsd	-8(%rbp), %xmm0, %xmm0  # 8-byte Folded Reload
.Ltmp25:
	#DEBUG_VALUE: value <- $xmm0
	.loc	1 106 1                 # backgroundfield/integratefunction.cpp:106:1
	movq	$_ZTV11T1DFunction+16, -40(%rbp)
	jmp	.LBB0_7
.Ltmp26:
.LBB0_6:                                # %L..inline.9786
	#DEBUG_VALUE: line <- $esi
	#DEBUG_VALUE: r1 <- $rdx
	#DEBUG_VALUE: lineAverage:r1 <- $rdx
	#DEBUG_VALUE: lineAverage:line <- $esi
	#DEBUG_VALUE: lineAverage:f1 <- $rdi
	.loc	1 570 1                 # backgroundfield/integratefunction.cpp:570:1
	movl	$_ZSt4cerr, %edi
.Ltmp27:
	movl	$.S08448, %esi
.Ltmp28:
	movl	$24, %edx
.Ltmp29:
	callq	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l
	vxorpd	%xmm0, %xmm0, %xmm0
.Ltmp30:
	#DEBUG_VALUE: value <- 0.000000e+00
.LBB0_7:                                # %L_T31759352_8084
	#DEBUG_VALUE: value <- $xmm0
	.loc	1 78 1                  # backgroundfield/integratefunction.cpp:78:1
	addq	$112, %rsp
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp31:
.Lfunc_end0:
	.size	_Z11lineAverageRK11T3DFunction10coordinatedPKdd, .Lfunc_end0-_Z11lineAverageRK11T3DFunction10coordinatedPKdd
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _Z14surfaceAverageRK11T3DFunction10coordinatedPKddd
.LCPI1_0:
	.quad	4607182418800017408     # double 1
	.text
	.globl	_Z14surfaceAverageRK11T3DFunction10coordinatedPKddd
	.p2align	4, 0x90
	.type	_Z14surfaceAverageRK11T3DFunction10coordinatedPKddd,@function
_Z14surfaceAverageRK11T3DFunction10coordinatedPKddd: # @_Z14surfaceAverageRK11T3DFunction10coordinatedPKddd
.Lfunc_begin1:
	.loc	1 88 0                  # backgroundfield/integratefunction.cpp:88:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: surfaceAverage:f1 <- $rdi
	#DEBUG_VALUE: surfaceAverage:face <- $esi
	#DEBUG_VALUE: surfaceAverage:accuracy <- $xmm0
	#DEBUG_VALUE: surfaceAverage:r1 <- $rdx
	#DEBUG_VALUE: surfaceAverage:L1 <- $xmm1
	#DEBUG_VALUE: surfaceAverage:L2 <- $xmm2
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$96, %rsp
	vmovapd	%xmm2, %xmm3
.Ltmp32:
	#DEBUG_VALUE: L2 <- $xmm3
	#DEBUG_VALUE: L1 <- $xmm1
	#DEBUG_VALUE: r1 <- undef
	#DEBUG_VALUE: accuracy <- $xmm0
	#DEBUG_VALUE: face <- $esi
	#DEBUG_VALUE: f1 <- undef
	#DEBUG_VALUE: surfaceAverage:L2 <- $xmm3
	#DEBUG_VALUE: surfaceAverage:L1 <- $xmm1
	#DEBUG_VALUE: surfaceAverage:r1 <- $rdx
	#DEBUG_VALUE: surfaceAverage:accuracy <- $xmm0
	#DEBUG_VALUE: surfaceAverage:face <- $esi
	#DEBUG_VALUE: surfaceAverage:f1 <- $rdi
	.loc	1 92 1 prologue_end     # backgroundfield/integratefunction.cpp:92:1
	vmulsd	%xmm1, %xmm0, %xmm0
.Ltmp33:
	vmulsd	%xmm2, %xmm0, %xmm4
.Ltmp34:
	#DEBUG_VALUE: norm <- undef
	#DEBUG_VALUE: acc <- $xmm4
	.loc	1 94 1                  # backgroundfield/integratefunction.cpp:94:1
	testl	%esi, %esi
	vmovsd	%xmm1, -16(%rbp)        # 8-byte Spill
	vmovsd	%xmm2, -8(%rbp)         # 8-byte Spill
	je	.LBB1_6
.Ltmp35:
# %bb.1:                                # %L.entry
	#DEBUG_VALUE: acc <- $xmm4
	#DEBUG_VALUE: surfaceAverage:L2 <- $xmm3
	#DEBUG_VALUE: face <- $esi
	#DEBUG_VALUE: L1 <- $xmm1
	#DEBUG_VALUE: L2 <- $xmm3
	#DEBUG_VALUE: surfaceAverage:L1 <- $xmm1
	#DEBUG_VALUE: surfaceAverage:r1 <- $rdx
	#DEBUG_VALUE: surfaceAverage:face <- $esi
	#DEBUG_VALUE: surfaceAverage:f1 <- $rdi
	cmpl	$1, %esi
	je	.LBB1_4
.Ltmp36:
# %bb.2:                                # %L.entry
	#DEBUG_VALUE: acc <- $xmm4
	#DEBUG_VALUE: surfaceAverage:L2 <- $xmm3
	#DEBUG_VALUE: face <- $esi
	#DEBUG_VALUE: L1 <- $xmm1
	#DEBUG_VALUE: L2 <- $xmm3
	#DEBUG_VALUE: surfaceAverage:L1 <- $xmm1
	#DEBUG_VALUE: surfaceAverage:r1 <- $rdx
	#DEBUG_VALUE: surfaceAverage:face <- $esi
	#DEBUG_VALUE: surfaceAverage:f1 <- $rdi
	cmpl	$2, %esi
	jne	.LBB1_5
.Ltmp37:
# %bb.3:                                # %L.B0023
	#DEBUG_VALUE: acc <- $xmm4
	#DEBUG_VALUE: surfaceAverage:L2 <- $xmm3
	#DEBUG_VALUE: face <- $esi
	#DEBUG_VALUE: L1 <- $xmm1
	#DEBUG_VALUE: L2 <- $xmm3
	#DEBUG_VALUE: surfaceAverage:L1 <- $xmm1
	#DEBUG_VALUE: surfaceAverage:r1 <- $rdx
	#DEBUG_VALUE: surfaceAverage:face <- $esi
	#DEBUG_VALUE: surfaceAverage:f1 <- $rdi
	.loc	1 109 1                 # backgroundfield/integratefunction.cpp:109:1
	movq	16(%rdx), %rax
.Ltmp38:
	.loc	1 82 1                  # backgroundfield/integratefunction.cpp:82:1
	movq	$_ZTV8T3D_fix3+16, -88(%rbp)
	movq	%rdi, -80(%rbp)
	movq	%rax, -72(%rbp)
.Ltmp39:
	.loc	1 110 1                 # backgroundfield/integratefunction.cpp:110:1
	vmovsd	(%rdx), %xmm0           # xmm0 = mem[0],zero
	vmovsd	8(%rdx), %xmm2          # xmm2 = mem[0],zero
	vaddsd	%xmm1, %xmm0, %xmm1
.Ltmp40:
	vaddsd	%xmm3, %xmm2, %xmm3
.Ltmp41:
	.loc	1 0 1 is_stmt 0         # backgroundfield/integratefunction.cpp:0:1
	leaq	-88(%rbp), %rdi
.Ltmp42:
	.loc	1 110 1                 # backgroundfield/integratefunction.cpp:110:1
	callq	_Z7RombergRK11T2DFunctionddddd
.Ltmp43:
	#DEBUG_VALUE: value <- undef
	.loc	1 84 1 is_stmt 1        # backgroundfield/integratefunction.cpp:84:1
	movq	$_ZTV11T2DFunction+16, -88(%rbp)
	jmp	.LBB1_7
.Ltmp44:
.LBB1_4:                                # %L.B0020
	#DEBUG_VALUE: acc <- $xmm4
	#DEBUG_VALUE: surfaceAverage:L2 <- $xmm3
	#DEBUG_VALUE: face <- $esi
	#DEBUG_VALUE: L1 <- $xmm1
	#DEBUG_VALUE: L2 <- $xmm3
	#DEBUG_VALUE: surfaceAverage:L1 <- $xmm1
	#DEBUG_VALUE: surfaceAverage:r1 <- $rdx
	#DEBUG_VALUE: surfaceAverage:face <- $esi
	#DEBUG_VALUE: surfaceAverage:f1 <- $rdi
	.loc	1 103 1                 # backgroundfield/integratefunction.cpp:103:1
	movq	8(%rdx), %rax
.Ltmp45:
	.loc	1 72 1                  # backgroundfield/integratefunction.cpp:72:1
	movq	$_ZTV8T3D_fix2+16, -40(%rbp)
	movq	%rdi, -32(%rbp)
	movq	%rax, -24(%rbp)
.Ltmp46:
	.loc	1 104 1                 # backgroundfield/integratefunction.cpp:104:1
	vmovsd	(%rdx), %xmm0           # xmm0 = mem[0],zero
	vmovsd	16(%rdx), %xmm2         # xmm2 = mem[0],zero
	vaddsd	%xmm1, %xmm0, %xmm1
.Ltmp47:
	vaddsd	%xmm3, %xmm2, %xmm3
.Ltmp48:
	.loc	1 0 1 is_stmt 0         # backgroundfield/integratefunction.cpp:0:1
	leaq	-40(%rbp), %rdi
.Ltmp49:
	.loc	1 104 1                 # backgroundfield/integratefunction.cpp:104:1
	callq	_Z7RombergRK11T2DFunctionddddd
.Ltmp50:
	#DEBUG_VALUE: value <- undef
	.loc	1 74 1 is_stmt 1        # backgroundfield/integratefunction.cpp:74:1
	movq	$_ZTV11T2DFunction+16, -40(%rbp)
	jmp	.LBB1_7
.Ltmp51:
.LBB1_6:                                # %L.B0017
	#DEBUG_VALUE: acc <- $xmm4
	#DEBUG_VALUE: surfaceAverage:L2 <- $xmm3
	#DEBUG_VALUE: face <- $esi
	#DEBUG_VALUE: L1 <- $xmm1
	#DEBUG_VALUE: L2 <- $xmm3
	#DEBUG_VALUE: surfaceAverage:L1 <- $xmm1
	#DEBUG_VALUE: surfaceAverage:r1 <- $rdx
	#DEBUG_VALUE: surfaceAverage:face <- $esi
	#DEBUG_VALUE: surfaceAverage:f1 <- $rdi
	.loc	1 97 1                  # backgroundfield/integratefunction.cpp:97:1
	movq	(%rdx), %rax
.Ltmp52:
	.loc	1 62 1                  # backgroundfield/integratefunction.cpp:62:1
	movq	$_ZTV8T3D_fix1+16, -64(%rbp)
	movq	%rdi, -56(%rbp)
	movq	%rax, -48(%rbp)
.Ltmp53:
	.loc	1 98 1                  # backgroundfield/integratefunction.cpp:98:1
	vmovsd	8(%rdx), %xmm0          # xmm0 = mem[0],zero
	vmovsd	16(%rdx), %xmm2         # xmm2 = mem[0],zero
	vaddsd	%xmm1, %xmm0, %xmm1
.Ltmp54:
	vaddsd	%xmm3, %xmm2, %xmm3
.Ltmp55:
	.loc	1 0 1 is_stmt 0         # backgroundfield/integratefunction.cpp:0:1
	leaq	-64(%rbp), %rdi
.Ltmp56:
	.loc	1 98 1                  # backgroundfield/integratefunction.cpp:98:1
	callq	_Z7RombergRK11T2DFunctionddddd
.Ltmp57:
	#DEBUG_VALUE: value <- undef
	.loc	1 64 1 is_stmt 1        # backgroundfield/integratefunction.cpp:64:1
	movq	$_ZTV11T2DFunction+16, -64(%rbp)
.Ltmp58:
.LBB1_7:                                # %L_T31592664_8085
	.loc	1 0 1 is_stmt 0         # backgroundfield/integratefunction.cpp:0:1
	vmovsd	-8(%rbp), %xmm1         # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vmulsd	-16(%rbp), %xmm1, %xmm1 # 8-byte Folded Reload
	vmovsd	.LCPI1_0(%rip), %xmm2   # xmm2 = mem[0],zero
.Ltmp59:
	vdivsd	%xmm1, %xmm2, %xmm1
.Ltmp60:
	#DEBUG_VALUE: norm <- $xmm1
	vmulsd	%xmm0, %xmm1, %xmm0
.Ltmp61:
	#DEBUG_VALUE: value <- $xmm0
	.loc	1 119 1 is_stmt 1       # backgroundfield/integratefunction.cpp:119:1
	addq	$96, %rsp
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp62:
.LBB1_5:                                # %L..inline.10054
	.cfi_def_cfa %rbp, 16
	#DEBUG_VALUE: acc <- $xmm4
	#DEBUG_VALUE: surfaceAverage:L2 <- $xmm3
	#DEBUG_VALUE: face <- $esi
	#DEBUG_VALUE: L1 <- $xmm1
	#DEBUG_VALUE: L2 <- $xmm3
	#DEBUG_VALUE: surfaceAverage:L1 <- $xmm1
	#DEBUG_VALUE: surfaceAverage:r1 <- $rdx
	#DEBUG_VALUE: surfaceAverage:face <- $esi
	#DEBUG_VALUE: surfaceAverage:f1 <- $rdi
	.loc	1 570 1                 # backgroundfield/integratefunction.cpp:570:1
	movl	$_ZSt4cerr, %edi
.Ltmp63:
	movl	$.S08496, %esi
.Ltmp64:
	movl	$27, %edx
.Ltmp65:
	callq	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l
.Ltmp66:
	.loc	1 115 1                 # backgroundfield/integratefunction.cpp:115:1
	movl	$1, %edi
	callq	exit
.Ltmp67:
.Lfunc_end1:
	.size	_Z14surfaceAverageRK11T3DFunction10coordinatedPKddd, .Lfunc_end1-_Z14surfaceAverageRK11T3DFunction10coordinatedPKddd
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _Z13volumeAverageRK11T3DFunctiondPKdS3_
.LCPI2_0:
	.quad	4607182418800017408     # double 1
	.text
	.globl	_Z13volumeAverageRK11T3DFunctiondPKdS3_
	.p2align	4, 0x90
	.type	_Z13volumeAverageRK11T3DFunctiondPKdS3_,@function
_Z13volumeAverageRK11T3DFunctiondPKdS3_: # @_Z13volumeAverageRK11T3DFunctiondPKdS3_
.Lfunc_begin2:
	.loc	1 128 0                 # backgroundfield/integratefunction.cpp:128:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: volumeAverage:f1 <- $rdi
	#DEBUG_VALUE: volumeAverage:f1 <- $rdi
	#DEBUG_VALUE: volumeAverage:accuracy <- $xmm0
	#DEBUG_VALUE: volumeAverage:accuracy <- $xmm0
	#DEBUG_VALUE: volumeAverage:r1 <- $rsi
	#DEBUG_VALUE: volumeAverage:r1 <- $rsi
	#DEBUG_VALUE: volumeAverage:r2 <- $rdx
	#DEBUG_VALUE: volumeAverage:r2 <- $rdx
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$16, %rsp
.Ltmp68:
	#DEBUG_VALUE: f1 <- $rdi
	#DEBUG_VALUE: f1 <- $rdi
	#DEBUG_VALUE: accuracy <- $xmm0
	#DEBUG_VALUE: accuracy <- $xmm0
	#DEBUG_VALUE: r1 <- $rsi
	#DEBUG_VALUE: r1 <- $rsi
	#DEBUG_VALUE: r2 <- $rdx
	#DEBUG_VALUE: r2 <- $rdx
	.loc	1 132 1 prologue_end    # backgroundfield/integratefunction.cpp:132:1
	vmovsd	16(%rdx), %xmm5         # xmm5 = mem[0],zero
	vmovsd	16(%rsi), %xmm4         # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm5, %xmm9
	vmovsd	(%rdx), %xmm1           # xmm1 = mem[0],zero
	vmovsd	8(%rdx), %xmm3          # xmm3 = mem[0],zero
	vmovsd	(%rsi), %xmm8           # xmm8 = mem[0],zero
	vmovsd	8(%rsi), %xmm2          # xmm2 = mem[0],zero
	vsubsd	%xmm2, %xmm3, %xmm10
	vsubsd	%xmm8, %xmm1, %xmm7
	vmulsd	%xmm0, %xmm7, %xmm0
.Ltmp69:
	vmulsd	%xmm0, %xmm10, %xmm0
	vmulsd	%xmm0, %xmm9, %xmm6
.Ltmp70:
	#DEBUG_VALUE: acc <- $xmm6
	.loc	1 133 1                 # backgroundfield/integratefunction.cpp:133:1
	vmulsd	%xmm7, %xmm10, %xmm0
	vmulsd	%xmm0, %xmm9, %xmm0
	vmovsd	.LCPI2_0(%rip), %xmm7   # xmm7 = mem[0],zero
	vdivsd	%xmm0, %xmm7, %xmm0
.Ltmp71:
	#DEBUG_VALUE: norm <- $xmm0
	.loc	1 0 1 is_stmt 0         # backgroundfield/integratefunction.cpp:0:1
	vmovsd	%xmm0, -8(%rbp)         # 8-byte Spill
.Ltmp72:
	#DEBUG_VALUE: norm <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: norm <- [DW_OP_constu 8, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f1 <- $rdi
	#DEBUG_VALUE: volumeAverage:f1 <- $rdi
	.loc	1 134 1 is_stmt 1       # backgroundfield/integratefunction.cpp:134:1
	vmovapd	%xmm8, %xmm0
	callq	_Z7RombergRK11T3DFunctionddddddd
.Ltmp73:
	vmulsd	-8(%rbp), %xmm0, %xmm0  # 8-byte Folded Reload
.Ltmp74:
	.loc	1 136 1                 # backgroundfield/integratefunction.cpp:136:1
	addq	$16, %rsp
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp75:
.Lfunc_end2:
	.size	_Z13volumeAverageRK11T3DFunctiondPKdS3_, .Lfunc_end2-_Z13volumeAverageRK11T3DFunctiondPKdS3_
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T1DFunctionD1Ev    # -- Begin function _ZN11T1DFunctionD1Ev
	.p2align	4, 0x90
	.type	_ZN11T1DFunctionD1Ev,@function
_ZN11T1DFunctionD1Ev:                   # @_ZN11T1DFunctionD1Ev
.Lfunc_begin3:
	.file	2 "/home/talgat/vlasiator/backgroundfield/functions.hpp"
	.loc	2 29 0                  # backgroundfield/functions.hpp:29:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T1DFunction: <- $rdi
	#DEBUG_VALUE: ~T1DFunction: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp76:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 29 1 prologue_end     # backgroundfield/functions.hpp:29:1
	movq	$_ZTV11T1DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp77:
.Lfunc_end3:
	.size	_ZN11T1DFunctionD1Ev, .Lfunc_end3-_ZN11T1DFunctionD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T1DFunctionD0Ev    # -- Begin function _ZN11T1DFunctionD0Ev
	.p2align	4, 0x90
	.type	_ZN11T1DFunctionD0Ev,@function
_ZN11T1DFunctionD0Ev:                   # @_ZN11T1DFunctionD0Ev
.Lfunc_begin4:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp78:
	.loc	1 29 1 prologue_end     # backgroundfield/integratefunction.cpp:29:1
	movl	$8, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp79:
.Lfunc_end4:
	.size	_ZN11T1DFunctionD0Ev, .Lfunc_end4-_ZN11T1DFunctionD0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T2DFunctionD1Ev    # -- Begin function _ZN11T2DFunctionD1Ev
	.p2align	4, 0x90
	.type	_ZN11T2DFunctionD1Ev,@function
_ZN11T2DFunctionD1Ev:                   # @_ZN11T2DFunctionD1Ev
.Lfunc_begin5:
	.loc	2 30 0                  # backgroundfield/functions.hpp:30:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T2DFunction: <- $rdi
	#DEBUG_VALUE: ~T2DFunction: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp80:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 30 1 prologue_end     # backgroundfield/functions.hpp:30:1
	movq	$_ZTV11T2DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp81:
.Lfunc_end5:
	.size	_ZN11T2DFunctionD1Ev, .Lfunc_end5-_ZN11T2DFunctionD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T2DFunctionD0Ev    # -- Begin function _ZN11T2DFunctionD0Ev
	.p2align	4, 0x90
	.type	_ZN11T2DFunctionD0Ev,@function
_ZN11T2DFunctionD0Ev:                   # @_ZN11T2DFunctionD0Ev
.Lfunc_begin6:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp82:
	.loc	1 30 1 prologue_end     # backgroundfield/integratefunction.cpp:30:1
	movl	$8, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp83:
.Lfunc_end6:
	.size	_ZN11T2DFunctionD0Ev, .Lfunc_end6-_ZN11T2DFunctionD0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK8T3D_fix14callEdd   # -- Begin function _ZNK8T3D_fix14callEdd
	.p2align	4, 0x90
	.type	_ZNK8T3D_fix14callEdd,@function
_ZNK8T3D_fix14callEdd:                  # @_ZNK8T3D_fix14callEdd
.Lfunc_begin7:
	.loc	2 63 0                  # backgroundfield/functions.hpp:63:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:y <- $xmm0
	#DEBUG_VALUE: call:z <- $xmm1
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp84:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: y <- $xmm0
	#DEBUG_VALUE: z <- $xmm1
	vmovaps	%xmm1, %xmm2
	vmovaps	%xmm0, %xmm1
.Ltmp85:
	#DEBUG_VALUE: z <- $xmm2
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	2 63 1 prologue_end     # backgroundfield/functions.hpp:63:1
	movq	8(%rdi), %rax
	vmovsd	16(%rdi), %xmm0         # xmm0 = mem[0],zero
	movq	(%rax), %rcx
	movq	(%rcx), %rcx
	movq	%rax, %rdi
.Ltmp86:
	#DEBUG_VALUE: z <- $xmm2
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmpq	*%rcx                   # TAILCALL
.Ltmp87:
.Lfunc_end7:
	.size	_ZNK8T3D_fix14callEdd, .Lfunc_end7-_ZNK8T3D_fix14callEdd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T3D_fix1D1Ev        # -- Begin function _ZN8T3D_fix1D1Ev
	.p2align	4, 0x90
	.type	_ZN8T3D_fix1D1Ev,@function
_ZN8T3D_fix1D1Ev:                       # @_ZN8T3D_fix1D1Ev
.Lfunc_begin8:
	.loc	2 64 0                  # backgroundfield/functions.hpp:64:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T3D_fix1: <- $rdi
	#DEBUG_VALUE: ~T3D_fix1: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp88:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 64 1 prologue_end     # backgroundfield/functions.hpp:64:1
	movq	$_ZTV11T2DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp89:
.Lfunc_end8:
	.size	_ZN8T3D_fix1D1Ev, .Lfunc_end8-_ZN8T3D_fix1D1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T3D_fix1D0Ev        # -- Begin function _ZN8T3D_fix1D0Ev
	.p2align	4, 0x90
	.type	_ZN8T3D_fix1D0Ev,@function
_ZN8T3D_fix1D0Ev:                       # @_ZN8T3D_fix1D0Ev
.Lfunc_begin9:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp90:
	.loc	1 64 1 prologue_end     # backgroundfield/integratefunction.cpp:64:1
	movl	$24, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp91:
.Lfunc_end9:
	.size	_ZN8T3D_fix1D0Ev, .Lfunc_end9-_ZN8T3D_fix1D0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK8T3D_fix24callEdd   # -- Begin function _ZNK8T3D_fix24callEdd
	.p2align	4, 0x90
	.type	_ZNK8T3D_fix24callEdd,@function
_ZNK8T3D_fix24callEdd:                  # @_ZNK8T3D_fix24callEdd
.Lfunc_begin10:
	.loc	2 73 0                  # backgroundfield/functions.hpp:73:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:x <- $xmm0
	#DEBUG_VALUE: call:z <- $xmm1
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp92:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: z <- $xmm1
	vmovaps	%xmm1, %xmm2
.Ltmp93:
	#DEBUG_VALUE: z <- $xmm2
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: call:x <- $xmm0
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	2 73 1 prologue_end     # backgroundfield/functions.hpp:73:1
	movq	8(%rdi), %rax
	vmovsd	16(%rdi), %xmm1         # xmm1 = mem[0],zero
	movq	(%rax), %rcx
	movq	(%rcx), %rcx
	movq	%rax, %rdi
.Ltmp94:
	#DEBUG_VALUE: z <- $xmm2
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: call:x <- $xmm0
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmpq	*%rcx                   # TAILCALL
.Ltmp95:
.Lfunc_end10:
	.size	_ZNK8T3D_fix24callEdd, .Lfunc_end10-_ZNK8T3D_fix24callEdd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T3D_fix2D1Ev        # -- Begin function _ZN8T3D_fix2D1Ev
	.p2align	4, 0x90
	.type	_ZN8T3D_fix2D1Ev,@function
_ZN8T3D_fix2D1Ev:                       # @_ZN8T3D_fix2D1Ev
.Lfunc_begin11:
	.loc	2 74 0                  # backgroundfield/functions.hpp:74:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T3D_fix2: <- $rdi
	#DEBUG_VALUE: ~T3D_fix2: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp96:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 74 1 prologue_end     # backgroundfield/functions.hpp:74:1
	movq	$_ZTV11T2DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp97:
.Lfunc_end11:
	.size	_ZN8T3D_fix2D1Ev, .Lfunc_end11-_ZN8T3D_fix2D1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T3D_fix2D0Ev        # -- Begin function _ZN8T3D_fix2D0Ev
	.p2align	4, 0x90
	.type	_ZN8T3D_fix2D0Ev,@function
_ZN8T3D_fix2D0Ev:                       # @_ZN8T3D_fix2D0Ev
.Lfunc_begin12:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp98:
	.loc	1 74 1 prologue_end     # backgroundfield/integratefunction.cpp:74:1
	movl	$24, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp99:
.Lfunc_end12:
	.size	_ZN8T3D_fix2D0Ev, .Lfunc_end12-_ZN8T3D_fix2D0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK8T3D_fix34callEdd   # -- Begin function _ZNK8T3D_fix34callEdd
	.p2align	4, 0x90
	.type	_ZNK8T3D_fix34callEdd,@function
_ZNK8T3D_fix34callEdd:                  # @_ZNK8T3D_fix34callEdd
.Lfunc_begin13:
	.loc	2 83 0                  # backgroundfield/functions.hpp:83:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:x <- $xmm0
	#DEBUG_VALUE: call:x <- $xmm0
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp100:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: y <- $xmm1
	.loc	2 83 1 prologue_end     # backgroundfield/functions.hpp:83:1
	movq	8(%rdi), %rax
	vmovsd	16(%rdi), %xmm2         # xmm2 = mem[0],zero
	movq	(%rax), %rcx
	movq	(%rcx), %rcx
	movq	%rax, %rdi
.Ltmp101:
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: call:x <- $xmm0
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmpq	*%rcx                   # TAILCALL
.Ltmp102:
.Lfunc_end13:
	.size	_ZNK8T3D_fix34callEdd, .Lfunc_end13-_ZNK8T3D_fix34callEdd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T3D_fix3D1Ev        # -- Begin function _ZN8T3D_fix3D1Ev
	.p2align	4, 0x90
	.type	_ZN8T3D_fix3D1Ev,@function
_ZN8T3D_fix3D1Ev:                       # @_ZN8T3D_fix3D1Ev
.Lfunc_begin14:
	.loc	2 84 0                  # backgroundfield/functions.hpp:84:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T3D_fix3: <- $rdi
	#DEBUG_VALUE: ~T3D_fix3: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp103:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 84 1 prologue_end     # backgroundfield/functions.hpp:84:1
	movq	$_ZTV11T2DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp104:
.Lfunc_end14:
	.size	_ZN8T3D_fix3D1Ev, .Lfunc_end14-_ZN8T3D_fix3D1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T3D_fix3D0Ev        # -- Begin function _ZN8T3D_fix3D0Ev
	.p2align	4, 0x90
	.type	_ZN8T3D_fix3D0Ev,@function
_ZN8T3D_fix3D0Ev:                       # @_ZN8T3D_fix3D0Ev
.Lfunc_begin15:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp105:
	.loc	1 84 1 prologue_end     # backgroundfield/integratefunction.cpp:84:1
	movl	$24, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp106:
.Lfunc_end15:
	.size	_ZN8T3D_fix3D0Ev, .Lfunc_end15-_ZN8T3D_fix3D0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK9T3D_fix124callEd   # -- Begin function _ZNK9T3D_fix124callEd
	.p2align	4, 0x90
	.type	_ZNK9T3D_fix124callEd,@function
_ZNK9T3D_fix124callEd:                  # @_ZNK9T3D_fix124callEd
.Lfunc_begin16:
	.loc	2 95 0                  # backgroundfield/functions.hpp:95:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:z <- $xmm0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp107:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: z <- $xmm0
	vmovaps	%xmm0, %xmm2
.Ltmp108:
	#DEBUG_VALUE: z <- $xmm2
	#DEBUG_VALUE: call:z <- $xmm2
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	2 95 1 prologue_end     # backgroundfield/functions.hpp:95:1
	movq	8(%rdi), %rax
	vmovsd	16(%rdi), %xmm0         # xmm0 = mem[0],zero
	vmovsd	24(%rdi), %xmm1         # xmm1 = mem[0],zero
	movq	(%rax), %rcx
	movq	(%rcx), %rcx
	movq	%rax, %rdi
.Ltmp109:
	#DEBUG_VALUE: z <- $xmm2
	#DEBUG_VALUE: call:z <- $xmm2
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmpq	*%rcx                   # TAILCALL
.Ltmp110:
.Lfunc_end16:
	.size	_ZNK9T3D_fix124callEd, .Lfunc_end16-_ZNK9T3D_fix124callEd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN9T3D_fix12D1Ev       # -- Begin function _ZN9T3D_fix12D1Ev
	.p2align	4, 0x90
	.type	_ZN9T3D_fix12D1Ev,@function
_ZN9T3D_fix12D1Ev:                      # @_ZN9T3D_fix12D1Ev
.Lfunc_begin17:
	.loc	2 96 0                  # backgroundfield/functions.hpp:96:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T3D_fix12: <- $rdi
	#DEBUG_VALUE: ~T3D_fix12: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp111:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 96 1 prologue_end     # backgroundfield/functions.hpp:96:1
	movq	$_ZTV11T1DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp112:
.Lfunc_end17:
	.size	_ZN9T3D_fix12D1Ev, .Lfunc_end17-_ZN9T3D_fix12D1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN9T3D_fix12D0Ev       # -- Begin function _ZN9T3D_fix12D0Ev
	.p2align	4, 0x90
	.type	_ZN9T3D_fix12D0Ev,@function
_ZN9T3D_fix12D0Ev:                      # @_ZN9T3D_fix12D0Ev
.Lfunc_begin18:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp113:
	.loc	1 96 1 prologue_end     # backgroundfield/integratefunction.cpp:96:1
	movl	$32, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp114:
.Lfunc_end18:
	.size	_ZN9T3D_fix12D0Ev, .Lfunc_end18-_ZN9T3D_fix12D0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK9T3D_fix134callEd   # -- Begin function _ZNK9T3D_fix134callEd
	.p2align	4, 0x90
	.type	_ZNK9T3D_fix134callEd,@function
_ZNK9T3D_fix134callEd:                  # @_ZNK9T3D_fix134callEd
.Lfunc_begin19:
	.loc	2 105 0                 # backgroundfield/functions.hpp:105:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:y <- $xmm0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp115:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: y <- $xmm0
	vmovaps	%xmm0, %xmm1
.Ltmp116:
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	2 105 1 prologue_end    # backgroundfield/functions.hpp:105:1
	movq	8(%rdi), %rax
	vmovsd	16(%rdi), %xmm0         # xmm0 = mem[0],zero
	vmovsd	24(%rdi), %xmm2         # xmm2 = mem[0],zero
	movq	(%rax), %rcx
	movq	(%rcx), %rcx
	movq	%rax, %rdi
.Ltmp117:
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmpq	*%rcx                   # TAILCALL
.Ltmp118:
.Lfunc_end19:
	.size	_ZNK9T3D_fix134callEd, .Lfunc_end19-_ZNK9T3D_fix134callEd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN9T3D_fix13D1Ev       # -- Begin function _ZN9T3D_fix13D1Ev
	.p2align	4, 0x90
	.type	_ZN9T3D_fix13D1Ev,@function
_ZN9T3D_fix13D1Ev:                      # @_ZN9T3D_fix13D1Ev
.Lfunc_begin20:
	.loc	2 106 0                 # backgroundfield/functions.hpp:106:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T3D_fix13: <- $rdi
	#DEBUG_VALUE: ~T3D_fix13: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp119:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 106 1 prologue_end    # backgroundfield/functions.hpp:106:1
	movq	$_ZTV11T1DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp120:
.Lfunc_end20:
	.size	_ZN9T3D_fix13D1Ev, .Lfunc_end20-_ZN9T3D_fix13D1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN9T3D_fix13D0Ev       # -- Begin function _ZN9T3D_fix13D0Ev
	.p2align	4, 0x90
	.type	_ZN9T3D_fix13D0Ev,@function
_ZN9T3D_fix13D0Ev:                      # @_ZN9T3D_fix13D0Ev
.Lfunc_begin21:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp121:
	.loc	1 106 1 prologue_end    # backgroundfield/integratefunction.cpp:106:1
	movl	$32, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp122:
.Lfunc_end21:
	.size	_ZN9T3D_fix13D0Ev, .Lfunc_end21-_ZN9T3D_fix13D0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK9T3D_fix234callEd   # -- Begin function _ZNK9T3D_fix234callEd
	.p2align	4, 0x90
	.type	_ZNK9T3D_fix234callEd,@function
_ZNK9T3D_fix234callEd:                  # @_ZNK9T3D_fix234callEd
.Lfunc_begin22:
	.loc	2 115 0                 # backgroundfield/functions.hpp:115:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:x <- $xmm0
	#DEBUG_VALUE: call:x <- $xmm0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp123:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: x <- $xmm0
	.loc	2 115 1 prologue_end    # backgroundfield/functions.hpp:115:1
	movq	8(%rdi), %rax
	vmovsd	16(%rdi), %xmm1         # xmm1 = mem[0],zero
	vmovsd	24(%rdi), %xmm2         # xmm2 = mem[0],zero
	movq	(%rax), %rcx
	movq	(%rcx), %rcx
	movq	%rax, %rdi
.Ltmp124:
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: call:x <- $xmm0
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmpq	*%rcx                   # TAILCALL
.Ltmp125:
.Lfunc_end22:
	.size	_ZNK9T3D_fix234callEd, .Lfunc_end22-_ZNK9T3D_fix234callEd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN9T3D_fix23D1Ev       # -- Begin function _ZN9T3D_fix23D1Ev
	.p2align	4, 0x90
	.type	_ZN9T3D_fix23D1Ev,@function
_ZN9T3D_fix23D1Ev:                      # @_ZN9T3D_fix23D1Ev
.Lfunc_begin23:
	.loc	2 116 0                 # backgroundfield/functions.hpp:116:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T3D_fix23: <- $rdi
	#DEBUG_VALUE: ~T3D_fix23: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp126:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 116 1 prologue_end    # backgroundfield/functions.hpp:116:1
	movq	$_ZTV11T1DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp127:
.Lfunc_end23:
	.size	_ZN9T3D_fix23D1Ev, .Lfunc_end23-_ZN9T3D_fix23D1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN9T3D_fix23D0Ev       # -- Begin function _ZN9T3D_fix23D0Ev
	.p2align	4, 0x90
	.type	_ZN9T3D_fix23D0Ev,@function
_ZN9T3D_fix23D0Ev:                      # @_ZN9T3D_fix23D0Ev
.Lfunc_begin24:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp128:
	.loc	1 116 1 prologue_end    # backgroundfield/integratefunction.cpp:116:1
	movl	$32, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp129:
.Lfunc_end24:
	.size	_ZN9T3D_fix23D0Ev, .Lfunc_end24-_ZN9T3D_fix23D0Ev
	.cfi_endproc
                                        # -- End function
	.globl	__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee # -- Begin function __sti___37_backgroundfield_integratefunction_cpp_f19cb1ee
	.p2align	4, 0x90
	.type	__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee,@function
__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee: # @__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee
.Lfunc_begin25:
	.loc	1 0 0                   # backgroundfield/integratefunction.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp130:
	.loc	1 74 1 prologue_end     # backgroundfield/integratefunction.cpp:74:1
	cmpl	$1, __I___37_backgroundfield_integratefunction_cpp_f19cb1ee(%rip)
	jne	.LBB25_2
# %bb.1:                                # %L.B0034
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB25_2:                               # %L.B0290
	.cfi_def_cfa %rbp, 16
	movl	$1, __I___37_backgroundfield_integratefunction_cpp_f19cb1ee(%rip)
	movl	$_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE, %edi
	callq	_ZNSt8ios_base4InitC1Ev
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	movl	$_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE, %esi
	movl	$__dso_handle, %edx
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	__cxa_atexit            # TAILCALL
.Ltmp131:
.Lfunc_end25:
	.size	__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee, .Lfunc_end25-__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee
	.cfi_endproc
                                        # -- End function
	.type	_ZTV11T1DFunction,@object # @_ZTV11T1DFunction
	.data
	.weak	_ZTV11T1DFunction
	.p2align	4
_ZTV11T1DFunction:
	.quad	0
	.quad	_ZTI11T1DFunction
	.quad	__cxa_pure_virtual
	.quad	_ZN11T1DFunctionD1Ev
	.quad	_ZN11T1DFunctionD0Ev
	.size	_ZTV11T1DFunction, 40

	.type	_ZTV11T2DFunction,@object # @_ZTV11T2DFunction
	.weak	_ZTV11T2DFunction
	.p2align	4
_ZTV11T2DFunction:
	.quad	0
	.quad	_ZTI11T2DFunction
	.quad	__cxa_pure_virtual
	.quad	_ZN11T2DFunctionD1Ev
	.quad	_ZN11T2DFunctionD0Ev
	.size	_ZTV11T2DFunction, 40

	.type	_ZTV8T3D_fix1,@object   # @_ZTV8T3D_fix1
	.weak	_ZTV8T3D_fix1
	.p2align	4
_ZTV8T3D_fix1:
	.quad	0
	.quad	_ZTI8T3D_fix1
	.quad	_ZNK8T3D_fix14callEdd
	.quad	_ZN8T3D_fix1D1Ev
	.quad	_ZN8T3D_fix1D0Ev
	.size	_ZTV8T3D_fix1, 40

	.type	_ZTV8T3D_fix2,@object   # @_ZTV8T3D_fix2
	.weak	_ZTV8T3D_fix2
	.p2align	4
_ZTV8T3D_fix2:
	.quad	0
	.quad	_ZTI8T3D_fix2
	.quad	_ZNK8T3D_fix24callEdd
	.quad	_ZN8T3D_fix2D1Ev
	.quad	_ZN8T3D_fix2D0Ev
	.size	_ZTV8T3D_fix2, 40

	.type	_ZTV8T3D_fix3,@object   # @_ZTV8T3D_fix3
	.weak	_ZTV8T3D_fix3
	.p2align	4
_ZTV8T3D_fix3:
	.quad	0
	.quad	_ZTI8T3D_fix3
	.quad	_ZNK8T3D_fix34callEdd
	.quad	_ZN8T3D_fix3D1Ev
	.quad	_ZN8T3D_fix3D0Ev
	.size	_ZTV8T3D_fix3, 40

	.type	_ZTV9T3D_fix12,@object  # @_ZTV9T3D_fix12
	.weak	_ZTV9T3D_fix12
	.p2align	4
_ZTV9T3D_fix12:
	.quad	0
	.quad	_ZTI9T3D_fix12
	.quad	_ZNK9T3D_fix124callEd
	.quad	_ZN9T3D_fix12D1Ev
	.quad	_ZN9T3D_fix12D0Ev
	.size	_ZTV9T3D_fix12, 40

	.type	_ZTV9T3D_fix13,@object  # @_ZTV9T3D_fix13
	.weak	_ZTV9T3D_fix13
	.p2align	4
_ZTV9T3D_fix13:
	.quad	0
	.quad	_ZTI9T3D_fix13
	.quad	_ZNK9T3D_fix134callEd
	.quad	_ZN9T3D_fix13D1Ev
	.quad	_ZN9T3D_fix13D0Ev
	.size	_ZTV9T3D_fix13, 40

	.type	_ZTV9T3D_fix23,@object  # @_ZTV9T3D_fix23
	.weak	_ZTV9T3D_fix23
	.p2align	4
_ZTV9T3D_fix23:
	.quad	0
	.quad	_ZTI9T3D_fix23
	.quad	_ZNK9T3D_fix234callEd
	.quad	_ZN9T3D_fix23D1Ev
	.quad	_ZN9T3D_fix23D0Ev
	.size	_ZTV9T3D_fix23, 40

	.type	_ZTI11T1DFunction,@object # @_ZTI11T1DFunction
	.weak	_ZTI11T1DFunction
	.p2align	4
_ZTI11T1DFunction:
	.quad	_ZTVN10__cxxabiv117__class_type_infoE+16
	.quad	_ZTS11T1DFunction
	.size	_ZTI11T1DFunction, 16

	.type	_ZTI11T2DFunction,@object # @_ZTI11T2DFunction
	.weak	_ZTI11T2DFunction
	.p2align	4
_ZTI11T2DFunction:
	.quad	_ZTVN10__cxxabiv117__class_type_infoE+16
	.quad	_ZTS11T2DFunction
	.size	_ZTI11T2DFunction, 16

	.type	_ZTI8T3D_fix1,@object   # @_ZTI8T3D_fix1
	.weak	_ZTI8T3D_fix1
	.p2align	4
_ZTI8T3D_fix1:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS8T3D_fix1
	.quad	_ZTI11T2DFunction
	.size	_ZTI8T3D_fix1, 24

	.type	_ZTI8T3D_fix2,@object   # @_ZTI8T3D_fix2
	.weak	_ZTI8T3D_fix2
	.p2align	4
_ZTI8T3D_fix2:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS8T3D_fix2
	.quad	_ZTI11T2DFunction
	.size	_ZTI8T3D_fix2, 24

	.type	_ZTI8T3D_fix3,@object   # @_ZTI8T3D_fix3
	.weak	_ZTI8T3D_fix3
	.p2align	4
_ZTI8T3D_fix3:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS8T3D_fix3
	.quad	_ZTI11T2DFunction
	.size	_ZTI8T3D_fix3, 24

	.type	_ZTI9T3D_fix12,@object  # @_ZTI9T3D_fix12
	.weak	_ZTI9T3D_fix12
	.p2align	4
_ZTI9T3D_fix12:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS9T3D_fix12
	.quad	_ZTI11T1DFunction
	.size	_ZTI9T3D_fix12, 24

	.type	_ZTI9T3D_fix13,@object  # @_ZTI9T3D_fix13
	.weak	_ZTI9T3D_fix13
	.p2align	4
_ZTI9T3D_fix13:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS9T3D_fix13
	.quad	_ZTI11T1DFunction
	.size	_ZTI9T3D_fix13, 24

	.type	_ZTI9T3D_fix23,@object  # @_ZTI9T3D_fix23
	.weak	_ZTI9T3D_fix23
	.p2align	4
_ZTI9T3D_fix23:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS9T3D_fix23
	.quad	_ZTI11T1DFunction
	.size	_ZTI9T3D_fix23, 24

	.type	_ZTS11T1DFunction,@object # @_ZTS11T1DFunction
	.weak	_ZTS11T1DFunction
	.p2align	3
_ZTS11T1DFunction:
	.asciz	"11T1DFunction"
	.size	_ZTS11T1DFunction, 14

	.type	_ZTS11T2DFunction,@object # @_ZTS11T2DFunction
	.weak	_ZTS11T2DFunction
	.p2align	3
_ZTS11T2DFunction:
	.asciz	"11T2DFunction"
	.size	_ZTS11T2DFunction, 14

	.type	_ZTS8T3D_fix1,@object   # @_ZTS8T3D_fix1
	.weak	_ZTS8T3D_fix1
	.p2align	3
_ZTS8T3D_fix1:
	.asciz	"8T3D_fix1"
	.size	_ZTS8T3D_fix1, 10

	.type	_ZTS8T3D_fix2,@object   # @_ZTS8T3D_fix2
	.weak	_ZTS8T3D_fix2
	.p2align	3
_ZTS8T3D_fix2:
	.asciz	"8T3D_fix2"
	.size	_ZTS8T3D_fix2, 10

	.type	_ZTS8T3D_fix3,@object   # @_ZTS8T3D_fix3
	.weak	_ZTS8T3D_fix3
	.p2align	3
_ZTS8T3D_fix3:
	.asciz	"8T3D_fix3"
	.size	_ZTS8T3D_fix3, 10

	.type	_ZTS9T3D_fix12,@object  # @_ZTS9T3D_fix12
	.weak	_ZTS9T3D_fix12
	.p2align	3
_ZTS9T3D_fix12:
	.asciz	"9T3D_fix12"
	.size	_ZTS9T3D_fix12, 11

	.type	_ZTS9T3D_fix13,@object  # @_ZTS9T3D_fix13
	.weak	_ZTS9T3D_fix13
	.p2align	3
_ZTS9T3D_fix13:
	.asciz	"9T3D_fix13"
	.size	_ZTS9T3D_fix13, 11

	.type	_ZTS9T3D_fix23,@object  # @_ZTS9T3D_fix23
	.weak	_ZTS9T3D_fix23
	.p2align	3
_ZTS9T3D_fix23:
	.asciz	"9T3D_fix23"
	.size	_ZTS9T3D_fix23, 11

	.type	__I___37_backgroundfield_integratefunction_cpp_f19cb1ee,@object # @__I___37_backgroundfield_integratefunction_cpp_f19cb1ee
	.bss
	.globl	__I___37_backgroundfield_integratefunction_cpp_f19cb1ee
	.p2align	2
__I___37_backgroundfield_integratefunction_cpp_f19cb1ee:
	.long	0                       # 0x0
	.size	__I___37_backgroundfield_integratefunction_cpp_f19cb1ee, 4

	.type	_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE,@object # @_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE
	.local	_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE
	.comm	_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE,1,1
	.type	.S08496,@object         # @.S08496
	.section	.rodata,"a",@progbits
.S08496:
	.asciz	"*** SurfaceAverage  is bad\n"
	.size	.S08496, 28

	.type	.S08448,@object         # @.S08448
.S08448:
	.asciz	"*** lineAverage  is bad\n"
	.size	.S08448, 25

	.section	.init_array,"aw",@init_array
	.p2align	3
	.quad	__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee
	.section	.debug_str,"MS",@progbits,1
.Linfo_string0:
	.asciz	" NVC++ 21.2-0"         # string offset=0
.Linfo_string1:
	.asciz	"backgroundfield/integratefunction.cpp" # string offset=14
.Linfo_string2:
	.asciz	"/home/talgat/vlasiator" # string offset=52
.Linfo_string3:
	.asciz	"_ZTV11T1DFunction"     # string offset=75
.Linfo_string4:
	.asciz	"int"                   # string offset=93
.Linfo_string5:
	.asciz	"__ARRAY_SIZE_TYPE__"   # string offset=97
.Linfo_string6:
	.asciz	"_ZTV9T3D_fix12"        # string offset=117
.Linfo_string7:
	.asciz	"std"                   # string offset=132
.Linfo_string8:
	.asciz	"_ZSt4cerr"             # string offset=136
.Linfo_string9:
	.asciz	"__vptr"                # string offset=146
.Linfo_string10:
	.asciz	"__v_St9basic_iosIcSt11char_traitsIcEE" # string offset=153
.Linfo_string11:
	.asciz	"__b_St8ios_base"       # string offset=191
.Linfo_string12:
	.asciz	"_M_precision"          # string offset=207
.Linfo_string13:
	.asciz	"long"                  # string offset=220
.Linfo_string14:
	.asciz	"_M_width"              # string offset=225
.Linfo_string15:
	.asciz	"_M_flags"              # string offset=234
.Linfo_string16:
	.asciz	"_ZSt12_S_boolalpha"    # string offset=243
.Linfo_string17:
	.asciz	"_ZSt6_S_dec"           # string offset=262
.Linfo_string18:
	.asciz	"_ZSt8_S_fixed"         # string offset=274
.Linfo_string19:
	.asciz	"_ZSt6_S_hex"           # string offset=288
.Linfo_string20:
	.asciz	"_ZSt11_S_internal"     # string offset=300
.Linfo_string21:
	.asciz	"_ZSt7_S_left"          # string offset=318
.Linfo_string22:
	.asciz	"_ZSt6_S_oct"           # string offset=331
.Linfo_string23:
	.asciz	"_ZSt8_S_right"         # string offset=343
.Linfo_string24:
	.asciz	"_ZSt13_S_scientific"   # string offset=357
.Linfo_string25:
	.asciz	"_ZSt11_S_showbase"     # string offset=377
.Linfo_string26:
	.asciz	"_ZSt12_S_showpoint"    # string offset=395
.Linfo_string27:
	.asciz	"_ZSt10_S_showpos"      # string offset=414
.Linfo_string28:
	.asciz	"_ZSt9_S_skipws"        # string offset=431
.Linfo_string29:
	.asciz	"_ZSt10_S_unitbuf"      # string offset=446
.Linfo_string30:
	.asciz	"_ZSt12_S_uppercase"    # string offset=463
.Linfo_string31:
	.asciz	"_ZSt14_S_adjustfield"  # string offset=482
.Linfo_string32:
	.asciz	"_ZSt12_S_basefield"    # string offset=503
.Linfo_string33:
	.asciz	"_ZSt13_S_floatfield"   # string offset=522
.Linfo_string34:
	.asciz	"_ZSt19_S_ios_fmtflags_end" # string offset=542
.Linfo_string35:
	.asciz	"_ZSt19_S_ios_fmtflags_max" # string offset=568
.Linfo_string36:
	.asciz	"_ZSt19_S_ios_fmtflags_min" # string offset=594
.Linfo_string37:
	.asciz	"_ZSt13_Ios_Fmtflags"   # string offset=620
.Linfo_string38:
	.asciz	"_M_exception"          # string offset=640
.Linfo_string39:
	.asciz	"_ZSt10_S_goodbit"      # string offset=653
.Linfo_string40:
	.asciz	"_ZSt9_S_badbit"        # string offset=670
.Linfo_string41:
	.asciz	"_ZSt9_S_eofbit"        # string offset=685
.Linfo_string42:
	.asciz	"_ZSt10_S_failbit"      # string offset=700
.Linfo_string43:
	.asciz	"_ZSt18_S_ios_iostate_end" # string offset=717
.Linfo_string44:
	.asciz	"_ZSt18_S_ios_iostate_max" # string offset=742
.Linfo_string45:
	.asciz	"_ZSt18_S_ios_iostate_min" # string offset=767
.Linfo_string46:
	.asciz	"_ZSt12_Ios_Iostate"    # string offset=792
.Linfo_string47:
	.asciz	"_M_streambuf_state"    # string offset=811
.Linfo_string48:
	.asciz	"_M_callbacks"          # string offset=830
.Linfo_string49:
	.asciz	"_M_next"               # string offset=843
.Linfo_string50:
	.asciz	"_M_fn"                 # string offset=851
.Linfo_string51:
	.asciz	"_ZNSt8ios_base11erase_eventE" # string offset=857
.Linfo_string52:
	.asciz	"_ZNSt8ios_base11imbue_eventE" # string offset=886
.Linfo_string53:
	.asciz	"_ZNSt8ios_base13copyfmt_eventE" # string offset=915
.Linfo_string54:
	.asciz	"_ZNSt8ios_base5eventE" # string offset=946
.Linfo_string55:
	.asciz	"_M_index"              # string offset=968
.Linfo_string56:
	.asciz	"_M_refcount"           # string offset=977
.Linfo_string57:
	.asciz	"_ZNSt8ios_base14_Callback_listE" # string offset=989
.Linfo_string58:
	.asciz	"_M_word_zero"          # string offset=1021
.Linfo_string59:
	.asciz	"_M_pword"              # string offset=1034
.Linfo_string60:
	.asciz	"void"                  # string offset=1043
.Linfo_string61:
	.asciz	"_M_iword"              # string offset=1048
.Linfo_string62:
	.asciz	"_ZNSt8ios_base6_WordsE" # string offset=1057
.Linfo_string63:
	.asciz	"_M_local_word"         # string offset=1080
.Linfo_string64:
	.asciz	"_M_word_size"          # string offset=1094
.Linfo_string65:
	.asciz	"_M_word"               # string offset=1107
.Linfo_string66:
	.asciz	"_M_ios_locale"         # string offset=1115
.Linfo_string67:
	.asciz	"_M_impl"               # string offset=1129
.Linfo_string68:
	.asciz	"_M_facets"             # string offset=1137
.Linfo_string69:
	.asciz	"_ZNSt6locale5facetE"   # string offset=1147
.Linfo_string70:
	.asciz	"_M_facets_size"        # string offset=1167
.Linfo_string71:
	.asciz	"unsigned long"         # string offset=1182
.Linfo_string72:
	.asciz	"_M_caches"             # string offset=1196
.Linfo_string73:
	.asciz	"_M_names"              # string offset=1206
.Linfo_string74:
	.asciz	"signed char"           # string offset=1215
.Linfo_string75:
	.asciz	"_ZNSt6locale5_ImplE"   # string offset=1227
.Linfo_string76:
	.asciz	"_ZSt6locale"           # string offset=1247
.Linfo_string77:
	.asciz	"_ZSt8ios_base"         # string offset=1259
.Linfo_string78:
	.asciz	"_M_tie"                # string offset=1273
.Linfo_string79:
	.asciz	"_M_fill"               # string offset=1280
.Linfo_string80:
	.asciz	"_M_fill_init"          # string offset=1288
.Linfo_string81:
	.asciz	"_M_streambuf"          # string offset=1301
.Linfo_string82:
	.asciz	"_M_in_beg"             # string offset=1314
.Linfo_string83:
	.asciz	"_M_in_cur"             # string offset=1324
.Linfo_string84:
	.asciz	"_M_in_end"             # string offset=1334
.Linfo_string85:
	.asciz	"_M_out_beg"            # string offset=1344
.Linfo_string86:
	.asciz	"_M_out_cur"            # string offset=1355
.Linfo_string87:
	.asciz	"_M_out_end"            # string offset=1366
.Linfo_string88:
	.asciz	"_M_buf_locale"         # string offset=1377
.Linfo_string89:
	.asciz	"_ZSt15basic_streambufIcSt11char_traitsIcEE" # string offset=1391
.Linfo_string90:
	.asciz	"_M_ctype"              # string offset=1434
.Linfo_string91:
	.asciz	"__b_NSt6locale5facetE" # string offset=1443
.Linfo_string92:
	.asciz	"__SO__NSt6locale5facetE" # string offset=1465
.Linfo_string93:
	.asciz	"_M_c_locale_ctype"     # string offset=1489
.Linfo_string94:
	.asciz	"__locales"             # string offset=1507
.Linfo_string95:
	.asciz	"__locale_data"         # string offset=1517
.Linfo_string96:
	.asciz	"__ctype_b"             # string offset=1531
.Linfo_string97:
	.asciz	"unsigned short"        # string offset=1541
.Linfo_string98:
	.asciz	"__ctype_tolower"       # string offset=1556
.Linfo_string99:
	.asciz	"__ctype_toupper"       # string offset=1572
.Linfo_string100:
	.asciz	"__names"               # string offset=1588
.Linfo_string101:
	.asciz	"__locale_struct"       # string offset=1596
.Linfo_string102:
	.asciz	"_M_del"                # string offset=1612
.Linfo_string103:
	.asciz	"_M_toupper"            # string offset=1619
.Linfo_string104:
	.asciz	"_M_tolower"            # string offset=1630
.Linfo_string105:
	.asciz	"_M_table"              # string offset=1641
.Linfo_string106:
	.asciz	"_M_widen_ok"           # string offset=1650
.Linfo_string107:
	.asciz	"_M_widen"              # string offset=1662
.Linfo_string108:
	.asciz	"_M_narrow"             # string offset=1671
.Linfo_string109:
	.asciz	"_M_narrow_ok"          # string offset=1681
.Linfo_string110:
	.asciz	"_ZSt5ctypeIcE"         # string offset=1694
.Linfo_string111:
	.asciz	"_M_num_put"            # string offset=1708
.Linfo_string112:
	.asciz	"_ZSt7num_putIcSt19ostreambuf_iteratorIcSt11char_traitsIcEEE" # string offset=1719
.Linfo_string113:
	.asciz	"_M_num_get"            # string offset=1779
.Linfo_string114:
	.asciz	"_ZSt7num_getIcSt19istreambuf_iteratorIcSt11char_traitsIcEEE" # string offset=1790
.Linfo_string115:
	.asciz	"_ZSt9basic_iosIcSt11char_traitsIcEE" # string offset=1850
.Linfo_string116:
	.asciz	"_ZSo"                  # string offset=1886
.Linfo_string117:
	.asciz	"ostream"               # string offset=1891
.Linfo_string118:
	.asciz	"_ZTV9T3D_fix23"        # string offset=1899
.Linfo_string119:
	.asciz	"_ZTV9T3D_fix13"        # string offset=1914
.Linfo_string120:
	.asciz	"_ZTV11T2DFunction"     # string offset=1929
.Linfo_string121:
	.asciz	"_ZTV8T3D_fix3"         # string offset=1947
.Linfo_string122:
	.asciz	"_ZTV8T3D_fix1"         # string offset=1961
.Linfo_string123:
	.asciz	"_ZTV8T3D_fix2"         # string offset=1975
.Linfo_string124:
	.asciz	"__I___37_backgroundfield_integratefunction_cpp_f19cb1ee" # string offset=1989
.Linfo_string125:
	.asciz	"_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE" # string offset=2045
.Linfo_string126:
	.asciz	"_ZNSt8ios_base4InitE"  # string offset=2122
.Linfo_string127:
	.asciz	"__dso_handle"          # string offset=2143
.Linfo_string128:
	.asciz	"WID"                   # string offset=2156
.Linfo_string129:
	.asciz	"WID2"                  # string offset=2160
.Linfo_string130:
	.asciz	"WID3"                  # string offset=2165
.Linfo_string131:
	.asciz	"_ZTI11T1DFunction"     # string offset=2170
.Linfo_string132:
	.asciz	"__class_type_info"     # string offset=2188
.Linfo_string133:
	.asciz	"_ZTI11T2DFunction"     # string offset=2206
.Linfo_string134:
	.asciz	"_ZTI8T3D_fix1"         # string offset=2224
.Linfo_string135:
	.asciz	"__si_class_type_info"  # string offset=2238
.Linfo_string136:
	.asciz	"_ZTI8T3D_fix2"         # string offset=2259
.Linfo_string137:
	.asciz	"_ZTI8T3D_fix3"         # string offset=2273
.Linfo_string138:
	.asciz	"_ZTI9T3D_fix12"        # string offset=2287
.Linfo_string139:
	.asciz	"_ZTI9T3D_fix13"        # string offset=2302
.Linfo_string140:
	.asciz	"_ZTI9T3D_fix23"        # string offset=2317
.Linfo_string141:
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # string offset=2332
.Linfo_string142:
	.asciz	"_ZTS11T1DFunction"     # string offset=2370
.Linfo_string143:
	.asciz	"_ZTS11T2DFunction"     # string offset=2388
.Linfo_string144:
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # string offset=2406
.Linfo_string145:
	.asciz	"_ZTS8T3D_fix1"         # string offset=2447
.Linfo_string146:
	.asciz	"_ZTS8T3D_fix2"         # string offset=2461
.Linfo_string147:
	.asciz	"_ZTS8T3D_fix3"         # string offset=2475
.Linfo_string148:
	.asciz	"_ZTS9T3D_fix12"        # string offset=2489
.Linfo_string149:
	.asciz	"_ZTS9T3D_fix13"        # string offset=2504
.Linfo_string150:
	.asciz	"_ZTS9T3D_fix23"        # string offset=2519
.Linfo_string151:
	.asciz	"_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1ee5vmesh15INVALID_LOCALIDE" # string offset=2534
.Linfo_string152:
	.asciz	"unsigned"              # string offset=2623
.Linfo_string153:
	.asciz	"T1DFunction"           # string offset=2632
.Linfo_string154:
	.asciz	"T2DFunction"           # string offset=2644
.Linfo_string155:
	.asciz	"f"                     # string offset=2656
.Linfo_string156:
	.asciz	"T3DFunction"           # string offset=2658
.Linfo_string157:
	.asciz	"x"                     # string offset=2670
.Linfo_string158:
	.asciz	"double"                # string offset=2672
.Linfo_string159:
	.asciz	"call"                  # string offset=2679
.Linfo_string160:
	.asciz	"T3D_fix1"              # string offset=2684
.Linfo_string161:
	.asciz	"y"                     # string offset=2693
.Linfo_string162:
	.asciz	"T3D_fix2"              # string offset=2695
.Linfo_string163:
	.asciz	"z"                     # string offset=2704
.Linfo_string164:
	.asciz	"T3D_fix3"              # string offset=2706
.Linfo_string165:
	.asciz	"T3D_fix12"             # string offset=2715
.Linfo_string166:
	.asciz	"T3D_fix13"             # string offset=2725
.Linfo_string167:
	.asciz	"T3D_fix23"             # string offset=2735
.Linfo_string168:
	.asciz	"lineAverage"           # string offset=2745
.Linfo_string169:
	.asciz	"surfaceAverage"        # string offset=2757
.Linfo_string170:
	.asciz	"volumeAverage"         # string offset=2772
.Linfo_string171:
	.asciz	"~T1DFunction"          # string offset=2786
.Linfo_string172:
	.asciz	"_ZN11T1DFunctionD0Ev"  # string offset=2799
.Linfo_string173:
	.asciz	"~T2DFunction"          # string offset=2820
.Linfo_string174:
	.asciz	"_ZN11T2DFunctionD0Ev"  # string offset=2833
.Linfo_string175:
	.asciz	"~T3D_fix1"             # string offset=2854
.Linfo_string176:
	.asciz	"_ZN8T3D_fix1D0Ev"      # string offset=2864
.Linfo_string177:
	.asciz	"~T3D_fix2"             # string offset=2881
.Linfo_string178:
	.asciz	"_ZN8T3D_fix2D0Ev"      # string offset=2891
.Linfo_string179:
	.asciz	"~T3D_fix3"             # string offset=2908
.Linfo_string180:
	.asciz	"_ZN8T3D_fix3D0Ev"      # string offset=2918
.Linfo_string181:
	.asciz	"~T3D_fix12"            # string offset=2935
.Linfo_string182:
	.asciz	"_ZN9T3D_fix12D0Ev"     # string offset=2946
.Linfo_string183:
	.asciz	"~T3D_fix13"            # string offset=2964
.Linfo_string184:
	.asciz	"_ZN9T3D_fix13D0Ev"     # string offset=2975
.Linfo_string185:
	.asciz	"~T3D_fix23"            # string offset=2993
.Linfo_string186:
	.asciz	"_ZN9T3D_fix23D0Ev"     # string offset=3004
.Linfo_string187:
	.asciz	"__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee" # string offset=3022
.Linfo_string188:
	.asciz	"f1"                    # string offset=3080
.Linfo_string189:
	.asciz	"line"                  # string offset=3083
.Linfo_string190:
	.asciz	"X"                     # string offset=3088
.Linfo_string191:
	.asciz	"Y"                     # string offset=3090
.Linfo_string192:
	.asciz	"Z"                     # string offset=3092
.Linfo_string193:
	.asciz	"coordinate"            # string offset=3094
.Linfo_string194:
	.asciz	"accuracy"              # string offset=3105
.Linfo_string195:
	.asciz	"r1"                    # string offset=3114
.Linfo_string196:
	.asciz	"L"                     # string offset=3117
.Linfo_string197:
	.asciz	"norm"                  # string offset=3119
.Linfo_string198:
	.asciz	"acc"                   # string offset=3124
.Linfo_string199:
	.asciz	"a"                     # string offset=3128
.Linfo_string200:
	.asciz	"b"                     # string offset=3130
.Linfo_string201:
	.asciz	"value"                 # string offset=3132
.Linfo_string202:
	.asciz	"face"                  # string offset=3138
.Linfo_string203:
	.asciz	"L1"                    # string offset=3143
.Linfo_string204:
	.asciz	"L2"                    # string offset=3146
.Linfo_string205:
	.asciz	"r2"                    # string offset=3149
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp17-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp22-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp26-.Lfunc_begin0
	.quad	.Ltmp27-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc1:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp26-.Lfunc_begin0
	.quad	.Ltmp28-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	0
	.quad	0
.Ldebug_loc2:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc3:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp26-.Lfunc_begin0
	.quad	.Ltmp29-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	0
	.quad	0
.Ldebug_loc4:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp5-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc5:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp5-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc6:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp26-.Lfunc_begin0
	.quad	.Ltmp29-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	0
	.quad	0
.Ldebug_loc7:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc8:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp26-.Lfunc_begin0
	.quad	.Ltmp28-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	0
	.quad	0
.Ldebug_loc9:
	.quad	.Ltmp2-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	120                     # -8
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Ltmp20-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	120                     # -8
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp23-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp23-.Lfunc_begin0
	.quad	.Ltmp26-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	120                     # -8
	.quad	0
	.quad	0
.Ldebug_loc10:
	.quad	.Ltmp3-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc11:
	.quad	.Ltmp4-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc12:
	.quad	.Ltmp5-.Lfunc_begin0
	.quad	.Ltmp12-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp14-.Lfunc_begin0
	.quad	.Ltmp18-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp24-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc13:
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	.Ltmp14-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp19-.Lfunc_begin0
	.quad	.Ltmp20-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp25-.Lfunc_begin0
	.quad	.Ltmp26-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp30-.Lfunc_begin0
	.quad	.Lfunc_end0-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc14:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp42-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp49-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp56-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp63-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc15:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp43-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp50-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp57-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp64-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	0
	.quad	0
.Ldebug_loc16:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp33-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc17:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp43-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp50-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp57-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp65-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	0
	.quad	0
.Ldebug_loc18:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp40-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp47-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp54-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp66-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc19:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp32-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp48-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp55-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp66-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc20:
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp48-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp55-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp66-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc21:
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp40-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp47-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp54-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp66-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc22:
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp33-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc23:
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp43-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp50-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp57-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp64-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # super-register DW_OP_reg4
	.quad	0
	.quad	0
.Ldebug_loc24:
	.quad	.Ltmp60-.Lfunc_begin0
	.quad	.Ltmp62-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc25:
	.quad	.Ltmp34-.Lfunc_begin0
	.quad	.Ltmp43-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp50-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp57-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp62-.Lfunc_begin0
	.quad	.Ltmp66-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	0
	.quad	0
.Ldebug_loc26:
	.quad	.Ltmp61-.Lfunc_begin0
	.quad	.Ltmp62-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc27:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp73-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc28:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp69-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc29:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp73-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # DW_OP_reg4
	.quad	0
	.quad	0
.Ldebug_loc30:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp73-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	0
	.quad	0
.Ldebug_loc31:
	.quad	.Ltmp68-.Lfunc_begin0
	.quad	.Ltmp73-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc32:
	.quad	.Ltmp68-.Lfunc_begin0
	.quad	.Ltmp69-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc33:
	.quad	.Ltmp68-.Lfunc_begin0
	.quad	.Ltmp73-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	84                      # DW_OP_reg4
	.quad	0
	.quad	0
.Ldebug_loc34:
	.quad	.Ltmp68-.Lfunc_begin0
	.quad	.Ltmp73-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	81                      # DW_OP_reg1
	.quad	0
	.quad	0
.Ldebug_loc35:
	.quad	.Ltmp70-.Lfunc_begin0
	.quad	.Ltmp73-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	0
	.quad	0
.Ldebug_loc36:
	.quad	.Ltmp71-.Lfunc_begin0
	.quad	.Ltmp72-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp72-.Lfunc_begin0
	.quad	.Lfunc_end2-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	120                     # -8
	.quad	0
	.quad	0
.Ldebug_loc37:
	.quad	.Lfunc_begin7-.Lfunc_begin0
	.quad	.Ltmp86-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc38:
	.quad	.Lfunc_begin7-.Lfunc_begin0
	.quad	.Ltmp85-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp85-.Lfunc_begin0
	.quad	.Ltmp87-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc39:
	.quad	.Lfunc_begin7-.Lfunc_begin0
	.quad	.Ltmp85-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp85-.Lfunc_begin0
	.quad	.Ltmp87-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc40:
	.quad	.Ltmp84-.Lfunc_begin0
	.quad	.Ltmp86-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc41:
	.quad	.Ltmp84-.Lfunc_begin0
	.quad	.Ltmp85-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp85-.Lfunc_begin0
	.quad	.Ltmp87-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc42:
	.quad	.Ltmp84-.Lfunc_begin0
	.quad	.Ltmp85-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp85-.Lfunc_begin0
	.quad	.Ltmp87-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc43:
	.quad	.Lfunc_begin10-.Lfunc_begin0
	.quad	.Ltmp94-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc44:
	.quad	.Lfunc_begin10-.Lfunc_begin0
	.quad	.Ltmp95-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc45:
	.quad	.Lfunc_begin10-.Lfunc_begin0
	.quad	.Ltmp93-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp93-.Lfunc_begin0
	.quad	.Ltmp95-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc46:
	.quad	.Ltmp92-.Lfunc_begin0
	.quad	.Ltmp94-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc47:
	.quad	.Ltmp92-.Lfunc_begin0
	.quad	.Ltmp95-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc48:
	.quad	.Ltmp92-.Lfunc_begin0
	.quad	.Ltmp93-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp93-.Lfunc_begin0
	.quad	.Ltmp95-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc49:
	.quad	.Lfunc_begin13-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc50:
	.quad	.Lfunc_begin13-.Lfunc_begin0
	.quad	.Ltmp102-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc51:
	.quad	.Lfunc_begin13-.Lfunc_begin0
	.quad	.Ltmp102-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc52:
	.quad	.Ltmp100-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc53:
	.quad	.Ltmp100-.Lfunc_begin0
	.quad	.Ltmp102-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc54:
	.quad	.Ltmp100-.Lfunc_begin0
	.quad	.Ltmp102-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc55:
	.quad	.Lfunc_begin16-.Lfunc_begin0
	.quad	.Ltmp109-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc56:
	.quad	.Lfunc_begin16-.Lfunc_begin0
	.quad	.Ltmp108-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp108-.Lfunc_begin0
	.quad	.Ltmp110-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc57:
	.quad	.Ltmp107-.Lfunc_begin0
	.quad	.Ltmp109-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc58:
	.quad	.Ltmp107-.Lfunc_begin0
	.quad	.Ltmp108-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp108-.Lfunc_begin0
	.quad	.Ltmp110-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc59:
	.quad	.Lfunc_begin19-.Lfunc_begin0
	.quad	.Ltmp117-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc60:
	.quad	.Lfunc_begin19-.Lfunc_begin0
	.quad	.Ltmp116-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp116-.Lfunc_begin0
	.quad	.Ltmp118-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc61:
	.quad	.Ltmp115-.Lfunc_begin0
	.quad	.Ltmp117-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc62:
	.quad	.Ltmp115-.Lfunc_begin0
	.quad	.Ltmp116-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp116-.Lfunc_begin0
	.quad	.Ltmp118-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc63:
	.quad	.Lfunc_begin22-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc64:
	.quad	.Lfunc_begin22-.Lfunc_begin0
	.quad	.Ltmp125-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc65:
	.quad	.Ltmp123-.Lfunc_begin0
	.quad	.Ltmp124-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc66:
	.quad	.Ltmp123-.Lfunc_begin0
	.quad	.Ltmp125-.Lfunc_begin0
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
	.byte	57                      # DW_TAG_namespace
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
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
	.byte	22                      # DW_TAG_typedef
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	11                      # DW_AT_byte_size
	.byte	5                       # DW_FORM_data2
	.byte	58                      # DW_AT_decl_file
	.byte	11                      # DW_FORM_data1
	.byte	59                      # DW_AT_decl_line
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	12                      # Abbreviation Code
	.byte	19                      # DW_TAG_structure_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	11                      # DW_AT_byte_size
	.byte	5                       # DW_FORM_data2
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	13                      # Abbreviation Code
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
	.byte	14                      # Abbreviation Code
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
	.byte	15                      # Abbreviation Code
	.byte	4                       # DW_TAG_enumeration_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	11                      # DW_AT_byte_size
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	16                      # Abbreviation Code
	.byte	40                      # DW_TAG_enumerator
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	28                      # DW_AT_const_value
	.byte	13                      # DW_FORM_sdata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	17                      # Abbreviation Code
	.byte	21                      # DW_TAG_subroutine_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	18                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	19                      # Abbreviation Code
	.byte	59                      # DW_TAG_unspecified_type
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	20                      # Abbreviation Code
	.byte	19                      # DW_TAG_structure_type
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	11                      # DW_AT_byte_size
	.byte	11                      # DW_FORM_data1
	.ascii	"\210\001"              # DW_AT_alignment
	.byte	15                      # DW_FORM_udata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	21                      # Abbreviation Code
	.byte	33                      # DW_TAG_subrange_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	55                      # DW_AT_count
	.byte	5                       # DW_FORM_data2
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	22                      # Abbreviation Code
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
	.byte	23                      # Abbreviation Code
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
	.byte	24                      # Abbreviation Code
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
	.byte	25                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	26                      # Abbreviation Code
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
	.byte	27                      # Abbreviation Code
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
	.byte	28                      # Abbreviation Code
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
	.byte	29                      # Abbreviation Code
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
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
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	85                      # DW_AT_ranges
	.byte	6                       # DW_FORM_data4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	32                      # Abbreviation Code
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
	.byte	33                      # Abbreviation Code
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
	.byte	34                      # Abbreviation Code
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
	.byte	35                      # Abbreviation Code
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
	.byte	36                      # Abbreviation Code
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
	.byte	37                      # Abbreviation Code
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
	.byte	41                      # Abbreviation Code
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
	.byte	42                      # Abbreviation Code
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
	.byte	56                      # DW_AT_data_member_location
	.byte	10                      # DW_FORM_block1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	43                      # Abbreviation Code
	.byte	2                       # DW_TAG_class_type
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
	.byte	44                      # Abbreviation Code
	.byte	21                      # DW_TAG_subroutine_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	45                      # Abbreviation Code
	.byte	16                      # DW_TAG_reference_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	46                      # Abbreviation Code
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
	.byte	0                       # EOM(3)
	.section	.debug_info,"",@progbits
.Lcu_begin0:
	.long	.Ldebug_info_end0-.Ldebug_info_start0 # Length of Unit
.Ldebug_info_start0:
	.short	2                       # DWARF version number
	.long	.debug_abbrev           # Offset Into Abbrev. Section
	.byte	8                       # Address Size (in bytes)
	.byte	1                       # Abbrev [1] 0xb:0x13c8 DW_TAG_compile_unit
	.long	.Linfo_string0          # DW_AT_producer
	.short	4                       # DW_AT_language
	.long	.Linfo_string1          # DW_AT_name
	.long	.Lline_table_start0     # DW_AT_stmt_list
	.long	.Linfo_string2          # DW_AT_comp_dir
	.byte	1                       # DW_AT_GNU_pubnames
	.quad	.Lfunc_begin0           # DW_AT_low_pc
	.quad	.Lfunc_end25            # DW_AT_high_pc
	.byte	2                       # Abbrev [2] 0x2f:0x14 DW_TAG_variable
	.long	.Linfo_string3          # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV11T1DFunction
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
	.quad	_ZTV9T3D_fix12
	.byte	9                       # Abbrev [9] 0x85:0x1d DW_TAG_namespace
	.long	.Linfo_string7          # DW_AT_name
	.byte	10                      # Abbrev [10] 0x8a:0xa DW_TAG_variable
	.long	.Linfo_string8          # DW_AT_name
	.long	148                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	11                      # Abbrev [11] 0x94:0xd DW_TAG_typedef
	.long	162                     # DW_AT_type
	.long	.Linfo_string117        # DW_AT_name
	.short	272                     # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	0                       # End Of Children Mark
	.byte	12                      # Abbrev [12] 0xa2:0x23 DW_TAG_structure_type
	.long	.Linfo_string116        # DW_AT_name
	.short	272                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0xaa:0xd DW_TAG_member
	.long	.Linfo_string9          # DW_AT_name
	.long	79                      # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0xb7:0xd DW_TAG_member
	.long	.Linfo_string10         # DW_AT_name
	.long	197                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	0                       # End Of Children Mark
	.byte	12                      # Abbrev [12] 0xc5:0x78 DW_TAG_structure_type
	.long	.Linfo_string115        # DW_AT_name
	.short	264                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0xcd:0xd DW_TAG_member
	.long	.Linfo_string11         # DW_AT_name
	.long	317                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0xda:0xe DW_TAG_member
	.long	.Linfo_string78         # DW_AT_name
	.long	1055                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\330\001"
	.byte	13                      # Abbrev [13] 0xe8:0xe DW_TAG_member
	.long	.Linfo_string79         # DW_AT_name
	.long	1048                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\340\001"
	.byte	13                      # Abbrev [13] 0xf6:0xe DW_TAG_member
	.long	.Linfo_string80         # DW_AT_name
	.long	1048                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\341\001"
	.byte	13                      # Abbrev [13] 0x104:0xe DW_TAG_member
	.long	.Linfo_string81         # DW_AT_name
	.long	1060                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\350\001"
	.byte	13                      # Abbrev [13] 0x112:0xe DW_TAG_member
	.long	.Linfo_string90         # DW_AT_name
	.long	1177                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\360\001"
	.byte	13                      # Abbrev [13] 0x120:0xe DW_TAG_member
	.long	.Linfo_string111        # DW_AT_name
	.long	1502                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\370\001"
	.byte	13                      # Abbrev [13] 0x12e:0xe DW_TAG_member
	.long	.Linfo_string113        # DW_AT_name
	.long	1528                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\200\002"
	.byte	0                       # End Of Children Mark
	.byte	14                      # Abbrev [14] 0x13d:0xa7 DW_TAG_structure_type
	.long	.Linfo_string77         # DW_AT_name
	.byte	216                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x144:0xd DW_TAG_member
	.long	.Linfo_string9          # DW_AT_name
	.long	79                      # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x151:0xd DW_TAG_member
	.long	.Linfo_string12         # DW_AT_name
	.long	484                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	13                      # Abbrev [13] 0x15e:0xd DW_TAG_member
	.long	.Linfo_string14         # DW_AT_name
	.long	484                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	13                      # Abbrev [13] 0x16b:0xd DW_TAG_member
	.long	.Linfo_string15         # DW_AT_name
	.long	491                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	13                      # Abbrev [13] 0x178:0xd DW_TAG_member
	.long	.Linfo_string38         # DW_AT_name
	.long	649                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	28
	.byte	13                      # Abbrev [13] 0x185:0xd DW_TAG_member
	.long	.Linfo_string47         # DW_AT_name
	.long	649                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	13                      # Abbrev [13] 0x192:0xd DW_TAG_member
	.long	.Linfo_string48         # DW_AT_name
	.long	709                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	13                      # Abbrev [13] 0x19f:0xd DW_TAG_member
	.long	.Linfo_string58         # DW_AT_name
	.long	827                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	13                      # Abbrev [13] 0x1ac:0xd DW_TAG_member
	.long	.Linfo_string63         # DW_AT_name
	.long	871                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	64
	.byte	13                      # Abbrev [13] 0x1b9:0xe DW_TAG_member
	.long	.Linfo_string64         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\300\001"
	.byte	13                      # Abbrev [13] 0x1c7:0xe DW_TAG_member
	.long	.Linfo_string65         # DW_AT_name
	.long	883                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\310\001"
	.byte	13                      # Abbrev [13] 0x1d5:0xe DW_TAG_member
	.long	.Linfo_string66         # DW_AT_name
	.long	888                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\320\001"
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x1e4:0x7 DW_TAG_base_type
	.long	.Linfo_string13         # DW_AT_name
	.byte	5                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	15                      # Abbrev [15] 0x1eb:0x9e DW_TAG_enumeration_type
	.long	.Linfo_string37         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	4                       # DW_AT_alignment
	.byte	16                      # Abbrev [16] 0x1f2:0x6 DW_TAG_enumerator
	.long	.Linfo_string16         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x1f8:0x6 DW_TAG_enumerator
	.long	.Linfo_string17         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x1fe:0x6 DW_TAG_enumerator
	.long	.Linfo_string18         # DW_AT_name
	.byte	4                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x204:0x6 DW_TAG_enumerator
	.long	.Linfo_string19         # DW_AT_name
	.byte	8                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x20a:0x6 DW_TAG_enumerator
	.long	.Linfo_string20         # DW_AT_name
	.byte	16                      # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x210:0x6 DW_TAG_enumerator
	.long	.Linfo_string21         # DW_AT_name
	.byte	32                      # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x216:0x7 DW_TAG_enumerator
	.long	.Linfo_string22         # DW_AT_name
	.asciz	"\300"                  # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x21d:0x7 DW_TAG_enumerator
	.long	.Linfo_string23         # DW_AT_name
	.ascii	"\200\001"              # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x224:0x7 DW_TAG_enumerator
	.long	.Linfo_string24         # DW_AT_name
	.ascii	"\200\002"              # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x22b:0x7 DW_TAG_enumerator
	.long	.Linfo_string25         # DW_AT_name
	.ascii	"\200\004"              # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x232:0x7 DW_TAG_enumerator
	.long	.Linfo_string26         # DW_AT_name
	.ascii	"\200\b"                # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x239:0x7 DW_TAG_enumerator
	.long	.Linfo_string27         # DW_AT_name
	.ascii	"\200\020"              # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x240:0x7 DW_TAG_enumerator
	.long	.Linfo_string28         # DW_AT_name
	.ascii	"\200 "                 # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x247:0x8 DW_TAG_enumerator
	.long	.Linfo_string29         # DW_AT_name
	.asciz	"\200\300"              # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x24f:0x8 DW_TAG_enumerator
	.long	.Linfo_string30         # DW_AT_name
	.ascii	"\200\200\001"          # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x257:0x7 DW_TAG_enumerator
	.long	.Linfo_string31         # DW_AT_name
	.ascii	"\260\001"              # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x25e:0x7 DW_TAG_enumerator
	.long	.Linfo_string32         # DW_AT_name
	.asciz	"\312"                  # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x265:0x7 DW_TAG_enumerator
	.long	.Linfo_string33         # DW_AT_name
	.ascii	"\204\002"              # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x26c:0x8 DW_TAG_enumerator
	.long	.Linfo_string34         # DW_AT_name
	.ascii	"\200\200\004"          # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x274:0xa DW_TAG_enumerator
	.long	.Linfo_string35         # DW_AT_name
	.ascii	"\377\377\377\377\007"  # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x27e:0xa DW_TAG_enumerator
	.long	.Linfo_string36         # DW_AT_name
	.ascii	"\200\200\200\200x"     # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	15                      # Abbrev [15] 0x289:0x3c DW_TAG_enumeration_type
	.long	.Linfo_string46         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	4                       # DW_AT_alignment
	.byte	16                      # Abbrev [16] 0x290:0x6 DW_TAG_enumerator
	.long	.Linfo_string39         # DW_AT_name
	.byte	0                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x296:0x6 DW_TAG_enumerator
	.long	.Linfo_string40         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x29c:0x6 DW_TAG_enumerator
	.long	.Linfo_string41         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x2a2:0x6 DW_TAG_enumerator
	.long	.Linfo_string42         # DW_AT_name
	.byte	4                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x2a8:0x8 DW_TAG_enumerator
	.long	.Linfo_string43         # DW_AT_name
	.ascii	"\200\200\004"          # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x2b0:0xa DW_TAG_enumerator
	.long	.Linfo_string44         # DW_AT_name
	.ascii	"\377\377\377\377\007"  # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x2ba:0xa DW_TAG_enumerator
	.long	.Linfo_string45         # DW_AT_name
	.ascii	"\200\200\200\200x"     # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x2c5:0x5 DW_TAG_pointer_type
	.long	714                     # DW_AT_type
	.byte	14                      # Abbrev [14] 0x2ca:0x3c DW_TAG_structure_type
	.long	.Linfo_string57         # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x2d1:0xd DW_TAG_member
	.long	.Linfo_string49         # DW_AT_name
	.long	709                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x2de:0xd DW_TAG_member
	.long	.Linfo_string50         # DW_AT_name
	.long	774                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	13                      # Abbrev [13] 0x2eb:0xd DW_TAG_member
	.long	.Linfo_string55         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	13                      # Abbrev [13] 0x2f8:0xd DW_TAG_member
	.long	.Linfo_string56         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	20
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x306:0x5 DW_TAG_pointer_type
	.long	779                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x30b:0x11 DW_TAG_subroutine_type
	.byte	18                      # Abbrev [18] 0x30c:0x5 DW_TAG_formal_parameter
	.long	796                     # DW_AT_type
	.byte	18                      # Abbrev [18] 0x311:0x5 DW_TAG_formal_parameter
	.long	822                     # DW_AT_type
	.byte	18                      # Abbrev [18] 0x316:0x5 DW_TAG_formal_parameter
	.long	99                      # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	15                      # Abbrev [15] 0x31c:0x1a DW_TAG_enumeration_type
	.long	.Linfo_string54         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	4                       # DW_AT_alignment
	.byte	16                      # Abbrev [16] 0x323:0x6 DW_TAG_enumerator
	.long	.Linfo_string51         # DW_AT_name
	.byte	0                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x329:0x6 DW_TAG_enumerator
	.long	.Linfo_string52         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x32f:0x6 DW_TAG_enumerator
	.long	.Linfo_string53         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x336:0x5 DW_TAG_pointer_type
	.long	317                     # DW_AT_type
	.byte	14                      # Abbrev [14] 0x33b:0x22 DW_TAG_structure_type
	.long	.Linfo_string62         # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x342:0xd DW_TAG_member
	.long	.Linfo_string59         # DW_AT_name
	.long	861                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x34f:0xd DW_TAG_member
	.long	.Linfo_string61         # DW_AT_name
	.long	484                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x35d:0x5 DW_TAG_pointer_type
	.long	866                     # DW_AT_type
	.byte	19                      # Abbrev [19] 0x362:0x5 DW_TAG_unspecified_type
	.long	.Linfo_string60         # DW_AT_name
	.byte	3                       # Abbrev [3] 0x367:0xc DW_TAG_array_type
	.long	827                     # DW_AT_type
	.byte	4                       # Abbrev [4] 0x36c:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	8                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x373:0x5 DW_TAG_pointer_type
	.long	827                     # DW_AT_type
	.byte	14                      # Abbrev [14] 0x378:0x15 DW_TAG_structure_type
	.long	.Linfo_string76         # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x37f:0xd DW_TAG_member
	.long	.Linfo_string67         # DW_AT_name
	.long	909                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x38d:0x5 DW_TAG_pointer_type
	.long	914                     # DW_AT_type
	.byte	14                      # Abbrev [14] 0x392:0x49 DW_TAG_structure_type
	.long	.Linfo_string75         # DW_AT_name
	.byte	40                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x399:0xd DW_TAG_member
	.long	.Linfo_string56         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x3a6:0xd DW_TAG_member
	.long	.Linfo_string68         # DW_AT_name
	.long	987                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	13                      # Abbrev [13] 0x3b3:0xd DW_TAG_member
	.long	.Linfo_string70         # DW_AT_name
	.long	1031                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	13                      # Abbrev [13] 0x3c0:0xd DW_TAG_member
	.long	.Linfo_string72         # DW_AT_name
	.long	987                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	13                      # Abbrev [13] 0x3cd:0xd DW_TAG_member
	.long	.Linfo_string73         # DW_AT_name
	.long	1038                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x3db:0x5 DW_TAG_pointer_type
	.long	992                     # DW_AT_type
	.byte	5                       # Abbrev [5] 0x3e0:0x5 DW_TAG_pointer_type
	.long	997                     # DW_AT_type
	.byte	14                      # Abbrev [14] 0x3e5:0x22 DW_TAG_structure_type
	.long	.Linfo_string69         # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x3ec:0xd DW_TAG_member
	.long	.Linfo_string9          # DW_AT_name
	.long	79                      # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x3f9:0xd DW_TAG_member
	.long	.Linfo_string56         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x407:0x7 DW_TAG_base_type
	.long	.Linfo_string71         # DW_AT_name
	.byte	7                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	5                       # Abbrev [5] 0x40e:0x5 DW_TAG_pointer_type
	.long	1043                    # DW_AT_type
	.byte	5                       # Abbrev [5] 0x413:0x5 DW_TAG_pointer_type
	.long	1048                    # DW_AT_type
	.byte	7                       # Abbrev [7] 0x418:0x7 DW_TAG_base_type
	.long	.Linfo_string74         # DW_AT_name
	.byte	6                       # DW_AT_encoding
	.byte	1                       # DW_AT_byte_size
	.byte	5                       # Abbrev [5] 0x41f:0x5 DW_TAG_pointer_type
	.long	162                     # DW_AT_type
	.byte	5                       # Abbrev [5] 0x424:0x5 DW_TAG_pointer_type
	.long	1065                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x429:0x70 DW_TAG_structure_type
	.long	.Linfo_string89         # DW_AT_name
	.byte	64                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x430:0xd DW_TAG_member
	.long	.Linfo_string9          # DW_AT_name
	.long	79                      # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x43d:0xd DW_TAG_member
	.long	.Linfo_string82         # DW_AT_name
	.long	1043                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	13                      # Abbrev [13] 0x44a:0xd DW_TAG_member
	.long	.Linfo_string83         # DW_AT_name
	.long	1043                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	13                      # Abbrev [13] 0x457:0xd DW_TAG_member
	.long	.Linfo_string84         # DW_AT_name
	.long	1043                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	13                      # Abbrev [13] 0x464:0xd DW_TAG_member
	.long	.Linfo_string85         # DW_AT_name
	.long	1043                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	13                      # Abbrev [13] 0x471:0xd DW_TAG_member
	.long	.Linfo_string86         # DW_AT_name
	.long	1043                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	13                      # Abbrev [13] 0x47e:0xd DW_TAG_member
	.long	.Linfo_string87         # DW_AT_name
	.long	1043                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	13                      # Abbrev [13] 0x48b:0xd DW_TAG_member
	.long	.Linfo_string88         # DW_AT_name
	.long	888                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	56
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x499:0x5 DW_TAG_pointer_type
	.long	1182                    # DW_AT_type
	.byte	12                      # Abbrev [12] 0x49e:0x8d DW_TAG_structure_type
	.long	.Linfo_string110        # DW_AT_name
	.short	576                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x4a6:0xd DW_TAG_member
	.long	.Linfo_string91         # DW_AT_name
	.long	1323                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x4b3:0xd DW_TAG_member
	.long	.Linfo_string93         # DW_AT_name
	.long	1357                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	13                      # Abbrev [13] 0x4c0:0xd DW_TAG_member
	.long	.Linfo_string102        # DW_AT_name
	.long	1048                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	13                      # Abbrev [13] 0x4cd:0xd DW_TAG_member
	.long	.Linfo_string103        # DW_AT_name
	.long	1472                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	13                      # Abbrev [13] 0x4da:0xd DW_TAG_member
	.long	.Linfo_string104        # DW_AT_name
	.long	1472                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	13                      # Abbrev [13] 0x4e7:0xd DW_TAG_member
	.long	.Linfo_string105        # DW_AT_name
	.long	1460                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	13                      # Abbrev [13] 0x4f4:0xd DW_TAG_member
	.long	.Linfo_string106        # DW_AT_name
	.long	1048                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	56
	.byte	13                      # Abbrev [13] 0x501:0xd DW_TAG_member
	.long	.Linfo_string107        # DW_AT_name
	.long	1489                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	57
	.byte	13                      # Abbrev [13] 0x50e:0xe DW_TAG_member
	.long	.Linfo_string108        # DW_AT_name
	.long	1489                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\271\002"
	.byte	13                      # Abbrev [13] 0x51c:0xe DW_TAG_member
	.long	.Linfo_string109        # DW_AT_name
	.long	1048                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\271\004"
	.byte	0                       # End Of Children Mark
	.byte	14                      # Abbrev [14] 0x52b:0x22 DW_TAG_structure_type
	.long	.Linfo_string92         # DW_AT_name
	.byte	12                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x532:0xd DW_TAG_member
	.long	.Linfo_string9          # DW_AT_name
	.long	79                      # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x53f:0xd DW_TAG_member
	.long	.Linfo_string56         # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x54d:0x5 DW_TAG_pointer_type
	.long	1362                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x552:0x4a DW_TAG_structure_type
	.long	.Linfo_string101        # DW_AT_name
	.byte	232                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x559:0xd DW_TAG_member
	.long	.Linfo_string94         # DW_AT_name
	.long	1436                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	13                      # Abbrev [13] 0x566:0xd DW_TAG_member
	.long	.Linfo_string96         # DW_AT_name
	.long	1460                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	104
	.byte	13                      # Abbrev [13] 0x573:0xd DW_TAG_member
	.long	.Linfo_string98         # DW_AT_name
	.long	1472                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	112
	.byte	13                      # Abbrev [13] 0x580:0xd DW_TAG_member
	.long	.Linfo_string99         # DW_AT_name
	.long	1472                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	120
	.byte	13                      # Abbrev [13] 0x58d:0xe DW_TAG_member
	.long	.Linfo_string100        # DW_AT_name
	.long	1477                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\200\001"
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x59c:0xc DW_TAG_array_type
	.long	1448                    # DW_AT_type
	.byte	4                       # Abbrev [4] 0x5a1:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	13                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x5a8:0x5 DW_TAG_pointer_type
	.long	1453                    # DW_AT_type
	.byte	20                      # Abbrev [20] 0x5ad:0x7 DW_TAG_structure_type
	.long	.Linfo_string95         # DW_AT_name
	.byte	0                       # DW_AT_byte_size
	.byte	1                       # DW_AT_alignment
	.byte	5                       # Abbrev [5] 0x5b4:0x5 DW_TAG_pointer_type
	.long	1465                    # DW_AT_type
	.byte	7                       # Abbrev [7] 0x5b9:0x7 DW_TAG_base_type
	.long	.Linfo_string97         # DW_AT_name
	.byte	7                       # DW_AT_encoding
	.byte	2                       # DW_AT_byte_size
	.byte	5                       # Abbrev [5] 0x5c0:0x5 DW_TAG_pointer_type
	.long	99                      # DW_AT_type
	.byte	3                       # Abbrev [3] 0x5c5:0xc DW_TAG_array_type
	.long	1043                    # DW_AT_type
	.byte	4                       # Abbrev [4] 0x5ca:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	13                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x5d1:0xd DW_TAG_array_type
	.long	1048                    # DW_AT_type
	.byte	21                      # Abbrev [21] 0x5d6:0x7 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.short	256                     # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x5de:0x5 DW_TAG_pointer_type
	.long	1507                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x5e3:0x15 DW_TAG_structure_type
	.long	.Linfo_string112        # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x5ea:0xd DW_TAG_member
	.long	.Linfo_string91         # DW_AT_name
	.long	1323                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x5f8:0x5 DW_TAG_pointer_type
	.long	1533                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x5fd:0x15 DW_TAG_structure_type
	.long	.Linfo_string114        # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x604:0xd DW_TAG_member
	.long	.Linfo_string91         # DW_AT_name
	.long	1323                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	2                       # Abbrev [2] 0x612:0x14 DW_TAG_variable
	.long	.Linfo_string118        # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV9T3D_fix23
	.byte	2                       # Abbrev [2] 0x626:0x14 DW_TAG_variable
	.long	.Linfo_string119        # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV9T3D_fix13
	.byte	2                       # Abbrev [2] 0x63a:0x14 DW_TAG_variable
	.long	.Linfo_string120        # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV11T2DFunction
	.byte	2                       # Abbrev [2] 0x64e:0x14 DW_TAG_variable
	.long	.Linfo_string121        # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV8T3D_fix3
	.byte	2                       # Abbrev [2] 0x662:0x14 DW_TAG_variable
	.long	.Linfo_string122        # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV8T3D_fix1
	.byte	2                       # Abbrev [2] 0x676:0x14 DW_TAG_variable
	.long	.Linfo_string123        # DW_AT_name
	.long	67                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV8T3D_fix2
	.byte	22                      # Abbrev [22] 0x68a:0x17 DW_TAG_variable
	.long	.Linfo_string124        # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	1                       # DW_AT_decl_file
	.short	8148                    # DW_AT_decl_line
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	__I___37_backgroundfield_integratefunction_cpp_f19cb1ee
	.byte	23                      # Abbrev [23] 0x6a1:0x13 DW_TAG_variable
	.long	.Linfo_string125        # DW_AT_name
	.long	1716                    # DW_AT_type
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE
	.byte	14                      # Abbrev [14] 0x6b4:0x11 DW_TAG_structure_type
	.long	.Linfo_string126        # DW_AT_name
	.byte	1                       # DW_AT_byte_size
	.byte	1                       # DW_AT_alignment
	.byte	24                      # Abbrev [24] 0x6bb:0x9 DW_TAG_member
	.long	1733                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x6c5:0xc DW_TAG_array_type
	.long	1048                    # DW_AT_type
	.byte	4                       # Abbrev [4] 0x6ca:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	0                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	10                      # Abbrev [10] 0x6d1:0xa DW_TAG_variable
	.long	.Linfo_string127        # DW_AT_name
	.long	861                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	25                      # Abbrev [25] 0x6db:0x9 DW_TAG_variable
	.long	.Linfo_string128        # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	25                      # Abbrev [25] 0x6e4:0x9 DW_TAG_variable
	.long	.Linfo_string129        # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	25                      # Abbrev [25] 0x6ed:0x9 DW_TAG_variable
	.long	.Linfo_string130        # DW_AT_name
	.long	99                      # DW_AT_type
	.byte	2                       # Abbrev [2] 0x6f6:0x14 DW_TAG_variable
	.long	.Linfo_string131        # DW_AT_name
	.long	1802                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI11T1DFunction
	.byte	26                      # Abbrev [26] 0x70a:0xa DW_TAG_structure_type
	.long	.Linfo_string132        # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	2                       # Abbrev [2] 0x714:0x14 DW_TAG_variable
	.long	.Linfo_string133        # DW_AT_name
	.long	1802                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI11T2DFunction
	.byte	2                       # Abbrev [2] 0x728:0x14 DW_TAG_variable
	.long	.Linfo_string134        # DW_AT_name
	.long	1852                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI8T3D_fix1
	.byte	26                      # Abbrev [26] 0x73c:0xa DW_TAG_structure_type
	.long	.Linfo_string135        # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	2                       # Abbrev [2] 0x746:0x14 DW_TAG_variable
	.long	.Linfo_string136        # DW_AT_name
	.long	1852                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI8T3D_fix2
	.byte	2                       # Abbrev [2] 0x75a:0x14 DW_TAG_variable
	.long	.Linfo_string137        # DW_AT_name
	.long	1852                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI8T3D_fix3
	.byte	2                       # Abbrev [2] 0x76e:0x14 DW_TAG_variable
	.long	.Linfo_string138        # DW_AT_name
	.long	1852                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI9T3D_fix12
	.byte	2                       # Abbrev [2] 0x782:0x14 DW_TAG_variable
	.long	.Linfo_string139        # DW_AT_name
	.long	1852                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI9T3D_fix13
	.byte	2                       # Abbrev [2] 0x796:0x14 DW_TAG_variable
	.long	.Linfo_string140        # DW_AT_name
	.long	1852                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI9T3D_fix23
	.byte	10                      # Abbrev [10] 0x7aa:0xa DW_TAG_variable
	.long	.Linfo_string141        # DW_AT_name
	.long	1972                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	3                       # Abbrev [3] 0x7b4:0xc DW_TAG_array_type
	.long	79                      # DW_AT_type
	.byte	4                       # Abbrev [4] 0x7b9:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	0                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	2                       # Abbrev [2] 0x7c0:0x14 DW_TAG_variable
	.long	.Linfo_string142        # DW_AT_name
	.long	2004                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS11T1DFunction
	.byte	3                       # Abbrev [3] 0x7d4:0xc DW_TAG_array_type
	.long	1048                    # DW_AT_type
	.byte	4                       # Abbrev [4] 0x7d9:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	14                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	2                       # Abbrev [2] 0x7e0:0x14 DW_TAG_variable
	.long	.Linfo_string143        # DW_AT_name
	.long	2004                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS11T2DFunction
	.byte	10                      # Abbrev [10] 0x7f4:0xa DW_TAG_variable
	.long	.Linfo_string144        # DW_AT_name
	.long	1972                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	2                       # Abbrev [2] 0x7fe:0x14 DW_TAG_variable
	.long	.Linfo_string145        # DW_AT_name
	.long	2066                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS8T3D_fix1
	.byte	3                       # Abbrev [3] 0x812:0xc DW_TAG_array_type
	.long	1048                    # DW_AT_type
	.byte	4                       # Abbrev [4] 0x817:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	10                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	2                       # Abbrev [2] 0x81e:0x14 DW_TAG_variable
	.long	.Linfo_string146        # DW_AT_name
	.long	2066                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS8T3D_fix2
	.byte	2                       # Abbrev [2] 0x832:0x14 DW_TAG_variable
	.long	.Linfo_string147        # DW_AT_name
	.long	2066                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS8T3D_fix3
	.byte	2                       # Abbrev [2] 0x846:0x14 DW_TAG_variable
	.long	.Linfo_string148        # DW_AT_name
	.long	2138                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS9T3D_fix12
	.byte	3                       # Abbrev [3] 0x85a:0xc DW_TAG_array_type
	.long	1048                    # DW_AT_type
	.byte	4                       # Abbrev [4] 0x85f:0x6 DW_TAG_subrange_type
	.long	106                     # DW_AT_type
	.byte	11                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	2                       # Abbrev [2] 0x866:0x14 DW_TAG_variable
	.long	.Linfo_string149        # DW_AT_name
	.long	2138                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS9T3D_fix13
	.byte	2                       # Abbrev [2] 0x87a:0x14 DW_TAG_variable
	.long	.Linfo_string150        # DW_AT_name
	.long	2138                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS9T3D_fix23
	.byte	25                      # Abbrev [25] 0x88e:0x9 DW_TAG_variable
	.long	.Linfo_string151        # DW_AT_name
	.long	2199                    # DW_AT_type
	.byte	7                       # Abbrev [7] 0x897:0x7 DW_TAG_base_type
	.long	.Linfo_string152        # DW_AT_name
	.byte	7                       # DW_AT_encoding
	.byte	4                       # DW_AT_byte_size
	.byte	27                      # Abbrev [27] 0x89e:0x12e DW_TAG_subprogram
	.quad	.Lfunc_begin0           # DW_AT_low_pc
	.quad	.Lfunc_end0             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string168        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	43                      # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	28                      # Abbrev [28] 0x8bc:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc0            # DW_AT_location
	.long	.Linfo_string188        # DW_AT_name
	.long	5004                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x8c9:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc1            # DW_AT_location
	.long	.Linfo_string189        # DW_AT_name
	.long	5030                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x8d6:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc2            # DW_AT_location
	.long	.Linfo_string194        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x8e3:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc3            # DW_AT_location
	.long	.Linfo_string195        # DW_AT_name
	.long	5059                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x8f0:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc4            # DW_AT_location
	.long	.Linfo_string196        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0x8fd:0xce DW_TAG_lexical_block
	.quad	.Ltmp1                  # DW_AT_low_pc
	.quad	.Ltmp31                 # DW_AT_high_pc
	.byte	25                      # Abbrev [25] 0x90e:0x9 DW_TAG_variable
	.long	.Linfo_string188        # DW_AT_name
	.long	5004                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x917:0xd DW_TAG_variable
	.long	.Ldebug_loc5            # DW_AT_location
	.long	.Linfo_string196        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x924:0xd DW_TAG_variable
	.long	.Ldebug_loc6            # DW_AT_location
	.long	.Linfo_string195        # DW_AT_name
	.long	5059                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x931:0xd DW_TAG_variable
	.long	.Ldebug_loc7            # DW_AT_location
	.long	.Linfo_string194        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x93e:0xd DW_TAG_variable
	.long	.Ldebug_loc8            # DW_AT_location
	.long	.Linfo_string189        # DW_AT_name
	.long	5030                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x94b:0xd DW_TAG_variable
	.long	.Ldebug_loc13           # DW_AT_location
	.long	.Linfo_string201        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x958:0x72 DW_TAG_lexical_block
	.long	.Ldebug_ranges3         # DW_AT_ranges
	.byte	30                      # Abbrev [30] 0x95d:0xd DW_TAG_variable
	.long	.Ldebug_loc9            # DW_AT_location
	.long	.Linfo_string197        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x96a:0xd DW_TAG_variable
	.long	.Ldebug_loc10           # DW_AT_location
	.long	.Linfo_string198        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x977:0xd DW_TAG_variable
	.long	.Ldebug_loc11           # DW_AT_location
	.long	.Linfo_string199        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x984:0xd DW_TAG_variable
	.long	.Ldebug_loc12           # DW_AT_location
	.long	.Linfo_string200        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x991:0x13 DW_TAG_lexical_block
	.long	.Ldebug_ranges0         # DW_AT_ranges
	.byte	32                      # Abbrev [32] 0x996:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\230\177"
	.long	.Linfo_string155        # DW_AT_name
	.long	4217                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	31                      # Abbrev [31] 0x9a4:0x13 DW_TAG_lexical_block
	.long	.Ldebug_ranges1         # DW_AT_ranges
	.byte	32                      # Abbrev [32] 0x9a9:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\270\177"
	.long	.Linfo_string155        # DW_AT_name
	.long	4789                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	31                      # Abbrev [31] 0x9b7:0x12 DW_TAG_lexical_block
	.long	.Ldebug_ranges2         # DW_AT_ranges
	.byte	32                      # Abbrev [32] 0x9bc:0xc DW_TAG_variable
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	88
	.long	.Linfo_string155        # DW_AT_name
	.long	4503                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	27                      # Abbrev [27] 0x9cc:0x129 DW_TAG_subprogram
	.quad	.Lfunc_begin1           # DW_AT_low_pc
	.quad	.Lfunc_end1             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string169        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	88                      # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	28                      # Abbrev [28] 0x9ea:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc14           # DW_AT_location
	.long	.Linfo_string188        # DW_AT_name
	.long	5004                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x9f7:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc15           # DW_AT_location
	.long	.Linfo_string202        # DW_AT_name
	.long	5030                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xa04:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc16           # DW_AT_location
	.long	.Linfo_string194        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xa11:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc17           # DW_AT_location
	.long	.Linfo_string195        # DW_AT_name
	.long	5059                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xa1e:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc18           # DW_AT_location
	.long	.Linfo_string203        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xa2b:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc19           # DW_AT_location
	.long	.Linfo_string204        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0xa38:0xbc DW_TAG_lexical_block
	.quad	.Ltmp32                 # DW_AT_low_pc
	.quad	.Ltmp67                 # DW_AT_high_pc
	.byte	30                      # Abbrev [30] 0xa49:0xd DW_TAG_variable
	.long	.Ldebug_loc20           # DW_AT_location
	.long	.Linfo_string204        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xa56:0xd DW_TAG_variable
	.long	.Ldebug_loc21           # DW_AT_location
	.long	.Linfo_string203        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	25                      # Abbrev [25] 0xa63:0x9 DW_TAG_variable
	.long	.Linfo_string195        # DW_AT_name
	.long	5059                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xa6c:0xd DW_TAG_variable
	.long	.Ldebug_loc22           # DW_AT_location
	.long	.Linfo_string194        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xa79:0xd DW_TAG_variable
	.long	.Ldebug_loc23           # DW_AT_location
	.long	.Linfo_string202        # DW_AT_name
	.long	5030                    # DW_AT_type
	.byte	25                      # Abbrev [25] 0xa86:0x9 DW_TAG_variable
	.long	.Linfo_string188        # DW_AT_name
	.long	5004                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xa8f:0xd DW_TAG_variable
	.long	.Ldebug_loc26           # DW_AT_location
	.long	.Linfo_string201        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xa9c:0x57 DW_TAG_lexical_block
	.long	.Ldebug_ranges7         # DW_AT_ranges
	.byte	30                      # Abbrev [30] 0xaa1:0xd DW_TAG_variable
	.long	.Ldebug_loc24           # DW_AT_location
	.long	.Linfo_string197        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xaae:0xd DW_TAG_variable
	.long	.Ldebug_loc25           # DW_AT_location
	.long	.Linfo_string198        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xabb:0x13 DW_TAG_lexical_block
	.long	.Ldebug_ranges4         # DW_AT_ranges
	.byte	32                      # Abbrev [32] 0xac0:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\250\177"
	.long	.Linfo_string155        # DW_AT_name
	.long	3941                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	31                      # Abbrev [31] 0xace:0x12 DW_TAG_lexical_block
	.long	.Ldebug_ranges5         # DW_AT_ranges
	.byte	32                      # Abbrev [32] 0xad3:0xc DW_TAG_variable
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	88
	.long	.Linfo_string155        # DW_AT_name
	.long	3639                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	31                      # Abbrev [31] 0xae0:0x12 DW_TAG_lexical_block
	.long	.Ldebug_ranges6         # DW_AT_ranges
	.byte	32                      # Abbrev [32] 0xae5:0xc DW_TAG_variable
	.byte	2                       # DW_AT_location
	.byte	145
	.byte	64
	.long	.Linfo_string155        # DW_AT_name
	.long	3315                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	27                      # Abbrev [27] 0xaf5:0xc5 DW_TAG_subprogram
	.quad	.Lfunc_begin2           # DW_AT_low_pc
	.quad	.Lfunc_end2             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string170        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	128                     # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	28                      # Abbrev [28] 0xb13:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc27           # DW_AT_location
	.long	.Linfo_string188        # DW_AT_name
	.long	5004                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xb20:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc28           # DW_AT_location
	.long	.Linfo_string194        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xb2d:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc29           # DW_AT_location
	.long	.Linfo_string195        # DW_AT_name
	.long	5059                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xb3a:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc30           # DW_AT_location
	.long	.Linfo_string205        # DW_AT_name
	.long	5059                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0xb47:0x72 DW_TAG_lexical_block
	.quad	.Ltmp68                 # DW_AT_low_pc
	.quad	.Ltmp75                 # DW_AT_high_pc
	.byte	30                      # Abbrev [30] 0xb58:0xd DW_TAG_variable
	.long	.Ldebug_loc31           # DW_AT_location
	.long	.Linfo_string188        # DW_AT_name
	.long	5004                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xb65:0xd DW_TAG_variable
	.long	.Ldebug_loc32           # DW_AT_location
	.long	.Linfo_string194        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xb72:0xd DW_TAG_variable
	.long	.Ldebug_loc33           # DW_AT_location
	.long	.Linfo_string195        # DW_AT_name
	.long	5059                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xb7f:0xd DW_TAG_variable
	.long	.Ldebug_loc34           # DW_AT_location
	.long	.Linfo_string205        # DW_AT_name
	.long	5059                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0xb8c:0x2c DW_TAG_lexical_block
	.quad	.Ltmp68                 # DW_AT_low_pc
	.quad	.Ltmp74                 # DW_AT_high_pc
	.byte	30                      # Abbrev [30] 0xb9d:0xd DW_TAG_variable
	.long	.Ldebug_loc35           # DW_AT_location
	.long	.Linfo_string198        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xbaa:0xd DW_TAG_variable
	.long	.Ldebug_loc36           # DW_AT_location
	.long	.Linfo_string197        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	33                      # Abbrev [33] 0xbba:0x48 DW_TAG_class_type
	.long	.Linfo_string153        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	34                      # Abbrev [34] 0xbc4:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin3           # DW_AT_low_pc
	.quad	.Lfunc_end3             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string171        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	29                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	35                      # Abbrev [35] 0xbde:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	5064                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0xbe6:0x1a DW_TAG_lexical_block
	.quad	.Ltmp76                 # DW_AT_low_pc
	.quad	.Ltmp77                 # DW_AT_high_pc
	.byte	36                      # Abbrev [36] 0xbf7:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	5064                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	37                      # Abbrev [37] 0xc02:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin4           # DW_AT_low_pc
	.quad	.Lfunc_end4             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string172        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	33                      # Abbrev [33] 0xc1a:0x48 DW_TAG_class_type
	.long	.Linfo_string154        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	34                      # Abbrev [34] 0xc24:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin5           # DW_AT_low_pc
	.quad	.Lfunc_end5             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string173        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	30                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	35                      # Abbrev [35] 0xc3e:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	5069                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0xc46:0x1a DW_TAG_lexical_block
	.quad	.Ltmp80                 # DW_AT_low_pc
	.quad	.Ltmp81                 # DW_AT_high_pc
	.byte	36                      # Abbrev [36] 0xc57:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	5069                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	37                      # Abbrev [37] 0xc62:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin6           # DW_AT_low_pc
	.quad	.Lfunc_end6             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string174        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0xc7a:0x79 DW_TAG_subprogram
	.quad	.Lfunc_begin7           # DW_AT_low_pc
	.quad	.Lfunc_end7             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string159        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	63                      # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	38                      # Abbrev [38] 0xc98:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc37           # DW_AT_location
	.long	3489                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	28                      # Abbrev [28] 0xca2:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc38           # DW_AT_location
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xcaf:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc39           # DW_AT_location
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0xcbc:0x36 DW_TAG_lexical_block
	.quad	.Ltmp85                 # DW_AT_low_pc
	.quad	.Ltmp87                 # DW_AT_high_pc
	.byte	39                      # Abbrev [39] 0xccd:0xa DW_TAG_variable
	.long	.Ldebug_loc40           # DW_AT_location
	.long	3489                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	30                      # Abbrev [30] 0xcd7:0xd DW_TAG_variable
	.long	.Ldebug_loc41           # DW_AT_location
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xce4:0xd DW_TAG_variable
	.long	.Ldebug_loc42           # DW_AT_location
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	33                      # Abbrev [33] 0xcf3:0x83 DW_TAG_class_type
	.long	.Linfo_string160        # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0xcfd:0xf DW_TAG_inheritance
	.long	.Linfo_string154        # DW_AT_name
	.long	3098                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0xd0c:0xf DW_TAG_member
	.long	.Linfo_string155        # DW_AT_name
	.long	3446                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0xd1b:0xf DW_TAG_member
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	42                      # Abbrev [42] 0xd2a:0xe DW_TAG_member
	.long	.Linfo_string159        # DW_AT_name
	.long	3468                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	63                      # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	34                      # Abbrev [34] 0xd38:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin8           # DW_AT_low_pc
	.quad	.Lfunc_end8             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string175        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	64                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	35                      # Abbrev [35] 0xd52:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	3489                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0xd5a:0x1a DW_TAG_lexical_block
	.quad	.Ltmp88                 # DW_AT_low_pc
	.quad	.Ltmp89                 # DW_AT_high_pc
	.byte	36                      # Abbrev [36] 0xd6b:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	3489                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0xd76:0x5 DW_TAG_pointer_type
	.long	3451                    # DW_AT_type
	.byte	43                      # Abbrev [43] 0xd7b:0xa DW_TAG_class_type
	.long	.Linfo_string156        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	7                       # Abbrev [7] 0xd85:0x7 DW_TAG_base_type
	.long	.Linfo_string158        # DW_AT_name
	.byte	4                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	44                      # Abbrev [44] 0xd8c:0x15 DW_TAG_subroutine_type
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xd91:0x5 DW_TAG_formal_parameter
	.long	3489                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xd96:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xd9b:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0xda1:0x5 DW_TAG_pointer_type
	.long	3315                    # DW_AT_type
	.byte	37                      # Abbrev [37] 0xda6:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin9           # DW_AT_low_pc
	.quad	.Lfunc_end9             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string176        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0xdbe:0x79 DW_TAG_subprogram
	.quad	.Lfunc_begin10          # DW_AT_low_pc
	.quad	.Lfunc_end10            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string159        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	73                      # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	38                      # Abbrev [38] 0xddc:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc43           # DW_AT_location
	.long	3791                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	28                      # Abbrev [28] 0xde6:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc44           # DW_AT_location
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xdf3:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc45           # DW_AT_location
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0xe00:0x36 DW_TAG_lexical_block
	.quad	.Ltmp93                 # DW_AT_low_pc
	.quad	.Ltmp95                 # DW_AT_high_pc
	.byte	39                      # Abbrev [39] 0xe11:0xa DW_TAG_variable
	.long	.Ldebug_loc46           # DW_AT_location
	.long	3791                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	30                      # Abbrev [30] 0xe1b:0xd DW_TAG_variable
	.long	.Ldebug_loc47           # DW_AT_location
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xe28:0xd DW_TAG_variable
	.long	.Ldebug_loc48           # DW_AT_location
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	33                      # Abbrev [33] 0xe37:0x83 DW_TAG_class_type
	.long	.Linfo_string162        # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0xe41:0xf DW_TAG_inheritance
	.long	.Linfo_string154        # DW_AT_name
	.long	3098                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0xe50:0xf DW_TAG_member
	.long	.Linfo_string155        # DW_AT_name
	.long	3446                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0xe5f:0xf DW_TAG_member
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	42                      # Abbrev [42] 0xe6e:0xe DW_TAG_member
	.long	.Linfo_string159        # DW_AT_name
	.long	3770                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	73                      # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	34                      # Abbrev [34] 0xe7c:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin11          # DW_AT_low_pc
	.quad	.Lfunc_end11            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string177        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	74                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	35                      # Abbrev [35] 0xe96:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	3791                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0xe9e:0x1a DW_TAG_lexical_block
	.quad	.Ltmp96                 # DW_AT_low_pc
	.quad	.Ltmp97                 # DW_AT_high_pc
	.byte	36                      # Abbrev [36] 0xeaf:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	3791                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	44                      # Abbrev [44] 0xeba:0x15 DW_TAG_subroutine_type
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xebf:0x5 DW_TAG_formal_parameter
	.long	3791                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xec4:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xec9:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0xecf:0x5 DW_TAG_pointer_type
	.long	3639                    # DW_AT_type
	.byte	37                      # Abbrev [37] 0xed4:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin12          # DW_AT_low_pc
	.quad	.Lfunc_end12            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string178        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0xeec:0x79 DW_TAG_subprogram
	.quad	.Lfunc_begin13          # DW_AT_low_pc
	.quad	.Lfunc_end13            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string159        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	83                      # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	38                      # Abbrev [38] 0xf0a:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc49           # DW_AT_location
	.long	4093                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	28                      # Abbrev [28] 0xf14:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc50           # DW_AT_location
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xf21:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc51           # DW_AT_location
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0xf2e:0x36 DW_TAG_lexical_block
	.quad	.Ltmp100                # DW_AT_low_pc
	.quad	.Ltmp102                # DW_AT_high_pc
	.byte	39                      # Abbrev [39] 0xf3f:0xa DW_TAG_variable
	.long	.Ldebug_loc52           # DW_AT_location
	.long	4093                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	30                      # Abbrev [30] 0xf49:0xd DW_TAG_variable
	.long	.Ldebug_loc53           # DW_AT_location
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0xf56:0xd DW_TAG_variable
	.long	.Ldebug_loc54           # DW_AT_location
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	33                      # Abbrev [33] 0xf65:0x83 DW_TAG_class_type
	.long	.Linfo_string164        # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0xf6f:0xf DW_TAG_inheritance
	.long	.Linfo_string154        # DW_AT_name
	.long	3098                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0xf7e:0xf DW_TAG_member
	.long	.Linfo_string155        # DW_AT_name
	.long	3446                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0xf8d:0xf DW_TAG_member
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	42                      # Abbrev [42] 0xf9c:0xe DW_TAG_member
	.long	.Linfo_string159        # DW_AT_name
	.long	4072                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	83                      # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	34                      # Abbrev [34] 0xfaa:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin14          # DW_AT_low_pc
	.quad	.Lfunc_end14            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string179        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	84                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	35                      # Abbrev [35] 0xfc4:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4093                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0xfcc:0x1a DW_TAG_lexical_block
	.quad	.Ltmp103                # DW_AT_low_pc
	.quad	.Ltmp104                # DW_AT_high_pc
	.byte	36                      # Abbrev [36] 0xfdd:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4093                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	44                      # Abbrev [44] 0xfe8:0x15 DW_TAG_subroutine_type
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xfed:0x5 DW_TAG_formal_parameter
	.long	4093                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xff2:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0xff7:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0xffd:0x5 DW_TAG_pointer_type
	.long	3941                    # DW_AT_type
	.byte	37                      # Abbrev [37] 0x1002:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin15          # DW_AT_low_pc
	.quad	.Lfunc_end15            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string180        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0x101a:0x5f DW_TAG_subprogram
	.quad	.Lfunc_begin16          # DW_AT_low_pc
	.quad	.Lfunc_end16            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string159        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	95                      # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	38                      # Abbrev [38] 0x1038:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc55           # DW_AT_location
	.long	4379                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	28                      # Abbrev [28] 0x1042:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc56           # DW_AT_location
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0x104f:0x29 DW_TAG_lexical_block
	.quad	.Ltmp108                # DW_AT_low_pc
	.quad	.Ltmp110                # DW_AT_high_pc
	.byte	39                      # Abbrev [39] 0x1060:0xa DW_TAG_variable
	.long	.Ldebug_loc57           # DW_AT_location
	.long	4379                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	30                      # Abbrev [30] 0x106a:0xd DW_TAG_variable
	.long	.Ldebug_loc58           # DW_AT_location
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	33                      # Abbrev [33] 0x1079:0x92 DW_TAG_class_type
	.long	.Linfo_string165        # DW_AT_name
	.byte	32                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0x1083:0xf DW_TAG_inheritance
	.long	.Linfo_string153        # DW_AT_name
	.long	3002                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0x1092:0xf DW_TAG_member
	.long	.Linfo_string155        # DW_AT_name
	.long	3446                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0x10a1:0xf DW_TAG_member
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	41                      # Abbrev [41] 0x10b0:0xf DW_TAG_member
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	42                      # Abbrev [42] 0x10bf:0xe DW_TAG_member
	.long	.Linfo_string159        # DW_AT_name
	.long	4363                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	95                      # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	34                      # Abbrev [34] 0x10cd:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin17          # DW_AT_low_pc
	.quad	.Lfunc_end17            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string181        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	96                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	35                      # Abbrev [35] 0x10e7:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4379                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0x10ef:0x1a DW_TAG_lexical_block
	.quad	.Ltmp111                # DW_AT_low_pc
	.quad	.Ltmp112                # DW_AT_high_pc
	.byte	36                      # Abbrev [36] 0x1100:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4379                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	44                      # Abbrev [44] 0x110b:0x10 DW_TAG_subroutine_type
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0x1110:0x5 DW_TAG_formal_parameter
	.long	4379                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0x1115:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x111b:0x5 DW_TAG_pointer_type
	.long	4217                    # DW_AT_type
	.byte	37                      # Abbrev [37] 0x1120:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin18          # DW_AT_low_pc
	.quad	.Lfunc_end18            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string182        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0x1138:0x5f DW_TAG_subprogram
	.quad	.Lfunc_begin19          # DW_AT_low_pc
	.quad	.Lfunc_end19            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string159        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	105                     # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	38                      # Abbrev [38] 0x1156:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc59           # DW_AT_location
	.long	4665                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	28                      # Abbrev [28] 0x1160:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc60           # DW_AT_location
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0x116d:0x29 DW_TAG_lexical_block
	.quad	.Ltmp116                # DW_AT_low_pc
	.quad	.Ltmp118                # DW_AT_high_pc
	.byte	39                      # Abbrev [39] 0x117e:0xa DW_TAG_variable
	.long	.Ldebug_loc61           # DW_AT_location
	.long	4665                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	30                      # Abbrev [30] 0x1188:0xd DW_TAG_variable
	.long	.Ldebug_loc62           # DW_AT_location
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	33                      # Abbrev [33] 0x1197:0x92 DW_TAG_class_type
	.long	.Linfo_string166        # DW_AT_name
	.byte	32                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0x11a1:0xf DW_TAG_inheritance
	.long	.Linfo_string153        # DW_AT_name
	.long	3002                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0x11b0:0xf DW_TAG_member
	.long	.Linfo_string155        # DW_AT_name
	.long	3446                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0x11bf:0xf DW_TAG_member
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	41                      # Abbrev [41] 0x11ce:0xf DW_TAG_member
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	42                      # Abbrev [42] 0x11dd:0xe DW_TAG_member
	.long	.Linfo_string159        # DW_AT_name
	.long	4649                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	105                     # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	34                      # Abbrev [34] 0x11eb:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin20          # DW_AT_low_pc
	.quad	.Lfunc_end20            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string183        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	106                     # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	35                      # Abbrev [35] 0x1205:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4665                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0x120d:0x1a DW_TAG_lexical_block
	.quad	.Ltmp119                # DW_AT_low_pc
	.quad	.Ltmp120                # DW_AT_high_pc
	.byte	36                      # Abbrev [36] 0x121e:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4665                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	44                      # Abbrev [44] 0x1229:0x10 DW_TAG_subroutine_type
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0x122e:0x5 DW_TAG_formal_parameter
	.long	4665                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0x1233:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x1239:0x5 DW_TAG_pointer_type
	.long	4503                    # DW_AT_type
	.byte	37                      # Abbrev [37] 0x123e:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin21          # DW_AT_low_pc
	.quad	.Lfunc_end21            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string184        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0x1256:0x5f DW_TAG_subprogram
	.quad	.Lfunc_begin22          # DW_AT_low_pc
	.quad	.Lfunc_end22            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string159        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	115                     # DW_AT_decl_line
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	38                      # Abbrev [38] 0x1274:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc63           # DW_AT_location
	.long	4951                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	28                      # Abbrev [28] 0x127e:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc64           # DW_AT_location
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0x128b:0x29 DW_TAG_lexical_block
	.quad	.Ltmp123                # DW_AT_low_pc
	.quad	.Ltmp125                # DW_AT_high_pc
	.byte	39                      # Abbrev [39] 0x129c:0xa DW_TAG_variable
	.long	.Ldebug_loc65           # DW_AT_location
	.long	4951                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	30                      # Abbrev [30] 0x12a6:0xd DW_TAG_variable
	.long	.Ldebug_loc66           # DW_AT_location
	.long	.Linfo_string157        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	33                      # Abbrev [33] 0x12b5:0x92 DW_TAG_class_type
	.long	.Linfo_string167        # DW_AT_name
	.byte	32                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0x12bf:0xf DW_TAG_inheritance
	.long	.Linfo_string153        # DW_AT_name
	.long	3002                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0x12ce:0xf DW_TAG_member
	.long	.Linfo_string155        # DW_AT_name
	.long	3446                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0x12dd:0xf DW_TAG_member
	.long	.Linfo_string161        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	41                      # Abbrev [41] 0x12ec:0xf DW_TAG_member
	.long	.Linfo_string163        # DW_AT_name
	.long	3461                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	42                      # Abbrev [42] 0x12fb:0xe DW_TAG_member
	.long	.Linfo_string159        # DW_AT_name
	.long	4935                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	115                     # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	34                      # Abbrev [34] 0x1309:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin23          # DW_AT_low_pc
	.quad	.Lfunc_end23            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string185        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	116                     # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	35                      # Abbrev [35] 0x1323:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4951                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	29                      # Abbrev [29] 0x132b:0x1a DW_TAG_lexical_block
	.quad	.Ltmp126                # DW_AT_low_pc
	.quad	.Ltmp127                # DW_AT_high_pc
	.byte	36                      # Abbrev [36] 0x133c:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4951                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	44                      # Abbrev [44] 0x1347:0x10 DW_TAG_subroutine_type
	.long	3461                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0x134c:0x5 DW_TAG_formal_parameter
	.long	4951                    # DW_AT_type
	.byte	18                      # Abbrev [18] 0x1351:0x5 DW_TAG_formal_parameter
	.long	3461                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x1357:0x5 DW_TAG_pointer_type
	.long	4789                    # DW_AT_type
	.byte	37                      # Abbrev [37] 0x135c:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin24          # DW_AT_low_pc
	.quad	.Lfunc_end24            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string186        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	37                      # Abbrev [37] 0x1374:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin25          # DW_AT_low_pc
	.quad	.Lfunc_end25            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string187        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	45                      # Abbrev [45] 0x138c:0x5 DW_TAG_reference_type
	.long	5009                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x1391:0x15 DW_TAG_structure_type
	.long	.Linfo_string156        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	13                      # Abbrev [13] 0x1398:0xd DW_TAG_member
	.long	.Linfo_string9          # DW_AT_name
	.long	79                      # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	46                      # Abbrev [46] 0x13a6:0x1d DW_TAG_enumeration_type
	.long	.Linfo_string193        # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	4                       # DW_AT_alignment
	.byte	16                      # Abbrev [16] 0x13b0:0x6 DW_TAG_enumerator
	.long	.Linfo_string190        # DW_AT_name
	.byte	0                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x13b6:0x6 DW_TAG_enumerator
	.long	.Linfo_string191        # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	16                      # Abbrev [16] 0x13bc:0x6 DW_TAG_enumerator
	.long	.Linfo_string192        # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x13c3:0x5 DW_TAG_pointer_type
	.long	3461                    # DW_AT_type
	.byte	5                       # Abbrev [5] 0x13c8:0x5 DW_TAG_pointer_type
	.long	3002                    # DW_AT_type
	.byte	5                       # Abbrev [5] 0x13cd:0x5 DW_TAG_pointer_type
	.long	3098                    # DW_AT_type
	.byte	0                       # End Of Children Mark
.Ldebug_info_end0:
	.section	.debug_ranges,"",@progbits
.Ldebug_ranges0:
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges1:
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Ltmp16-.Lfunc_begin0
	.quad	.Ltmp17-.Lfunc_begin0
	.quad	.Ltmp19-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges2:
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
	.quad	.Ltmp23-.Lfunc_begin0
	.quad	.Ltmp25-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges3:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp10-.Lfunc_begin0
	.quad	.Ltmp11-.Lfunc_begin0
	.quad	.Ltmp13-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Ltmp16-.Lfunc_begin0
	.quad	.Ltmp17-.Lfunc_begin0
	.quad	.Ltmp19-.Lfunc_begin0
	.quad	.Ltmp20-.Lfunc_begin0
	.quad	.Ltmp21-.Lfunc_begin0
	.quad	.Ltmp23-.Lfunc_begin0
	.quad	.Ltmp25-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges4:
	.quad	.Ltmp37-.Lfunc_begin0
	.quad	.Ltmp38-.Lfunc_begin0
	.quad	.Ltmp39-.Lfunc_begin0
	.quad	.Ltmp43-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges5:
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp45-.Lfunc_begin0
	.quad	.Ltmp46-.Lfunc_begin0
	.quad	.Ltmp50-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges6:
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp52-.Lfunc_begin0
	.quad	.Ltmp53-.Lfunc_begin0
	.quad	.Ltmp57-.Lfunc_begin0
	.quad	0
	.quad	0
.Ldebug_ranges7:
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp38-.Lfunc_begin0
	.quad	.Ltmp39-.Lfunc_begin0
	.quad	.Ltmp43-.Lfunc_begin0
	.quad	.Ltmp44-.Lfunc_begin0
	.quad	.Ltmp45-.Lfunc_begin0
	.quad	.Ltmp46-.Lfunc_begin0
	.quad	.Ltmp50-.Lfunc_begin0
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp52-.Lfunc_begin0
	.quad	.Ltmp53-.Lfunc_begin0
	.quad	.Ltmp57-.Lfunc_begin0
	.quad	.Ltmp59-.Lfunc_begin0
	.quad	.Ltmp61-.Lfunc_begin0
	.quad	.Ltmp66-.Lfunc_begin0
	.quad	.Ltmp67-.Lfunc_begin0
	.quad	0
	.quad	0
	.section	.debug_pubnames,"",@progbits
	.long	.LpubNames_end0-.LpubNames_begin0 # Length of Public Names Info
.LpubNames_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	5075                    # Compilation Unit Length
	.long	562                     # DIE offset
	.asciz	"_ZSt12_S_showpoint"    # External Name
	.long	2190                    # DIE offset
	.asciz	"_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1ee5vmesh15INVALID_LOCALIDE" # External Name
	.long	4694                    # DIE offset
	.asciz	"call"                  # External Name
	.long	1697                    # DIE offset
	.asciz	"_ZN59_INTERNAL_37_backgroundfield_integratefunction_cpp_f19cb1eeSt8__ioinitE" # External Name
	.long	1782                    # DIE offset
	.asciz	"_ZTI11T1DFunction"     # External Name
	.long	3074                    # DIE offset
	.asciz	"_ZN11T1DFunctionD0Ev"  # External Name
	.long	133                     # DIE offset
	.asciz	"std"                   # External Name
	.long	576                     # DIE offset
	.asciz	"_ZSt9_S_skipws"        # External Name
	.long	522                     # DIE offset
	.asciz	"_ZSt11_S_internal"     # External Name
	.long	662                     # DIE offset
	.asciz	"_ZSt9_S_badbit"        # External Name
	.long	1554                    # DIE offset
	.asciz	"_ZTV9T3D_fix23"        # External Name
	.long	504                     # DIE offset
	.asciz	"_ZSt6_S_dec"           # External Name
	.long	2036                    # DIE offset
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # External Name
	.long	510                     # DIE offset
	.asciz	"_ZSt8_S_fixed"         # External Name
	.long	2046                    # DIE offset
	.asciz	"_ZTS8T3D_fix1"         # External Name
	.long	2078                    # DIE offset
	.asciz	"_ZTS8T3D_fix2"         # External Name
	.long	2098                    # DIE offset
	.asciz	"_ZTS8T3D_fix3"         # External Name
	.long	4301                    # DIE offset
	.asciz	"T3D_fix12::~T3D_fix12" # External Name
	.long	3708                    # DIE offset
	.asciz	"T3D_fix2::~T3D_fix2"   # External Name
	.long	4587                    # DIE offset
	.asciz	"T3D_fix13::~T3D_fix13" # External Name
	.long	1902                    # DIE offset
	.asciz	"_ZTI9T3D_fix12"        # External Name
	.long	1922                    # DIE offset
	.asciz	"_ZTI9T3D_fix13"        # External Name
	.long	569                     # DIE offset
	.asciz	"_ZSt10_S_showpos"      # External Name
	.long	809                     # DIE offset
	.asciz	"_ZNSt8ios_base11imbue_eventE" # External Name
	.long	2508                    # DIE offset
	.asciz	"surfaceAverage"        # External Name
	.long	1812                    # DIE offset
	.asciz	"_ZTI11T2DFunction"     # External Name
	.long	1962                    # DIE offset
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # External Name
	.long	1674                    # DIE offset
	.asciz	"__I___37_backgroundfield_integratefunction_cpp_f19cb1ee" # External Name
	.long	3170                    # DIE offset
	.asciz	"_ZN11T2DFunctionD0Ev"  # External Name
	.long	516                     # DIE offset
	.asciz	"_ZSt6_S_hex"           # External Name
	.long	534                     # DIE offset
	.asciz	"_ZSt6_S_oct"           # External Name
	.long	2118                    # DIE offset
	.asciz	"_ZTS9T3D_fix12"        # External Name
	.long	555                     # DIE offset
	.asciz	"_ZSt11_S_showbase"     # External Name
	.long	2150                    # DIE offset
	.asciz	"_ZTS9T3D_fix13"        # External Name
	.long	528                     # DIE offset
	.asciz	"_ZSt7_S_left"          # External Name
	.long	47                      # DIE offset
	.asciz	"_ZTV11T1DFunction"     # External Name
	.long	3384                    # DIE offset
	.asciz	"T3D_fix1::~T3D_fix1"   # External Name
	.long	638                     # DIE offset
	.asciz	"_ZSt19_S_ios_fmtflags_min" # External Name
	.long	628                     # DIE offset
	.asciz	"_ZSt19_S_ios_fmtflags_max" # External Name
	.long	698                     # DIE offset
	.asciz	"_ZSt18_S_ios_iostate_min" # External Name
	.long	688                     # DIE offset
	.asciz	"_ZSt18_S_ios_iostate_max" # External Name
	.long	668                     # DIE offset
	.asciz	"_ZSt9_S_eofbit"        # External Name
	.long	803                     # DIE offset
	.asciz	"_ZNSt8ios_base11erase_eventE" # External Name
	.long	1755                    # DIE offset
	.asciz	"WID"                   # External Name
	.long	1942                    # DIE offset
	.asciz	"_ZTI9T3D_fix23"        # External Name
	.long	541                     # DIE offset
	.asciz	"_ZSt8_S_right"         # External Name
	.long	2206                    # DIE offset
	.asciz	"lineAverage"           # External Name
	.long	4670                    # DIE offset
	.asciz	"_ZN9T3D_fix13D0Ev"     # External Name
	.long	4384                    # DIE offset
	.asciz	"_ZN9T3D_fix12D0Ev"     # External Name
	.long	613                     # DIE offset
	.asciz	"_ZSt13_S_floatfield"   # External Name
	.long	620                     # DIE offset
	.asciz	"_ZSt19_S_ios_fmtflags_end" # External Name
	.long	1832                    # DIE offset
	.asciz	"_ZTI8T3D_fix1"         # External Name
	.long	1745                    # DIE offset
	.asciz	"__dso_handle"          # External Name
	.long	606                     # DIE offset
	.asciz	"_ZSt12_S_basefield"    # External Name
	.long	680                     # DIE offset
	.asciz	"_ZSt18_S_ios_iostate_end" # External Name
	.long	1594                    # DIE offset
	.asciz	"_ZTV11T2DFunction"     # External Name
	.long	548                     # DIE offset
	.asciz	"_ZSt13_S_scientific"   # External Name
	.long	1862                    # DIE offset
	.asciz	"_ZTI8T3D_fix2"         # External Name
	.long	1882                    # DIE offset
	.asciz	"_ZTI8T3D_fix3"         # External Name
	.long	2170                    # DIE offset
	.asciz	"_ZTS9T3D_fix23"        # External Name
	.long	3012                    # DIE offset
	.asciz	"T1DFunction::~T1DFunction" # External Name
	.long	1984                    # DIE offset
	.asciz	"_ZTS11T1DFunction"     # External Name
	.long	3108                    # DIE offset
	.asciz	"T2DFunction::~T2DFunction" # External Name
	.long	2805                    # DIE offset
	.asciz	"volumeAverage"         # External Name
	.long	4873                    # DIE offset
	.asciz	"T3D_fix23::~T3D_fix23" # External Name
	.long	138                     # DIE offset
	.asciz	"std::_ZSt4cerr"        # External Name
	.long	3796                    # DIE offset
	.asciz	"_ZN8T3D_fix2D0Ev"      # External Name
	.long	4098                    # DIE offset
	.asciz	"_ZN8T3D_fix3D0Ev"      # External Name
	.long	5040                    # DIE offset
	.asciz	"X"                     # External Name
	.long	5046                    # DIE offset
	.asciz	"Y"                     # External Name
	.long	5052                    # DIE offset
	.asciz	"Z"                     # External Name
	.long	599                     # DIE offset
	.asciz	"_ZSt14_S_adjustfield"  # External Name
	.long	498                     # DIE offset
	.asciz	"_ZSt12_S_boolalpha"    # External Name
	.long	4956                    # DIE offset
	.asciz	"_ZN9T3D_fix23D0Ev"     # External Name
	.long	3494                    # DIE offset
	.asciz	"_ZN8T3D_fix1D0Ev"      # External Name
	.long	674                     # DIE offset
	.asciz	"_ZSt10_S_failbit"      # External Name
	.long	583                     # DIE offset
	.asciz	"_ZSt10_S_unitbuf"      # External Name
	.long	113                     # DIE offset
	.asciz	"_ZTV9T3D_fix12"        # External Name
	.long	1574                    # DIE offset
	.asciz	"_ZTV9T3D_fix13"        # External Name
	.long	2016                    # DIE offset
	.asciz	"_ZTS11T2DFunction"     # External Name
	.long	4980                    # DIE offset
	.asciz	"__sti___37_backgroundfield_integratefunction_cpp_f19cb1ee" # External Name
	.long	1764                    # DIE offset
	.asciz	"WID2"                  # External Name
	.long	815                     # DIE offset
	.asciz	"_ZNSt8ios_base13copyfmt_eventE" # External Name
	.long	1773                    # DIE offset
	.asciz	"WID3"                  # External Name
	.long	656                     # DIE offset
	.asciz	"_ZSt10_S_goodbit"      # External Name
	.long	4010                    # DIE offset
	.asciz	"T3D_fix3::~T3D_fix3"   # External Name
	.long	591                     # DIE offset
	.asciz	"_ZSt12_S_uppercase"    # External Name
	.long	1634                    # DIE offset
	.asciz	"_ZTV8T3D_fix1"         # External Name
	.long	1654                    # DIE offset
	.asciz	"_ZTV8T3D_fix2"         # External Name
	.long	1614                    # DIE offset
	.asciz	"_ZTV8T3D_fix3"         # External Name
	.long	0                       # End Mark
.LpubNames_end0:
	.section	.debug_pubtypes,"",@progbits
	.long	.LpubTypes_end0-.LpubTypes_begin0 # Length of Public Types Info
.LpubTypes_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	5075                    # Compilation Unit Length
	.long	1802                    # DIE offset
	.asciz	"__class_type_info"     # External Name
	.long	491                     # DIE offset
	.asciz	"_ZSt13_Ios_Fmtflags"   # External Name
	.long	649                     # DIE offset
	.asciz	"_ZSt12_Ios_Iostate"    # External Name
	.long	888                     # DIE offset
	.asciz	"_ZSt6locale"           # External Name
	.long	796                     # DIE offset
	.asciz	"_ZNSt8ios_base5eventE" # External Name
	.long	99                      # DIE offset
	.asciz	"int"                   # External Name
	.long	1182                    # DIE offset
	.asciz	"_ZSt5ctypeIcE"         # External Name
	.long	1031                    # DIE offset
	.asciz	"unsigned long"         # External Name
	.long	3002                    # DIE offset
	.asciz	"T1DFunction"           # External Name
	.long	997                     # DIE offset
	.asciz	"_ZNSt6locale5facetE"   # External Name
	.long	866                     # DIE offset
	.asciz	"void"                  # External Name
	.long	5009                    # DIE offset
	.asciz	"T3DFunction"           # External Name
	.long	4217                    # DIE offset
	.asciz	"T3D_fix12"             # External Name
	.long	1323                    # DIE offset
	.asciz	"__SO__NSt6locale5facetE" # External Name
	.long	1048                    # DIE offset
	.asciz	"signed char"           # External Name
	.long	1716                    # DIE offset
	.asciz	"_ZNSt8ios_base4InitE"  # External Name
	.long	148                     # DIE offset
	.asciz	"std::ostream"          # External Name
	.long	3461                    # DIE offset
	.asciz	"double"                # External Name
	.long	714                     # DIE offset
	.asciz	"_ZNSt8ios_base14_Callback_listE" # External Name
	.long	317                     # DIE offset
	.asciz	"_ZSt8ios_base"         # External Name
	.long	2199                    # DIE offset
	.asciz	"unsigned"              # External Name
	.long	914                     # DIE offset
	.asciz	"_ZNSt6locale5_ImplE"   # External Name
	.long	1507                    # DIE offset
	.asciz	"_ZSt7num_putIcSt19ostreambuf_iteratorIcSt11char_traitsIcEEE" # External Name
	.long	1533                    # DIE offset
	.asciz	"_ZSt7num_getIcSt19istreambuf_iteratorIcSt11char_traitsIcEEE" # External Name
	.long	3315                    # DIE offset
	.asciz	"T3D_fix1"              # External Name
	.long	197                     # DIE offset
	.asciz	"_ZSt9basic_iosIcSt11char_traitsIcEE" # External Name
	.long	3941                    # DIE offset
	.asciz	"T3D_fix3"              # External Name
	.long	3639                    # DIE offset
	.asciz	"T3D_fix2"              # External Name
	.long	1453                    # DIE offset
	.asciz	"__locale_data"         # External Name
	.long	5030                    # DIE offset
	.asciz	"coordinate"            # External Name
	.long	4503                    # DIE offset
	.asciz	"T3D_fix13"             # External Name
	.long	1465                    # DIE offset
	.asciz	"unsigned short"        # External Name
	.long	1065                    # DIE offset
	.asciz	"_ZSt15basic_streambufIcSt11char_traitsIcEE" # External Name
	.long	484                     # DIE offset
	.asciz	"long"                  # External Name
	.long	3098                    # DIE offset
	.asciz	"T2DFunction"           # External Name
	.long	1362                    # DIE offset
	.asciz	"__locale_struct"       # External Name
	.long	827                     # DIE offset
	.asciz	"_ZNSt8ios_base6_WordsE" # External Name
	.long	4789                    # DIE offset
	.asciz	"T3D_fix23"             # External Name
	.long	162                     # DIE offset
	.asciz	"_ZSo"                  # External Name
	.long	1852                    # DIE offset
	.asciz	"__si_class_type_info"  # External Name
	.long	0                       # End Mark
.LpubTypes_end0:
	.section	".note.GNU-stack","",@progbits
	.section	.debug_line,"",@progbits
.Lline_table_start0:
