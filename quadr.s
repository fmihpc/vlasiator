	.text
	.file	"quadr.ll"
	.file	1 "/home/talgat/vlasiator/backgroundfield/quadr.cpp"
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _Z14Romberg_simpleRK11T1DFunctionddd
.LCPI0_0:
	.quad	4602678819172646912     # double 0.5
.LCPI0_3:
	.quad	4598175219545276416     # double 0.25
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4
.LCPI0_1:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
.LCPI0_2:
	.quad	9223372036854775807     # double NaN
	.quad	9223372036854775807     # double NaN
	.text
	.globl	_Z14Romberg_simpleRK11T1DFunctionddd
	.p2align	4, 0x90
	.type	_Z14Romberg_simpleRK11T1DFunctionddd,@function
_Z14Romberg_simpleRK11T1DFunctionddd:   # @_Z14Romberg_simpleRK11T1DFunctionddd
.Lfunc_begin0:
	.loc	1 156 0                 # backgroundfield/quadr.cpp:156:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: Romberg_simple:func <- $rdi
	#DEBUG_VALUE: Romberg_simple:a <- $xmm0
	#DEBUG_VALUE: Romberg_simple:b <- $xmm1
	#DEBUG_VALUE: Romberg_simple:absacc <- $xmm2
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$600, %rsp              # imm = 0x258
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	vmovsd	%xmm2, -160(%rbp)       # 8-byte Spill
	movq	%rdi, %r13
.Ltmp0:
	#DEBUG_VALUE: absacc <- undef
	#DEBUG_VALUE: func <- undef
	movabsq	$4607182418800017408, %rax # imm = 0x3FF0000000000000
.Ltmp1:
	#DEBUG_VALUE: it <- 0
	#DEBUG_VALUE: b <- $xmm1
	#DEBUG_VALUE: a <- $xmm0
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- $xmm1
	#DEBUG_VALUE: Romberg_simple:a <- $xmm0
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	.loc	1 162 1 prologue_end    # backgroundfield/quadr.cpp:162:1
	movq	%rax, -240(%rbp)
.Ltmp2:
	#DEBUG_VALUE: j <- 1
	#DEBUG_VALUE: dresult <- 0.000000e+00
	.loc	1 164 1                 # backgroundfield/quadr.cpp:164:1
	leaq	-304(%rbp), %rcx
	leaq	-232(%rbp), %rax
	movq	%rax, -96(%rbp)         # 8-byte Spill
	vmovsd	%xmm1, -144(%rbp)       # 8-byte Spill
.Ltmp3:
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vmovsd	%xmm0, -88(%rbp)        # 8-byte Spill
.Ltmp4:
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	vsubsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -80(%rbp)        # 8-byte Spill
.Ltmp5:
	vmulsd	.LCPI0_0(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, -136(%rbp)       # 8-byte Spill
	xorl	%ebx, %ebx
	vxorpd	%xmm0, %xmm0, %xmm0
	vmovapd	%xmm0, -128(%rbp)       # 16-byte Spill
	movl	$1, %r14d
	xorl	%r15d, %r15d
.Ltmp6:
	.p2align	4, 0x90
.LBB0_1:                                # %L.B0033
                                        # =>This Loop Header: Depth=1
                                        #     Child Loop BB0_6 Depth 2
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r15d
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: j <- $r14
	.loc	1 165 1 is_stmt 1       # backgroundfield/quadr.cpp:165:1
	leaq	1(%rbx), %rax
	cmpl	$3, %eax
	movq	%rax, -168(%rbp)        # 8-byte Spill
	movl	%eax, %r12d
	movl	$3, %eax
	cmovgeq	%rax, %r12
	cmpq	$1, %r14
	movq	%rcx, -56(%rbp)         # 8-byte Spill
	jne	.LBB0_3
.Ltmp7:
# %bb.2:                                # %L..inline.10302.thread
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r15d
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 49 1                  # backgroundfield/quadr.cpp:49:1
	movq	(%r13), %rax
	movq	%r13, %rdi
	vmovsd	-144(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	callq	*(%rax)
	vmovsd	%xmm0, -48(%rbp)        # 8-byte Spill
	movq	(%r13), %rax
	movq	%r13, %rdi
	vmovsd	-88(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	callq	*(%rax)
	vaddsd	-48(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmulsd	-136(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmovsd	%xmm0, -320(%rbp,%r14,8)
	movl	$1, -48(%rbp)           # 4-byte Folded Spill
.Ltmp8:
	#DEBUG_VALUE: it <- 1
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	jmp	.LBB0_8
.Ltmp9:
	.p2align	4, 0x90
.LBB0_3:                                # %L..inline.10296
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r15d
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 52 1 is_stmt 1        # backgroundfield/quadr.cpp:52:1
	vcvtsi2sd	%r15d, %xmm6, %xmm1
	.loc	1 56 1                  # backgroundfield/quadr.cpp:56:1
	testl	%r15d, %r15d
	movq	%rbx, -112(%rbp)        # 8-byte Spill
	vmovsd	%xmm1, -152(%rbp)       # 8-byte Spill
	jle	.LBB0_4
.Ltmp10:
# %bb.5:                                # %L.B0204
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r15d
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vmovsd	-80(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vdivsd	%xmm1, %xmm0, %xmm1
	vmovsd	-88(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	%xmm1, -64(%rbp)        # 8-byte Spill
	.loc	1 53 1 is_stmt 1        # backgroundfield/quadr.cpp:53:1
	vfmadd231sd	.LCPI0_0(%rip), %xmm1, %xmm0 # xmm0 = (xmm1 * mem) + xmm0
	vxorpd	%xmm1, %xmm1, %xmm1
	.loc	1 56 1                  # backgroundfield/quadr.cpp:56:1
	movl	%r15d, %ebx
.Ltmp11:
	.p2align	4, 0x90
.LBB0_6:                                # %L..inline.10310
                                        #   Parent Loop BB0_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r15d
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vmovsd	%xmm0, -72(%rbp)        # 8-byte Spill
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	.loc	1 57 1 is_stmt 1        # backgroundfield/quadr.cpp:57:1
	movq	(%r13), %rax
	movq	%r13, %rdi
	callq	*(%rax)
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	%xmm0, %xmm1, %xmm1
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	vmovsd	-72(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	-64(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	.loc	1 59 1                  # backgroundfield/quadr.cpp:59:1
	decl	%ebx
	jg	.LBB0_6
	jmp	.LBB0_7
.Ltmp12:
.LBB0_4:                                #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r15d
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vxorpd	%xmm1, %xmm1, %xmm1
.Ltmp13:
.LBB0_7:                                # %L..inline.10302
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r15d
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 59 1                  # backgroundfield/quadr.cpp:59:1
	vmulsd	-80(%rbp), %xmm1, %xmm0 # 8-byte Folded Reload
	vdivsd	-152(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vaddsd	-320(%rbp,%r14,8), %xmm0, %xmm0
	vmulsd	.LCPI0_0(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, -320(%rbp,%r14,8)
	.loc	1 60 1 is_stmt 1        # backgroundfield/quadr.cpp:60:1
	addl	%r15d, %r15d
.Ltmp14:
	#DEBUG_VALUE: it <- $r15d
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movl	%r15d, -48(%rbp)        # 4-byte Spill
.Ltmp15:
	#DEBUG_VALUE: it <- [DW_OP_constu 48, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- [DW_OP_constu 48, DW_OP_minus] [$rbp+0]
	movq	-112(%rbp), %rbx        # 8-byte Reload
.Ltmp16:
.LBB0_8:                                # %L.B0205
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 168 1 is_stmt 1       # backgroundfield/quadr.cpp:168:1
	cmpq	$3, %r14
	movl	$3, %r15d
	cmovbq	%r14, %r15
	.loc	1 169 1                 # backgroundfield/quadr.cpp:169:1
	movl	%r14d, %eax
	subl	%r15d, %eax
	cltq
	leaq	-312(,%rax,8), %rcx
	addq	%rbp, %rcx
.Ltmp17:
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	testl	%r15d, %r15d
	je	.LBB0_9
.Ltmp18:
# %bb.10:                               # %L.B0206
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%rcx, -72(%rbp)         # 8-byte Spill
	leaq	-304(,%rbx,8), %rbx
	addq	%rbp, %rbx
	shlq	$3, %r12
	subq	%r12, %rbx
	movq	%rax, -64(%rbp)         # 8-byte Spill
.Ltmp19:
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	vmovsd	-240(%rbp,%rax,8), %xmm0 # xmm0 = mem[0],zero
	vmovapd	%xmm0, -112(%rbp)       # 16-byte Spill
	leaq	-632(%rbp), %rdi
	movq	%rbx, %rsi
	movq	%r12, %rdx
	callq	memcpy
	leaq	-472(%rbp), %rdi
	movq	%rbx, %rsi
	movq	%r12, %rdx
	callq	memcpy
	xorl	%eax, %eax
	cmpq	$2, %r15
	vmovapd	.LCPI0_2(%rip), %xmm5   # xmm5 = [NaN,NaN]
	jb	.LBB0_12
.Ltmp20:
# %bb.11:                               # %L.B0197.L.B0197_crit_edge.lr.ph
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovapd	-112(%rbp), %xmm0       # 16-byte Reload
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	vxorpd	.LCPI0_1(%rip), %xmm0, %xmm0
	vmovddup	%xmm0, %xmm0    # xmm0 = xmm0[0,0]
	vandpd	%xmm5, %xmm0, %xmm0
	vpermilpd	$1, %xmm0, %xmm1 # xmm1 = xmm0[1,0]
	vminsd	%xmm1, %xmm0, %xmm0
	movq	-64(%rbp), %rcx         # 8-byte Reload
	vmovsd	-232(%rbp,%rcx,8), %xmm1 # xmm1 = mem[0],zero
	vandpd	%xmm5, %xmm1, %xmm1
	xorl	%eax, %eax
	vucomisd	%xmm1, %xmm0
	seta	%al
	cmpq	$2, %r15
	jbe	.LBB0_12
.Ltmp21:
# %bb.19:                               # %L.B0197.L.B0197_crit_edge.1
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	-224(%rbp,%rcx,8), %xmm1 # xmm1 = mem[0],zero
	vandpd	%xmm5, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$2, %ecx
	cmoval	%ecx, %eax
.Ltmp22:
.LBB0_12:                               # %L.B0198
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 83 1 is_stmt 1        # backgroundfield/quadr.cpp:83:1
	movl	%eax, %ecx
	movq	-72(%rbp), %rdx         # 8-byte Reload
	vmovsd	(%rdx,%rcx,8), %xmm0    # xmm0 = mem[0],zero
.Ltmp23:
	#DEBUG_VALUE: result <- $xmm0
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	cmpl	$2, %r15d
	jb	.LBB0_15
.Ltmp24:
# %bb.13:                               # %L.B0208
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	vmovsd	-464(%rbp), %xmm1       # xmm1 = mem[0],zero
	vsubsd	-632(%rbp), %xmm1, %xmm1
	#DEBUG_VALUE: result <- $xmm0
	movq	-64(%rbp), %rbx         # 8-byte Reload
	vmovsd	-240(%rbp,%rbx,8), %xmm2 # xmm2 = mem[0],zero
	vmovsd	-232(%rbp,%rbx,8), %xmm3 # xmm3 = mem[0],zero
	vsubsd	%xmm3, %xmm2, %xmm4
	vdivsd	%xmm4, %xmm1, %xmm1
	vmulsd	%xmm1, %xmm3, %xmm3
	vmovsd	%xmm3, -632(%rbp)
.Ltmp25:
	.loc	1 86 1                  # backgroundfield/quadr.cpp:86:1
	decl	%eax
.Ltmp26:
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	vmulsd	%xmm1, %xmm2, %xmm1
	vmovsd	%xmm1, -472(%rbp)
.Ltmp27:
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	leaq	-1(%r15), %r8
.Ltmp28:
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	cmpq	$2, %r15
	jle	.LBB0_14
.Ltmp29:
# %bb.20:                               # %L.B0201.1
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	vmovsd	-456(%rbp), %xmm1       # xmm1 = mem[0],zero
	vsubsd	-624(%rbp), %xmm1, %xmm1
	vmovsd	-232(%rbp,%rbx,8), %xmm2 # xmm2 = mem[0],zero
	vmovsd	-224(%rbp,%rbx,8), %xmm3 # xmm3 = mem[0],zero
	vsubsd	%xmm3, %xmm2, %xmm4
	vdivsd	%xmm4, %xmm1, %xmm1
	vmulsd	%xmm1, %xmm3, %xmm3
	vmovsd	%xmm3, -624(%rbp)
	vmulsd	%xmm1, %xmm2, %xmm1
	vmovsd	%xmm1, -464(%rbp)
.Ltmp30:
.LBB0_14:                               # %L.B0202
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	leal	(%rax,%rax), %esi
	addl	$2, %esi
	leal	-1(%r15), %edi
	xorl	%edx, %edx
	cmpl	%edi, %esi
	setge	%r9b
	movslq	%eax, %rsi
	leaq	-632(%rbp,%rsi,8), %rdi
	leaq	-464(%rbp), %rcx
	leaq	(%rcx,%rsi,8), %rsi
	cmovgeq	%rdi, %rsi
	vmovsd	(%rsi), %xmm1           # xmm1 = mem[0],zero
.Ltmp31:
	#DEBUG_VALUE: dresult <- $xmm1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovapd	%xmm1, -128(%rbp)       # 16-byte Spill
.Ltmp32:
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	vaddsd	%xmm1, %xmm0, %xmm0
.Ltmp33:
	#DEBUG_VALUE: result <- $xmm0
	cmpq	$1, %r8
	jbe	.LBB0_15
.Ltmp34:
# %bb.21:                               # %L.B0199.1
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	movb	%r9b, %dl
	subl	%edx, %eax
	#DEBUG_VALUE: result <- $xmm0
.Ltmp35:
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	cmpl	$3, %r15d
	jne	.LBB0_23
.Ltmp36:
# %bb.22:                               # %L.B0201.159
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	vmovsd	-464(%rbp), %xmm1       # xmm1 = mem[0],zero
	vsubsd	-632(%rbp), %xmm1, %xmm1
	vmovsd	-240(%rbp,%rbx,8), %xmm2 # xmm2 = mem[0],zero
	vmovsd	-224(%rbp,%rbx,8), %xmm3 # xmm3 = mem[0],zero
	vsubsd	%xmm3, %xmm2, %xmm4
	vdivsd	%xmm4, %xmm1, %xmm1
	vmulsd	%xmm1, %xmm3, %xmm3
	vmovsd	%xmm3, -632(%rbp)
	vmulsd	%xmm1, %xmm2, %xmm1
	vmovsd	%xmm1, -472(%rbp)
.Ltmp37:
.LBB0_23:                               # %L.B0202.1
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	leal	(%rax,%rax), %ecx
	addl	$2, %ecx
	leal	-2(%r15), %edx
	cmpl	%edx, %ecx
	cltq
	leaq	-632(%rbp,%rax,8), %rcx
	leaq	-464(%rbp), %rdx
	leaq	(%rdx,%rax,8), %rax
	cmovgeq	%rcx, %rax
	vmovsd	(%rax), %xmm1           # xmm1 = mem[0],zero
.Ltmp38:
	#DEBUG_VALUE: dresult <- $xmm1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovapd	%xmm1, -128(%rbp)       # 16-byte Spill
.Ltmp39:
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	.loc	1 87 1                  # backgroundfield/quadr.cpp:87:1
	vaddsd	%xmm1, %xmm0, %xmm0
.Ltmp40:
	#DEBUG_VALUE: result <- $xmm0
.LBB0_15:                               #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-56(%rbp), %rcx         # 8-byte Reload
	jmp	.LBB0_16
.Ltmp41:
	.p2align	4, 0x90
.LBB0_9:                                # %L.B0198.thread
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 83 1 is_stmt 1        # backgroundfield/quadr.cpp:83:1
	vmovsd	(%rcx), %xmm0           # xmm0 = mem[0],zero
	movq	-56(%rbp), %rcx         # 8-byte Reload
	vmovapd	.LCPI0_2(%rip), %xmm5   # xmm5 = [NaN,NaN]
.Ltmp42:
.LBB0_16:                               # %L.B0200
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movl	-48(%rbp), %r15d        # 4-byte Reload
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
.Ltmp43:
	#DEBUG_VALUE: result <- $xmm0
	.loc	1 170 1 is_stmt 1       # backgroundfield/quadr.cpp:170:1
	vandpd	-128(%rbp), %xmm5, %xmm1 # 16-byte Folded Reload
	vmovsd	-160(%rbp), %xmm2       # 8-byte Reload
                                        # xmm2 = mem[0],zero
	vucomisd	%xmm1, %xmm2
	ja	.LBB0_18
.Ltmp44:
# %bb.17:                               # %L.B0035
                                        #   in Loop: Header=BB0_1 Depth=1
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- $xmm0
	.loc	1 172 1                 # backgroundfield/quadr.cpp:172:1
	movq	-8(%rcx), %rax
	movq	%rax, (%rcx)
	movq	-96(%rbp), %rax         # 8-byte Reload
	.loc	1 173 1                 # backgroundfield/quadr.cpp:173:1
	vmovsd	.LCPI0_3(%rip), %xmm1   # xmm1 = mem[0],zero
	vmulsd	-8(%rax), %xmm1, %xmm1
	vmovsd	%xmm1, (%rax)
	.loc	1 174 1                 # backgroundfield/quadr.cpp:174:1
	incq	%r14
.Ltmp45:
	#DEBUG_VALUE: j <- $r14
	.loc	1 164 1                 # backgroundfield/quadr.cpp:164:1
	addq	$8, %rcx
	addq	$8, %rax
	movq	%rax, -96(%rbp)         # 8-byte Spill
	movq	-168(%rbp), %rax        # 8-byte Reload
	movq	%rax, %rbx
	.loc	1 174 1                 # backgroundfield/quadr.cpp:174:1
	cmpq	$8, %rax
	jne	.LBB0_1
.Ltmp46:
.LBB0_18:                               # %L.N0005
	#DEBUG_VALUE: result <- $xmm0
	#DEBUG_VALUE: j <- $r14
	#DEBUG_VALUE: dresult <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 88, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg_simple:func <- $r13
	#DEBUG_VALUE: Romberg_simple:absacc <- [DW_OP_constu 160, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- $xmm0
	.loc	1 176 1                 # backgroundfield/quadr.cpp:176:1
	addq	$600, %rsp              # imm = 0x258
	popq	%rbx
	popq	%r12
	popq	%r13
.Ltmp47:
	popq	%r14
.Ltmp48:
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp49:
.Lfunc_end0:
	.size	_Z14Romberg_simpleRK11T1DFunctionddd, .Lfunc_end0-_Z14Romberg_simpleRK11T1DFunctionddd
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _Z7RombergRK11T1DFunctionddd
.LCPI1_0:
	.quad	4602678819172646912     # double 0.5
.LCPI1_3:
	.quad	4158027847206421152     # double 1.0000000000000001E-30
.LCPI1_4:
	.quad	4598175219545276416     # double 0.25
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4
.LCPI1_1:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
.LCPI1_2:
	.quad	9223372036854775807     # double NaN
	.quad	9223372036854775807     # double NaN
	.text
	.globl	_Z7RombergRK11T1DFunctionddd
	.p2align	4, 0x90
	.type	_Z7RombergRK11T1DFunctionddd,@function
_Z7RombergRK11T1DFunctionddd:           # @_Z7RombergRK11T1DFunctionddd
.Lfunc_begin1:
	.loc	1 180 0                 # backgroundfield/quadr.cpp:180:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: Romberg:func <- $rdi
	#DEBUG_VALUE: Romberg:a <- $xmm0
	#DEBUG_VALUE: Romberg:b <- $xmm1
	#DEBUG_VALUE: Romberg:absacc <- $xmm2
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$984, %rsp              # imm = 0x3D8
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	vmovsd	%xmm2, -200(%rbp)       # 8-byte Spill
	movq	%rdi, %r12
.Ltmp50:
	#DEBUG_VALUE: absacc <- undef
	#DEBUG_VALUE: func <- undef
	movabsq	$4607182418800017408, %rax # imm = 0x3FF0000000000000
.Ltmp51:
	#DEBUG_VALUE: it <- 0
	#DEBUG_VALUE: b <- $xmm1
	#DEBUG_VALUE: a <- $xmm0
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- $xmm1
	#DEBUG_VALUE: Romberg:a <- $xmm0
	#DEBUG_VALUE: Romberg:func <- $r12
	.loc	1 187 1 prologue_end    # backgroundfield/quadr.cpp:187:1
	movq	%rax, -696(%rbp)
.Ltmp52:
	#DEBUG_VALUE: j <- 1
	#DEBUG_VALUE: dresult_pol <- 0.000000e+00
	.loc	1 190 1                 # backgroundfield/quadr.cpp:190:1
	leaq	-616(%rbp), %r14
	leaq	-688(%rbp), %rax
	vmovsd	%xmm1, -224(%rbp)       # 8-byte Spill
.Ltmp53:
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vmovsd	%xmm0, -144(%rbp)       # 8-byte Spill
.Ltmp54:
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	vsubsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -136(%rbp)       # 8-byte Spill
.Ltmp55:
	vmulsd	.LCPI1_0(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, -216(%rbp)       # 8-byte Spill
	xorl	%r13d, %r13d
	movl	$1, %ecx
                                        # implicit-def: $xmm2
	.loc	1 189 1 is_stmt 1       # backgroundfield/quadr.cpp:189:1
	vxorpd	%xmm9, %xmm9, %xmm9
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	xorl	%r8d, %r8d
	jmp	.LBB1_1
.Ltmp56:
	.p2align	4, 0x90
.LBB1_2:                                # %L..inline.10796.thread
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: j <- $rcx
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 49 1 is_stmt 1        # backgroundfield/quadr.cpp:49:1
	movq	(%r12), %rax
	movq	%r12, %rdi
	vmovsd	-224(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	movq	%rcx, %rbx
.Ltmp57:
	#DEBUG_VALUE: j <- $rbx
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%rdx, %r13
	.loc	1 49 1                  # backgroundfield/quadr.cpp:49:1
	callq	*(%rax)
.Ltmp58:
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovsd	%xmm0, -48(%rbp)        # 8-byte Spill
	.loc	1 49 1                  # backgroundfield/quadr.cpp:49:1
	movq	(%r12), %rax
	movq	%r12, %rdi
	vmovsd	-144(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	callq	*(%rax)
	vmovapd	-176(%rbp), %xmm2       # 16-byte Reload
.Ltmp59:
	#DEBUG_VALUE: dresult_rat <- $xmm2
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovapd	-128(%rbp), %xmm9       # 16-byte Reload
.Ltmp60:
	#DEBUG_VALUE: dresult_pol <- $xmm9
	movq	%rbx, %rcx
	.loc	1 49 1                  # backgroundfield/quadr.cpp:49:1
	vaddsd	-48(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmulsd	-216(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmovsd	%xmm0, -632(%rbp,%rbx,8)
.Ltmp61:
	#DEBUG_VALUE: it <- 1
	.loc	1 192 1 is_stmt 1       # backgroundfield/quadr.cpp:192:1
	leaq	-8(%r14), %rdx
	vmovsd	-8(%r14), %xmm8         # xmm8 = mem[0],zero
	movl	$1, %r8d
.Ltmp62:
.LBB1_3:                                # %L.B0041
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: dresult_rat <- $xmm2
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: dresult_rat <- $xmm2
	.loc	1 210 1                 # backgroundfield/quadr.cpp:210:1
	movq	(%rdx), %rax
	movq	%rax, (%r14)
	.loc	1 211 1                 # backgroundfield/quadr.cpp:211:1
	vmovsd	.LCPI1_4(%rip), %xmm0   # xmm0 = mem[0],zero
	movq	-232(%rbp), %rax        # 8-byte Reload
	vmulsd	-8(%rax), %xmm0, %xmm0
	vmovsd	%xmm0, (%rax)
	.loc	1 212 1                 # backgroundfield/quadr.cpp:212:1
	incq	%rcx
.Ltmp63:
	#DEBUG_VALUE: j <- $rcx
	.loc	1 190 1                 # backgroundfield/quadr.cpp:190:1
	addq	$8, %r14
	addq	$8, %rax
	.loc	1 212 1                 # backgroundfield/quadr.cpp:212:1
	cmpq	$8, %r13
	je	.LBB1_38
.Ltmp64:
.LBB1_1:                                # %L.B0039
                                        # =>This Loop Header: Depth=1
                                        #     Child Loop BB1_7 Depth 2
                                        #     Child Loop BB1_17 Depth 2
                                        #     Child Loop BB1_32 Depth 2
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%rax, -232(%rbp)        # 8-byte Spill
.Ltmp65:
	#DEBUG_VALUE: dresult_rat <- $xmm2
	#DEBUG_VALUE: j <- $rcx
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: it <- $r8d
	.loc	1 191 1 is_stmt 1       # backgroundfield/quadr.cpp:191:1
	cmpq	$5, %rcx
	movl	$5, %r15d
	cmovbq	%rcx, %r15
	leaq	1(%r13), %rdx
	cmpl	$5, %edx
	movl	%edx, %ebx
	movl	$5, %eax
	cmovgeq	%rax, %rbx
	cmpq	$1, %rcx
	vmovapd	%xmm9, -128(%rbp)       # 16-byte Spill
.Ltmp66:
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vmovapd	%xmm2, -176(%rbp)       # 16-byte Spill
.Ltmp67:
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	je	.LBB1_2
.Ltmp68:
# %bb.4:                                # %L..inline.10790
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: j <- $rcx
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	%rdx, -72(%rbp)         # 8-byte Spill
	movq	%rcx, -80(%rbp)         # 8-byte Spill
.Ltmp69:
	#DEBUG_VALUE: j <- [DW_OP_constu 80, DW_OP_minus] [$rbp+0]
	movq	%r14, -152(%rbp)        # 8-byte Spill
.Ltmp70:
	.loc	1 52 1 is_stmt 1        # backgroundfield/quadr.cpp:52:1
	vcvtsi2sd	%r8d, %xmm10, %xmm1
	.loc	1 56 1                  # backgroundfield/quadr.cpp:56:1
	testl	%r8d, %r8d
	movl	%r8d, -52(%rbp)         # 4-byte Spill
	vmovsd	%xmm1, -64(%rbp)        # 8-byte Spill
	jle	.LBB1_5
.Ltmp71:
# %bb.6:                                # %L.B0233
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: j <- [DW_OP_constu 80, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vmovsd	-136(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vdivsd	%xmm1, %xmm0, %xmm1
	vmovsd	-144(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	%xmm1, -88(%rbp)        # 8-byte Spill
	.loc	1 53 1 is_stmt 1        # backgroundfield/quadr.cpp:53:1
	vfmadd231sd	.LCPI1_0(%rip), %xmm1, %xmm0 # xmm0 = (xmm1 * mem) + xmm0
	vxorpd	%xmm1, %xmm1, %xmm1
	.loc	1 56 1                  # backgroundfield/quadr.cpp:56:1
	movl	%r8d, %r14d
.Ltmp72:
	#DEBUG_VALUE: it <- $r14d
	.p2align	4, 0x90
.LBB1_7:                                # %L..inline.10804
                                        #   Parent Loop BB1_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: j <- [DW_OP_constu 80, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vmovsd	%xmm0, -96(%rbp)        # 8-byte Spill
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	.loc	1 57 1 is_stmt 1        # backgroundfield/quadr.cpp:57:1
	movq	(%r12), %rax
	movq	%r12, %rdi
	callq	*(%rax)
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	%xmm0, %xmm1, %xmm1
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	vmovsd	-96(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	-88(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	.loc	1 59 1                  # backgroundfield/quadr.cpp:59:1
	decl	%r14d
	jg	.LBB1_7
	jmp	.LBB1_8
.Ltmp73:
.LBB1_5:                                #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: j <- [DW_OP_constu 80, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vxorpd	%xmm1, %xmm1, %xmm1
.Ltmp74:
.LBB1_8:                                # %L..inline.10796
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: j <- [DW_OP_constu 80, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 59 1                  # backgroundfield/quadr.cpp:59:1
	vmulsd	-136(%rbp), %xmm1, %xmm0 # 8-byte Folded Reload
	vdivsd	-64(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	movq	-80(%rbp), %rcx         # 8-byte Reload
.Ltmp75:
	#DEBUG_VALUE: j <- $rcx
	vaddsd	-632(%rbp,%rcx,8), %xmm0, %xmm0
	vmulsd	.LCPI1_0(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, -632(%rbp,%rcx,8)
	movl	-52(%rbp), %r8d         # 4-byte Reload
	.loc	1 60 1 is_stmt 1        # backgroundfield/quadr.cpp:60:1
	addl	%r8d, %r8d
.Ltmp76:
	#DEBUG_VALUE: it <- $r8d
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	-152(%rbp), %r14        # 8-byte Reload
.Ltmp77:
	.loc	1 192 1 is_stmt 1       # backgroundfield/quadr.cpp:192:1
	leaq	-8(%r14), %rdx
	vmovsd	-8(%r14), %xmm8         # xmm8 = mem[0],zero
.Ltmp78:
	#DEBUG_VALUE: result <- $xmm8
	.loc	1 193 1                 # backgroundfield/quadr.cpp:193:1
	cmpq	$2, %rcx
	jae	.LBB1_10
.Ltmp79:
# %bb.9:                                #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result <- $xmm8
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: j <- $rcx
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	-72(%rbp), %r13         # 8-byte Reload
	vmovapd	-128(%rbp), %xmm9       # 16-byte Reload
.Ltmp80:
	#DEBUG_VALUE: dresult_pol <- $xmm9
	vmovapd	-176(%rbp), %xmm2       # 16-byte Reload
.Ltmp81:
	#DEBUG_VALUE: dresult_rat <- $xmm2
	jmp	.LBB1_3
.Ltmp82:
	.p2align	4, 0x90
.LBB1_10:                               # %L.B0234
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result <- $xmm8
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: j <- $rcx
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 194 1 is_stmt 1       # backgroundfield/quadr.cpp:194:1
	cmpq	$5, %rcx
	movl	$5, %edi
	cmovbq	%rcx, %rdi
.Ltmp83:
	#DEBUG_VALUE: k1 <- $edi
	.loc	1 197 1                 # backgroundfield/quadr.cpp:197:1
	movl	%ecx, %eax
	subl	%edi, %eax
	cltq
	leaq	-624(,%rax,8), %rsi
	addq	%rbp, %rsi
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	testl	%edi, %edi
	je	.LBB1_11
.Ltmp84:
# %bb.12:                               # %L.B0235
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: k1 <- $edi
	#DEBUG_VALUE: result <- $xmm8
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: j <- $rcx
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%rdx, -104(%rbp)        # 8-byte Spill
	vmovsd	%xmm8, -64(%rbp)        # 8-byte Spill
.Ltmp85:
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	movl	%r8d, -52(%rbp)         # 4-byte Spill
.Ltmp86:
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
                                        # kill: def $ecx killed $ecx killed $rcx
	subl	%r15d, %ecx
.Ltmp87:
	movl	%ecx, -48(%rbp)         # 4-byte Spill
	leaq	-616(,%r13,8), %r13
	addq	%rbp, %r13
	shlq	$3, %rbx
	subq	%rbx, %r13
	leaq	-696(,%rax,8), %rax
	addq	%rbp, %rax
	movq	%rax, -208(%rbp)        # 8-byte Spill
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vmovsd	(%rax), %xmm0           # xmm0 = mem[0],zero
	vmovapd	%xmm0, -192(%rbp)       # 16-byte Spill
	movq	%rdi, -88(%rbp)         # 8-byte Spill
	leaq	-1016(%rbp), %rdi
.Ltmp88:
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	%rsi, -96(%rbp)         # 8-byte Spill
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	movq	%r13, %rsi
	movq	%rbx, %rdx
	callq	memcpy
	leaq	-856(%rbp), %rdi
	movq	%r13, %rsi
	movq	%rbx, %rdx
	callq	memcpy
	movq	-88(%rbp), %rsi         # 8-byte Reload
	movq	-96(%rbp), %rcx         # 8-byte Reload
	movq	-208(%rbp), %r10        # 8-byte Reload
	xorl	%edx, %edx
	cmpq	$2, %rsi
	vmovapd	.LCPI1_2(%rip), %xmm6   # xmm6 = [NaN,NaN]
	jb	.LBB1_14
.Ltmp89:
# %bb.13:                               # %L.B0226.L.B0226_crit_edge.lr.ph
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovapd	-192(%rbp), %xmm0       # 16-byte Reload
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vxorpd	.LCPI1_1(%rip), %xmm0, %xmm0
	vmovddup	%xmm0, %xmm0    # xmm0 = xmm0[0,0]
	vandpd	%xmm6, %xmm0, %xmm0
	vpermilpd	$1, %xmm0, %xmm1 # xmm1 = xmm0[1,0]
	vminsd	%xmm1, %xmm0, %xmm0
	vmovsd	8(%r10), %xmm1          # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	xorl	%edx, %edx
	vucomisd	%xmm1, %xmm0
	seta	%dl
	cmpq	$2, %rsi
	jbe	.LBB1_14
.Ltmp90:
# %bb.39:                               # %L.B0226.L.B0226_crit_edge.1
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	16(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$2, %eax
	cmoval	%eax, %edx
	cmpq	$3, %rsi
	je	.LBB1_14
.Ltmp91:
# %bb.40:                               # %L.B0226.L.B0226_crit_edge.2
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	24(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$3, %eax
	cmoval	%eax, %edx
	cmpq	$5, %rsi
	jb	.LBB1_14
.Ltmp92:
# %bb.41:                               # %L.B0226.L.B0226_crit_edge.3
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	32(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$4, %eax
	cmoval	%eax, %edx
.Ltmp93:
.LBB1_14:                               # %L.B0227
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	movslq	-48(%rbp), %rdi         # 4-byte Folded Reload
	.loc	1 83 1 is_stmt 1        # backgroundfield/quadr.cpp:83:1
	movl	%edx, %eax
	vmovsd	(%rcx,%rax,8), %xmm0    # xmm0 = mem[0],zero
.Ltmp94:
	#DEBUG_VALUE: result_pol <- $xmm0
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	cmpl	$2, %esi
	setb	%al
	vxorpd	%xmm8, %xmm8, %xmm8
	jae	.LBB1_16
.Ltmp95:
# %bb.15:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movb	$1, %r9b
	movl	-52(%rbp), %r8d         # 4-byte Reload
.Ltmp96:
	#DEBUG_VALUE: it <- $r8d
	vmovapd	-128(%rbp), %xmm9       # 16-byte Reload
.Ltmp97:
	#DEBUG_VALUE: dresult_pol <- $xmm9
	jmp	.LBB1_21
.Ltmp98:
.LBB1_11:                               # %L.B0229.thread
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: k1 <- $edi
	#DEBUG_VALUE: result <- $xmm8
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: j <- $rcx
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 83 1 is_stmt 1        # backgroundfield/quadr.cpp:83:1
	vmovsd	(%rsi), %xmm0           # xmm0 = mem[0],zero
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vmovapd	%xmm0, %xmm1
	movq	-72(%rbp), %r13         # 8-byte Reload
	vmovapd	-128(%rbp), %xmm9       # 16-byte Reload
.Ltmp99:
	#DEBUG_VALUE: dresult_pol <- $xmm9
	vmovapd	-176(%rbp), %xmm6       # 16-byte Reload
.Ltmp100:
	#DEBUG_VALUE: dresult_rat <- $xmm6
	jmp	.LBB1_36
.Ltmp101:
.LBB1_16:                               # %L.B0237
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- [DW_OP_constu 128, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	movb	%al, -192(%rbp)         # 1-byte Spill
	leaq	-688(%rbp), %rax
	movq	%rdi, -128(%rbp)        # 8-byte Spill
.Ltmp102:
	#DEBUG_VALUE: dresult_pol <- undef
	leaq	(%rax,%rdi,8), %rax
	movq	%rax, -48(%rbp)         # 8-byte Spill
	leaq	-1(%rsi), %r11
	.loc	1 86 1 is_stmt 1        # backgroundfield/quadr.cpp:86:1
	decl	%edx
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vmovsd	-824(%rbp), %xmm1       # xmm1 = mem[0],zero
	.loc	1 136 1 is_stmt 1       # backgroundfield/quadr.cpp:136:1
	movq	%r15, %r13
	movabsq	$4294967296, %rax       # imm = 0x100000000
	movabsq	$8589934592, %r9        # imm = 0x200000000
	movabsq	$12884901888, %r14      # imm = 0x300000000
	xorl	%edi, %edi
	jmp	.LBB1_17
.Ltmp103:
	.p2align	4, 0x90
.LBB1_19:                               # %L.B0231
                                        #   in Loop: Header=BB1_17 Depth=2
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	leal	(%rdx,%rdx), %ebx
	addl	$2, %ebx
	xorl	%esi, %esi
	cmpl	%r13d, %ebx
	setge	%sil
	movslq	%edx, %rdx
	leaq	-1016(%rbp,%rdx,8), %rbx
	leaq	-848(%rbp), %rcx
	leaq	(%rcx,%rdx,8), %r8
	cmovgeq	%rbx, %r8
	subl	%esi, %edx
	vmovsd	(%r8), %xmm9            # xmm9 = mem[0],zero
.Ltmp104:
	#DEBUG_VALUE: dresult_pol <- $xmm9
	vaddsd	%xmm9, %xmm0, %xmm0
.Ltmp105:
	#DEBUG_VALUE: result_pol <- $xmm0
	incq	%rdi
	movabsq	$4294967296, %rsi       # imm = 0x100000000
	addq	%rsi, %r14
	addq	%rsi, %r9
	addq	%rsi, %rax
	cmpq	%r11, %rdi
	jae	.LBB1_20
.Ltmp106:
.LBB1_17:                               # %L.B0228
                                        #   Parent Loop BB1_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	movq	%r13, %rbx
	#DEBUG_VALUE: result_pol <- $xmm0
	decq	%r13
	cmpl	$2, %ebx
	jl	.LBB1_19
.Ltmp107:
# %bb.18:                               # %L.B0230
                                        #   in Loop: Header=BB1_17 Depth=2
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	vmovsd	(%r10), %xmm2           # xmm2 = mem[0],zero
	vmovsd	-848(%rbp), %xmm3       # xmm3 = mem[0],zero
	vsubsd	-1016(%rbp), %xmm3, %xmm3
	movq	-48(%rbp), %rcx         # 8-byte Reload
	vmovsd	(%rcx,%rdi,8), %xmm4    # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm2, %xmm5
	vdivsd	%xmm5, %xmm3, %xmm3
	vmulsd	%xmm3, %xmm4, %xmm4
	vmovsd	%xmm4, -1016(%rbp)
	vmulsd	%xmm3, %xmm2, %xmm2
	vmovsd	%xmm2, -856(%rbp)
	cmpq	$1, %r13
	jle	.LBB1_19
.Ltmp108:
# %bb.42:                               # %L.B0230.1
                                        #   in Loop: Header=BB1_17 Depth=2
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	vmovsd	-840(%rbp), %xmm2       # xmm2 = mem[0],zero
	vsubsd	-1008(%rbp), %xmm2, %xmm2
	movq	%r10, %rcx
	vmovsd	8(%r10), %xmm3          # xmm3 = mem[0],zero
	movq	%rax, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm3, %xmm5
	vdivsd	%xmm5, %xmm2, %xmm2
	vmulsd	%xmm2, %xmm4, %xmm4
	vmovsd	%xmm4, -1008(%rbp)
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovsd	%xmm2, -848(%rbp)
	cmpq	$3, %rbx
	je	.LBB1_19
.Ltmp109:
# %bb.43:                               # %L.B0230.2
                                        #   in Loop: Header=BB1_17 Depth=2
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	vmovsd	-832(%rbp), %xmm2       # xmm2 = mem[0],zero
	vsubsd	-1000(%rbp), %xmm2, %xmm2
	movq	%r10, %rcx
	vmovsd	16(%r10), %xmm3         # xmm3 = mem[0],zero
	movq	%r9, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm3, %xmm5
	vdivsd	%xmm5, %xmm2, %xmm2
	vmulsd	%xmm2, %xmm4, %xmm4
	vmovsd	%xmm4, -1000(%rbp)
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovsd	%xmm2, -840(%rbp)
	cmpq	$4, %r13
	jl	.LBB1_19
.Ltmp110:
# %bb.44:                               # %L.B0230.3
                                        #   in Loop: Header=BB1_17 Depth=2
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r10, %rcx
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vmovsd	24(%r10), %xmm2         # xmm2 = mem[0],zero
	vsubsd	-992(%rbp), %xmm1, %xmm3
	movq	%r14, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm2, %xmm5
	vdivsd	%xmm5, %xmm3, %xmm3
	vmulsd	%xmm3, %xmm4, %xmm4
	vmovsd	%xmm4, -992(%rbp)
	vmulsd	%xmm3, %xmm2, %xmm2
	vmovsd	%xmm2, -832(%rbp)
	jmp	.LBB1_19
.Ltmp111:
.LBB1_20:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: it <- [DW_OP_constu 52, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-152(%rbp), %r14        # 8-byte Reload
	movl	-52(%rbp), %r8d         # 4-byte Reload
.Ltmp112:
	#DEBUG_VALUE: it <- $r8d
	movq	-96(%rbp), %rcx         # 8-byte Reload
	movq	-88(%rbp), %rsi         # 8-byte Reload
	movq	-128(%rbp), %rdi        # 8-byte Reload
	movb	-192(%rbp), %r9b        # 1-byte Reload
.Ltmp113:
.LBB1_21:                               # %L.B0229
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: result_pol <- $xmm0
	xorl	%edx, %edx
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	testl	%esi, %esi
	.loc	1 119 1 is_stmt 1       # backgroundfield/quadr.cpp:119:1
	je	.LBB1_29
.Ltmp114:
# %bb.22:                               # %L.B0240
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 116 1                 # backgroundfield/quadr.cpp:116:1
	vmovsd	(%r10), %xmm1           # xmm1 = mem[0],zero
	.loc	1 120 1                 # backgroundfield/quadr.cpp:120:1
	vucomisd	%xmm8, %xmm1
	jne	.LBB1_25
	jnp	.LBB1_23
.Ltmp115:
.LBB1_25:                               # %L..inline.10861.lr.ph
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 129 1                 # backgroundfield/quadr.cpp:129:1
	movq	(%rcx), %rax
	movq	%rax, -552(%rbp)
	vmovq	%rax, %xmm2
	vaddsd	.LCPI1_3(%rip), %xmm2, %xmm2
	vmovsd	%xmm2, -392(%rbp)
	.loc	1 131 1                 # backgroundfield/quadr.cpp:131:1
	cmpq	$1, %rsi
	jbe	.LBB1_29
.Ltmp116:
# %bb.26:                               # %L..inline.10861.L..inline.10857_crit_edge
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 119 1                 # backgroundfield/quadr.cpp:119:1
	vmovsd	8(%r10), %xmm2          # xmm2 = mem[0],zero
	.loc	1 120 1                 # backgroundfield/quadr.cpp:120:1
	vucomisd	%xmm8, %xmm2
	jne	.LBB1_45
	jnp	.LBB1_27
.Ltmp117:
.LBB1_45:                               # %L..inline.10861.1
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vandpd	%xmm6, %xmm1, %xmm1
	vminsd	%xmm1, %xmm1, %xmm1
	vandpd	%xmm6, %xmm2, %xmm2
	.loc	1 125 1 is_stmt 1       # backgroundfield/quadr.cpp:125:1
	xorl	%edx, %edx
	vucomisd	%xmm2, %xmm1
	seta	%dl
	.loc	1 129 1                 # backgroundfield/quadr.cpp:129:1
	movq	8(%rcx), %rax
	movq	%rax, -544(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI1_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -384(%rbp)
	.loc	1 131 1                 # backgroundfield/quadr.cpp:131:1
	cmpq	$3, %rsi
	jb	.LBB1_29
.Ltmp118:
# %bb.46:                               # %L..inline.10861.L..inline.10857_crit_edge.1
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 119 1                 # backgroundfield/quadr.cpp:119:1
	vmovsd	16(%r10), %xmm3         # xmm3 = mem[0],zero
	.loc	1 120 1                 # backgroundfield/quadr.cpp:120:1
	vucomisd	%xmm8, %xmm3
	jne	.LBB1_48
	jnp	.LBB1_47
.Ltmp119:
.LBB1_48:                               # %L..inline.10861.2
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vminsd	%xmm1, %xmm2, %xmm1
	vandpd	%xmm6, %xmm3, %xmm2
	.loc	1 125 1 is_stmt 1       # backgroundfield/quadr.cpp:125:1
	vucomisd	%xmm2, %xmm1
	movl	$2, %eax
	cmoval	%eax, %edx
	.loc	1 129 1                 # backgroundfield/quadr.cpp:129:1
	movq	16(%rcx), %rax
	movq	%rax, -536(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI1_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -376(%rbp)
	.loc	1 131 1                 # backgroundfield/quadr.cpp:131:1
	cmpq	$4, %rsi
	jb	.LBB1_29
.Ltmp120:
# %bb.49:                               # %L..inline.10861.L..inline.10857_crit_edge.2
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 119 1                 # backgroundfield/quadr.cpp:119:1
	vmovsd	24(%r10), %xmm3         # xmm3 = mem[0],zero
	.loc	1 120 1                 # backgroundfield/quadr.cpp:120:1
	vucomisd	%xmm8, %xmm3
	jne	.LBB1_51
	jnp	.LBB1_50
.Ltmp121:
.LBB1_51:                               # %L..inline.10861.3
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vminsd	%xmm1, %xmm2, %xmm1
	vandpd	%xmm6, %xmm3, %xmm2
	.loc	1 125 1 is_stmt 1       # backgroundfield/quadr.cpp:125:1
	vucomisd	%xmm2, %xmm1
	movl	$3, %eax
	cmoval	%eax, %edx
	.loc	1 129 1                 # backgroundfield/quadr.cpp:129:1
	movq	24(%rcx), %rax
	movq	%rax, -528(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI1_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -368(%rbp)
	.loc	1 131 1                 # backgroundfield/quadr.cpp:131:1
	cmpq	$5, %rsi
	jb	.LBB1_29
.Ltmp122:
# %bb.52:                               # %L..inline.10861.L..inline.10857_crit_edge.3
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 119 1                 # backgroundfield/quadr.cpp:119:1
	vmovsd	32(%r10), %xmm3         # xmm3 = mem[0],zero
	.loc	1 120 1                 # backgroundfield/quadr.cpp:120:1
	vucomisd	%xmm8, %xmm3
	jne	.LBB1_28
	jnp	.LBB1_53
.Ltmp123:
.LBB1_28:                               # %L..inline.10861.4
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 119 1                 # backgroundfield/quadr.cpp:119:1
	vandpd	%xmm6, %xmm3, %xmm3
	.loc	1 125 1                 # backgroundfield/quadr.cpp:125:1
	vminsd	%xmm1, %xmm2, %xmm1
	vucomisd	%xmm3, %xmm1
	movl	$4, %eax
	cmoval	%eax, %edx
	.loc	1 129 1                 # backgroundfield/quadr.cpp:129:1
	movq	32(%rcx), %rax
	movq	%rax, -520(%rbp)
	vmovq	%rax, %xmm1
	vaddsd	.LCPI1_3(%rip), %xmm1, %xmm1
	vmovsd	%xmm1, -360(%rbp)
.Ltmp124:
.LBB1_29:                               # %L..inline.10859
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 133 1                 # backgroundfield/quadr.cpp:133:1
	movl	%edx, %eax
	vmovsd	(%rcx,%rax,8), %xmm1    # xmm1 = mem[0],zero
.Ltmp125:
	#DEBUG_VALUE: result_rat <- $xmm1
	.loc	1 135 1                 # backgroundfield/quadr.cpp:135:1
	testb	%r9b, %r9b
	je	.LBB1_31
.Ltmp126:
# %bb.30:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	-80(%rbp), %rcx         # 8-byte Reload
	movq	-72(%rbp), %r13         # 8-byte Reload
	vmovsd	-64(%rbp), %xmm8        # 8-byte Reload
                                        # xmm8 = mem[0],zero
.Ltmp127:
	#DEBUG_VALUE: result <- $xmm8
	movq	-104(%rbp), %rdx        # 8-byte Reload
	vmovapd	-176(%rbp), %xmm6       # 16-byte Reload
.Ltmp128:
	#DEBUG_VALUE: dresult_rat <- $xmm6
	jmp	.LBB1_36
.Ltmp129:
.LBB1_31:                               # %L.B0243
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	leaq	-688(%rbp), %rax
	leaq	(%rax,%rdi,8), %rax
	addq	$24, %rax
	decq	%r15
	.loc	1 133 1 is_stmt 1       # backgroundfield/quadr.cpp:133:1
	decl	%edx
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vmovsd	(%r10), %xmm2           # xmm2 = mem[0],zero
	vmovsd	-520(%rbp), %xmm3       # xmm3 = mem[0],zero
	jmp	.LBB1_32
.Ltmp130:
.LBB1_59:                               # %L..inline.10877.3
                                        #   in Loop: Header=BB1_32 Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 136 1 is_stmt 1       # backgroundfield/quadr.cpp:136:1
	vsubsd	%xmm5, %xmm3, %xmm5
	.loc	1 146 1                 # backgroundfield/quadr.cpp:146:1
	vdivsd	%xmm6, %xmm5, %xmm5
	vmulsd	%xmm5, %xmm3, %xmm6
	vmovsd	%xmm6, -368(%rbp)
	.loc	1 147 1                 # backgroundfield/quadr.cpp:147:1
	vmulsd	%xmm5, %xmm4, %xmm4
	vmovsd	%xmm4, -528(%rbp)
.Ltmp131:
	.p2align	4, 0x90
.LBB1_34:                               # %L..inline.10873
                                        #   in Loop: Header=BB1_32 Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 149 1                 # backgroundfield/quadr.cpp:149:1
	leal	(%rdx,%rdx), %ecx
	addl	$2, %ecx
	movslq	%ecx, %rcx
	xorl	%esi, %esi
	cmpq	%rcx, %r15
	setle	%sil
	movslq	%edx, %rdx
	leaq	-392(%rbp,%rdx,8), %rcx
	leaq	-544(%rbp), %rdi
	leaq	(%rdi,%rdx,8), %rdi
	cmovleq	%rcx, %rdi
	subl	%esi, %edx
	.loc	1 151 1                 # backgroundfield/quadr.cpp:151:1
	vmovsd	(%rdi), %xmm6           # xmm6 = mem[0],zero
.Ltmp132:
	#DEBUG_VALUE: dresult_rat <- $xmm6
	vaddsd	%xmm6, %xmm1, %xmm1
.Ltmp133:
	#DEBUG_VALUE: result_rat <- $xmm1
	.loc	1 152 1                 # backgroundfield/quadr.cpp:152:1
	addq	$8, %rax
	decq	%r15
	testl	%r15d, %r15d
	jle	.LBB1_35
.Ltmp134:
.LBB1_32:                               # %L..inline.10872.preheader
                                        #   Parent Loop BB1_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result_rat <- $xmm1
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vmovsd	-544(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-392(%rbp), %xmm6       # xmm6 = mem[0],zero
	.loc	1 138 1                 # backgroundfield/quadr.cpp:138:1
	vmulsd	%xmm2, %xmm6, %xmm5
	vdivsd	-24(%rax), %xmm5, %xmm5
	.loc	1 139 1                 # backgroundfield/quadr.cpp:139:1
	vsubsd	%xmm4, %xmm5, %xmm7
	.loc	1 140 1                 # backgroundfield/quadr.cpp:140:1
	vucomisd	%xmm8, %xmm7
	jne	.LBB1_33
	jnp	.LBB1_60
.Ltmp135:
.LBB1_33:                               # %L..inline.10877
                                        #   in Loop: Header=BB1_32 Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vsubsd	%xmm6, %xmm4, %xmm6
	.loc	1 146 1                 # backgroundfield/quadr.cpp:146:1
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -392(%rbp)
	.loc	1 147 1                 # backgroundfield/quadr.cpp:147:1
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -552(%rbp)
	.loc	1 149 1                 # backgroundfield/quadr.cpp:149:1
	cmpq	$1, %r15
	jle	.LBB1_34
.Ltmp136:
# %bb.54:                               # %L..inline.10872.1
                                        #   in Loop: Header=BB1_32 Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vmovsd	-536(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-384(%rbp), %xmm6       # xmm6 = mem[0],zero
	.loc	1 138 1                 # backgroundfield/quadr.cpp:138:1
	vmulsd	8(%r10), %xmm6, %xmm5
	vdivsd	-16(%rax), %xmm5, %xmm5
	.loc	1 139 1                 # backgroundfield/quadr.cpp:139:1
	vsubsd	%xmm4, %xmm5, %xmm7
	.loc	1 140 1                 # backgroundfield/quadr.cpp:140:1
	vucomisd	%xmm8, %xmm7
	jne	.LBB1_55
	jnp	.LBB1_60
.Ltmp137:
.LBB1_55:                               # %L..inline.10877.1
                                        #   in Loop: Header=BB1_32 Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vsubsd	%xmm6, %xmm4, %xmm6
	.loc	1 146 1                 # backgroundfield/quadr.cpp:146:1
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -384(%rbp)
	.loc	1 147 1                 # backgroundfield/quadr.cpp:147:1
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -544(%rbp)
	.loc	1 149 1                 # backgroundfield/quadr.cpp:149:1
	cmpq	$3, %r15
	jl	.LBB1_34
.Ltmp138:
# %bb.56:                               # %L..inline.10872.2
                                        #   in Loop: Header=BB1_32 Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vmovsd	-528(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-376(%rbp), %xmm6       # xmm6 = mem[0],zero
	.loc	1 138 1                 # backgroundfield/quadr.cpp:138:1
	vmulsd	16(%r10), %xmm6, %xmm5
	vdivsd	-8(%rax), %xmm5, %xmm5
	.loc	1 139 1                 # backgroundfield/quadr.cpp:139:1
	vsubsd	%xmm4, %xmm5, %xmm7
	.loc	1 140 1                 # backgroundfield/quadr.cpp:140:1
	vucomisd	%xmm8, %xmm7
	jne	.LBB1_57
	jnp	.LBB1_60
.Ltmp139:
.LBB1_57:                               # %L..inline.10877.2
                                        #   in Loop: Header=BB1_32 Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vsubsd	%xmm6, %xmm4, %xmm6
	.loc	1 146 1                 # backgroundfield/quadr.cpp:146:1
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -376(%rbp)
	.loc	1 147 1                 # backgroundfield/quadr.cpp:147:1
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -536(%rbp)
	.loc	1 149 1                 # backgroundfield/quadr.cpp:149:1
	cmpq	$4, %r15
	jl	.LBB1_34
.Ltmp140:
# %bb.58:                               # %L..inline.10872.3
                                        #   in Loop: Header=BB1_32 Depth=2
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 136 1                 # backgroundfield/quadr.cpp:136:1
	vmovsd	-368(%rbp), %xmm5       # xmm5 = mem[0],zero
	.loc	1 138 1                 # backgroundfield/quadr.cpp:138:1
	vmulsd	24(%r10), %xmm5, %xmm4
	vdivsd	(%rax), %xmm4, %xmm4
	.loc	1 139 1                 # backgroundfield/quadr.cpp:139:1
	vsubsd	%xmm3, %xmm4, %xmm6
	.loc	1 140 1                 # backgroundfield/quadr.cpp:140:1
	vucomisd	%xmm8, %xmm6
	jne	.LBB1_59
	jp	.LBB1_59
	jmp	.LBB1_60
.Ltmp141:
.LBB1_23:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	xorl	%eax, %eax
	jmp	.LBB1_24
.Ltmp142:
.LBB1_27:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movl	$1, %eax
	jmp	.LBB1_24
.Ltmp143:
.LBB1_47:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	movl	$2, %eax
	jmp	.LBB1_24
.Ltmp144:
.LBB1_50:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	movl	$3, %eax
	jmp	.LBB1_24
.Ltmp145:
.LBB1_53:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	movl	$4, %eax
.Ltmp146:
.LBB1_24:                               # %L.B0241
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- [DW_OP_constu 176, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 121 1 is_stmt 1       # backgroundfield/quadr.cpp:121:1
	vmovsd	(%rcx,%rax,8), %xmm1    # xmm1 = mem[0],zero
.Ltmp147:
	#DEBUG_VALUE: result_rat <- $xmm1
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	vxorpd	%xmm6, %xmm6, %xmm6
.Ltmp148:
	#DEBUG_VALUE: dresult_rat <- 0.000000e+00
.LBB1_35:                               #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	movq	-80(%rbp), %rcx         # 8-byte Reload
	movq	-72(%rbp), %r13         # 8-byte Reload
	vmovsd	-64(%rbp), %xmm8        # 8-byte Reload
                                        # xmm8 = mem[0],zero
.Ltmp149:
	#DEBUG_VALUE: result <- $xmm8
	movq	-104(%rbp), %rdx        # 8-byte Reload
.Ltmp150:
.LBB1_36:                               # %L..inline.10916
                                        #   in Loop: Header=BB1_1 Depth=1
	#DEBUG_VALUE: result <- $xmm8
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: dresult_rat <- $xmm6
	#DEBUG_VALUE: result_rat <- $xmm1
	vmovapd	.LCPI1_2(%rip), %xmm5   # xmm5 = [NaN,NaN]
	.loc	1 199 1 is_stmt 1       # backgroundfield/quadr.cpp:199:1
	vandpd	%xmm5, %xmm9, %xmm9
.Ltmp151:
	#DEBUG_VALUE: dresult_pol <- $xmm9
	.loc	1 200 1                 # backgroundfield/quadr.cpp:200:1
	vandpd	%xmm5, %xmm6, %xmm2
.Ltmp152:
	#DEBUG_VALUE: dresult_rat <- $xmm2
	.loc	1 201 1                 # backgroundfield/quadr.cpp:201:1
	vmaxsd	%xmm9, %xmm2, %xmm3
	vcmpunordsd	%xmm9, %xmm9, %xmm4
	vblendvpd	%xmm4, %xmm2, %xmm3, %xmm3
	vcmpunordsd	%xmm6, %xmm6, %xmm4
	vblendvpd	%xmm4, %xmm9, %xmm3, %xmm3
	.loc	1 202 1                 # backgroundfield/quadr.cpp:202:1
	vsubsd	%xmm1, %xmm0, %xmm4
	vandpd	%xmm5, %xmm4, %xmm5
	.loc	1 203 1                 # backgroundfield/quadr.cpp:203:1
	vcmpunordsd	%xmm3, %xmm3, %xmm6
	vmaxsd	%xmm3, %xmm5, %xmm7
	vblendvpd	%xmm6, %xmm5, %xmm7, %xmm5
	vcmpunordsd	%xmm4, %xmm4, %xmm4
	vblendvpd	%xmm4, %xmm3, %xmm5, %xmm3
	vmovsd	-200(%rbp), %xmm4       # 8-byte Reload
                                        # xmm4 = mem[0],zero
	.loc	1 204 1                 # backgroundfield/quadr.cpp:204:1
	vucomisd	%xmm3, %xmm4
	jbe	.LBB1_3
.Ltmp153:
# %bb.37:                               # %L.B0248
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result <- $xmm8
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: dresult_rat <- $xmm2
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result_rat <- $xmm1
	.loc	1 206 1                 # backgroundfield/quadr.cpp:206:1
	vaddsd	%xmm1, %xmm0, %xmm0
	vmulsd	.LCPI1_0(%rip), %xmm0, %xmm8
.Ltmp154:
	#DEBUG_VALUE: result <- $xmm8
.LBB1_38:                               # %L_T22936352_7633
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: dresult_rat <- $xmm2
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: result <- $xmm8
	.loc	1 214 1                 # backgroundfield/quadr.cpp:214:1
	vmovapd	%xmm8, %xmm0
	addq	$984, %rsp              # imm = 0x3D8
	popq	%rbx
	popq	%r12
.Ltmp155:
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp156:
.LBB1_60:                               # %L..inline.10882
	.cfi_def_cfa %rbp, 16
	#DEBUG_VALUE: result_rat <- $xmm1
	#DEBUG_VALUE: result_pol <- $xmm0
	#DEBUG_VALUE: result <- [DW_OP_constu 64, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: it <- $r8d
	#DEBUG_VALUE: dresult_pol <- $xmm9
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 144, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 224, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $r12
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 200, DW_OP_minus] [$rbp+0]
	.loc	1 142 1                 # backgroundfield/quadr.cpp:142:1
	movl	$_ZSt4cerr, %edi
	movl	$.S08003, %esi
	movl	$20, %edx
	callq	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l
.Ltmp157:
	.loc	1 143 1                 # backgroundfield/quadr.cpp:143:1
	movl	$111, %edi
	callq	exit
.Ltmp158:
.Lfunc_end1:
	.size	_Z7RombergRK11T1DFunctionddd, .Lfunc_end1-_Z7RombergRK11T1DFunctionddd
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _Z7RombergRK11T2DFunctionddddd
.LCPI2_0:
	.quad	4602678819172646912     # double 0.5
.LCPI2_3:
	.quad	4158027847206421152     # double 1.0000000000000001E-30
.LCPI2_4:
	.quad	4598175219545276416     # double 0.25
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4
.LCPI2_1:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
.LCPI2_2:
	.quad	9223372036854775807     # double NaN
	.quad	9223372036854775807     # double NaN
	.text
	.globl	_Z7RombergRK11T2DFunctionddddd
	.p2align	4, 0x90
	.type	_Z7RombergRK11T2DFunctionddddd,@function
_Z7RombergRK11T2DFunctionddddd:         # @_Z7RombergRK11T2DFunctionddddd
.Lfunc_begin2:
	.loc	1 230 0                 # backgroundfield/quadr.cpp:230:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: Romberg:func <- $rdi
	#DEBUG_VALUE: Romberg:a <- $xmm0
	#DEBUG_VALUE: Romberg:b <- $xmm1
	#DEBUG_VALUE: Romberg:c <- $xmm2
	#DEBUG_VALUE: Romberg:d <- $xmm3
	#DEBUG_VALUE: Romberg:absacc <- $xmm4
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$1016, %rsp             # imm = 0x3F8
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	vmovsd	%xmm1, -256(%rbp)       # 8-byte Spill
	vmovsd	%xmm0, -136(%rbp)       # 8-byte Spill
.Ltmp159:
	#DEBUG_VALUE: absacc <- $xmm4
	#DEBUG_VALUE: d <- $xmm3
	#DEBUG_VALUE: c <- $xmm2
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: func <- $rdi
	#DEBUG_VALUE: Romberg:absacc <- $xmm4
	#DEBUG_VALUE: Romberg:d <- $xmm3
	#DEBUG_VALUE: Romberg:c <- $xmm2
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:func <- $rdi
	.loc	1 232 1 prologue_end    # backgroundfield/quadr.cpp:232:1
	vsubsd	%xmm0, %xmm1, %xmm1
	vmovsd	%xmm4, -232(%rbp)       # 8-byte Spill
.Ltmp160:
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	vdivsd	%xmm1, %xmm4, %xmm0
	.loc	1 225 1                 # backgroundfield/quadr.cpp:225:1
	movq	$_ZTV9Tinty_f2D+16, -184(%rbp)
	movq	%rdi, -176(%rbp)
	vmovsd	%xmm2, -168(%rbp)
	vmovsd	%xmm3, -160(%rbp)
	vmovsd	%xmm0, -152(%rbp)
	movabsq	$4607182418800017408, %rax # imm = 0x3FF0000000000000
	.loc	1 184 1                 # backgroundfield/quadr.cpp:184:1
	movq	%rax, -736(%rbp)
	.loc	1 190 1                 # backgroundfield/quadr.cpp:190:1
	leaq	-656(%rbp), %r15
	leaq	-728(%rbp), %rax
	vmovsd	%xmm1, -128(%rbp)       # 8-byte Spill
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vmulsd	.LCPI2_0(%rip), %xmm1, %xmm0
	vmovsd	%xmm0, -248(%rbp)       # 8-byte Spill
	xorl	%r13d, %r13d
	movl	$1, %r14d
	leaq	-184(%rbp), %rdi
.Ltmp161:
                                        # implicit-def: $xmm2
	.loc	1 188 1 is_stmt 1       # backgroundfield/quadr.cpp:188:1
	vxorpd	%xmm9, %xmm9, %xmm9
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	xorl	%edx, %edx
	jmp	.LBB2_1
.Ltmp162:
	.p2align	4, 0x90
.LBB2_2:                                # %L.B0272
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 191 1 is_stmt 1       # backgroundfield/quadr.cpp:191:1
	movq	-184(%rbp), %rax
	vmovsd	-256(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	movq	%rcx, %r13
	callq	*(%rax)
	vmovsd	%xmm0, -48(%rbp)        # 8-byte Spill
	movq	-184(%rbp), %rax
	leaq	-184(%rbp), %rdi
	vmovsd	-136(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	callq	*(%rax)
	vmovapd	-208(%rbp), %xmm2       # 16-byte Reload
	vmovapd	-112(%rbp), %xmm9       # 16-byte Reload
	leaq	-184(%rbp), %rdi
	vaddsd	-48(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmulsd	-248(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmovsd	%xmm0, -672(%rbp,%r14,8)
	leaq	-8(%r15), %r8
	vmovsd	-8(%r15), %xmm8         # xmm8 = mem[0],zero
	movl	$1, %edx
.Ltmp163:
.LBB2_3:                                # %L..inline.11466
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 210 1                 # backgroundfield/quadr.cpp:210:1
	movq	(%r8), %rax
	movq	%rax, (%r15)
	vmovsd	.LCPI2_4(%rip), %xmm0   # xmm0 = mem[0],zero
	movq	-272(%rbp), %rax        # 8-byte Reload
	vmulsd	-8(%rax), %xmm0, %xmm0
	vmovsd	%xmm0, (%rax)
	.loc	1 212 1                 # backgroundfield/quadr.cpp:212:1
	incq	%r14
	.loc	1 190 1                 # backgroundfield/quadr.cpp:190:1
	addq	$8, %r15
	addq	$8, %rax
	.loc	1 212 1                 # backgroundfield/quadr.cpp:212:1
	cmpq	$8, %r13
	je	.LBB2_38
.Ltmp164:
.LBB2_1:                                # %L..inline.11439
                                        # =>This Loop Header: Depth=1
                                        #     Child Loop BB2_7 Depth 2
                                        #     Child Loop BB2_17 Depth 2
                                        #     Child Loop BB2_32 Depth 2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%rax, -272(%rbp)        # 8-byte Spill
	.loc	1 191 1 is_stmt 1       # backgroundfield/quadr.cpp:191:1
	cmpq	$5, %r14
	movl	$5, %r12d
	cmovbq	%r14, %r12
	leaq	1(%r13), %rcx
	cmpl	$5, %ecx
	movl	%ecx, %ebx
	movl	$5, %eax
	cmovgeq	%rax, %rbx
	cmpq	$1, %r14
	vmovapd	%xmm9, -112(%rbp)       # 16-byte Spill
	vmovapd	%xmm2, -208(%rbp)       # 16-byte Spill
	je	.LBB2_2
.Ltmp165:
# %bb.4:                                # %L..inline.11453
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r14, -264(%rbp)        # 8-byte Spill
	movq	%rcx, -64(%rbp)         # 8-byte Spill
	movq	%r15, -144(%rbp)        # 8-byte Spill
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vcvtsi2sd	%edx, %xmm10, %xmm1
	testl	%edx, %edx
	vmovsd	%xmm1, -56(%rbp)        # 8-byte Spill
	jle	.LBB2_5
.Ltmp166:
# %bb.6:                                # %L.B0275
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovsd	-128(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vdivsd	%xmm1, %xmm0, %xmm1
	vmovsd	-136(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	%xmm1, -72(%rbp)        # 8-byte Spill
	vfmadd231sd	.LCPI2_0(%rip), %xmm1, %xmm0 # xmm0 = (xmm1 * mem) + xmm0
	vxorpd	%xmm1, %xmm1, %xmm1
	movl	%edx, %r14d
	movl	%edx, %r15d
.Ltmp167:
	.p2align	4, 0x90
.LBB2_7:                                # %L..inline.11463
                                        #   Parent Loop BB2_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovsd	%xmm0, -80(%rbp)        # 8-byte Spill
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	movq	-184(%rbp), %rax
	callq	*(%rax)
	leaq	-184(%rbp), %rdi
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	%xmm0, %xmm1, %xmm1
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	vmovsd	-80(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	-72(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	decl	%r15d
	jg	.LBB2_7
	jmp	.LBB2_8
.Ltmp168:
.LBB2_5:                                #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movl	%edx, %r14d
	vxorpd	%xmm1, %xmm1, %xmm1
.Ltmp169:
.LBB2_8:                                # %L..inline.11455
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vmulsd	-128(%rbp), %xmm1, %xmm0 # 8-byte Folded Reload
	vdivsd	-56(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	movq	-264(%rbp), %rax        # 8-byte Reload
	vaddsd	-672(%rbp,%rax,8), %xmm0, %xmm0
	vmulsd	.LCPI2_0(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, -672(%rbp,%rax,8)
	movl	%r14d, %edx
	movq	%rax, %r14
	addl	%edx, %edx
	movq	-144(%rbp), %r15        # 8-byte Reload
	leaq	-8(%r15), %r8
	vmovsd	-8(%r15), %xmm8         # xmm8 = mem[0],zero
	.loc	1 192 1 is_stmt 1       # backgroundfield/quadr.cpp:192:1
	cmpq	$2, %rax
	jae	.LBB2_10
.Ltmp170:
# %bb.9:                                #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	vmovapd	-112(%rbp), %xmm9       # 16-byte Reload
	vmovapd	-208(%rbp), %xmm2       # 16-byte Reload
	jmp	.LBB2_3
.Ltmp171:
	.p2align	4, 0x90
.LBB2_10:                               # %L.B0277
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 194 1 is_stmt 1       # backgroundfield/quadr.cpp:194:1
	cmpq	$5, %r14
	movl	$5, %esi
	cmovbq	%r14, %rsi
	.loc	1 197 1                 # backgroundfield/quadr.cpp:197:1
	movl	%r14d, %eax
	subl	%esi, %eax
	cltq
	leaq	-664(,%rax,8), %rcx
	addq	%rbp, %rcx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	testl	%esi, %esi
	je	.LBB2_11
.Ltmp172:
# %bb.12:                               # %L.B0278
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r8, -120(%rbp)         # 8-byte Spill
	vmovsd	%xmm8, -88(%rbp)        # 8-byte Spill
	movl	%edx, -56(%rbp)         # 4-byte Spill
	movl	%r14d, %edx
	subl	%r12d, %edx
	movl	%edx, -48(%rbp)         # 4-byte Spill
	leaq	-656(,%r13,8), %r13
	addq	%rbp, %r13
	shlq	$3, %rbx
	subq	%rbx, %r13
	leaq	-736(,%rax,8), %rax
	addq	%rbp, %rax
	movq	%rax, -240(%rbp)        # 8-byte Spill
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	(%rax), %xmm0           # xmm0 = mem[0],zero
	vmovapd	%xmm0, -224(%rbp)       # 16-byte Spill
	leaq	-1056(%rbp), %rdi
	movq	%rsi, -72(%rbp)         # 8-byte Spill
	movq	%r13, %rsi
	movq	%rbx, %rdx
	movq	%rcx, -80(%rbp)         # 8-byte Spill
	callq	memcpy
	leaq	-896(%rbp), %rdi
	movq	%r13, %rsi
	movq	%rbx, %rdx
	callq	memcpy
	movq	-72(%rbp), %rsi         # 8-byte Reload
	movq	-80(%rbp), %rcx         # 8-byte Reload
	movq	-240(%rbp), %r10        # 8-byte Reload
	xorl	%edx, %edx
	cmpq	$2, %rsi
	vmovapd	.LCPI2_2(%rip), %xmm6   # xmm6 = [NaN,NaN]
	jb	.LBB2_14
.Ltmp173:
# %bb.13:                               # %L.B0266.L.B0266_crit_edge.lr.ph
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovapd	-224(%rbp), %xmm0       # 16-byte Reload
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vxorpd	.LCPI2_1(%rip), %xmm0, %xmm0
	vmovddup	%xmm0, %xmm0    # xmm0 = xmm0[0,0]
	vandpd	%xmm6, %xmm0, %xmm0
	vpermilpd	$1, %xmm0, %xmm1 # xmm1 = xmm0[1,0]
	vminsd	%xmm1, %xmm0, %xmm0
	vmovsd	8(%r10), %xmm1          # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	xorl	%edx, %edx
	vucomisd	%xmm1, %xmm0
	seta	%dl
	cmpq	$2, %rsi
	jbe	.LBB2_14
.Ltmp174:
# %bb.39:                               # %L.B0266.L.B0266_crit_edge.1
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	16(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$2, %eax
	cmoval	%eax, %edx
	cmpq	$3, %rsi
	je	.LBB2_14
.Ltmp175:
# %bb.40:                               # %L.B0266.L.B0266_crit_edge.2
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	24(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$3, %eax
	cmoval	%eax, %edx
	cmpq	$5, %rsi
	jb	.LBB2_14
.Ltmp176:
# %bb.41:                               # %L.B0266.L.B0266_crit_edge.3
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	32(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$4, %eax
	cmoval	%eax, %edx
.Ltmp177:
.LBB2_14:                               # %L.B0267
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	movslq	-48(%rbp), %rdi         # 4-byte Folded Reload
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	movl	%edx, %eax
	vmovsd	(%rcx,%rax,8), %xmm0    # xmm0 = mem[0],zero
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	cmpl	$2, %esi
	setb	%al
	vxorpd	%xmm8, %xmm8, %xmm8
	jae	.LBB2_16
.Ltmp178:
# %bb.15:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movb	$1, %r9b
	vmovapd	-112(%rbp), %xmm9       # 16-byte Reload
	movq	-120(%rbp), %r8         # 8-byte Reload
	jmp	.LBB2_21
.Ltmp179:
.LBB2_11:                               # %L.B0269.thread
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	vmovsd	(%rcx), %xmm0           # xmm0 = mem[0],zero
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovapd	%xmm0, %xmm1
	movq	-64(%rbp), %r13         # 8-byte Reload
	vmovapd	-112(%rbp), %xmm9       # 16-byte Reload
	vmovapd	-208(%rbp), %xmm6       # 16-byte Reload
	jmp	.LBB2_36
.Ltmp180:
.LBB2_16:                               # %L.B0280
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movb	%al, -224(%rbp)         # 1-byte Spill
	leaq	-728(%rbp), %rax
	movq	%rdi, -112(%rbp)        # 8-byte Spill
	leaq	(%rax,%rdi,8), %rax
	movq	%rax, -48(%rbp)         # 8-byte Spill
	leaq	-1(%rsi), %r11
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	decl	%edx
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vmovsd	-864(%rbp), %xmm1       # xmm1 = mem[0],zero
	.loc	1 198 1 is_stmt 1       # backgroundfield/quadr.cpp:198:1
	movq	%r12, %r13
	movabsq	$4294967296, %rax       # imm = 0x100000000
	movabsq	$8589934592, %r9        # imm = 0x200000000
	movabsq	$12884901888, %r15      # imm = 0x300000000
	xorl	%edi, %edi
	jmp	.LBB2_17
.Ltmp181:
	.p2align	4, 0x90
.LBB2_19:                               # %L.B0271
                                        #   in Loop: Header=BB2_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	leal	(%rdx,%rdx), %ebx
	addl	$2, %ebx
	xorl	%esi, %esi
	cmpl	%r13d, %ebx
	setge	%sil
	movslq	%edx, %rdx
	leaq	-1056(%rbp,%rdx,8), %rbx
	leaq	-888(%rbp), %rcx
	leaq	(%rcx,%rdx,8), %r8
	cmovgeq	%rbx, %r8
	subl	%esi, %edx
	vmovsd	(%r8), %xmm9            # xmm9 = mem[0],zero
	vaddsd	%xmm9, %xmm0, %xmm0
	incq	%rdi
	movabsq	$4294967296, %rsi       # imm = 0x100000000
	addq	%rsi, %r15
	addq	%rsi, %r9
	addq	%rsi, %rax
	cmpq	%r11, %rdi
	jae	.LBB2_20
.Ltmp182:
.LBB2_17:                               # %L.B0268
                                        #   Parent Loop BB2_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	movq	%r13, %rbx
	decq	%r13
	cmpl	$2, %ebx
	jl	.LBB2_19
.Ltmp183:
# %bb.18:                               # %L.B0270
                                        #   in Loop: Header=BB2_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	(%r10), %xmm2           # xmm2 = mem[0],zero
	vmovsd	-888(%rbp), %xmm3       # xmm3 = mem[0],zero
	vsubsd	-1056(%rbp), %xmm3, %xmm3
	movq	-48(%rbp), %rcx         # 8-byte Reload
	vmovsd	(%rcx,%rdi,8), %xmm4    # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm2, %xmm5
	vdivsd	%xmm5, %xmm3, %xmm3
	vmulsd	%xmm3, %xmm4, %xmm4
	vmovsd	%xmm4, -1056(%rbp)
	vmulsd	%xmm3, %xmm2, %xmm2
	vmovsd	%xmm2, -896(%rbp)
	cmpq	$1, %r13
	jle	.LBB2_19
.Ltmp184:
# %bb.42:                               # %L.B0270.1
                                        #   in Loop: Header=BB2_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	-880(%rbp), %xmm2       # xmm2 = mem[0],zero
	vsubsd	-1048(%rbp), %xmm2, %xmm2
	movq	%r10, %rcx
	vmovsd	8(%r10), %xmm3          # xmm3 = mem[0],zero
	movq	%rax, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm3, %xmm5
	vdivsd	%xmm5, %xmm2, %xmm2
	vmulsd	%xmm2, %xmm4, %xmm4
	vmovsd	%xmm4, -1048(%rbp)
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovsd	%xmm2, -888(%rbp)
	cmpq	$3, %rbx
	je	.LBB2_19
.Ltmp185:
# %bb.43:                               # %L.B0270.2
                                        #   in Loop: Header=BB2_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	-872(%rbp), %xmm2       # xmm2 = mem[0],zero
	vsubsd	-1040(%rbp), %xmm2, %xmm2
	movq	%r10, %rcx
	vmovsd	16(%r10), %xmm3         # xmm3 = mem[0],zero
	movq	%r9, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm3, %xmm5
	vdivsd	%xmm5, %xmm2, %xmm2
	vmulsd	%xmm2, %xmm4, %xmm4
	vmovsd	%xmm4, -1040(%rbp)
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovsd	%xmm2, -880(%rbp)
	cmpq	$4, %r13
	jl	.LBB2_19
.Ltmp186:
# %bb.44:                               # %L.B0270.3
                                        #   in Loop: Header=BB2_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r10, %rcx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	24(%r10), %xmm2         # xmm2 = mem[0],zero
	vsubsd	-1032(%rbp), %xmm1, %xmm3
	movq	%r15, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm2, %xmm5
	vdivsd	%xmm5, %xmm3, %xmm3
	vmulsd	%xmm3, %xmm4, %xmm4
	vmovsd	%xmm4, -1032(%rbp)
	vmulsd	%xmm3, %xmm2, %xmm2
	vmovsd	%xmm2, -872(%rbp)
	jmp	.LBB2_19
.Ltmp187:
.LBB2_20:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-144(%rbp), %r15        # 8-byte Reload
	movq	-120(%rbp), %r8         # 8-byte Reload
	movq	-80(%rbp), %rcx         # 8-byte Reload
	movq	-72(%rbp), %rsi         # 8-byte Reload
	movq	-112(%rbp), %rdi        # 8-byte Reload
	movb	-224(%rbp), %r9b        # 1-byte Reload
.Ltmp188:
.LBB2_21:                               # %L.B0269
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	xorl	%edx, %edx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	testl	%esi, %esi
	je	.LBB2_29
.Ltmp189:
# %bb.22:                               # %L.B0283
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	(%r10), %xmm1           # xmm1 = mem[0],zero
	vucomisd	%xmm8, %xmm1
	jne	.LBB2_25
	jnp	.LBB2_23
.Ltmp190:
.LBB2_25:                               # %L..inline.11527.lr.ph
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	movq	(%rcx), %rax
	movq	%rax, -592(%rbp)
	vmovq	%rax, %xmm2
	vaddsd	.LCPI2_3(%rip), %xmm2, %xmm2
	vmovsd	%xmm2, -432(%rbp)
	cmpq	$1, %rsi
	jbe	.LBB2_29
.Ltmp191:
# %bb.26:                               # %L..inline.11527.L..inline.11524_crit_edge
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	8(%r10), %xmm2          # xmm2 = mem[0],zero
	vucomisd	%xmm8, %xmm2
	jne	.LBB2_45
	jnp	.LBB2_27
.Ltmp192:
.LBB2_45:                               # %L..inline.11527.1
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vandpd	%xmm6, %xmm1, %xmm1
	vminsd	%xmm1, %xmm1, %xmm1
	vandpd	%xmm6, %xmm2, %xmm2
	xorl	%edx, %edx
	vucomisd	%xmm2, %xmm1
	seta	%dl
	movq	8(%rcx), %rax
	movq	%rax, -584(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI2_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -424(%rbp)
	cmpq	$3, %rsi
	jb	.LBB2_29
.Ltmp193:
# %bb.46:                               # %L..inline.11527.L..inline.11524_crit_edge.1
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	16(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB2_48
	jnp	.LBB2_47
.Ltmp194:
.LBB2_48:                               # %L..inline.11527.2
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vminsd	%xmm1, %xmm2, %xmm1
	vandpd	%xmm6, %xmm3, %xmm2
	vucomisd	%xmm2, %xmm1
	movl	$2, %eax
	cmoval	%eax, %edx
	movq	16(%rcx), %rax
	movq	%rax, -576(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI2_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -416(%rbp)
	cmpq	$4, %rsi
	jb	.LBB2_29
.Ltmp195:
# %bb.49:                               # %L..inline.11527.L..inline.11524_crit_edge.2
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	24(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB2_51
	jnp	.LBB2_50
.Ltmp196:
.LBB2_51:                               # %L..inline.11527.3
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vminsd	%xmm1, %xmm2, %xmm1
	vandpd	%xmm6, %xmm3, %xmm2
	vucomisd	%xmm2, %xmm1
	movl	$3, %eax
	cmoval	%eax, %edx
	movq	24(%rcx), %rax
	movq	%rax, -568(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI2_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -408(%rbp)
	cmpq	$5, %rsi
	jb	.LBB2_29
.Ltmp197:
# %bb.52:                               # %L..inline.11527.L..inline.11524_crit_edge.3
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	32(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB2_28
	jnp	.LBB2_53
.Ltmp198:
.LBB2_28:                               # %L..inline.11527.4
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vandpd	%xmm6, %xmm3, %xmm3
	vminsd	%xmm1, %xmm2, %xmm1
	vucomisd	%xmm3, %xmm1
	movl	$4, %eax
	cmoval	%eax, %edx
	movq	32(%rcx), %rax
	movq	%rax, -560(%rbp)
	vmovq	%rax, %xmm1
	vaddsd	.LCPI2_3(%rip), %xmm1, %xmm1
	vmovsd	%xmm1, -400(%rbp)
.Ltmp199:
.LBB2_29:                               # %L..inline.11525
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	movl	%edx, %eax
	vmovsd	(%rcx,%rax,8), %xmm1    # xmm1 = mem[0],zero
	testb	%r9b, %r9b
	je	.LBB2_31
.Ltmp200:
# %bb.30:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	leaq	-184(%rbp), %rdi
	movl	-56(%rbp), %edx         # 4-byte Reload
	vmovsd	-88(%rbp), %xmm8        # 8-byte Reload
                                        # xmm8 = mem[0],zero
	vmovapd	-208(%rbp), %xmm6       # 16-byte Reload
	jmp	.LBB2_36
.Ltmp201:
.LBB2_31:                               # %L.B0286
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	leaq	-728(%rbp), %rax
	leaq	(%rax,%rdi,8), %rax
	addq	$24, %rax
	decq	%r12
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	decl	%edx
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	vmovsd	(%r10), %xmm2           # xmm2 = mem[0],zero
	vmovsd	-560(%rbp), %xmm3       # xmm3 = mem[0],zero
	jmp	.LBB2_32
.Ltmp202:
.LBB2_59:                               # %L..inline.11543.3
                                        #   in Loop: Header=BB2_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vsubsd	%xmm5, %xmm3, %xmm5
	vdivsd	%xmm6, %xmm5, %xmm5
	vmulsd	%xmm5, %xmm3, %xmm6
	vmovsd	%xmm6, -408(%rbp)
	vmulsd	%xmm5, %xmm4, %xmm4
	vmovsd	%xmm4, -568(%rbp)
.Ltmp203:
	.p2align	4, 0x90
.LBB2_34:                               # %L..inline.11539
                                        #   in Loop: Header=BB2_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	leal	(%rdx,%rdx), %ecx
	addl	$2, %ecx
	movslq	%ecx, %rcx
	xorl	%esi, %esi
	cmpq	%rcx, %r12
	setle	%sil
	movslq	%edx, %rdx
	leaq	-432(%rbp,%rdx,8), %rcx
	leaq	-584(%rbp), %rdi
	leaq	(%rdi,%rdx,8), %rdi
	cmovleq	%rcx, %rdi
	subl	%esi, %edx
	vmovsd	(%rdi), %xmm6           # xmm6 = mem[0],zero
	vaddsd	%xmm6, %xmm1, %xmm1
	addq	$8, %rax
	decq	%r12
	testl	%r12d, %r12d
	jle	.LBB2_35
.Ltmp204:
.LBB2_32:                               # %L..inline.11538.preheader
                                        #   Parent Loop BB2_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	-584(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-432(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	%xmm2, %xmm6, %xmm5
	vdivsd	-24(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB2_33
	jnp	.LBB2_60
.Ltmp205:
.LBB2_33:                               # %L..inline.11543
                                        #   in Loop: Header=BB2_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -432(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -592(%rbp)
	cmpq	$1, %r12
	jle	.LBB2_34
.Ltmp206:
# %bb.54:                               # %L..inline.11538.1
                                        #   in Loop: Header=BB2_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	-576(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-424(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	8(%r10), %xmm6, %xmm5
	vdivsd	-16(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB2_55
	jnp	.LBB2_60
.Ltmp207:
.LBB2_55:                               # %L..inline.11543.1
                                        #   in Loop: Header=BB2_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -424(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -584(%rbp)
	cmpq	$3, %r12
	jl	.LBB2_34
.Ltmp208:
# %bb.56:                               # %L..inline.11538.2
                                        #   in Loop: Header=BB2_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	-568(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-416(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	16(%r10), %xmm6, %xmm5
	vdivsd	-8(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB2_57
	jnp	.LBB2_60
.Ltmp209:
.LBB2_57:                               # %L..inline.11543.2
                                        #   in Loop: Header=BB2_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -416(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -576(%rbp)
	cmpq	$4, %r12
	jl	.LBB2_34
.Ltmp210:
# %bb.58:                               # %L..inline.11538.3
                                        #   in Loop: Header=BB2_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovsd	-408(%rbp), %xmm5       # xmm5 = mem[0],zero
	vmulsd	24(%r10), %xmm5, %xmm4
	vdivsd	(%rax), %xmm4, %xmm4
	vsubsd	%xmm3, %xmm4, %xmm6
	vucomisd	%xmm8, %xmm6
	jne	.LBB2_59
	jp	.LBB2_59
	jmp	.LBB2_60
.Ltmp211:
.LBB2_23:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	xorl	%eax, %eax
	jmp	.LBB2_24
.Ltmp212:
.LBB2_27:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movl	$1, %eax
	jmp	.LBB2_24
.Ltmp213:
.LBB2_47:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	movl	$2, %eax
	jmp	.LBB2_24
.Ltmp214:
.LBB2_50:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	movl	$3, %eax
	jmp	.LBB2_24
.Ltmp215:
.LBB2_53:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	movl	$4, %eax
.Ltmp216:
.LBB2_24:                               # %L.B0284
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	(%rcx,%rax,8), %xmm1    # xmm1 = mem[0],zero
	vxorpd	%xmm6, %xmm6, %xmm6
.Ltmp217:
.LBB2_35:                               #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	leaq	-184(%rbp), %rdi
	movl	-56(%rbp), %edx         # 4-byte Reload
	vmovsd	-88(%rbp), %xmm8        # 8-byte Reload
                                        # xmm8 = mem[0],zero
.Ltmp218:
.LBB2_36:                               # %L..inline.11530
                                        #   in Loop: Header=BB2_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vmovapd	.LCPI2_2(%rip), %xmm5   # xmm5 = [NaN,NaN]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vandpd	%xmm5, %xmm9, %xmm9
	.loc	1 199 1 is_stmt 1       # backgroundfield/quadr.cpp:199:1
	vandpd	%xmm5, %xmm6, %xmm2
	.loc	1 200 1                 # backgroundfield/quadr.cpp:200:1
	vmaxsd	%xmm9, %xmm2, %xmm3
	vcmpunordsd	%xmm9, %xmm9, %xmm4
	vblendvpd	%xmm4, %xmm2, %xmm3, %xmm3
	vcmpunordsd	%xmm6, %xmm6, %xmm4
	vblendvpd	%xmm4, %xmm9, %xmm3, %xmm3
	.loc	1 201 1                 # backgroundfield/quadr.cpp:201:1
	vsubsd	%xmm1, %xmm0, %xmm4
	vandpd	%xmm5, %xmm4, %xmm5
	.loc	1 202 1                 # backgroundfield/quadr.cpp:202:1
	vcmpunordsd	%xmm3, %xmm3, %xmm6
	vmaxsd	%xmm3, %xmm5, %xmm7
	vblendvpd	%xmm6, %xmm5, %xmm7, %xmm5
	vcmpunordsd	%xmm4, %xmm4, %xmm4
	vblendvpd	%xmm4, %xmm3, %xmm5, %xmm3
	vmovsd	-232(%rbp), %xmm4       # 8-byte Reload
                                        # xmm4 = mem[0],zero
	.loc	1 204 1                 # backgroundfield/quadr.cpp:204:1
	vucomisd	%xmm3, %xmm4
	jbe	.LBB2_3
.Ltmp219:
# %bb.37:                               # %L.B0293
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	vaddsd	%xmm1, %xmm0, %xmm0
	vmulsd	.LCPI2_0(%rip), %xmm0, %xmm8
.Ltmp220:
.LBB2_38:                               # %L..inline.11589
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 232 1                 # backgroundfield/quadr.cpp:232:1
	vmovapd	%xmm8, %xmm0
	addq	$1016, %rsp             # imm = 0x3F8
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp221:
.LBB2_60:                               # %L..inline.11548
	.cfi_def_cfa %rbp, 16
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 232, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: a <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: b <- [DW_OP_constu 256, DW_OP_minus] [$rbp+0]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	movl	$_ZSt4cerr, %edi
	movl	$.S08003, %esi
	movl	$20, %edx
	callq	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l
	movl	$111, %edi
	callq	exit
.Ltmp222:
.Lfunc_end2:
	.size	_Z7RombergRK11T2DFunctionddddd, .Lfunc_end2-_Z7RombergRK11T2DFunctionddddd
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _Z7RombergRK11T3DFunctionddddddd
.LCPI3_0:
	.quad	4602678819172646912     # double 0.5
.LCPI3_3:
	.quad	4158027847206421152     # double 1.0000000000000001E-30
.LCPI3_4:
	.quad	4598175219545276416     # double 0.25
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4
.LCPI3_1:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
.LCPI3_2:
	.quad	9223372036854775807     # double NaN
	.quad	9223372036854775807     # double NaN
	.text
	.globl	_Z7RombergRK11T3DFunctionddddddd
	.p2align	4, 0x90
	.type	_Z7RombergRK11T3DFunctionddddddd,@function
_Z7RombergRK11T3DFunctionddddddd:       # @_Z7RombergRK11T3DFunctionddddddd
.Lfunc_begin3:
	.loc	1 248 0                 # backgroundfield/quadr.cpp:248:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: Romberg:func <- $rdi
	#DEBUG_VALUE: Romberg:a <- $xmm0
	#DEBUG_VALUE: Romberg:b <- $xmm1
	#DEBUG_VALUE: Romberg:c <- $xmm2
	#DEBUG_VALUE: Romberg:d <- $xmm3
	#DEBUG_VALUE: Romberg:e <- $xmm4
	#DEBUG_VALUE: Romberg:f <- $xmm5
	#DEBUG_VALUE: Romberg:absacc <- $xmm6
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$1032, %rsp             # imm = 0x408
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	vmovsd	%xmm5, -272(%rbp)       # 8-byte Spill
	vmovsd	%xmm4, -136(%rbp)       # 8-byte Spill
.Ltmp223:
	#DEBUG_VALUE: absacc <- $xmm6
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: d <- $xmm3
	#DEBUG_VALUE: c <- $xmm2
	#DEBUG_VALUE: b <- $xmm1
	#DEBUG_VALUE: a <- $xmm0
	#DEBUG_VALUE: func <- $rdi
	#DEBUG_VALUE: Romberg:absacc <- $xmm6
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:d <- $xmm3
	#DEBUG_VALUE: Romberg:c <- $xmm2
	#DEBUG_VALUE: Romberg:b <- $xmm1
	#DEBUG_VALUE: Romberg:a <- $xmm0
	#DEBUG_VALUE: Romberg:func <- $rdi
	.loc	1 249 1 prologue_end    # backgroundfield/quadr.cpp:249:1
	vsubsd	%xmm4, %xmm5, %xmm5
	vmovsd	%xmm6, -248(%rbp)       # 8-byte Spill
.Ltmp224:
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	vdivsd	%xmm5, %xmm6, %xmm4
	.loc	1 243 1                 # backgroundfield/quadr.cpp:243:1
	movq	$_ZTV10Tintxy_f3D+16, -216(%rbp)
	movq	%rdi, -208(%rbp)
	vmovsd	%xmm0, -200(%rbp)
	vmovsd	%xmm1, -192(%rbp)
	vmovsd	%xmm2, -184(%rbp)
	vmovsd	%xmm3, -176(%rbp)
	vmovsd	%xmm4, -168(%rbp)
	movabsq	$4607182418800017408, %rax # imm = 0x3FF0000000000000
	.loc	1 184 1                 # backgroundfield/quadr.cpp:184:1
	movq	%rax, -752(%rbp)
	.loc	1 190 1                 # backgroundfield/quadr.cpp:190:1
	leaq	-672(%rbp), %r15
	leaq	-744(%rbp), %rax
	vmovsd	%xmm5, -128(%rbp)       # 8-byte Spill
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vmulsd	.LCPI3_0(%rip), %xmm5, %xmm0
.Ltmp225:
	vmovsd	%xmm0, -264(%rbp)       # 8-byte Spill
	xorl	%r13d, %r13d
	movl	$1, %r14d
	leaq	-216(%rbp), %rdi
.Ltmp226:
                                        # implicit-def: $xmm2
	.loc	1 188 1 is_stmt 1       # backgroundfield/quadr.cpp:188:1
	vxorpd	%xmm9, %xmm9, %xmm9
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	xorl	%edx, %edx
	jmp	.LBB3_1
.Ltmp227:
	.p2align	4, 0x90
.LBB3_2:                                # %L.B0318
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 191 1 is_stmt 1       # backgroundfield/quadr.cpp:191:1
	movq	-216(%rbp), %rax
	vmovsd	-272(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	movq	%rcx, %r13
	callq	*(%rax)
	vmovsd	%xmm0, -48(%rbp)        # 8-byte Spill
	movq	-216(%rbp), %rax
	leaq	-216(%rbp), %rdi
	vmovsd	-136(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	callq	*(%rax)
	vmovapd	-160(%rbp), %xmm2       # 16-byte Reload
	vmovapd	-112(%rbp), %xmm9       # 16-byte Reload
	leaq	-216(%rbp), %rdi
	vaddsd	-48(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmulsd	-264(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmovsd	%xmm0, -688(%rbp,%r14,8)
	leaq	-8(%r15), %r8
	vmovsd	-8(%r15), %xmm8         # xmm8 = mem[0],zero
	movl	$1, %edx
.Ltmp228:
.LBB3_3:                                # %L..inline.12205
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 210 1                 # backgroundfield/quadr.cpp:210:1
	movq	(%r8), %rax
	movq	%rax, (%r15)
	vmovsd	.LCPI3_4(%rip), %xmm0   # xmm0 = mem[0],zero
	movq	-288(%rbp), %rax        # 8-byte Reload
	vmulsd	-8(%rax), %xmm0, %xmm0
	vmovsd	%xmm0, (%rax)
	.loc	1 212 1                 # backgroundfield/quadr.cpp:212:1
	incq	%r14
	.loc	1 190 1                 # backgroundfield/quadr.cpp:190:1
	addq	$8, %r15
	addq	$8, %rax
	.loc	1 212 1                 # backgroundfield/quadr.cpp:212:1
	cmpq	$8, %r13
	je	.LBB3_38
.Ltmp229:
.LBB3_1:                                # %L..inline.12178
                                        # =>This Loop Header: Depth=1
                                        #     Child Loop BB3_7 Depth 2
                                        #     Child Loop BB3_17 Depth 2
                                        #     Child Loop BB3_32 Depth 2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%rax, -288(%rbp)        # 8-byte Spill
	.loc	1 191 1 is_stmt 1       # backgroundfield/quadr.cpp:191:1
	cmpq	$5, %r14
	movl	$5, %r12d
	cmovbq	%r14, %r12
	leaq	1(%r13), %rcx
	cmpl	$5, %ecx
	movl	%ecx, %ebx
	movl	$5, %eax
	cmovgeq	%rax, %rbx
	cmpq	$1, %r14
	vmovapd	%xmm9, -112(%rbp)       # 16-byte Spill
	vmovapd	%xmm2, -160(%rbp)       # 16-byte Spill
	je	.LBB3_2
.Ltmp230:
# %bb.4:                                # %L..inline.12192
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r14, -280(%rbp)        # 8-byte Spill
	movq	%rcx, -64(%rbp)         # 8-byte Spill
	movq	%r15, -144(%rbp)        # 8-byte Spill
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vcvtsi2sd	%edx, %xmm10, %xmm1
	testl	%edx, %edx
	vmovsd	%xmm1, -56(%rbp)        # 8-byte Spill
	jle	.LBB3_5
.Ltmp231:
# %bb.6:                                # %L.B0321
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovsd	-128(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vdivsd	%xmm1, %xmm0, %xmm1
	vmovsd	-136(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	%xmm1, -72(%rbp)        # 8-byte Spill
	vfmadd231sd	.LCPI3_0(%rip), %xmm1, %xmm0 # xmm0 = (xmm1 * mem) + xmm0
	vxorpd	%xmm1, %xmm1, %xmm1
	movl	%edx, %r14d
	movl	%edx, %r15d
.Ltmp232:
	.p2align	4, 0x90
.LBB3_7:                                # %L..inline.12202
                                        #   Parent Loop BB3_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovsd	%xmm0, -80(%rbp)        # 8-byte Spill
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	movq	-216(%rbp), %rax
	callq	*(%rax)
	leaq	-216(%rbp), %rdi
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	%xmm0, %xmm1, %xmm1
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	vmovsd	-80(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	-72(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	decl	%r15d
	jg	.LBB3_7
	jmp	.LBB3_8
.Ltmp233:
.LBB3_5:                                #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movl	%edx, %r14d
	vxorpd	%xmm1, %xmm1, %xmm1
.Ltmp234:
.LBB3_8:                                # %L..inline.12194
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vmulsd	-128(%rbp), %xmm1, %xmm0 # 8-byte Folded Reload
	vdivsd	-56(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	movq	-280(%rbp), %rax        # 8-byte Reload
	vaddsd	-688(%rbp,%rax,8), %xmm0, %xmm0
	vmulsd	.LCPI3_0(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, -688(%rbp,%rax,8)
	movl	%r14d, %edx
	movq	%rax, %r14
	addl	%edx, %edx
	movq	-144(%rbp), %r15        # 8-byte Reload
	leaq	-8(%r15), %r8
	vmovsd	-8(%r15), %xmm8         # xmm8 = mem[0],zero
	.loc	1 192 1 is_stmt 1       # backgroundfield/quadr.cpp:192:1
	cmpq	$2, %rax
	jae	.LBB3_10
.Ltmp235:
# %bb.9:                                #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	vmovapd	-112(%rbp), %xmm9       # 16-byte Reload
	vmovapd	-160(%rbp), %xmm2       # 16-byte Reload
	jmp	.LBB3_3
.Ltmp236:
	.p2align	4, 0x90
.LBB3_10:                               # %L.B0323
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 194 1 is_stmt 1       # backgroundfield/quadr.cpp:194:1
	cmpq	$5, %r14
	movl	$5, %esi
	cmovbq	%r14, %rsi
	.loc	1 197 1                 # backgroundfield/quadr.cpp:197:1
	movl	%r14d, %eax
	subl	%esi, %eax
	cltq
	leaq	-680(,%rax,8), %rcx
	addq	%rbp, %rcx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	testl	%esi, %esi
	je	.LBB3_11
.Ltmp237:
# %bb.12:                               # %L.B0324
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r8, -120(%rbp)         # 8-byte Spill
	vmovsd	%xmm8, -88(%rbp)        # 8-byte Spill
	movl	%edx, -56(%rbp)         # 4-byte Spill
	movl	%r14d, %edx
	subl	%r12d, %edx
	movl	%edx, -48(%rbp)         # 4-byte Spill
	leaq	-672(,%r13,8), %r13
	addq	%rbp, %r13
	shlq	$3, %rbx
	subq	%rbx, %r13
	leaq	-752(,%rax,8), %rax
	addq	%rbp, %rax
	movq	%rax, -256(%rbp)        # 8-byte Spill
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	(%rax), %xmm0           # xmm0 = mem[0],zero
	vmovapd	%xmm0, -240(%rbp)       # 16-byte Spill
	leaq	-1072(%rbp), %rdi
	movq	%rsi, -72(%rbp)         # 8-byte Spill
	movq	%r13, %rsi
	movq	%rbx, %rdx
	movq	%rcx, -80(%rbp)         # 8-byte Spill
	callq	memcpy
	leaq	-912(%rbp), %rdi
	movq	%r13, %rsi
	movq	%rbx, %rdx
	callq	memcpy
	movq	-72(%rbp), %rsi         # 8-byte Reload
	movq	-80(%rbp), %rcx         # 8-byte Reload
	movq	-256(%rbp), %r10        # 8-byte Reload
	xorl	%edx, %edx
	cmpq	$2, %rsi
	vmovapd	.LCPI3_2(%rip), %xmm6   # xmm6 = [NaN,NaN]
	jb	.LBB3_14
.Ltmp238:
# %bb.13:                               # %L.B0312.L.B0312_crit_edge.lr.ph
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovapd	-240(%rbp), %xmm0       # 16-byte Reload
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vxorpd	.LCPI3_1(%rip), %xmm0, %xmm0
	vmovddup	%xmm0, %xmm0    # xmm0 = xmm0[0,0]
	vandpd	%xmm6, %xmm0, %xmm0
	vpermilpd	$1, %xmm0, %xmm1 # xmm1 = xmm0[1,0]
	vminsd	%xmm1, %xmm0, %xmm0
	vmovsd	8(%r10), %xmm1          # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	xorl	%edx, %edx
	vucomisd	%xmm1, %xmm0
	seta	%dl
	cmpq	$2, %rsi
	jbe	.LBB3_14
.Ltmp239:
# %bb.39:                               # %L.B0312.L.B0312_crit_edge.1
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	16(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$2, %eax
	cmoval	%eax, %edx
	cmpq	$3, %rsi
	je	.LBB3_14
.Ltmp240:
# %bb.40:                               # %L.B0312.L.B0312_crit_edge.2
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	24(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$3, %eax
	cmoval	%eax, %edx
	cmpq	$5, %rsi
	jb	.LBB3_14
.Ltmp241:
# %bb.41:                               # %L.B0312.L.B0312_crit_edge.3
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	32(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$4, %eax
	cmoval	%eax, %edx
.Ltmp242:
.LBB3_14:                               # %L.B0313
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	movslq	-48(%rbp), %rdi         # 4-byte Folded Reload
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	movl	%edx, %eax
	vmovsd	(%rcx,%rax,8), %xmm0    # xmm0 = mem[0],zero
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	cmpl	$2, %esi
	setb	%al
	vxorpd	%xmm8, %xmm8, %xmm8
	jae	.LBB3_16
.Ltmp243:
# %bb.15:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movb	$1, %r9b
	vmovapd	-112(%rbp), %xmm9       # 16-byte Reload
	movq	-120(%rbp), %r8         # 8-byte Reload
	jmp	.LBB3_21
.Ltmp244:
.LBB3_11:                               # %L.B0315.thread
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	vmovsd	(%rcx), %xmm0           # xmm0 = mem[0],zero
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovapd	%xmm0, %xmm1
	movq	-64(%rbp), %r13         # 8-byte Reload
	vmovapd	-112(%rbp), %xmm9       # 16-byte Reload
	vmovapd	-160(%rbp), %xmm6       # 16-byte Reload
	jmp	.LBB3_36
.Ltmp245:
.LBB3_16:                               # %L.B0326
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movb	%al, -240(%rbp)         # 1-byte Spill
	leaq	-744(%rbp), %rax
	movq	%rdi, -112(%rbp)        # 8-byte Spill
	leaq	(%rax,%rdi,8), %rax
	movq	%rax, -48(%rbp)         # 8-byte Spill
	leaq	-1(%rsi), %r11
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	decl	%edx
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vmovsd	-880(%rbp), %xmm1       # xmm1 = mem[0],zero
	.loc	1 198 1 is_stmt 1       # backgroundfield/quadr.cpp:198:1
	movq	%r12, %r13
	movabsq	$4294967296, %rax       # imm = 0x100000000
	movabsq	$8589934592, %r9        # imm = 0x200000000
	movabsq	$12884901888, %r15      # imm = 0x300000000
	xorl	%edi, %edi
	jmp	.LBB3_17
.Ltmp246:
	.p2align	4, 0x90
.LBB3_19:                               # %L.B0317
                                        #   in Loop: Header=BB3_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	leal	(%rdx,%rdx), %ebx
	addl	$2, %ebx
	xorl	%esi, %esi
	cmpl	%r13d, %ebx
	setge	%sil
	movslq	%edx, %rdx
	leaq	-1072(%rbp,%rdx,8), %rbx
	leaq	-904(%rbp), %rcx
	leaq	(%rcx,%rdx,8), %r8
	cmovgeq	%rbx, %r8
	subl	%esi, %edx
	vmovsd	(%r8), %xmm9            # xmm9 = mem[0],zero
	vaddsd	%xmm9, %xmm0, %xmm0
	incq	%rdi
	movabsq	$4294967296, %rsi       # imm = 0x100000000
	addq	%rsi, %r15
	addq	%rsi, %r9
	addq	%rsi, %rax
	cmpq	%r11, %rdi
	jae	.LBB3_20
.Ltmp247:
.LBB3_17:                               # %L.B0314
                                        #   Parent Loop BB3_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	movq	%r13, %rbx
	decq	%r13
	cmpl	$2, %ebx
	jl	.LBB3_19
.Ltmp248:
# %bb.18:                               # %L.B0316
                                        #   in Loop: Header=BB3_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	(%r10), %xmm2           # xmm2 = mem[0],zero
	vmovsd	-904(%rbp), %xmm3       # xmm3 = mem[0],zero
	vsubsd	-1072(%rbp), %xmm3, %xmm3
	movq	-48(%rbp), %rcx         # 8-byte Reload
	vmovsd	(%rcx,%rdi,8), %xmm4    # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm2, %xmm5
	vdivsd	%xmm5, %xmm3, %xmm3
	vmulsd	%xmm3, %xmm4, %xmm4
	vmovsd	%xmm4, -1072(%rbp)
	vmulsd	%xmm3, %xmm2, %xmm2
	vmovsd	%xmm2, -912(%rbp)
	cmpq	$1, %r13
	jle	.LBB3_19
.Ltmp249:
# %bb.42:                               # %L.B0316.1
                                        #   in Loop: Header=BB3_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	-896(%rbp), %xmm2       # xmm2 = mem[0],zero
	vsubsd	-1064(%rbp), %xmm2, %xmm2
	movq	%r10, %rcx
	vmovsd	8(%r10), %xmm3          # xmm3 = mem[0],zero
	movq	%rax, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm3, %xmm5
	vdivsd	%xmm5, %xmm2, %xmm2
	vmulsd	%xmm2, %xmm4, %xmm4
	vmovsd	%xmm4, -1064(%rbp)
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovsd	%xmm2, -904(%rbp)
	cmpq	$3, %rbx
	je	.LBB3_19
.Ltmp250:
# %bb.43:                               # %L.B0316.2
                                        #   in Loop: Header=BB3_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	-888(%rbp), %xmm2       # xmm2 = mem[0],zero
	vsubsd	-1056(%rbp), %xmm2, %xmm2
	movq	%r10, %rcx
	vmovsd	16(%r10), %xmm3         # xmm3 = mem[0],zero
	movq	%r9, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm3, %xmm5
	vdivsd	%xmm5, %xmm2, %xmm2
	vmulsd	%xmm2, %xmm4, %xmm4
	vmovsd	%xmm4, -1056(%rbp)
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovsd	%xmm2, -896(%rbp)
	cmpq	$4, %r13
	jl	.LBB3_19
.Ltmp251:
# %bb.44:                               # %L.B0316.3
                                        #   in Loop: Header=BB3_17 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r10, %rcx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	24(%r10), %xmm2         # xmm2 = mem[0],zero
	vsubsd	-1048(%rbp), %xmm1, %xmm3
	movq	%r15, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm2, %xmm5
	vdivsd	%xmm5, %xmm3, %xmm3
	vmulsd	%xmm3, %xmm4, %xmm4
	vmovsd	%xmm4, -1048(%rbp)
	vmulsd	%xmm3, %xmm2, %xmm2
	vmovsd	%xmm2, -888(%rbp)
	jmp	.LBB3_19
.Ltmp252:
.LBB3_20:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-144(%rbp), %r15        # 8-byte Reload
	movq	-120(%rbp), %r8         # 8-byte Reload
	movq	-80(%rbp), %rcx         # 8-byte Reload
	movq	-72(%rbp), %rsi         # 8-byte Reload
	movq	-112(%rbp), %rdi        # 8-byte Reload
	movb	-240(%rbp), %r9b        # 1-byte Reload
.Ltmp253:
.LBB3_21:                               # %L.B0315
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	xorl	%edx, %edx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	testl	%esi, %esi
	je	.LBB3_29
.Ltmp254:
# %bb.22:                               # %L.B0329
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	(%r10), %xmm1           # xmm1 = mem[0],zero
	vucomisd	%xmm8, %xmm1
	jne	.LBB3_25
	jnp	.LBB3_23
.Ltmp255:
.LBB3_25:                               # %L..inline.12266.lr.ph
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	movq	(%rcx), %rax
	movq	%rax, -608(%rbp)
	vmovq	%rax, %xmm2
	vaddsd	.LCPI3_3(%rip), %xmm2, %xmm2
	vmovsd	%xmm2, -448(%rbp)
	cmpq	$1, %rsi
	jbe	.LBB3_29
.Ltmp256:
# %bb.26:                               # %L..inline.12266.L..inline.12263_crit_edge
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	8(%r10), %xmm2          # xmm2 = mem[0],zero
	vucomisd	%xmm8, %xmm2
	jne	.LBB3_45
	jnp	.LBB3_27
.Ltmp257:
.LBB3_45:                               # %L..inline.12266.1
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vandpd	%xmm6, %xmm1, %xmm1
	vminsd	%xmm1, %xmm1, %xmm1
	vandpd	%xmm6, %xmm2, %xmm2
	xorl	%edx, %edx
	vucomisd	%xmm2, %xmm1
	seta	%dl
	movq	8(%rcx), %rax
	movq	%rax, -600(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI3_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -440(%rbp)
	cmpq	$3, %rsi
	jb	.LBB3_29
.Ltmp258:
# %bb.46:                               # %L..inline.12266.L..inline.12263_crit_edge.1
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	16(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB3_48
	jnp	.LBB3_47
.Ltmp259:
.LBB3_48:                               # %L..inline.12266.2
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vminsd	%xmm1, %xmm2, %xmm1
	vandpd	%xmm6, %xmm3, %xmm2
	vucomisd	%xmm2, %xmm1
	movl	$2, %eax
	cmoval	%eax, %edx
	movq	16(%rcx), %rax
	movq	%rax, -592(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI3_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -432(%rbp)
	cmpq	$4, %rsi
	jb	.LBB3_29
.Ltmp260:
# %bb.49:                               # %L..inline.12266.L..inline.12263_crit_edge.2
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	24(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB3_51
	jnp	.LBB3_50
.Ltmp261:
.LBB3_51:                               # %L..inline.12266.3
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vminsd	%xmm1, %xmm2, %xmm1
	vandpd	%xmm6, %xmm3, %xmm2
	vucomisd	%xmm2, %xmm1
	movl	$3, %eax
	cmoval	%eax, %edx
	movq	24(%rcx), %rax
	movq	%rax, -584(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI3_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -424(%rbp)
	cmpq	$5, %rsi
	jb	.LBB3_29
.Ltmp262:
# %bb.52:                               # %L..inline.12266.L..inline.12263_crit_edge.3
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	32(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB3_28
	jnp	.LBB3_53
.Ltmp263:
.LBB3_28:                               # %L..inline.12266.4
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vandpd	%xmm6, %xmm3, %xmm3
	vminsd	%xmm1, %xmm2, %xmm1
	vucomisd	%xmm3, %xmm1
	movl	$4, %eax
	cmoval	%eax, %edx
	movq	32(%rcx), %rax
	movq	%rax, -576(%rbp)
	vmovq	%rax, %xmm1
	vaddsd	.LCPI3_3(%rip), %xmm1, %xmm1
	vmovsd	%xmm1, -416(%rbp)
.Ltmp264:
.LBB3_29:                               # %L..inline.12264
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	movl	%edx, %eax
	vmovsd	(%rcx,%rax,8), %xmm1    # xmm1 = mem[0],zero
	testb	%r9b, %r9b
	je	.LBB3_31
.Ltmp265:
# %bb.30:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	leaq	-216(%rbp), %rdi
	movl	-56(%rbp), %edx         # 4-byte Reload
	vmovsd	-88(%rbp), %xmm8        # 8-byte Reload
                                        # xmm8 = mem[0],zero
	vmovapd	-160(%rbp), %xmm6       # 16-byte Reload
	jmp	.LBB3_36
.Ltmp266:
.LBB3_31:                               # %L.B0332
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	leaq	-744(%rbp), %rax
	leaq	(%rax,%rdi,8), %rax
	addq	$24, %rax
	decq	%r12
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	decl	%edx
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	vmovsd	(%r10), %xmm2           # xmm2 = mem[0],zero
	vmovsd	-576(%rbp), %xmm3       # xmm3 = mem[0],zero
	jmp	.LBB3_32
.Ltmp267:
.LBB3_59:                               # %L..inline.12282.3
                                        #   in Loop: Header=BB3_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vsubsd	%xmm5, %xmm3, %xmm5
	vdivsd	%xmm6, %xmm5, %xmm5
	vmulsd	%xmm5, %xmm3, %xmm6
	vmovsd	%xmm6, -424(%rbp)
	vmulsd	%xmm5, %xmm4, %xmm4
	vmovsd	%xmm4, -584(%rbp)
.Ltmp268:
	.p2align	4, 0x90
.LBB3_34:                               # %L..inline.12278
                                        #   in Loop: Header=BB3_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	leal	(%rdx,%rdx), %ecx
	addl	$2, %ecx
	movslq	%ecx, %rcx
	xorl	%esi, %esi
	cmpq	%rcx, %r12
	setle	%sil
	movslq	%edx, %rdx
	leaq	-448(%rbp,%rdx,8), %rcx
	leaq	-600(%rbp), %rdi
	leaq	(%rdi,%rdx,8), %rdi
	cmovleq	%rcx, %rdi
	subl	%esi, %edx
	vmovsd	(%rdi), %xmm6           # xmm6 = mem[0],zero
	vaddsd	%xmm6, %xmm1, %xmm1
	addq	$8, %rax
	decq	%r12
	testl	%r12d, %r12d
	jle	.LBB3_35
.Ltmp269:
.LBB3_32:                               # %L..inline.12277.preheader
                                        #   Parent Loop BB3_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	-600(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-448(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	%xmm2, %xmm6, %xmm5
	vdivsd	-24(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB3_33
	jnp	.LBB3_60
.Ltmp270:
.LBB3_33:                               # %L..inline.12282
                                        #   in Loop: Header=BB3_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -448(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -608(%rbp)
	cmpq	$1, %r12
	jle	.LBB3_34
.Ltmp271:
# %bb.54:                               # %L..inline.12277.1
                                        #   in Loop: Header=BB3_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	-592(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-440(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	8(%r10), %xmm6, %xmm5
	vdivsd	-16(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB3_55
	jnp	.LBB3_60
.Ltmp272:
.LBB3_55:                               # %L..inline.12282.1
                                        #   in Loop: Header=BB3_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -440(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -600(%rbp)
	cmpq	$3, %r12
	jl	.LBB3_34
.Ltmp273:
# %bb.56:                               # %L..inline.12277.2
                                        #   in Loop: Header=BB3_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	-584(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-432(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	16(%r10), %xmm6, %xmm5
	vdivsd	-8(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB3_57
	jnp	.LBB3_60
.Ltmp274:
.LBB3_57:                               # %L..inline.12282.2
                                        #   in Loop: Header=BB3_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -432(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -592(%rbp)
	cmpq	$4, %r12
	jl	.LBB3_34
.Ltmp275:
# %bb.58:                               # %L..inline.12277.3
                                        #   in Loop: Header=BB3_32 Depth=2
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovsd	-424(%rbp), %xmm5       # xmm5 = mem[0],zero
	vmulsd	24(%r10), %xmm5, %xmm4
	vdivsd	(%rax), %xmm4, %xmm4
	vsubsd	%xmm3, %xmm4, %xmm6
	vucomisd	%xmm8, %xmm6
	jne	.LBB3_59
	jp	.LBB3_59
	jmp	.LBB3_60
.Ltmp276:
.LBB3_23:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	xorl	%eax, %eax
	jmp	.LBB3_24
.Ltmp277:
.LBB3_27:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movl	$1, %eax
	jmp	.LBB3_24
.Ltmp278:
.LBB3_47:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	movl	$2, %eax
	jmp	.LBB3_24
.Ltmp279:
.LBB3_50:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	movl	$3, %eax
	jmp	.LBB3_24
.Ltmp280:
.LBB3_53:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	movl	$4, %eax
.Ltmp281:
.LBB3_24:                               # %L.B0330
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	(%rcx,%rax,8), %xmm1    # xmm1 = mem[0],zero
	vxorpd	%xmm6, %xmm6, %xmm6
.Ltmp282:
.LBB3_35:                               #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	leaq	-216(%rbp), %rdi
	movl	-56(%rbp), %edx         # 4-byte Reload
	vmovsd	-88(%rbp), %xmm8        # 8-byte Reload
                                        # xmm8 = mem[0],zero
.Ltmp283:
.LBB3_36:                               # %L..inline.12269
                                        #   in Loop: Header=BB3_1 Depth=1
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vmovapd	.LCPI3_2(%rip), %xmm5   # xmm5 = [NaN,NaN]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vandpd	%xmm5, %xmm9, %xmm9
	.loc	1 199 1 is_stmt 1       # backgroundfield/quadr.cpp:199:1
	vandpd	%xmm5, %xmm6, %xmm2
	.loc	1 200 1                 # backgroundfield/quadr.cpp:200:1
	vmaxsd	%xmm9, %xmm2, %xmm3
	vcmpunordsd	%xmm9, %xmm9, %xmm4
	vblendvpd	%xmm4, %xmm2, %xmm3, %xmm3
	vcmpunordsd	%xmm6, %xmm6, %xmm4
	vblendvpd	%xmm4, %xmm9, %xmm3, %xmm3
	.loc	1 201 1                 # backgroundfield/quadr.cpp:201:1
	vsubsd	%xmm1, %xmm0, %xmm4
	vandpd	%xmm5, %xmm4, %xmm5
	.loc	1 202 1                 # backgroundfield/quadr.cpp:202:1
	vcmpunordsd	%xmm3, %xmm3, %xmm6
	vmaxsd	%xmm3, %xmm5, %xmm7
	vblendvpd	%xmm6, %xmm5, %xmm7, %xmm5
	vcmpunordsd	%xmm4, %xmm4, %xmm4
	vblendvpd	%xmm4, %xmm3, %xmm5, %xmm3
	vmovsd	-248(%rbp), %xmm4       # 8-byte Reload
                                        # xmm4 = mem[0],zero
	.loc	1 204 1                 # backgroundfield/quadr.cpp:204:1
	vucomisd	%xmm3, %xmm4
	jbe	.LBB3_3
.Ltmp284:
# %bb.37:                               # %L.B0339
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	vaddsd	%xmm1, %xmm0, %xmm0
	vmulsd	.LCPI3_0(%rip), %xmm0, %xmm8
.Ltmp285:
.LBB3_38:                               # %L..inline.12328
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 249 1                 # backgroundfield/quadr.cpp:249:1
	vmovapd	%xmm8, %xmm0
	addq	$1032, %rsp             # imm = 0x408
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp286:
.LBB3_60:                               # %L..inline.12287
	.cfi_def_cfa %rbp, 16
	#DEBUG_VALUE: Romberg:absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: absacc <- [DW_OP_constu 248, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: Romberg:f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: e <- [DW_OP_constu 136, DW_OP_minus] [$rbp+0]
	#DEBUG_VALUE: f <- [DW_OP_constu 272, DW_OP_minus] [$rbp+0]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	movl	$_ZSt4cerr, %edi
	movl	$.S08003, %esi
	movl	$20, %edx
	callq	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l
	movl	$111, %edi
	callq	exit
.Ltmp287:
.Lfunc_end3:
	.size	_Z7RombergRK11T3DFunctionddddddd, .Lfunc_end3-_Z7RombergRK11T3DFunctionddddddd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T1DFunctionD1Ev    # -- Begin function _ZN11T1DFunctionD1Ev
	.p2align	4, 0x90
	.type	_ZN11T1DFunctionD1Ev,@function
_ZN11T1DFunctionD1Ev:                   # @_ZN11T1DFunctionD1Ev
.Lfunc_begin4:
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
.Ltmp288:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 29 1 prologue_end     # backgroundfield/functions.hpp:29:1
	movq	$_ZTV11T1DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp289:
.Lfunc_end4:
	.size	_ZN11T1DFunctionD1Ev, .Lfunc_end4-_ZN11T1DFunctionD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T1DFunctionD0Ev    # -- Begin function _ZN11T1DFunctionD0Ev
	.p2align	4, 0x90
	.type	_ZN11T1DFunctionD0Ev,@function
_ZN11T1DFunctionD0Ev:                   # @_ZN11T1DFunctionD0Ev
.Lfunc_begin5:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp290:
	.loc	1 29 1 prologue_end     # backgroundfield/quadr.cpp:29:1
	movl	$8, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp291:
.Lfunc_end5:
	.size	_ZN11T1DFunctionD0Ev, .Lfunc_end5-_ZN11T1DFunctionD0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T2DFunctionD1Ev    # -- Begin function _ZN11T2DFunctionD1Ev
	.p2align	4, 0x90
	.type	_ZN11T2DFunctionD1Ev,@function
_ZN11T2DFunctionD1Ev:                   # @_ZN11T2DFunctionD1Ev
.Lfunc_begin6:
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
.Ltmp292:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 30 1 prologue_end     # backgroundfield/functions.hpp:30:1
	movq	$_ZTV11T2DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp293:
.Lfunc_end6:
	.size	_ZN11T2DFunctionD1Ev, .Lfunc_end6-_ZN11T2DFunctionD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN11T2DFunctionD0Ev    # -- Begin function _ZN11T2DFunctionD0Ev
	.p2align	4, 0x90
	.type	_ZN11T2DFunctionD0Ev,@function
_ZN11T2DFunctionD0Ev:                   # @_ZN11T2DFunctionD0Ev
.Lfunc_begin7:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp294:
	.loc	1 30 1 prologue_end     # backgroundfield/quadr.cpp:30:1
	movl	$8, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp295:
.Lfunc_end7:
	.size	_ZN11T2DFunctionD0Ev, .Lfunc_end7-_ZN11T2DFunctionD0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK8T2D_fix14callEd    # -- Begin function _ZNK8T2D_fix14callEd
	.p2align	4, 0x90
	.type	_ZNK8T2D_fix14callEd,@function
_ZNK8T2D_fix14callEd:                   # @_ZNK8T2D_fix14callEd
.Lfunc_begin8:
	.loc	2 41 0                  # backgroundfield/functions.hpp:41:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:y <- $xmm0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp296:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: y <- $xmm0
	vmovaps	%xmm0, %xmm1
.Ltmp297:
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	.loc	2 41 1 prologue_end     # backgroundfield/functions.hpp:41:1
	movq	8(%rdi), %rax
	vmovsd	16(%rdi), %xmm0         # xmm0 = mem[0],zero
	movq	(%rax), %rcx
	movq	(%rcx), %rcx
	movq	%rax, %rdi
.Ltmp298:
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmpq	*%rcx                   # TAILCALL
.Ltmp299:
.Lfunc_end8:
	.size	_ZNK8T2D_fix14callEd, .Lfunc_end8-_ZNK8T2D_fix14callEd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T2D_fix1D1Ev        # -- Begin function _ZN8T2D_fix1D1Ev
	.p2align	4, 0x90
	.type	_ZN8T2D_fix1D1Ev,@function
_ZN8T2D_fix1D1Ev:                       # @_ZN8T2D_fix1D1Ev
.Lfunc_begin9:
	.loc	2 42 0                  # backgroundfield/functions.hpp:42:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~T2D_fix1: <- $rdi
	#DEBUG_VALUE: ~T2D_fix1: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp300:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 42 1 prologue_end     # backgroundfield/functions.hpp:42:1
	movq	$_ZTV11T1DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp301:
.Lfunc_end9:
	.size	_ZN8T2D_fix1D1Ev, .Lfunc_end9-_ZN8T2D_fix1D1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T2D_fix1D0Ev        # -- Begin function _ZN8T2D_fix1D0Ev
	.p2align	4, 0x90
	.type	_ZN8T2D_fix1D0Ev,@function
_ZN8T2D_fix1D0Ev:                       # @_ZN8T2D_fix1D0Ev
.Lfunc_begin10:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp302:
	.loc	1 42 1 prologue_end     # backgroundfield/quadr.cpp:42:1
	movl	$24, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp303:
.Lfunc_end10:
	.size	_ZN8T2D_fix1D0Ev, .Lfunc_end10-_ZN8T2D_fix1D0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK8T3D_fix34callEdd   # -- Begin function _ZNK8T3D_fix34callEdd
	.p2align	4, 0x90
	.type	_ZNK8T3D_fix34callEdd,@function
_ZNK8T3D_fix34callEdd:                  # @_ZNK8T3D_fix34callEdd
.Lfunc_begin11:
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
.Ltmp304:
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
.Ltmp305:
	#DEBUG_VALUE: y <- $xmm1
	#DEBUG_VALUE: call:y <- $xmm1
	#DEBUG_VALUE: x <- $xmm0
	#DEBUG_VALUE: call:x <- $xmm0
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmpq	*%rcx                   # TAILCALL
.Ltmp306:
.Lfunc_end11:
	.size	_ZNK8T3D_fix34callEdd, .Lfunc_end11-_ZNK8T3D_fix34callEdd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T3D_fix3D1Ev        # -- Begin function _ZN8T3D_fix3D1Ev
	.p2align	4, 0x90
	.type	_ZN8T3D_fix3D1Ev,@function
_ZN8T3D_fix3D1Ev:                       # @_ZN8T3D_fix3D1Ev
.Lfunc_begin12:
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
.Ltmp307:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	2 84 1 prologue_end     # backgroundfield/functions.hpp:84:1
	movq	$_ZTV11T2DFunction+16, (%rdi)
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp308:
.Lfunc_end12:
	.size	_ZN8T3D_fix3D1Ev, .Lfunc_end12-_ZN8T3D_fix3D1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN8T3D_fix3D0Ev        # -- Begin function _ZN8T3D_fix3D0Ev
	.p2align	4, 0x90
	.type	_ZN8T3D_fix3D0Ev,@function
_ZN8T3D_fix3D0Ev:                       # @_ZN8T3D_fix3D0Ev
.Lfunc_begin13:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp309:
	.loc	1 84 1 prologue_end     # backgroundfield/quadr.cpp:84:1
	movl	$24, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp310:
.Lfunc_end13:
	.size	_ZN8T3D_fix3D0Ev, .Lfunc_end13-_ZN8T3D_fix3D0Ev
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               # -- Begin function _ZNK9Tinty_f2D4callEd
.LCPI14_0:
	.quad	4602678819172646912     # double 0.5
.LCPI14_3:
	.quad	4158027847206421152     # double 1.0000000000000001E-30
.LCPI14_4:
	.quad	4598175219545276416     # double 0.25
	.section	.rodata.cst16,"aM",@progbits,16
	.p2align	4
.LCPI14_1:
	.quad	-9223372036854775808    # double -0
	.quad	-9223372036854775808    # double -0
.LCPI14_2:
	.quad	9223372036854775807     # double NaN
	.quad	9223372036854775807     # double NaN
	.text
	.weak	_ZNK9Tinty_f2D4callEd
	.p2align	4, 0x90
	.type	_ZNK9Tinty_f2D4callEd,@function
_ZNK9Tinty_f2D4callEd:                  # @_ZNK9Tinty_f2D4callEd
.Lfunc_begin14:
	.loc	1 226 0                 # backgroundfield/quadr.cpp:226:0
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
.Ltmp311:
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$984, %rsp              # imm = 0x3D8
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
.Ltmp312:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: x <- $xmm0
	.loc	1 226 1 prologue_end    # backgroundfield/quadr.cpp:226:1
	movq	8(%rdi), %rax
	.loc	1 40 1                  # backgroundfield/quadr.cpp:40:1
	movq	$_ZTV8T2D_fix1+16, -104(%rbp)
	movq	%rax, -96(%rbp)
	vmovsd	%xmm0, -88(%rbp)
	.loc	1 226 1                 # backgroundfield/quadr.cpp:226:1
	vmovsd	16(%rdi), %xmm0         # xmm0 = mem[0],zero
.Ltmp313:
	vmovsd	24(%rdi), %xmm1         # xmm1 = mem[0],zero
	vmovsd	32(%rdi), %xmm2         # xmm2 = mem[0],zero
	vmovsd	%xmm2, -200(%rbp)       # 8-byte Spill
	movabsq	$4607182418800017408, %rax # imm = 0x3FF0000000000000
	.loc	1 184 1                 # backgroundfield/quadr.cpp:184:1
	movq	%rax, -704(%rbp)
	.loc	1 190 1                 # backgroundfield/quadr.cpp:190:1
	leaq	-624(%rbp), %r15
	leaq	-696(%rbp), %rax
	vmovsd	%xmm0, -152(%rbp)       # 8-byte Spill
	vmovsd	%xmm1, -224(%rbp)       # 8-byte Spill
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vsubsd	%xmm0, %xmm1, %xmm0
	vmovsd	%xmm0, -144(%rbp)       # 8-byte Spill
	vmulsd	.LCPI14_0(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, -216(%rbp)       # 8-byte Spill
	xorl	%r13d, %r13d
	movl	$1, %r14d
	leaq	-104(%rbp), %rdi
.Ltmp314:
                                        # implicit-def: $xmm2
	.loc	1 188 1 is_stmt 1       # backgroundfield/quadr.cpp:188:1
	vxorpd	%xmm9, %xmm9, %xmm9
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	xorl	%edx, %edx
	jmp	.LBB14_1
	.p2align	4, 0x90
.LBB14_2:                               # %L.B0400
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 191 1 is_stmt 1       # backgroundfield/quadr.cpp:191:1
	movq	-104(%rbp), %rax
	vmovsd	-224(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	movq	%rcx, %r13
	callq	*(%rax)
	vmovsd	%xmm0, -48(%rbp)        # 8-byte Spill
	movq	-104(%rbp), %rax
	leaq	-104(%rbp), %rdi
	vmovsd	-152(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	callq	*(%rax)
	vmovapd	-176(%rbp), %xmm2       # 16-byte Reload
	vmovapd	-128(%rbp), %xmm9       # 16-byte Reload
	leaq	-104(%rbp), %rdi
	vaddsd	-48(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmulsd	-216(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	vmovsd	%xmm0, -640(%rbp,%r14,8)
	leaq	-8(%r15), %r8
	vmovsd	-8(%r15), %xmm8         # xmm8 = mem[0],zero
	movl	$1, %edx
.LBB14_3:                               # %L..inline.13280
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 210 1                 # backgroundfield/quadr.cpp:210:1
	movq	(%r8), %rax
	movq	%rax, (%r15)
	vmovsd	.LCPI14_4(%rip), %xmm0  # xmm0 = mem[0],zero
	movq	-240(%rbp), %rax        # 8-byte Reload
	vmulsd	-8(%rax), %xmm0, %xmm0
	vmovsd	%xmm0, (%rax)
	.loc	1 212 1                 # backgroundfield/quadr.cpp:212:1
	incq	%r14
	.loc	1 190 1                 # backgroundfield/quadr.cpp:190:1
	addq	$8, %r15
	addq	$8, %rax
	.loc	1 212 1                 # backgroundfield/quadr.cpp:212:1
	cmpq	$8, %r13
	je	.LBB14_38
.LBB14_1:                               # %L..inline.13253
                                        # =>This Loop Header: Depth=1
                                        #     Child Loop BB14_7 Depth 2
                                        #     Child Loop BB14_17 Depth 2
                                        #     Child Loop BB14_32 Depth 2
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%rax, -240(%rbp)        # 8-byte Spill
	.loc	1 191 1 is_stmt 1       # backgroundfield/quadr.cpp:191:1
	cmpq	$5, %r14
	movl	$5, %r12d
	cmovbq	%r14, %r12
	leaq	1(%r13), %rcx
	cmpl	$5, %ecx
	movl	%ecx, %ebx
	movl	$5, %eax
	cmovgeq	%rax, %rbx
	cmpq	$1, %r14
	vmovapd	%xmm9, -128(%rbp)       # 16-byte Spill
	vmovapd	%xmm2, -176(%rbp)       # 16-byte Spill
	je	.LBB14_2
# %bb.4:                                # %L..inline.13267
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r14, -232(%rbp)        # 8-byte Spill
	movq	%rcx, -64(%rbp)         # 8-byte Spill
	movq	%r15, -160(%rbp)        # 8-byte Spill
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vcvtsi2sd	%edx, %xmm10, %xmm1
	testl	%edx, %edx
	vmovsd	%xmm1, -56(%rbp)        # 8-byte Spill
	jle	.LBB14_5
# %bb.6:                                # %L.B0403
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovsd	-144(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vdivsd	%xmm1, %xmm0, %xmm1
	vmovsd	-152(%rbp), %xmm0       # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	%xmm1, -72(%rbp)        # 8-byte Spill
	vfmadd231sd	.LCPI14_0(%rip), %xmm1, %xmm0 # xmm0 = (xmm1 * mem) + xmm0
	vxorpd	%xmm1, %xmm1, %xmm1
	movl	%edx, %r14d
	movl	%edx, %r15d
	.p2align	4, 0x90
.LBB14_7:                               # %L..inline.13277
                                        #   Parent Loop BB14_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovsd	%xmm0, -80(%rbp)        # 8-byte Spill
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	movq	-104(%rbp), %rax
	callq	*(%rax)
	leaq	-104(%rbp), %rdi
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	%xmm0, %xmm1, %xmm1
	vmovsd	%xmm1, -48(%rbp)        # 8-byte Spill
	vmovsd	-80(%rbp), %xmm0        # 8-byte Reload
                                        # xmm0 = mem[0],zero
	vmovsd	-48(%rbp), %xmm1        # 8-byte Reload
                                        # xmm1 = mem[0],zero
	vaddsd	-72(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	decl	%r15d
	jg	.LBB14_7
	jmp	.LBB14_8
.LBB14_5:                               #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movl	%edx, %r14d
	vxorpd	%xmm1, %xmm1, %xmm1
.LBB14_8:                               # %L..inline.13269
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 191 1                 # backgroundfield/quadr.cpp:191:1
	vmulsd	-144(%rbp), %xmm1, %xmm0 # 8-byte Folded Reload
	vdivsd	-56(%rbp), %xmm0, %xmm0 # 8-byte Folded Reload
	movq	-232(%rbp), %rax        # 8-byte Reload
	vaddsd	-640(%rbp,%rax,8), %xmm0, %xmm0
	vmulsd	.LCPI14_0(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, -640(%rbp,%rax,8)
	movl	%r14d, %edx
	movq	%rax, %r14
	addl	%edx, %edx
	movq	-160(%rbp), %r15        # 8-byte Reload
	leaq	-8(%r15), %r8
	vmovsd	-8(%r15), %xmm8         # xmm8 = mem[0],zero
	.loc	1 192 1 is_stmt 1       # backgroundfield/quadr.cpp:192:1
	cmpq	$2, %rax
	jae	.LBB14_10
# %bb.9:                                #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	vmovapd	-128(%rbp), %xmm9       # 16-byte Reload
	vmovapd	-176(%rbp), %xmm2       # 16-byte Reload
	jmp	.LBB14_3
	.p2align	4, 0x90
.LBB14_10:                              # %L.B0405
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 194 1 is_stmt 1       # backgroundfield/quadr.cpp:194:1
	cmpq	$5, %r14
	movl	$5, %esi
	cmovbq	%r14, %rsi
	.loc	1 197 1                 # backgroundfield/quadr.cpp:197:1
	movl	%r14d, %eax
	subl	%esi, %eax
	cltq
	leaq	-632(,%rax,8), %rcx
	addq	%rbp, %rcx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	testl	%esi, %esi
	je	.LBB14_11
# %bb.12:                               # %L.B0406
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r8, -136(%rbp)         # 8-byte Spill
	vmovsd	%xmm8, -112(%rbp)       # 8-byte Spill
	movl	%edx, -56(%rbp)         # 4-byte Spill
	movl	%r14d, %edx
	subl	%r12d, %edx
	movl	%edx, -48(%rbp)         # 4-byte Spill
	leaq	-624(,%r13,8), %r13
	addq	%rbp, %r13
	shlq	$3, %rbx
	subq	%rbx, %r13
	leaq	-704(,%rax,8), %rax
	addq	%rbp, %rax
	movq	%rax, -208(%rbp)        # 8-byte Spill
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	(%rax), %xmm0           # xmm0 = mem[0],zero
	vmovapd	%xmm0, -192(%rbp)       # 16-byte Spill
	leaq	-1024(%rbp), %rdi
	movq	%rsi, -72(%rbp)         # 8-byte Spill
	movq	%r13, %rsi
	movq	%rbx, %rdx
	movq	%rcx, -80(%rbp)         # 8-byte Spill
	callq	memcpy
	leaq	-864(%rbp), %rdi
	movq	%r13, %rsi
	movq	%rbx, %rdx
	callq	memcpy
	movq	-72(%rbp), %rsi         # 8-byte Reload
	movq	-80(%rbp), %rcx         # 8-byte Reload
	movq	-208(%rbp), %r10        # 8-byte Reload
	xorl	%edx, %edx
	cmpq	$2, %rsi
	vmovapd	.LCPI14_2(%rip), %xmm6  # xmm6 = [NaN,NaN]
	jb	.LBB14_14
# %bb.13:                               # %L.B0394.L.B0394_crit_edge.lr.ph
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	vmovapd	-192(%rbp), %xmm0       # 16-byte Reload
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vxorpd	.LCPI14_1(%rip), %xmm0, %xmm0
	vmovddup	%xmm0, %xmm0    # xmm0 = xmm0[0,0]
	vandpd	%xmm6, %xmm0, %xmm0
	vpermilpd	$1, %xmm0, %xmm1 # xmm1 = xmm0[1,0]
	vminsd	%xmm1, %xmm0, %xmm0
	vmovsd	8(%r10), %xmm1          # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	xorl	%edx, %edx
	vucomisd	%xmm1, %xmm0
	seta	%dl
	cmpq	$2, %rsi
	jbe	.LBB14_14
# %bb.39:                               # %L.B0394.L.B0394_crit_edge.1
                                        #   in Loop: Header=BB14_1 Depth=1
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	16(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$2, %eax
	cmoval	%eax, %edx
	cmpq	$3, %rsi
	je	.LBB14_14
# %bb.40:                               # %L.B0394.L.B0394_crit_edge.2
                                        #   in Loop: Header=BB14_1 Depth=1
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	24(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$3, %eax
	cmoval	%eax, %edx
	cmpq	$5, %rsi
	jb	.LBB14_14
# %bb.41:                               # %L.B0394.L.B0394_crit_edge.3
                                        #   in Loop: Header=BB14_1 Depth=1
	vminsd	%xmm0, %xmm1, %xmm0
	vmovsd	32(%r10), %xmm1         # xmm1 = mem[0],zero
	vandpd	%xmm6, %xmm1, %xmm1
	vucomisd	%xmm1, %xmm0
	movl	$4, %eax
	cmoval	%eax, %edx
.LBB14_14:                              # %L.B0395
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	movslq	-48(%rbp), %rdi         # 4-byte Folded Reload
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	movl	%edx, %eax
	vmovsd	(%rcx,%rax,8), %xmm0    # xmm0 = mem[0],zero
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	cmpl	$2, %esi
	setb	%al
	vxorpd	%xmm8, %xmm8, %xmm8
	jae	.LBB14_16
# %bb.15:                               #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movb	$1, %r9b
	vmovapd	-128(%rbp), %xmm9       # 16-byte Reload
	movq	-136(%rbp), %r8         # 8-byte Reload
	jmp	.LBB14_21
.LBB14_11:                              # %L.B0397.thread
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	vmovsd	(%rcx), %xmm0           # xmm0 = mem[0],zero
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovapd	%xmm0, %xmm1
	movq	-64(%rbp), %r13         # 8-byte Reload
	vmovapd	-128(%rbp), %xmm9       # 16-byte Reload
	vmovapd	-176(%rbp), %xmm6       # 16-byte Reload
	jmp	.LBB14_36
.LBB14_16:                              # %L.B0408
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movb	%al, -192(%rbp)         # 1-byte Spill
	leaq	-696(%rbp), %rax
	movq	%rdi, -128(%rbp)        # 8-byte Spill
	leaq	(%rax,%rdi,8), %rax
	movq	%rax, -48(%rbp)         # 8-byte Spill
	leaq	-1(%rsi), %r11
	.loc	1 197 1 is_stmt 1       # backgroundfield/quadr.cpp:197:1
	decl	%edx
	.loc	1 0 0 is_stmt 0         # backgroundfield/quadr.cpp:0:0
	vmovsd	-832(%rbp), %xmm1       # xmm1 = mem[0],zero
	.loc	1 198 1 is_stmt 1       # backgroundfield/quadr.cpp:198:1
	movq	%r12, %r13
	movabsq	$4294967296, %rax       # imm = 0x100000000
	movabsq	$8589934592, %r9        # imm = 0x200000000
	movabsq	$12884901888, %r15      # imm = 0x300000000
	xorl	%edi, %edi
	jmp	.LBB14_17
	.p2align	4, 0x90
.LBB14_19:                              # %L.B0399
                                        #   in Loop: Header=BB14_17 Depth=2
	leal	(%rdx,%rdx), %ebx
	addl	$2, %ebx
	xorl	%esi, %esi
	cmpl	%r13d, %ebx
	setge	%sil
	movslq	%edx, %rdx
	leaq	-1024(%rbp,%rdx,8), %rbx
	leaq	-856(%rbp), %rcx
	leaq	(%rcx,%rdx,8), %r8
	cmovgeq	%rbx, %r8
	subl	%esi, %edx
	vmovsd	(%r8), %xmm9            # xmm9 = mem[0],zero
	vaddsd	%xmm9, %xmm0, %xmm0
	incq	%rdi
	movabsq	$4294967296, %rsi       # imm = 0x100000000
	addq	%rsi, %r15
	addq	%rsi, %r9
	addq	%rsi, %rax
	cmpq	%r11, %rdi
	jae	.LBB14_20
.LBB14_17:                              # %L.B0396
                                        #   Parent Loop BB14_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	movq	%r13, %rbx
	decq	%r13
	cmpl	$2, %ebx
	jl	.LBB14_19
# %bb.18:                               # %L.B0398
                                        #   in Loop: Header=BB14_17 Depth=2
	vmovsd	(%r10), %xmm2           # xmm2 = mem[0],zero
	vmovsd	-856(%rbp), %xmm3       # xmm3 = mem[0],zero
	vsubsd	-1024(%rbp), %xmm3, %xmm3
	movq	-48(%rbp), %rcx         # 8-byte Reload
	vmovsd	(%rcx,%rdi,8), %xmm4    # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm2, %xmm5
	vdivsd	%xmm5, %xmm3, %xmm3
	vmulsd	%xmm3, %xmm4, %xmm4
	vmovsd	%xmm4, -1024(%rbp)
	vmulsd	%xmm3, %xmm2, %xmm2
	vmovsd	%xmm2, -864(%rbp)
	cmpq	$1, %r13
	jle	.LBB14_19
# %bb.42:                               # %L.B0398.1
                                        #   in Loop: Header=BB14_17 Depth=2
	vmovsd	-848(%rbp), %xmm2       # xmm2 = mem[0],zero
	vsubsd	-1016(%rbp), %xmm2, %xmm2
	movq	%r10, %rcx
	vmovsd	8(%r10), %xmm3          # xmm3 = mem[0],zero
	movq	%rax, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm3, %xmm5
	vdivsd	%xmm5, %xmm2, %xmm2
	vmulsd	%xmm2, %xmm4, %xmm4
	vmovsd	%xmm4, -1016(%rbp)
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovsd	%xmm2, -856(%rbp)
	cmpq	$3, %rbx
	je	.LBB14_19
# %bb.43:                               # %L.B0398.2
                                        #   in Loop: Header=BB14_17 Depth=2
	vmovsd	-840(%rbp), %xmm2       # xmm2 = mem[0],zero
	vsubsd	-1008(%rbp), %xmm2, %xmm2
	movq	%r10, %rcx
	vmovsd	16(%r10), %xmm3         # xmm3 = mem[0],zero
	movq	%r9, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm3, %xmm5
	vdivsd	%xmm5, %xmm2, %xmm2
	vmulsd	%xmm2, %xmm4, %xmm4
	vmovsd	%xmm4, -1008(%rbp)
	vmulsd	%xmm2, %xmm3, %xmm2
	vmovsd	%xmm2, -848(%rbp)
	cmpq	$4, %r13
	jl	.LBB14_19
# %bb.44:                               # %L.B0398.3
                                        #   in Loop: Header=BB14_17 Depth=2
	.loc	1 0 1 is_stmt 0         # backgroundfield/quadr.cpp:0:1
	movq	%r10, %rcx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	24(%r10), %xmm2         # xmm2 = mem[0],zero
	vsubsd	-1000(%rbp), %xmm1, %xmm3
	movq	%r15, %rsi
	sarq	$29, %rsi
	vmovsd	8(%r10,%rsi), %xmm4     # xmm4 = mem[0],zero
	vsubsd	%xmm4, %xmm2, %xmm5
	vdivsd	%xmm5, %xmm3, %xmm3
	vmulsd	%xmm3, %xmm4, %xmm4
	vmovsd	%xmm4, -1000(%rbp)
	vmulsd	%xmm3, %xmm2, %xmm2
	vmovsd	%xmm2, -840(%rbp)
	jmp	.LBB14_19
.LBB14_20:                              #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-160(%rbp), %r15        # 8-byte Reload
	movq	-136(%rbp), %r8         # 8-byte Reload
	movq	-80(%rbp), %rcx         # 8-byte Reload
	movq	-72(%rbp), %rsi         # 8-byte Reload
	movq	-128(%rbp), %rdi        # 8-byte Reload
	movb	-192(%rbp), %r9b        # 1-byte Reload
.LBB14_21:                              # %L.B0397
                                        #   in Loop: Header=BB14_1 Depth=1
	xorl	%edx, %edx
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	testl	%esi, %esi
	je	.LBB14_29
# %bb.22:                               # %L.B0411
                                        #   in Loop: Header=BB14_1 Depth=1
	vmovsd	(%r10), %xmm1           # xmm1 = mem[0],zero
	vucomisd	%xmm8, %xmm1
	jne	.LBB14_25
	jnp	.LBB14_23
.LBB14_25:                              # %L..inline.13341.lr.ph
                                        #   in Loop: Header=BB14_1 Depth=1
	movq	(%rcx), %rax
	movq	%rax, -560(%rbp)
	vmovq	%rax, %xmm2
	vaddsd	.LCPI14_3(%rip), %xmm2, %xmm2
	vmovsd	%xmm2, -400(%rbp)
	cmpq	$1, %rsi
	jbe	.LBB14_29
# %bb.26:                               # %L..inline.13341.L..inline.13338_crit_edge
                                        #   in Loop: Header=BB14_1 Depth=1
	vmovsd	8(%r10), %xmm2          # xmm2 = mem[0],zero
	vucomisd	%xmm8, %xmm2
	jne	.LBB14_45
	jnp	.LBB14_27
.LBB14_45:                              # %L..inline.13341.1
                                        #   in Loop: Header=BB14_1 Depth=1
	vandpd	%xmm6, %xmm1, %xmm1
	vminsd	%xmm1, %xmm1, %xmm1
	vandpd	%xmm6, %xmm2, %xmm2
	xorl	%edx, %edx
	vucomisd	%xmm2, %xmm1
	seta	%dl
	movq	8(%rcx), %rax
	movq	%rax, -552(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI14_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -392(%rbp)
	cmpq	$3, %rsi
	jb	.LBB14_29
# %bb.46:                               # %L..inline.13341.L..inline.13338_crit_edge.1
                                        #   in Loop: Header=BB14_1 Depth=1
	vmovsd	16(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB14_48
	jnp	.LBB14_47
.LBB14_48:                              # %L..inline.13341.2
                                        #   in Loop: Header=BB14_1 Depth=1
	vminsd	%xmm1, %xmm2, %xmm1
	vandpd	%xmm6, %xmm3, %xmm2
	vucomisd	%xmm2, %xmm1
	movl	$2, %eax
	cmoval	%eax, %edx
	movq	16(%rcx), %rax
	movq	%rax, -544(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI14_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -384(%rbp)
	cmpq	$4, %rsi
	jb	.LBB14_29
# %bb.49:                               # %L..inline.13341.L..inline.13338_crit_edge.2
                                        #   in Loop: Header=BB14_1 Depth=1
	vmovsd	24(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB14_51
	jnp	.LBB14_50
.LBB14_51:                              # %L..inline.13341.3
                                        #   in Loop: Header=BB14_1 Depth=1
	vminsd	%xmm1, %xmm2, %xmm1
	vandpd	%xmm6, %xmm3, %xmm2
	vucomisd	%xmm2, %xmm1
	movl	$3, %eax
	cmoval	%eax, %edx
	movq	24(%rcx), %rax
	movq	%rax, -536(%rbp)
	vmovq	%rax, %xmm3
	vaddsd	.LCPI14_3(%rip), %xmm3, %xmm3
	vmovsd	%xmm3, -376(%rbp)
	cmpq	$5, %rsi
	jb	.LBB14_29
# %bb.52:                               # %L..inline.13341.L..inline.13338_crit_edge.3
                                        #   in Loop: Header=BB14_1 Depth=1
	vmovsd	32(%r10), %xmm3         # xmm3 = mem[0],zero
	vucomisd	%xmm8, %xmm3
	jne	.LBB14_28
	jnp	.LBB14_53
.LBB14_28:                              # %L..inline.13341.4
                                        #   in Loop: Header=BB14_1 Depth=1
	vandpd	%xmm6, %xmm3, %xmm3
	vminsd	%xmm1, %xmm2, %xmm1
	vucomisd	%xmm3, %xmm1
	movl	$4, %eax
	cmoval	%eax, %edx
	movq	32(%rcx), %rax
	movq	%rax, -528(%rbp)
	vmovq	%rax, %xmm1
	vaddsd	.LCPI14_3(%rip), %xmm1, %xmm1
	vmovsd	%xmm1, -368(%rbp)
.LBB14_29:                              # %L..inline.13339
                                        #   in Loop: Header=BB14_1 Depth=1
	movl	%edx, %eax
	vmovsd	(%rcx,%rax,8), %xmm1    # xmm1 = mem[0],zero
	testb	%r9b, %r9b
	je	.LBB14_31
# %bb.30:                               #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	leaq	-104(%rbp), %rdi
	movl	-56(%rbp), %edx         # 4-byte Reload
	vmovsd	-112(%rbp), %xmm8       # 8-byte Reload
                                        # xmm8 = mem[0],zero
	vmovapd	-176(%rbp), %xmm6       # 16-byte Reload
	jmp	.LBB14_36
.LBB14_31:                              # %L.B0414
                                        #   in Loop: Header=BB14_1 Depth=1
	leaq	-696(%rbp), %rax
	leaq	(%rax,%rdi,8), %rax
	addq	$24, %rax
	decq	%r12
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	decl	%edx
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	vmovsd	(%r10), %xmm2           # xmm2 = mem[0],zero
	vmovsd	-528(%rbp), %xmm3       # xmm3 = mem[0],zero
	jmp	.LBB14_32
.LBB14_59:                              # %L..inline.13357.3
                                        #   in Loop: Header=BB14_32 Depth=2
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vsubsd	%xmm5, %xmm3, %xmm5
	vdivsd	%xmm6, %xmm5, %xmm5
	vmulsd	%xmm5, %xmm3, %xmm6
	vmovsd	%xmm6, -376(%rbp)
	vmulsd	%xmm5, %xmm4, %xmm4
	vmovsd	%xmm4, -536(%rbp)
	.p2align	4, 0x90
.LBB14_34:                              # %L..inline.13353
                                        #   in Loop: Header=BB14_32 Depth=2
	leal	(%rdx,%rdx), %ecx
	addl	$2, %ecx
	movslq	%ecx, %rcx
	xorl	%esi, %esi
	cmpq	%rcx, %r12
	setle	%sil
	movslq	%edx, %rdx
	leaq	-400(%rbp,%rdx,8), %rcx
	leaq	-552(%rbp), %rdi
	leaq	(%rdi,%rdx,8), %rdi
	cmovleq	%rcx, %rdi
	subl	%esi, %edx
	vmovsd	(%rdi), %xmm6           # xmm6 = mem[0],zero
	vaddsd	%xmm6, %xmm1, %xmm1
	addq	$8, %rax
	decq	%r12
	testl	%r12d, %r12d
	jle	.LBB14_35
.LBB14_32:                              # %L..inline.13352.preheader
                                        #   Parent Loop BB14_1 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	vmovsd	-552(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-400(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	%xmm2, %xmm6, %xmm5
	vdivsd	-24(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB14_33
	jnp	.LBB14_60
.LBB14_33:                              # %L..inline.13357
                                        #   in Loop: Header=BB14_32 Depth=2
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -400(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -560(%rbp)
	cmpq	$1, %r12
	jle	.LBB14_34
# %bb.54:                               # %L..inline.13352.1
                                        #   in Loop: Header=BB14_32 Depth=2
	vmovsd	-544(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-392(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	8(%r10), %xmm6, %xmm5
	vdivsd	-16(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB14_55
	jnp	.LBB14_60
.LBB14_55:                              # %L..inline.13357.1
                                        #   in Loop: Header=BB14_32 Depth=2
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -392(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -552(%rbp)
	cmpq	$3, %r12
	jl	.LBB14_34
# %bb.56:                               # %L..inline.13352.2
                                        #   in Loop: Header=BB14_32 Depth=2
	vmovsd	-536(%rbp), %xmm4       # xmm4 = mem[0],zero
	vmovsd	-384(%rbp), %xmm6       # xmm6 = mem[0],zero
	vmulsd	16(%r10), %xmm6, %xmm5
	vdivsd	-8(%rax), %xmm5, %xmm5
	vsubsd	%xmm4, %xmm5, %xmm7
	vucomisd	%xmm8, %xmm7
	jne	.LBB14_57
	jnp	.LBB14_60
.LBB14_57:                              # %L..inline.13357.2
                                        #   in Loop: Header=BB14_32 Depth=2
	vsubsd	%xmm6, %xmm4, %xmm6
	vdivsd	%xmm7, %xmm6, %xmm6
	vmulsd	%xmm6, %xmm4, %xmm4
	vmovsd	%xmm4, -384(%rbp)
	vmulsd	%xmm6, %xmm5, %xmm4
	vmovsd	%xmm4, -544(%rbp)
	cmpq	$4, %r12
	jl	.LBB14_34
# %bb.58:                               # %L..inline.13352.3
                                        #   in Loop: Header=BB14_32 Depth=2
	vmovsd	-376(%rbp), %xmm5       # xmm5 = mem[0],zero
	vmulsd	24(%r10), %xmm5, %xmm4
	vdivsd	(%rax), %xmm4, %xmm4
	vsubsd	%xmm3, %xmm4, %xmm6
	vucomisd	%xmm8, %xmm6
	jne	.LBB14_59
	jp	.LBB14_59
	jmp	.LBB14_60
.LBB14_23:                              #   in Loop: Header=BB14_1 Depth=1
	xorl	%eax, %eax
	jmp	.LBB14_24
.LBB14_27:                              #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movl	$1, %eax
	jmp	.LBB14_24
.LBB14_47:                              #   in Loop: Header=BB14_1 Depth=1
	movl	$2, %eax
	jmp	.LBB14_24
.LBB14_50:                              #   in Loop: Header=BB14_1 Depth=1
	movl	$3, %eax
	jmp	.LBB14_24
.LBB14_53:                              #   in Loop: Header=BB14_1 Depth=1
	movl	$4, %eax
.LBB14_24:                              # %L.B0412
                                        #   in Loop: Header=BB14_1 Depth=1
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vmovsd	(%rcx,%rax,8), %xmm1    # xmm1 = mem[0],zero
	vxorpd	%xmm6, %xmm6, %xmm6
.LBB14_35:                              #   in Loop: Header=BB14_1 Depth=1
	.loc	1 0 1                   # backgroundfield/quadr.cpp:0:1
	movq	-64(%rbp), %r13         # 8-byte Reload
	leaq	-104(%rbp), %rdi
	movl	-56(%rbp), %edx         # 4-byte Reload
	vmovsd	-112(%rbp), %xmm8       # 8-byte Reload
                                        # xmm8 = mem[0],zero
.LBB14_36:                              # %L..inline.13344
                                        #   in Loop: Header=BB14_1 Depth=1
	vmovapd	.LCPI14_2(%rip), %xmm5  # xmm5 = [NaN,NaN]
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	vandpd	%xmm5, %xmm9, %xmm9
	.loc	1 199 1 is_stmt 1       # backgroundfield/quadr.cpp:199:1
	vandpd	%xmm5, %xmm6, %xmm2
	.loc	1 200 1                 # backgroundfield/quadr.cpp:200:1
	vmaxsd	%xmm9, %xmm2, %xmm3
	vcmpunordsd	%xmm9, %xmm9, %xmm4
	vblendvpd	%xmm4, %xmm2, %xmm3, %xmm3
	vcmpunordsd	%xmm6, %xmm6, %xmm4
	vblendvpd	%xmm4, %xmm9, %xmm3, %xmm3
	.loc	1 201 1                 # backgroundfield/quadr.cpp:201:1
	vsubsd	%xmm1, %xmm0, %xmm4
	vandpd	%xmm5, %xmm4, %xmm5
	.loc	1 202 1                 # backgroundfield/quadr.cpp:202:1
	vcmpunordsd	%xmm3, %xmm3, %xmm6
	vmaxsd	%xmm3, %xmm5, %xmm7
	vblendvpd	%xmm6, %xmm5, %xmm7, %xmm5
	vcmpunordsd	%xmm4, %xmm4, %xmm4
	vblendvpd	%xmm4, %xmm3, %xmm5, %xmm3
	vmovsd	-200(%rbp), %xmm4       # 8-byte Reload
                                        # xmm4 = mem[0],zero
	.loc	1 204 1                 # backgroundfield/quadr.cpp:204:1
	vucomisd	%xmm3, %xmm4
	jbe	.LBB14_3
# %bb.37:                               # %L.B0421
	vaddsd	%xmm1, %xmm0, %xmm0
	vmulsd	.LCPI14_0(%rip), %xmm0, %xmm8
.LBB14_38:                              # %L..inline.13403
	.loc	1 226 1                 # backgroundfield/quadr.cpp:226:1
	vmovapd	%xmm8, %xmm0
	addq	$984, %rsp              # imm = 0x3D8
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB14_60:                              # %L..inline.13362
	.cfi_def_cfa %rbp, 16
	.loc	1 198 1                 # backgroundfield/quadr.cpp:198:1
	movl	$_ZSt4cerr, %edi
	movl	$.S08003, %esi
	movl	$20, %edx
	callq	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l
	movl	$111, %edi
	callq	exit
.Ltmp315:
.Lfunc_end14:
	.size	_ZNK9Tinty_f2D4callEd, .Lfunc_end14-_ZNK9Tinty_f2D4callEd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN9Tinty_f2DD1Ev       # -- Begin function _ZN9Tinty_f2DD1Ev
	.p2align	4, 0x90
	.type	_ZN9Tinty_f2DD1Ev,@function
_ZN9Tinty_f2DD1Ev:                      # @_ZN9Tinty_f2DD1Ev
.Lfunc_begin15:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~Tinty_f2D: <- $rdi
	#DEBUG_VALUE: ~Tinty_f2D: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp316:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	1 198 1 prologue_end    # backgroundfield/quadr.cpp:198:1
	movq	$_ZTV11T1DFunction+16, (%rdi)
	.loc	1 29 1                  # backgroundfield/quadr.cpp:29:1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp317:
.Lfunc_end15:
	.size	_ZN9Tinty_f2DD1Ev, .Lfunc_end15-_ZN9Tinty_f2DD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN9Tinty_f2DD0Ev       # -- Begin function _ZN9Tinty_f2DD0Ev
	.p2align	4, 0x90
	.type	_ZN9Tinty_f2DD0Ev,@function
_ZN9Tinty_f2DD0Ev:                      # @_ZN9Tinty_f2DD0Ev
.Lfunc_begin16:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp318:
	.loc	1 29 1 prologue_end     # backgroundfield/quadr.cpp:29:1
	movl	$40, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp319:
.Lfunc_end16:
	.size	_ZN9Tinty_f2DD0Ev, .Lfunc_end16-_ZN9Tinty_f2DD0Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZNK10Tintxy_f3D4callEd # -- Begin function _ZNK10Tintxy_f3D4callEd
	.p2align	4, 0x90
	.type	_ZNK10Tintxy_f3D4callEd,@function
_ZNK10Tintxy_f3D4callEd:                # @_ZNK10Tintxy_f3D4callEd
.Lfunc_begin17:
	.loc	1 244 0                 # backgroundfield/quadr.cpp:244:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call: <- $rdi
	#DEBUG_VALUE: call:z <- $xmm0
	#DEBUG_VALUE: call:z <- $xmm0
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	subq	$32, %rsp
.Ltmp320:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE: z <- $xmm0
	.loc	1 244 1 prologue_end    # backgroundfield/quadr.cpp:244:1
	movq	8(%rdi), %rax
	.loc	1 82 1                  # backgroundfield/quadr.cpp:82:1
	movq	$_ZTV8T3D_fix3+16, -24(%rbp)
	movq	%rax, -16(%rbp)
	vmovsd	%xmm0, -8(%rbp)
	.loc	1 244 1                 # backgroundfield/quadr.cpp:244:1
	vmovsd	16(%rdi), %xmm0         # xmm0 = mem[0],zero
.Ltmp321:
	vmovsd	24(%rdi), %xmm1         # xmm1 = mem[0],zero
	vmovsd	32(%rdi), %xmm2         # xmm2 = mem[0],zero
	vmovsd	40(%rdi), %xmm3         # xmm3 = mem[0],zero
	vmovsd	48(%rdi), %xmm4         # xmm4 = mem[0],zero
	leaq	-24(%rbp), %rdi
.Ltmp322:
	callq	_Z7RombergRK11T2DFunctionddddd
	addq	$32, %rsp
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp323:
.Lfunc_end17:
	.size	_ZNK10Tintxy_f3D4callEd, .Lfunc_end17-_ZNK10Tintxy_f3D4callEd
	.cfi_endproc
                                        # -- End function
	.weak	_ZN10Tintxy_f3DD1Ev     # -- Begin function _ZN10Tintxy_f3DD1Ev
	.p2align	4, 0x90
	.type	_ZN10Tintxy_f3DD1Ev,@function
_ZN10Tintxy_f3DD1Ev:                    # @_ZN10Tintxy_f3DD1Ev
.Lfunc_begin18:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	#DEBUG_VALUE: ~Tintxy_f3D: <- $rdi
	#DEBUG_VALUE: ~Tintxy_f3D: <- $rdi
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp324:
	#DEBUG_VALUE:  <- $rdi
	#DEBUG_VALUE:  <- $rdi
	.loc	1 244 1 prologue_end    # backgroundfield/quadr.cpp:244:1
	movq	$_ZTV11T1DFunction+16, (%rdi)
	.loc	1 29 1                  # backgroundfield/quadr.cpp:29:1
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.Ltmp325:
.Lfunc_end18:
	.size	_ZN10Tintxy_f3DD1Ev, .Lfunc_end18-_ZN10Tintxy_f3DD1Ev
	.cfi_endproc
                                        # -- End function
	.weak	_ZN10Tintxy_f3DD0Ev     # -- Begin function _ZN10Tintxy_f3DD0Ev
	.p2align	4, 0x90
	.type	_ZN10Tintxy_f3DD0Ev,@function
_ZN10Tintxy_f3DD0Ev:                    # @_ZN10Tintxy_f3DD0Ev
.Lfunc_begin19:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp326:
	.loc	1 29 1 prologue_end     # backgroundfield/quadr.cpp:29:1
	movl	$56, %esi
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	_ZdlPvm                 # TAILCALL
.Ltmp327:
.Lfunc_end19:
	.size	_ZN10Tintxy_f3DD0Ev, .Lfunc_end19-_ZN10Tintxy_f3DD0Ev
	.cfi_endproc
                                        # -- End function
	.globl	__sti___25_backgroundfield_quadr_cpp_6e0a9365 # -- Begin function __sti___25_backgroundfield_quadr_cpp_6e0a9365
	.p2align	4, 0x90
	.type	__sti___25_backgroundfield_quadr_cpp_6e0a9365,@function
__sti___25_backgroundfield_quadr_cpp_6e0a9365: # @__sti___25_backgroundfield_quadr_cpp_6e0a9365
.Lfunc_begin20:
	.loc	1 0 0                   # backgroundfield/quadr.cpp:0:0
	.cfi_startproc
# %bb.0:                                # %L.entry
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
.Ltmp328:
	.loc	1 74 1 prologue_end     # backgroundfield/quadr.cpp:74:1
	cmpl	$1, __I___25_backgroundfield_quadr_cpp_6e0a9365(%rip)
	jne	.LBB20_2
# %bb.1:                                # %L.B0050
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	retq
.LBB20_2:                               # %L.B0456
	.cfi_def_cfa %rbp, 16
	movl	$1, __I___25_backgroundfield_quadr_cpp_6e0a9365(%rip)
	movl	$_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE, %edi
	callq	_ZNSt8ios_base4InitC1Ev
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	movl	$_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE, %esi
	movl	$__dso_handle, %edx
	popq	%rbp
	.cfi_def_cfa %rsp, 8
	jmp	__cxa_atexit            # TAILCALL
.Ltmp329:
.Lfunc_end20:
	.size	__sti___25_backgroundfield_quadr_cpp_6e0a9365, .Lfunc_end20-__sti___25_backgroundfield_quadr_cpp_6e0a9365
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

	.type	_ZTV8T2D_fix1,@object   # @_ZTV8T2D_fix1
	.weak	_ZTV8T2D_fix1
	.p2align	4
_ZTV8T2D_fix1:
	.quad	0
	.quad	_ZTI8T2D_fix1
	.quad	_ZNK8T2D_fix14callEd
	.quad	_ZN8T2D_fix1D1Ev
	.quad	_ZN8T2D_fix1D0Ev
	.size	_ZTV8T2D_fix1, 40

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

	.type	_ZTV9Tinty_f2D,@object  # @_ZTV9Tinty_f2D
	.weak	_ZTV9Tinty_f2D
	.p2align	4
_ZTV9Tinty_f2D:
	.quad	0
	.quad	_ZTI9Tinty_f2D
	.quad	_ZNK9Tinty_f2D4callEd
	.quad	_ZN9Tinty_f2DD1Ev
	.quad	_ZN9Tinty_f2DD0Ev
	.size	_ZTV9Tinty_f2D, 40

	.type	_ZTV10Tintxy_f3D,@object # @_ZTV10Tintxy_f3D
	.weak	_ZTV10Tintxy_f3D
	.p2align	4
_ZTV10Tintxy_f3D:
	.quad	0
	.quad	_ZTI10Tintxy_f3D
	.quad	_ZNK10Tintxy_f3D4callEd
	.quad	_ZN10Tintxy_f3DD1Ev
	.quad	_ZN10Tintxy_f3DD0Ev
	.size	_ZTV10Tintxy_f3D, 40

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

	.type	_ZTI8T2D_fix1,@object   # @_ZTI8T2D_fix1
	.weak	_ZTI8T2D_fix1
	.p2align	4
_ZTI8T2D_fix1:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS8T2D_fix1
	.quad	_ZTI11T1DFunction
	.size	_ZTI8T2D_fix1, 24

	.type	_ZTI8T3D_fix3,@object   # @_ZTI8T3D_fix3
	.weak	_ZTI8T3D_fix3
	.p2align	4
_ZTI8T3D_fix3:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS8T3D_fix3
	.quad	_ZTI11T2DFunction
	.size	_ZTI8T3D_fix3, 24

	.type	_ZTI9Tinty_f2D,@object  # @_ZTI9Tinty_f2D
	.weak	_ZTI9Tinty_f2D
	.p2align	4
_ZTI9Tinty_f2D:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS9Tinty_f2D
	.quad	_ZTI11T1DFunction
	.size	_ZTI9Tinty_f2D, 24

	.type	_ZTI10Tintxy_f3D,@object # @_ZTI10Tintxy_f3D
	.weak	_ZTI10Tintxy_f3D
	.p2align	4
_ZTI10Tintxy_f3D:
	.quad	_ZTVN10__cxxabiv120__si_class_type_infoE+16
	.quad	_ZTS10Tintxy_f3D
	.quad	_ZTI11T1DFunction
	.size	_ZTI10Tintxy_f3D, 24

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

	.type	_ZTS8T2D_fix1,@object   # @_ZTS8T2D_fix1
	.weak	_ZTS8T2D_fix1
	.p2align	3
_ZTS8T2D_fix1:
	.asciz	"8T2D_fix1"
	.size	_ZTS8T2D_fix1, 10

	.type	_ZTS8T3D_fix3,@object   # @_ZTS8T3D_fix3
	.weak	_ZTS8T3D_fix3
	.p2align	3
_ZTS8T3D_fix3:
	.asciz	"8T3D_fix3"
	.size	_ZTS8T3D_fix3, 10

	.type	_ZTS9Tinty_f2D,@object  # @_ZTS9Tinty_f2D
	.weak	_ZTS9Tinty_f2D
	.p2align	3
_ZTS9Tinty_f2D:
	.asciz	"9Tinty_f2D"
	.size	_ZTS9Tinty_f2D, 11

	.type	_ZTS10Tintxy_f3D,@object # @_ZTS10Tintxy_f3D
	.weak	_ZTS10Tintxy_f3D
	.p2align	3
_ZTS10Tintxy_f3D:
	.asciz	"10Tintxy_f3D"
	.size	_ZTS10Tintxy_f3D, 13

	.type	__I___25_backgroundfield_quadr_cpp_6e0a9365,@object # @__I___25_backgroundfield_quadr_cpp_6e0a9365
	.bss
	.globl	__I___25_backgroundfield_quadr_cpp_6e0a9365
	.p2align	2
__I___25_backgroundfield_quadr_cpp_6e0a9365:
	.long	0                       # 0x0
	.size	__I___25_backgroundfield_quadr_cpp_6e0a9365, 4

	.type	_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE,@object # @_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE
	.local	_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE
	.comm	_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE,1,1
	.type	.S08003,@object         # @.S08003
	.section	.rodata,"a",@progbits
.S08003:
	.asciz	"*** Error in ratint\n"
	.size	.S08003, 21

	.section	.init_array,"aw",@init_array
	.p2align	3
	.quad	__sti___25_backgroundfield_quadr_cpp_6e0a9365
	.section	.debug_str,"MS",@progbits,1
.Linfo_string0:
	.asciz	" NVC++ 21.2-0"         # string offset=0
.Linfo_string1:
	.asciz	"backgroundfield/quadr.cpp" # string offset=14
.Linfo_string2:
	.asciz	"/home/talgat/vlasiator" # string offset=40
.Linfo_string3:
	.asciz	"std"                   # string offset=63
.Linfo_string4:
	.asciz	"_ZSt4cerr"             # string offset=67
.Linfo_string5:
	.asciz	"__vptr"                # string offset=77
.Linfo_string6:
	.asciz	"int"                   # string offset=84
.Linfo_string7:
	.asciz	"__v_St9basic_iosIcSt11char_traitsIcEE" # string offset=88
.Linfo_string8:
	.asciz	"__b_St8ios_base"       # string offset=126
.Linfo_string9:
	.asciz	"_M_precision"          # string offset=142
.Linfo_string10:
	.asciz	"long"                  # string offset=155
.Linfo_string11:
	.asciz	"_M_width"              # string offset=160
.Linfo_string12:
	.asciz	"_M_flags"              # string offset=169
.Linfo_string13:
	.asciz	"_ZSt12_S_boolalpha"    # string offset=178
.Linfo_string14:
	.asciz	"_ZSt6_S_dec"           # string offset=197
.Linfo_string15:
	.asciz	"_ZSt8_S_fixed"         # string offset=209
.Linfo_string16:
	.asciz	"_ZSt6_S_hex"           # string offset=223
.Linfo_string17:
	.asciz	"_ZSt11_S_internal"     # string offset=235
.Linfo_string18:
	.asciz	"_ZSt7_S_left"          # string offset=253
.Linfo_string19:
	.asciz	"_ZSt6_S_oct"           # string offset=266
.Linfo_string20:
	.asciz	"_ZSt8_S_right"         # string offset=278
.Linfo_string21:
	.asciz	"_ZSt13_S_scientific"   # string offset=292
.Linfo_string22:
	.asciz	"_ZSt11_S_showbase"     # string offset=312
.Linfo_string23:
	.asciz	"_ZSt12_S_showpoint"    # string offset=330
.Linfo_string24:
	.asciz	"_ZSt10_S_showpos"      # string offset=349
.Linfo_string25:
	.asciz	"_ZSt9_S_skipws"        # string offset=366
.Linfo_string26:
	.asciz	"_ZSt10_S_unitbuf"      # string offset=381
.Linfo_string27:
	.asciz	"_ZSt12_S_uppercase"    # string offset=398
.Linfo_string28:
	.asciz	"_ZSt14_S_adjustfield"  # string offset=417
.Linfo_string29:
	.asciz	"_ZSt12_S_basefield"    # string offset=438
.Linfo_string30:
	.asciz	"_ZSt13_S_floatfield"   # string offset=457
.Linfo_string31:
	.asciz	"_ZSt19_S_ios_fmtflags_end" # string offset=477
.Linfo_string32:
	.asciz	"_ZSt19_S_ios_fmtflags_max" # string offset=503
.Linfo_string33:
	.asciz	"_ZSt19_S_ios_fmtflags_min" # string offset=529
.Linfo_string34:
	.asciz	"_ZSt13_Ios_Fmtflags"   # string offset=555
.Linfo_string35:
	.asciz	"_M_exception"          # string offset=575
.Linfo_string36:
	.asciz	"_ZSt10_S_goodbit"      # string offset=588
.Linfo_string37:
	.asciz	"_ZSt9_S_badbit"        # string offset=605
.Linfo_string38:
	.asciz	"_ZSt9_S_eofbit"        # string offset=620
.Linfo_string39:
	.asciz	"_ZSt10_S_failbit"      # string offset=635
.Linfo_string40:
	.asciz	"_ZSt18_S_ios_iostate_end" # string offset=652
.Linfo_string41:
	.asciz	"_ZSt18_S_ios_iostate_max" # string offset=677
.Linfo_string42:
	.asciz	"_ZSt18_S_ios_iostate_min" # string offset=702
.Linfo_string43:
	.asciz	"_ZSt12_Ios_Iostate"    # string offset=727
.Linfo_string44:
	.asciz	"_M_streambuf_state"    # string offset=746
.Linfo_string45:
	.asciz	"_M_callbacks"          # string offset=765
.Linfo_string46:
	.asciz	"_M_next"               # string offset=778
.Linfo_string47:
	.asciz	"_M_fn"                 # string offset=786
.Linfo_string48:
	.asciz	"_ZNSt8ios_base11erase_eventE" # string offset=792
.Linfo_string49:
	.asciz	"_ZNSt8ios_base11imbue_eventE" # string offset=821
.Linfo_string50:
	.asciz	"_ZNSt8ios_base13copyfmt_eventE" # string offset=850
.Linfo_string51:
	.asciz	"_ZNSt8ios_base5eventE" # string offset=881
.Linfo_string52:
	.asciz	"_M_index"              # string offset=903
.Linfo_string53:
	.asciz	"_M_refcount"           # string offset=912
.Linfo_string54:
	.asciz	"_ZNSt8ios_base14_Callback_listE" # string offset=924
.Linfo_string55:
	.asciz	"_M_word_zero"          # string offset=956
.Linfo_string56:
	.asciz	"_M_pword"              # string offset=969
.Linfo_string57:
	.asciz	"void"                  # string offset=978
.Linfo_string58:
	.asciz	"_M_iword"              # string offset=983
.Linfo_string59:
	.asciz	"_ZNSt8ios_base6_WordsE" # string offset=992
.Linfo_string60:
	.asciz	"_M_local_word"         # string offset=1015
.Linfo_string61:
	.asciz	"__ARRAY_SIZE_TYPE__"   # string offset=1029
.Linfo_string62:
	.asciz	"_M_word_size"          # string offset=1049
.Linfo_string63:
	.asciz	"_M_word"               # string offset=1062
.Linfo_string64:
	.asciz	"_M_ios_locale"         # string offset=1070
.Linfo_string65:
	.asciz	"_M_impl"               # string offset=1084
.Linfo_string66:
	.asciz	"_M_facets"             # string offset=1092
.Linfo_string67:
	.asciz	"_ZNSt6locale5facetE"   # string offset=1102
.Linfo_string68:
	.asciz	"_M_facets_size"        # string offset=1122
.Linfo_string69:
	.asciz	"unsigned long"         # string offset=1137
.Linfo_string70:
	.asciz	"_M_caches"             # string offset=1151
.Linfo_string71:
	.asciz	"_M_names"              # string offset=1161
.Linfo_string72:
	.asciz	"signed char"           # string offset=1170
.Linfo_string73:
	.asciz	"_ZNSt6locale5_ImplE"   # string offset=1182
.Linfo_string74:
	.asciz	"_ZSt6locale"           # string offset=1202
.Linfo_string75:
	.asciz	"_ZSt8ios_base"         # string offset=1214
.Linfo_string76:
	.asciz	"_M_tie"                # string offset=1228
.Linfo_string77:
	.asciz	"_M_fill"               # string offset=1235
.Linfo_string78:
	.asciz	"_M_fill_init"          # string offset=1243
.Linfo_string79:
	.asciz	"_M_streambuf"          # string offset=1256
.Linfo_string80:
	.asciz	"_M_in_beg"             # string offset=1269
.Linfo_string81:
	.asciz	"_M_in_cur"             # string offset=1279
.Linfo_string82:
	.asciz	"_M_in_end"             # string offset=1289
.Linfo_string83:
	.asciz	"_M_out_beg"            # string offset=1299
.Linfo_string84:
	.asciz	"_M_out_cur"            # string offset=1310
.Linfo_string85:
	.asciz	"_M_out_end"            # string offset=1321
.Linfo_string86:
	.asciz	"_M_buf_locale"         # string offset=1332
.Linfo_string87:
	.asciz	"_ZSt15basic_streambufIcSt11char_traitsIcEE" # string offset=1346
.Linfo_string88:
	.asciz	"_M_ctype"              # string offset=1389
.Linfo_string89:
	.asciz	"__b_NSt6locale5facetE" # string offset=1398
.Linfo_string90:
	.asciz	"__SO__NSt6locale5facetE" # string offset=1420
.Linfo_string91:
	.asciz	"_M_c_locale_ctype"     # string offset=1444
.Linfo_string92:
	.asciz	"__locales"             # string offset=1462
.Linfo_string93:
	.asciz	"__locale_data"         # string offset=1472
.Linfo_string94:
	.asciz	"__ctype_b"             # string offset=1486
.Linfo_string95:
	.asciz	"unsigned short"        # string offset=1496
.Linfo_string96:
	.asciz	"__ctype_tolower"       # string offset=1511
.Linfo_string97:
	.asciz	"__ctype_toupper"       # string offset=1527
.Linfo_string98:
	.asciz	"__names"               # string offset=1543
.Linfo_string99:
	.asciz	"__locale_struct"       # string offset=1551
.Linfo_string100:
	.asciz	"_M_del"                # string offset=1567
.Linfo_string101:
	.asciz	"_M_toupper"            # string offset=1574
.Linfo_string102:
	.asciz	"_M_tolower"            # string offset=1585
.Linfo_string103:
	.asciz	"_M_table"              # string offset=1596
.Linfo_string104:
	.asciz	"_M_widen_ok"           # string offset=1605
.Linfo_string105:
	.asciz	"_M_widen"              # string offset=1617
.Linfo_string106:
	.asciz	"_M_narrow"             # string offset=1626
.Linfo_string107:
	.asciz	"_M_narrow_ok"          # string offset=1636
.Linfo_string108:
	.asciz	"_ZSt5ctypeIcE"         # string offset=1649
.Linfo_string109:
	.asciz	"_M_num_put"            # string offset=1663
.Linfo_string110:
	.asciz	"_ZSt7num_putIcSt19ostreambuf_iteratorIcSt11char_traitsIcEEE" # string offset=1674
.Linfo_string111:
	.asciz	"_M_num_get"            # string offset=1734
.Linfo_string112:
	.asciz	"_ZSt7num_getIcSt19istreambuf_iteratorIcSt11char_traitsIcEEE" # string offset=1745
.Linfo_string113:
	.asciz	"_ZSt9basic_iosIcSt11char_traitsIcEE" # string offset=1805
.Linfo_string114:
	.asciz	"_ZSo"                  # string offset=1841
.Linfo_string115:
	.asciz	"ostream"               # string offset=1846
.Linfo_string116:
	.asciz	"_ZTV11T1DFunction"     # string offset=1854
.Linfo_string117:
	.asciz	"_ZTV9Tinty_f2D"        # string offset=1872
.Linfo_string118:
	.asciz	"_ZTV10Tintxy_f3D"      # string offset=1887
.Linfo_string119:
	.asciz	"_ZTV11T2DFunction"     # string offset=1904
.Linfo_string120:
	.asciz	"_ZTV8T2D_fix1"         # string offset=1922
.Linfo_string121:
	.asciz	"_ZTV8T3D_fix3"         # string offset=1936
.Linfo_string122:
	.asciz	"__I___25_backgroundfield_quadr_cpp_6e0a9365" # string offset=1950
.Linfo_string123:
	.asciz	"_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE" # string offset=1994
.Linfo_string124:
	.asciz	"_ZNSt8ios_base4InitE"  # string offset=2059
.Linfo_string125:
	.asciz	"__dso_handle"          # string offset=2080
.Linfo_string126:
	.asciz	"_ZTI11T1DFunction"     # string offset=2093
.Linfo_string127:
	.asciz	"__class_type_info"     # string offset=2111
.Linfo_string128:
	.asciz	"_ZTI11T2DFunction"     # string offset=2129
.Linfo_string129:
	.asciz	"_ZTI8T2D_fix1"         # string offset=2147
.Linfo_string130:
	.asciz	"__si_class_type_info"  # string offset=2161
.Linfo_string131:
	.asciz	"_ZTI8T3D_fix3"         # string offset=2182
.Linfo_string132:
	.asciz	"_ZTI9Tinty_f2D"        # string offset=2196
.Linfo_string133:
	.asciz	"_ZTI10Tintxy_f3D"      # string offset=2211
.Linfo_string134:
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # string offset=2228
.Linfo_string135:
	.asciz	"_ZTS11T1DFunction"     # string offset=2266
.Linfo_string136:
	.asciz	"_ZTS11T2DFunction"     # string offset=2284
.Linfo_string137:
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # string offset=2302
.Linfo_string138:
	.asciz	"_ZTS8T2D_fix1"         # string offset=2343
.Linfo_string139:
	.asciz	"_ZTS8T3D_fix3"         # string offset=2357
.Linfo_string140:
	.asciz	"_ZTS9Tinty_f2D"        # string offset=2371
.Linfo_string141:
	.asciz	"_ZTS10Tintxy_f3D"      # string offset=2386
.Linfo_string142:
	.asciz	"T1DFunction"           # string offset=2403
.Linfo_string143:
	.asciz	"T2DFunction"           # string offset=2415
.Linfo_string144:
	.asciz	"f"                     # string offset=2427
.Linfo_string145:
	.asciz	"x"                     # string offset=2429
.Linfo_string146:
	.asciz	"double"                # string offset=2431
.Linfo_string147:
	.asciz	"call"                  # string offset=2438
.Linfo_string148:
	.asciz	"T2D_fix1"              # string offset=2443
.Linfo_string149:
	.asciz	"T3DFunction"           # string offset=2452
.Linfo_string150:
	.asciz	"z"                     # string offset=2464
.Linfo_string151:
	.asciz	"T3D_fix3"              # string offset=2466
.Linfo_string152:
	.asciz	"Romberg_simple"        # string offset=2475
.Linfo_string153:
	.asciz	"Romberg"               # string offset=2490
.Linfo_string154:
	.asciz	"~T1DFunction"          # string offset=2498
.Linfo_string155:
	.asciz	"_ZN11T1DFunctionD0Ev"  # string offset=2511
.Linfo_string156:
	.asciz	"~T2DFunction"          # string offset=2532
.Linfo_string157:
	.asciz	"_ZN11T2DFunctionD0Ev"  # string offset=2545
.Linfo_string158:
	.asciz	"~T2D_fix1"             # string offset=2566
.Linfo_string159:
	.asciz	"_ZN8T2D_fix1D0Ev"      # string offset=2576
.Linfo_string160:
	.asciz	"~T3D_fix3"             # string offset=2593
.Linfo_string161:
	.asciz	"_ZN8T3D_fix3D0Ev"      # string offset=2603
.Linfo_string162:
	.asciz	"~Tinty_f2D"            # string offset=2620
.Linfo_string163:
	.asciz	"_ZN9Tinty_f2DD0Ev"     # string offset=2631
.Linfo_string164:
	.asciz	"~Tintxy_f3D"           # string offset=2649
.Linfo_string165:
	.asciz	"_ZN10Tintxy_f3DD0Ev"   # string offset=2661
.Linfo_string166:
	.asciz	"__sti___25_backgroundfield_quadr_cpp_6e0a9365" # string offset=2681
.Linfo_string167:
	.asciz	"H"                     # string offset=2727
.Linfo_string168:
	.asciz	"S"                     # string offset=2729
.Linfo_string169:
	.asciz	"func"                  # string offset=2731
.Linfo_string170:
	.asciz	"a"                     # string offset=2736
.Linfo_string171:
	.asciz	"b"                     # string offset=2738
.Linfo_string172:
	.asciz	"absacc"                # string offset=2740
.Linfo_string173:
	.asciz	"it"                    # string offset=2747
.Linfo_string174:
	.asciz	"j"                     # string offset=2750
.Linfo_string175:
	.asciz	"dresult"               # string offset=2752
.Linfo_string176:
	.asciz	"result"                # string offset=2760
.Linfo_string177:
	.asciz	"dresult_pol"           # string offset=2767
.Linfo_string178:
	.asciz	"dresult_rat"           # string offset=2779
.Linfo_string179:
	.asciz	"k1"                    # string offset=2791
.Linfo_string180:
	.asciz	"result_pol"            # string offset=2794
.Linfo_string181:
	.asciz	"result_rat"            # string offset=2805
.Linfo_string182:
	.asciz	"c"                     # string offset=2816
.Linfo_string183:
	.asciz	"d"                     # string offset=2818
.Linfo_string184:
	.asciz	"e"                     # string offset=2820
.Linfo_string185:
	.asciz	"y"                     # string offset=2822
.Linfo_string186:
	.asciz	"ymin2D"                # string offset=2824
.Linfo_string187:
	.asciz	"ymax2D"                # string offset=2831
.Linfo_string188:
	.asciz	"the_absacc2D"          # string offset=2838
.Linfo_string189:
	.asciz	"Tinty_f2D"             # string offset=2851
.Linfo_string190:
	.asciz	"xmin3D"                # string offset=2861
.Linfo_string191:
	.asciz	"xmax3D"                # string offset=2868
.Linfo_string192:
	.asciz	"ymin3D"                # string offset=2875
.Linfo_string193:
	.asciz	"ymax3D"                # string offset=2882
.Linfo_string194:
	.asciz	"the_absacc3D"          # string offset=2889
.Linfo_string195:
	.asciz	"Tintxy_f3D"            # string offset=2902
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp1-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp47-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	93                      # DW_OP_reg13
	.quad	0
	.quad	0
.Ldebug_loc1:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp4-.Lfunc_begin0
	.quad	.Lfunc_end0-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	168                     # -88
	.byte	127                     # 
	.quad	0
	.quad	0
.Ldebug_loc2:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp3-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp3-.Lfunc_begin0
	.quad	.Lfunc_end0-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	240                     # -144
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc3:
	.quad	.Lfunc_begin0-.Lfunc_begin0
	.quad	.Ltmp1-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Lfunc_end0-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	224                     # -160
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc4:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp6-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	17                      # DW_OP_consts
	.byte	0                       # 0
	.quad	.Ltmp6-.Lfunc_begin0
	.quad	.Ltmp8-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	95                      # super-register DW_OP_reg15
	.quad	.Ltmp8-.Lfunc_begin0
	.quad	.Ltmp9-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	17                      # DW_OP_consts
	.byte	1                       # 1
	.quad	.Ltmp9-.Lfunc_begin0
	.quad	.Ltmp15-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	95                      # super-register DW_OP_reg15
	.quad	.Ltmp15-.Lfunc_begin0
	.quad	.Ltmp16-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	80                      # -48
	.quad	0
	.quad	0
.Ldebug_loc5:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp3-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp3-.Lfunc_begin0
	.quad	.Lfunc_end0-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	240                     # -144
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc6:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp4-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp4-.Lfunc_begin0
	.quad	.Lfunc_end0-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	168                     # -88
	.byte	127                     # 
	.quad	0
	.quad	0
.Ldebug_loc7:
	.quad	.Ltmp2-.Lfunc_begin0
	.quad	.Ltmp6-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	17                      # DW_OP_consts
	.byte	1                       # 1
	.quad	.Ltmp6-.Lfunc_begin0
	.quad	.Ltmp48-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	94                      # DW_OP_reg14
	.quad	0
	.quad	0
.Ldebug_loc8:
	.quad	.Ltmp2-.Lfunc_begin0
	.quad	.Ltmp6-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	48                      # DW_OP_lit0
	.quad	.Ltmp6-.Lfunc_begin0
	.quad	.Ltmp31-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -128
	.byte	127                     # 
	.quad	.Ltmp31-.Lfunc_begin0
	.quad	.Ltmp32-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp32-.Lfunc_begin0
	.quad	.Ltmp38-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -128
	.byte	127                     # 
	.quad	.Ltmp38-.Lfunc_begin0
	.quad	.Ltmp39-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp39-.Lfunc_begin0
	.quad	.Lfunc_end0-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -128
	.byte	127                     # 
	.quad	0
	.quad	0
.Ldebug_loc9:
	.quad	.Ltmp23-.Lfunc_begin0
	.quad	.Ltmp41-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp43-.Lfunc_begin0
	.quad	.Lfunc_end0-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc10:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp51-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp155-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	92                      # DW_OP_reg12
	.quad	.Ltmp156-.Lfunc_begin0
	.quad	.Lfunc_end1-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	92                      # DW_OP_reg12
	.quad	0
	.quad	0
.Ldebug_loc11:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp54-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp54-.Lfunc_begin0
	.quad	.Lfunc_end1-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	240                     # -144
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc12:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp53-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp53-.Lfunc_begin0
	.quad	.Lfunc_end1-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	160                     # -224
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc13:
	.quad	.Lfunc_begin1-.Lfunc_begin0
	.quad	.Ltmp51-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Lfunc_end1-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	184                     # -200
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc14:
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp56-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	17                      # DW_OP_consts
	.byte	0                       # 0
	.quad	.Ltmp56-.Lfunc_begin0
	.quad	.Ltmp58-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	88                      # super-register DW_OP_reg8
	.quad	.Ltmp61-.Lfunc_begin0
	.quad	.Ltmp62-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	17                      # DW_OP_consts
	.byte	1                       # 1
	.quad	.Ltmp65-.Lfunc_begin0
	.quad	.Ltmp72-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	88                      # super-register DW_OP_reg8
	.quad	.Ltmp73-.Lfunc_begin0
	.quad	.Ltmp74-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	88                      # super-register DW_OP_reg8
	.quad	.Ltmp76-.Lfunc_begin0
	.quad	.Ltmp86-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	88                      # super-register DW_OP_reg8
	.quad	.Ltmp86-.Lfunc_begin0
	.quad	.Ltmp96-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	76                      # -52
	.quad	.Ltmp96-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	88                      # super-register DW_OP_reg8
	.quad	.Ltmp101-.Lfunc_begin0
	.quad	.Ltmp112-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	76                      # -52
	.quad	.Ltmp112-.Lfunc_begin0
	.quad	.Ltmp154-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	88                      # super-register DW_OP_reg8
	.quad	.Ltmp156-.Lfunc_begin0
	.quad	.Ltmp157-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	88                      # super-register DW_OP_reg8
	.quad	0
	.quad	0
.Ldebug_loc15:
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp53-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp53-.Lfunc_begin0
	.quad	.Lfunc_end1-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	160                     # -224
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc16:
	.quad	.Ltmp51-.Lfunc_begin0
	.quad	.Ltmp54-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp54-.Lfunc_begin0
	.quad	.Lfunc_end1-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	240                     # -144
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc17:
	.quad	.Ltmp52-.Lfunc_begin0
	.quad	.Ltmp56-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	17                      # DW_OP_consts
	.byte	1                       # 1
	.quad	.Ltmp56-.Lfunc_begin0
	.quad	.Ltmp57-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	82                      # DW_OP_reg2
	.quad	.Ltmp57-.Lfunc_begin0
	.quad	.Ltmp62-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	83                      # DW_OP_reg3
	.quad	.Ltmp63-.Lfunc_begin0
	.quad	.Ltmp64-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	82                      # DW_OP_reg2
	.quad	.Ltmp65-.Lfunc_begin0
	.quad	.Ltmp69-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	82                      # DW_OP_reg2
	.quad	.Ltmp69-.Lfunc_begin0
	.quad	.Ltmp75-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	176                     # -80
	.byte	127                     # 
	.quad	.Ltmp75-.Lfunc_begin0
	.quad	.Ltmp87-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	82                      # DW_OP_reg2
	.quad	.Ltmp98-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	82                      # DW_OP_reg2
	.quad	0
	.quad	0
.Ldebug_loc18:
	.quad	.Ltmp52-.Lfunc_begin0
	.quad	.Ltmp56-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	48                      # DW_OP_lit0
	.quad	.Ltmp56-.Lfunc_begin0
	.quad	.Ltmp60-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -128
	.byte	127                     # 
	.quad	.Ltmp60-.Lfunc_begin0
	.quad	.Ltmp64-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp65-.Lfunc_begin0
	.quad	.Ltmp66-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp66-.Lfunc_begin0
	.quad	.Ltmp80-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -128
	.byte	127                     # 
	.quad	.Ltmp80-.Lfunc_begin0
	.quad	.Ltmp82-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp82-.Lfunc_begin0
	.quad	.Ltmp97-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -128
	.byte	127                     # 
	.quad	.Ltmp97-.Lfunc_begin0
	.quad	.Ltmp98-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp98-.Lfunc_begin0
	.quad	.Ltmp99-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -128
	.byte	127                     # 
	.quad	.Ltmp99-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp101-.Lfunc_begin0
	.quad	.Ltmp102-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -128
	.byte	127                     # 
	.quad	.Ltmp104-.Lfunc_begin0
	.quad	.Ltmp106-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	.Ltmp111-.Lfunc_begin0
	.quad	.Ltmp157-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	106                     # DW_OP_reg26
	.quad	0
	.quad	0
.Ldebug_loc19:
	.quad	.Ltmp56-.Lfunc_begin0
	.quad	.Ltmp59-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	208                     # -176
	.byte	126                     # 
	.quad	.Ltmp59-.Lfunc_begin0
	.quad	.Ltmp64-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp65-.Lfunc_begin0
	.quad	.Ltmp67-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp67-.Lfunc_begin0
	.quad	.Ltmp81-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	208                     # -176
	.byte	126                     # 
	.quad	.Ltmp81-.Lfunc_begin0
	.quad	.Ltmp82-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	.Ltmp82-.Lfunc_begin0
	.quad	.Ltmp100-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	208                     # -176
	.byte	126                     # 
	.quad	.Ltmp100-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp101-.Lfunc_begin0
	.quad	.Ltmp128-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	208                     # -176
	.byte	126                     # 
	.quad	.Ltmp128-.Lfunc_begin0
	.quad	.Ltmp129-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp129-.Lfunc_begin0
	.quad	.Ltmp130-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	208                     # -176
	.byte	126                     # 
	.quad	.Ltmp132-.Lfunc_begin0
	.quad	.Ltmp134-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp141-.Lfunc_begin0
	.quad	.Ltmp148-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	208                     # -176
	.byte	126                     # 
	.quad	.Ltmp150-.Lfunc_begin0
	.quad	.Ltmp152-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp152-.Lfunc_begin0
	.quad	.Ltmp156-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc20:
	.quad	.Ltmp78-.Lfunc_begin0
	.quad	.Ltmp85-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	105                     # DW_OP_reg25
	.quad	.Ltmp85-.Lfunc_begin0
	.quad	.Ltmp98-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	64                      # -64
	.quad	.Ltmp98-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	105                     # DW_OP_reg25
	.quad	.Ltmp101-.Lfunc_begin0
	.quad	.Ltmp127-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	64                      # -64
	.quad	.Ltmp127-.Lfunc_begin0
	.quad	.Ltmp129-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	105                     # DW_OP_reg25
	.quad	.Ltmp129-.Lfunc_begin0
	.quad	.Ltmp149-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	64                      # -64
	.quad	.Ltmp149-.Lfunc_begin0
	.quad	.Ltmp156-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	105                     # DW_OP_reg25
	.quad	.Ltmp156-.Lfunc_begin0
	.quad	.Lfunc_end1-.Lfunc_begin0
	.short	2                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	64                      # -64
	.quad	0
	.quad	0
.Ldebug_loc21:
	.quad	.Ltmp83-.Lfunc_begin0
	.quad	.Ltmp88-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # super-register DW_OP_reg5
	.quad	.Ltmp98-.Lfunc_begin0
	.quad	.Ltmp101-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # super-register DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc22:
	.quad	.Ltmp94-.Lfunc_begin0
	.quad	.Ltmp98-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp101-.Lfunc_begin0
	.quad	.Ltmp150-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp156-.Lfunc_begin0
	.quad	.Ltmp157-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc23:
	.quad	.Ltmp125-.Lfunc_begin0
	.quad	.Ltmp141-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp147-.Lfunc_begin0
	.quad	.Ltmp154-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp156-.Lfunc_begin0
	.quad	.Ltmp157-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc24:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp161-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc25:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp159-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp159-.Lfunc_begin0
	.quad	.Lfunc_end2-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	248                     # -136
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc26:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp159-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	.Ltmp159-.Lfunc_begin0
	.quad	.Lfunc_end2-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	128                     # -256
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc27:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp162-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc28:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp162-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc29:
	.quad	.Lfunc_begin2-.Lfunc_begin0
	.quad	.Ltmp160-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp160-.Lfunc_begin0
	.quad	.Lfunc_end2-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	152                     # -232
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc30:
	.quad	.Ltmp159-.Lfunc_begin0
	.quad	.Ltmp160-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp160-.Lfunc_begin0
	.quad	.Lfunc_end2-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	152                     # -232
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc31:
	.quad	.Ltmp159-.Lfunc_begin0
	.quad	.Ltmp162-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc32:
	.quad	.Ltmp159-.Lfunc_begin0
	.quad	.Ltmp162-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc33:
	.quad	.Ltmp159-.Lfunc_begin0
	.quad	.Ltmp161-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc34:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp226-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc35:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp225-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc36:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp227-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc37:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp227-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc38:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp227-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc39:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp223-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	101                     # DW_OP_reg21
	.quad	.Ltmp223-.Lfunc_begin0
	.quad	.Lfunc_end3-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	248                     # -136
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc40:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp223-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	102                     # DW_OP_reg22
	.quad	.Ltmp223-.Lfunc_begin0
	.quad	.Lfunc_end3-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	240                     # -272
	.byte	125                     # 
	.quad	0
	.quad	0
.Ldebug_loc41:
	.quad	.Lfunc_begin3-.Lfunc_begin0
	.quad	.Ltmp224-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp224-.Lfunc_begin0
	.quad	.Lfunc_end3-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	136                     # -248
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc42:
	.quad	.Ltmp223-.Lfunc_begin0
	.quad	.Ltmp224-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	103                     # DW_OP_reg23
	.quad	.Ltmp224-.Lfunc_begin0
	.quad	.Lfunc_end3-.Lfunc_begin0
	.short	3                       # Loc expr size
	.byte	118                     # DW_OP_breg6
	.byte	136                     # -248
	.byte	126                     # 
	.quad	0
	.quad	0
.Ldebug_loc43:
	.quad	.Ltmp223-.Lfunc_begin0
	.quad	.Ltmp227-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	100                     # DW_OP_reg20
	.quad	0
	.quad	0
.Ldebug_loc44:
	.quad	.Ltmp223-.Lfunc_begin0
	.quad	.Ltmp227-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	99                      # DW_OP_reg19
	.quad	0
	.quad	0
.Ldebug_loc45:
	.quad	.Ltmp223-.Lfunc_begin0
	.quad	.Ltmp227-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc46:
	.quad	.Ltmp223-.Lfunc_begin0
	.quad	.Ltmp225-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc47:
	.quad	.Ltmp223-.Lfunc_begin0
	.quad	.Ltmp226-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc48:
	.quad	.Lfunc_begin8-.Lfunc_begin0
	.quad	.Ltmp298-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc49:
	.quad	.Lfunc_begin8-.Lfunc_begin0
	.quad	.Ltmp297-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp297-.Lfunc_begin0
	.quad	.Ltmp299-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc50:
	.quad	.Ltmp296-.Lfunc_begin0
	.quad	.Ltmp298-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc51:
	.quad	.Ltmp296-.Lfunc_begin0
	.quad	.Ltmp297-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	.Ltmp297-.Lfunc_begin0
	.quad	.Ltmp299-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc52:
	.quad	.Lfunc_begin11-.Lfunc_begin0
	.quad	.Ltmp305-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc53:
	.quad	.Lfunc_begin11-.Lfunc_begin0
	.quad	.Ltmp306-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc54:
	.quad	.Lfunc_begin11-.Lfunc_begin0
	.quad	.Ltmp306-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc55:
	.quad	.Ltmp304-.Lfunc_begin0
	.quad	.Ltmp305-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc56:
	.quad	.Ltmp304-.Lfunc_begin0
	.quad	.Ltmp306-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc57:
	.quad	.Ltmp304-.Lfunc_begin0
	.quad	.Ltmp306-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	98                      # DW_OP_reg18
	.quad	0
	.quad	0
.Ldebug_loc58:
	.quad	.Lfunc_begin14-.Lfunc_begin0
	.quad	.Ltmp314-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc59:
	.quad	.Lfunc_begin14-.Lfunc_begin0
	.quad	.Ltmp313-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc60:
	.quad	.Ltmp312-.Lfunc_begin0
	.quad	.Ltmp314-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc61:
	.quad	.Ltmp312-.Lfunc_begin0
	.quad	.Ltmp313-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc62:
	.quad	.Lfunc_begin17-.Lfunc_begin0
	.quad	.Ltmp322-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc63:
	.quad	.Lfunc_begin17-.Lfunc_begin0
	.quad	.Ltmp321-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	97                      # DW_OP_reg17
	.quad	0
	.quad	0
.Ldebug_loc64:
	.quad	.Ltmp320-.Lfunc_begin0
	.quad	.Ltmp322-.Lfunc_begin0
	.short	1                       # Loc expr size
	.byte	85                      # DW_OP_reg5
	.quad	0
	.quad	0
.Ldebug_loc65:
	.quad	.Ltmp320-.Lfunc_begin0
	.quad	.Ltmp321-.Lfunc_begin0
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
	.byte	57                      # DW_TAG_namespace
	.byte	1                       # DW_CHILDREN_yes
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	3                       # Abbreviation Code
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
	.byte	4                       # Abbreviation Code
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
	.byte	5                       # Abbreviation Code
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
	.byte	6                       # Abbreviation Code
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
	.byte	7                       # Abbreviation Code
	.byte	15                      # DW_TAG_pointer_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	8                       # Abbreviation Code
	.byte	21                      # DW_TAG_subroutine_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	9                       # Abbreviation Code
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
	.byte	10                      # Abbreviation Code
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
	.byte	11                      # Abbreviation Code
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
	.byte	12                      # Abbreviation Code
	.byte	40                      # DW_TAG_enumerator
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	28                      # DW_AT_const_value
	.byte	13                      # DW_FORM_sdata
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	13                      # Abbreviation Code
	.byte	21                      # DW_TAG_subroutine_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	14                      # Abbreviation Code
	.byte	5                       # DW_TAG_formal_parameter
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	15                      # Abbreviation Code
	.byte	59                      # DW_TAG_unspecified_type
	.byte	0                       # DW_CHILDREN_no
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	16                      # Abbreviation Code
	.byte	1                       # DW_TAG_array_type
	.byte	1                       # DW_CHILDREN_yes
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	17                      # Abbreviation Code
	.byte	33                      # DW_TAG_subrange_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	55                      # DW_AT_count
	.byte	11                      # DW_FORM_data1
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	18                      # Abbreviation Code
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
	.byte	19                      # Abbreviation Code
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
	.byte	20                      # Abbreviation Code
	.byte	33                      # DW_TAG_subrange_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	55                      # DW_AT_count
	.byte	5                       # DW_FORM_data2
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	21                      # Abbreviation Code
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
	.byte	26                      # Abbreviation Code
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
	.byte	10                      # DW_FORM_block1
	.byte	3                       # DW_AT_name
	.byte	14                      # DW_FORM_strp
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	29                      # Abbreviation Code
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	85                      # DW_AT_ranges
	.byte	6                       # DW_FORM_data4
	.byte	0                       # EOM(1)
	.byte	0                       # EOM(2)
	.byte	30                      # Abbreviation Code
	.byte	52                      # DW_TAG_variable
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
	.byte	11                      # DW_TAG_lexical_block
	.byte	1                       # DW_CHILDREN_yes
	.byte	17                      # DW_AT_low_pc
	.byte	1                       # DW_FORM_addr
	.byte	18                      # DW_AT_high_pc
	.byte	1                       # DW_FORM_addr
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
	.byte	34                      # Abbreviation Code
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
	.byte	35                      # Abbreviation Code
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
	.byte	36                      # Abbreviation Code
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
	.byte	37                      # Abbreviation Code
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
	.byte	38                      # Abbreviation Code
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
	.byte	39                      # Abbreviation Code
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
	.byte	46                      # Abbreviation Code
	.byte	16                      # DW_TAG_reference_type
	.byte	0                       # DW_CHILDREN_no
	.byte	73                      # DW_AT_type
	.byte	19                      # DW_FORM_ref4
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
	.byte	1                       # Abbrev [1] 0xb:0x11b7 DW_TAG_compile_unit
	.long	.Linfo_string0          # DW_AT_producer
	.short	4                       # DW_AT_language
	.long	.Linfo_string1          # DW_AT_name
	.long	.Lline_table_start0     # DW_AT_stmt_list
	.long	.Linfo_string2          # DW_AT_comp_dir
	.byte	1                       # DW_AT_GNU_pubnames
	.quad	.Lfunc_begin0           # DW_AT_low_pc
	.quad	.Lfunc_end20            # DW_AT_high_pc
	.byte	2                       # Abbrev [2] 0x2f:0x1d DW_TAG_namespace
	.long	.Linfo_string3          # DW_AT_name
	.byte	3                       # Abbrev [3] 0x34:0xa DW_TAG_variable
	.long	.Linfo_string4          # DW_AT_name
	.long	62                      # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	4                       # Abbrev [4] 0x3e:0xd DW_TAG_typedef
	.long	76                      # DW_AT_type
	.long	.Linfo_string115        # DW_AT_name
	.short	272                     # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	0                       # End Of Children Mark
	.byte	5                       # Abbrev [5] 0x4c:0x23 DW_TAG_structure_type
	.long	.Linfo_string114        # DW_AT_name
	.short	272                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x54:0xd DW_TAG_member
	.long	.Linfo_string5          # DW_AT_name
	.long	111                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x61:0xd DW_TAG_member
	.long	.Linfo_string7          # DW_AT_name
	.long	138                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x6f:0x5 DW_TAG_pointer_type
	.long	116                     # DW_AT_type
	.byte	8                       # Abbrev [8] 0x74:0x5 DW_TAG_subroutine_type
	.long	121                     # DW_AT_type
	.byte	7                       # Abbrev [7] 0x79:0x5 DW_TAG_pointer_type
	.long	126                     # DW_AT_type
	.byte	8                       # Abbrev [8] 0x7e:0x5 DW_TAG_subroutine_type
	.long	131                     # DW_AT_type
	.byte	9                       # Abbrev [9] 0x83:0x7 DW_TAG_base_type
	.long	.Linfo_string6          # DW_AT_name
	.byte	5                       # DW_AT_encoding
	.byte	4                       # DW_AT_byte_size
	.byte	5                       # Abbrev [5] 0x8a:0x78 DW_TAG_structure_type
	.long	.Linfo_string113        # DW_AT_name
	.short	264                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x92:0xd DW_TAG_member
	.long	.Linfo_string8          # DW_AT_name
	.long	258                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x9f:0xe DW_TAG_member
	.long	.Linfo_string76         # DW_AT_name
	.long	1003                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\330\001"
	.byte	6                       # Abbrev [6] 0xad:0xe DW_TAG_member
	.long	.Linfo_string77         # DW_AT_name
	.long	996                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\340\001"
	.byte	6                       # Abbrev [6] 0xbb:0xe DW_TAG_member
	.long	.Linfo_string78         # DW_AT_name
	.long	996                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\341\001"
	.byte	6                       # Abbrev [6] 0xc9:0xe DW_TAG_member
	.long	.Linfo_string79         # DW_AT_name
	.long	1008                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\350\001"
	.byte	6                       # Abbrev [6] 0xd7:0xe DW_TAG_member
	.long	.Linfo_string88         # DW_AT_name
	.long	1125                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\360\001"
	.byte	6                       # Abbrev [6] 0xe5:0xe DW_TAG_member
	.long	.Linfo_string109        # DW_AT_name
	.long	1450                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\370\001"
	.byte	6                       # Abbrev [6] 0xf3:0xe DW_TAG_member
	.long	.Linfo_string111        # DW_AT_name
	.long	1476                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\200\002"
	.byte	0                       # End Of Children Mark
	.byte	10                      # Abbrev [10] 0x102:0xa7 DW_TAG_structure_type
	.long	.Linfo_string75         # DW_AT_name
	.byte	216                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x109:0xd DW_TAG_member
	.long	.Linfo_string5          # DW_AT_name
	.long	111                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x116:0xd DW_TAG_member
	.long	.Linfo_string9          # DW_AT_name
	.long	425                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	6                       # Abbrev [6] 0x123:0xd DW_TAG_member
	.long	.Linfo_string11         # DW_AT_name
	.long	425                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	6                       # Abbrev [6] 0x130:0xd DW_TAG_member
	.long	.Linfo_string12         # DW_AT_name
	.long	432                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	6                       # Abbrev [6] 0x13d:0xd DW_TAG_member
	.long	.Linfo_string35         # DW_AT_name
	.long	590                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	28
	.byte	6                       # Abbrev [6] 0x14a:0xd DW_TAG_member
	.long	.Linfo_string44         # DW_AT_name
	.long	590                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	6                       # Abbrev [6] 0x157:0xd DW_TAG_member
	.long	.Linfo_string45         # DW_AT_name
	.long	650                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	6                       # Abbrev [6] 0x164:0xd DW_TAG_member
	.long	.Linfo_string55         # DW_AT_name
	.long	768                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	6                       # Abbrev [6] 0x171:0xd DW_TAG_member
	.long	.Linfo_string60         # DW_AT_name
	.long	812                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	64
	.byte	6                       # Abbrev [6] 0x17e:0xe DW_TAG_member
	.long	.Linfo_string62         # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\300\001"
	.byte	6                       # Abbrev [6] 0x18c:0xe DW_TAG_member
	.long	.Linfo_string63         # DW_AT_name
	.long	831                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\310\001"
	.byte	6                       # Abbrev [6] 0x19a:0xe DW_TAG_member
	.long	.Linfo_string64         # DW_AT_name
	.long	836                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\320\001"
	.byte	0                       # End Of Children Mark
	.byte	9                       # Abbrev [9] 0x1a9:0x7 DW_TAG_base_type
	.long	.Linfo_string10         # DW_AT_name
	.byte	5                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	11                      # Abbrev [11] 0x1b0:0x9e DW_TAG_enumeration_type
	.long	.Linfo_string34         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	4                       # DW_AT_alignment
	.byte	12                      # Abbrev [12] 0x1b7:0x6 DW_TAG_enumerator
	.long	.Linfo_string13         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1bd:0x6 DW_TAG_enumerator
	.long	.Linfo_string14         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1c3:0x6 DW_TAG_enumerator
	.long	.Linfo_string15         # DW_AT_name
	.byte	4                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1c9:0x6 DW_TAG_enumerator
	.long	.Linfo_string16         # DW_AT_name
	.byte	8                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1cf:0x6 DW_TAG_enumerator
	.long	.Linfo_string17         # DW_AT_name
	.byte	16                      # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1d5:0x6 DW_TAG_enumerator
	.long	.Linfo_string18         # DW_AT_name
	.byte	32                      # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1db:0x7 DW_TAG_enumerator
	.long	.Linfo_string19         # DW_AT_name
	.asciz	"\300"                  # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1e2:0x7 DW_TAG_enumerator
	.long	.Linfo_string20         # DW_AT_name
	.ascii	"\200\001"              # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1e9:0x7 DW_TAG_enumerator
	.long	.Linfo_string21         # DW_AT_name
	.ascii	"\200\002"              # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1f0:0x7 DW_TAG_enumerator
	.long	.Linfo_string22         # DW_AT_name
	.ascii	"\200\004"              # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1f7:0x7 DW_TAG_enumerator
	.long	.Linfo_string23         # DW_AT_name
	.ascii	"\200\b"                # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x1fe:0x7 DW_TAG_enumerator
	.long	.Linfo_string24         # DW_AT_name
	.ascii	"\200\020"              # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x205:0x7 DW_TAG_enumerator
	.long	.Linfo_string25         # DW_AT_name
	.ascii	"\200 "                 # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x20c:0x8 DW_TAG_enumerator
	.long	.Linfo_string26         # DW_AT_name
	.asciz	"\200\300"              # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x214:0x8 DW_TAG_enumerator
	.long	.Linfo_string27         # DW_AT_name
	.ascii	"\200\200\001"          # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x21c:0x7 DW_TAG_enumerator
	.long	.Linfo_string28         # DW_AT_name
	.ascii	"\260\001"              # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x223:0x7 DW_TAG_enumerator
	.long	.Linfo_string29         # DW_AT_name
	.asciz	"\312"                  # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x22a:0x7 DW_TAG_enumerator
	.long	.Linfo_string30         # DW_AT_name
	.ascii	"\204\002"              # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x231:0x8 DW_TAG_enumerator
	.long	.Linfo_string31         # DW_AT_name
	.ascii	"\200\200\004"          # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x239:0xa DW_TAG_enumerator
	.long	.Linfo_string32         # DW_AT_name
	.ascii	"\377\377\377\377\007"  # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x243:0xa DW_TAG_enumerator
	.long	.Linfo_string33         # DW_AT_name
	.ascii	"\200\200\200\200x"     # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	11                      # Abbrev [11] 0x24e:0x3c DW_TAG_enumeration_type
	.long	.Linfo_string43         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	4                       # DW_AT_alignment
	.byte	12                      # Abbrev [12] 0x255:0x6 DW_TAG_enumerator
	.long	.Linfo_string36         # DW_AT_name
	.byte	0                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x25b:0x6 DW_TAG_enumerator
	.long	.Linfo_string37         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x261:0x6 DW_TAG_enumerator
	.long	.Linfo_string38         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x267:0x6 DW_TAG_enumerator
	.long	.Linfo_string39         # DW_AT_name
	.byte	4                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x26d:0x8 DW_TAG_enumerator
	.long	.Linfo_string40         # DW_AT_name
	.ascii	"\200\200\004"          # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x275:0xa DW_TAG_enumerator
	.long	.Linfo_string41         # DW_AT_name
	.ascii	"\377\377\377\377\007"  # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x27f:0xa DW_TAG_enumerator
	.long	.Linfo_string42         # DW_AT_name
	.ascii	"\200\200\200\200x"     # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x28a:0x5 DW_TAG_pointer_type
	.long	655                     # DW_AT_type
	.byte	10                      # Abbrev [10] 0x28f:0x3c DW_TAG_structure_type
	.long	.Linfo_string54         # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x296:0xd DW_TAG_member
	.long	.Linfo_string46         # DW_AT_name
	.long	650                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x2a3:0xd DW_TAG_member
	.long	.Linfo_string47         # DW_AT_name
	.long	715                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	6                       # Abbrev [6] 0x2b0:0xd DW_TAG_member
	.long	.Linfo_string52         # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	6                       # Abbrev [6] 0x2bd:0xd DW_TAG_member
	.long	.Linfo_string53         # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	20
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x2cb:0x5 DW_TAG_pointer_type
	.long	720                     # DW_AT_type
	.byte	13                      # Abbrev [13] 0x2d0:0x11 DW_TAG_subroutine_type
	.byte	14                      # Abbrev [14] 0x2d1:0x5 DW_TAG_formal_parameter
	.long	737                     # DW_AT_type
	.byte	14                      # Abbrev [14] 0x2d6:0x5 DW_TAG_formal_parameter
	.long	763                     # DW_AT_type
	.byte	14                      # Abbrev [14] 0x2db:0x5 DW_TAG_formal_parameter
	.long	131                     # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	11                      # Abbrev [11] 0x2e1:0x1a DW_TAG_enumeration_type
	.long	.Linfo_string51         # DW_AT_name
	.byte	4                       # DW_AT_byte_size
	.byte	4                       # DW_AT_alignment
	.byte	12                      # Abbrev [12] 0x2e8:0x6 DW_TAG_enumerator
	.long	.Linfo_string48         # DW_AT_name
	.byte	0                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x2ee:0x6 DW_TAG_enumerator
	.long	.Linfo_string49         # DW_AT_name
	.byte	1                       # DW_AT_const_value
	.byte	12                      # Abbrev [12] 0x2f4:0x6 DW_TAG_enumerator
	.long	.Linfo_string50         # DW_AT_name
	.byte	2                       # DW_AT_const_value
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x2fb:0x5 DW_TAG_pointer_type
	.long	258                     # DW_AT_type
	.byte	10                      # Abbrev [10] 0x300:0x22 DW_TAG_structure_type
	.long	.Linfo_string59         # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x307:0xd DW_TAG_member
	.long	.Linfo_string56         # DW_AT_name
	.long	802                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x314:0xd DW_TAG_member
	.long	.Linfo_string58         # DW_AT_name
	.long	425                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x322:0x5 DW_TAG_pointer_type
	.long	807                     # DW_AT_type
	.byte	15                      # Abbrev [15] 0x327:0x5 DW_TAG_unspecified_type
	.long	.Linfo_string57         # DW_AT_name
	.byte	16                      # Abbrev [16] 0x32c:0xc DW_TAG_array_type
	.long	768                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x331:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	8                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	18                      # Abbrev [18] 0x338:0x7 DW_TAG_base_type
	.long	.Linfo_string61         # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	7                       # DW_AT_encoding
	.byte	7                       # Abbrev [7] 0x33f:0x5 DW_TAG_pointer_type
	.long	768                     # DW_AT_type
	.byte	10                      # Abbrev [10] 0x344:0x15 DW_TAG_structure_type
	.long	.Linfo_string74         # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x34b:0xd DW_TAG_member
	.long	.Linfo_string65         # DW_AT_name
	.long	857                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x359:0x5 DW_TAG_pointer_type
	.long	862                     # DW_AT_type
	.byte	10                      # Abbrev [10] 0x35e:0x49 DW_TAG_structure_type
	.long	.Linfo_string73         # DW_AT_name
	.byte	40                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x365:0xd DW_TAG_member
	.long	.Linfo_string53         # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x372:0xd DW_TAG_member
	.long	.Linfo_string66         # DW_AT_name
	.long	935                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	6                       # Abbrev [6] 0x37f:0xd DW_TAG_member
	.long	.Linfo_string68         # DW_AT_name
	.long	979                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	6                       # Abbrev [6] 0x38c:0xd DW_TAG_member
	.long	.Linfo_string70         # DW_AT_name
	.long	935                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	6                       # Abbrev [6] 0x399:0xd DW_TAG_member
	.long	.Linfo_string71         # DW_AT_name
	.long	986                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x3a7:0x5 DW_TAG_pointer_type
	.long	940                     # DW_AT_type
	.byte	7                       # Abbrev [7] 0x3ac:0x5 DW_TAG_pointer_type
	.long	945                     # DW_AT_type
	.byte	10                      # Abbrev [10] 0x3b1:0x22 DW_TAG_structure_type
	.long	.Linfo_string67         # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x3b8:0xd DW_TAG_member
	.long	.Linfo_string5          # DW_AT_name
	.long	111                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x3c5:0xd DW_TAG_member
	.long	.Linfo_string53         # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	0                       # End Of Children Mark
	.byte	9                       # Abbrev [9] 0x3d3:0x7 DW_TAG_base_type
	.long	.Linfo_string69         # DW_AT_name
	.byte	7                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	7                       # Abbrev [7] 0x3da:0x5 DW_TAG_pointer_type
	.long	991                     # DW_AT_type
	.byte	7                       # Abbrev [7] 0x3df:0x5 DW_TAG_pointer_type
	.long	996                     # DW_AT_type
	.byte	9                       # Abbrev [9] 0x3e4:0x7 DW_TAG_base_type
	.long	.Linfo_string72         # DW_AT_name
	.byte	6                       # DW_AT_encoding
	.byte	1                       # DW_AT_byte_size
	.byte	7                       # Abbrev [7] 0x3eb:0x5 DW_TAG_pointer_type
	.long	76                      # DW_AT_type
	.byte	7                       # Abbrev [7] 0x3f0:0x5 DW_TAG_pointer_type
	.long	1013                    # DW_AT_type
	.byte	10                      # Abbrev [10] 0x3f5:0x70 DW_TAG_structure_type
	.long	.Linfo_string87         # DW_AT_name
	.byte	64                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x3fc:0xd DW_TAG_member
	.long	.Linfo_string5          # DW_AT_name
	.long	111                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x409:0xd DW_TAG_member
	.long	.Linfo_string80         # DW_AT_name
	.long	991                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	6                       # Abbrev [6] 0x416:0xd DW_TAG_member
	.long	.Linfo_string81         # DW_AT_name
	.long	991                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	6                       # Abbrev [6] 0x423:0xd DW_TAG_member
	.long	.Linfo_string82         # DW_AT_name
	.long	991                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	6                       # Abbrev [6] 0x430:0xd DW_TAG_member
	.long	.Linfo_string83         # DW_AT_name
	.long	991                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	6                       # Abbrev [6] 0x43d:0xd DW_TAG_member
	.long	.Linfo_string84         # DW_AT_name
	.long	991                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	6                       # Abbrev [6] 0x44a:0xd DW_TAG_member
	.long	.Linfo_string85         # DW_AT_name
	.long	991                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	6                       # Abbrev [6] 0x457:0xd DW_TAG_member
	.long	.Linfo_string86         # DW_AT_name
	.long	836                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	56
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x465:0x5 DW_TAG_pointer_type
	.long	1130                    # DW_AT_type
	.byte	5                       # Abbrev [5] 0x46a:0x8d DW_TAG_structure_type
	.long	.Linfo_string108        # DW_AT_name
	.short	576                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x472:0xd DW_TAG_member
	.long	.Linfo_string89         # DW_AT_name
	.long	1271                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x47f:0xd DW_TAG_member
	.long	.Linfo_string91         # DW_AT_name
	.long	1305                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	6                       # Abbrev [6] 0x48c:0xd DW_TAG_member
	.long	.Linfo_string100        # DW_AT_name
	.long	996                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	6                       # Abbrev [6] 0x499:0xd DW_TAG_member
	.long	.Linfo_string101        # DW_AT_name
	.long	1420                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	6                       # Abbrev [6] 0x4a6:0xd DW_TAG_member
	.long	.Linfo_string102        # DW_AT_name
	.long	1420                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	6                       # Abbrev [6] 0x4b3:0xd DW_TAG_member
	.long	.Linfo_string103        # DW_AT_name
	.long	1408                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	6                       # Abbrev [6] 0x4c0:0xd DW_TAG_member
	.long	.Linfo_string104        # DW_AT_name
	.long	996                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	56
	.byte	6                       # Abbrev [6] 0x4cd:0xd DW_TAG_member
	.long	.Linfo_string105        # DW_AT_name
	.long	1437                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	57
	.byte	6                       # Abbrev [6] 0x4da:0xe DW_TAG_member
	.long	.Linfo_string106        # DW_AT_name
	.long	1437                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\271\002"
	.byte	6                       # Abbrev [6] 0x4e8:0xe DW_TAG_member
	.long	.Linfo_string107        # DW_AT_name
	.long	996                     # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\271\004"
	.byte	0                       # End Of Children Mark
	.byte	10                      # Abbrev [10] 0x4f7:0x22 DW_TAG_structure_type
	.long	.Linfo_string90         # DW_AT_name
	.byte	12                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x4fe:0xd DW_TAG_member
	.long	.Linfo_string5          # DW_AT_name
	.long	111                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x50b:0xd DW_TAG_member
	.long	.Linfo_string53         # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	4                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x519:0x5 DW_TAG_pointer_type
	.long	1310                    # DW_AT_type
	.byte	10                      # Abbrev [10] 0x51e:0x4a DW_TAG_structure_type
	.long	.Linfo_string99         # DW_AT_name
	.byte	232                     # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x525:0xd DW_TAG_member
	.long	.Linfo_string92         # DW_AT_name
	.long	1384                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	6                       # Abbrev [6] 0x532:0xd DW_TAG_member
	.long	.Linfo_string94         # DW_AT_name
	.long	1408                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	104
	.byte	6                       # Abbrev [6] 0x53f:0xd DW_TAG_member
	.long	.Linfo_string96         # DW_AT_name
	.long	1420                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	112
	.byte	6                       # Abbrev [6] 0x54c:0xd DW_TAG_member
	.long	.Linfo_string97         # DW_AT_name
	.long	1420                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	120
	.byte	6                       # Abbrev [6] 0x559:0xe DW_TAG_member
	.long	.Linfo_string98         # DW_AT_name
	.long	1425                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	3                       # DW_AT_data_member_location
	.byte	35
	.ascii	"\200\001"
	.byte	0                       # End Of Children Mark
	.byte	16                      # Abbrev [16] 0x568:0xc DW_TAG_array_type
	.long	1396                    # DW_AT_type
	.byte	17                      # Abbrev [17] 0x56d:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	13                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x574:0x5 DW_TAG_pointer_type
	.long	1401                    # DW_AT_type
	.byte	19                      # Abbrev [19] 0x579:0x7 DW_TAG_structure_type
	.long	.Linfo_string93         # DW_AT_name
	.byte	0                       # DW_AT_byte_size
	.byte	1                       # DW_AT_alignment
	.byte	7                       # Abbrev [7] 0x580:0x5 DW_TAG_pointer_type
	.long	1413                    # DW_AT_type
	.byte	9                       # Abbrev [9] 0x585:0x7 DW_TAG_base_type
	.long	.Linfo_string95         # DW_AT_name
	.byte	7                       # DW_AT_encoding
	.byte	2                       # DW_AT_byte_size
	.byte	7                       # Abbrev [7] 0x58c:0x5 DW_TAG_pointer_type
	.long	131                     # DW_AT_type
	.byte	16                      # Abbrev [16] 0x591:0xc DW_TAG_array_type
	.long	991                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x596:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	13                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	16                      # Abbrev [16] 0x59d:0xd DW_TAG_array_type
	.long	996                     # DW_AT_type
	.byte	20                      # Abbrev [20] 0x5a2:0x7 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.short	256                     # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x5aa:0x5 DW_TAG_pointer_type
	.long	1455                    # DW_AT_type
	.byte	10                      # Abbrev [10] 0x5af:0x15 DW_TAG_structure_type
	.long	.Linfo_string110        # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x5b6:0xd DW_TAG_member
	.long	.Linfo_string89         # DW_AT_name
	.long	1271                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x5c4:0x5 DW_TAG_pointer_type
	.long	1481                    # DW_AT_type
	.byte	10                      # Abbrev [10] 0x5c9:0x15 DW_TAG_structure_type
	.long	.Linfo_string112        # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x5d0:0xd DW_TAG_member
	.long	.Linfo_string89         # DW_AT_name
	.long	1271                    # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	21                      # Abbrev [21] 0x5de:0x14 DW_TAG_variable
	.long	.Linfo_string116        # DW_AT_name
	.long	1522                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV11T1DFunction
	.byte	16                      # Abbrev [16] 0x5f2:0xc DW_TAG_array_type
	.long	111                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x5f7:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	5                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	21                      # Abbrev [21] 0x5fe:0x14 DW_TAG_variable
	.long	.Linfo_string117        # DW_AT_name
	.long	1522                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV9Tinty_f2D
	.byte	21                      # Abbrev [21] 0x612:0x14 DW_TAG_variable
	.long	.Linfo_string118        # DW_AT_name
	.long	1522                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV10Tintxy_f3D
	.byte	21                      # Abbrev [21] 0x626:0x14 DW_TAG_variable
	.long	.Linfo_string119        # DW_AT_name
	.long	1522                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV11T2DFunction
	.byte	21                      # Abbrev [21] 0x63a:0x14 DW_TAG_variable
	.long	.Linfo_string120        # DW_AT_name
	.long	1522                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV8T2D_fix1
	.byte	21                      # Abbrev [21] 0x64e:0x14 DW_TAG_variable
	.long	.Linfo_string121        # DW_AT_name
	.long	1522                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTV8T3D_fix3
	.byte	22                      # Abbrev [22] 0x662:0x17 DW_TAG_variable
	.long	.Linfo_string122        # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	1                       # DW_AT_decl_file
	.short	7685                    # DW_AT_decl_line
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	__I___25_backgroundfield_quadr_cpp_6e0a9365
	.byte	23                      # Abbrev [23] 0x679:0x13 DW_TAG_variable
	.long	.Linfo_string123        # DW_AT_name
	.long	1676                    # DW_AT_type
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE
	.byte	10                      # Abbrev [10] 0x68c:0x11 DW_TAG_structure_type
	.long	.Linfo_string124        # DW_AT_name
	.byte	1                       # DW_AT_byte_size
	.byte	1                       # DW_AT_alignment
	.byte	24                      # Abbrev [24] 0x693:0x9 DW_TAG_member
	.long	1693                    # DW_AT_type
	.byte	1                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	16                      # Abbrev [16] 0x69d:0xc DW_TAG_array_type
	.long	996                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x6a2:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	0                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	3                       # Abbrev [3] 0x6a9:0xa DW_TAG_variable
	.long	.Linfo_string125        # DW_AT_name
	.long	802                     # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	21                      # Abbrev [21] 0x6b3:0x14 DW_TAG_variable
	.long	.Linfo_string126        # DW_AT_name
	.long	1735                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI11T1DFunction
	.byte	25                      # Abbrev [25] 0x6c7:0xa DW_TAG_structure_type
	.long	.Linfo_string127        # DW_AT_name
	.byte	16                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	21                      # Abbrev [21] 0x6d1:0x14 DW_TAG_variable
	.long	.Linfo_string128        # DW_AT_name
	.long	1735                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI11T2DFunction
	.byte	21                      # Abbrev [21] 0x6e5:0x14 DW_TAG_variable
	.long	.Linfo_string129        # DW_AT_name
	.long	1785                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI8T2D_fix1
	.byte	25                      # Abbrev [25] 0x6f9:0xa DW_TAG_structure_type
	.long	.Linfo_string130        # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	21                      # Abbrev [21] 0x703:0x14 DW_TAG_variable
	.long	.Linfo_string131        # DW_AT_name
	.long	1785                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI8T3D_fix3
	.byte	21                      # Abbrev [21] 0x717:0x14 DW_TAG_variable
	.long	.Linfo_string132        # DW_AT_name
	.long	1785                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI9Tinty_f2D
	.byte	21                      # Abbrev [21] 0x72b:0x14 DW_TAG_variable
	.long	.Linfo_string133        # DW_AT_name
	.long	1785                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTI10Tintxy_f3D
	.byte	3                       # Abbrev [3] 0x73f:0xa DW_TAG_variable
	.long	.Linfo_string134        # DW_AT_name
	.long	1865                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	16                      # Abbrev [16] 0x749:0xc DW_TAG_array_type
	.long	111                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x74e:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	0                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	21                      # Abbrev [21] 0x755:0x14 DW_TAG_variable
	.long	.Linfo_string135        # DW_AT_name
	.long	1897                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS11T1DFunction
	.byte	16                      # Abbrev [16] 0x769:0xc DW_TAG_array_type
	.long	996                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x76e:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	14                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	21                      # Abbrev [21] 0x775:0x14 DW_TAG_variable
	.long	.Linfo_string136        # DW_AT_name
	.long	1897                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS11T2DFunction
	.byte	3                       # Abbrev [3] 0x789:0xa DW_TAG_variable
	.long	.Linfo_string137        # DW_AT_name
	.long	1865                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	21                      # Abbrev [21] 0x793:0x14 DW_TAG_variable
	.long	.Linfo_string138        # DW_AT_name
	.long	1959                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS8T2D_fix1
	.byte	16                      # Abbrev [16] 0x7a7:0xc DW_TAG_array_type
	.long	996                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x7ac:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	10                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	21                      # Abbrev [21] 0x7b3:0x14 DW_TAG_variable
	.long	.Linfo_string139        # DW_AT_name
	.long	1959                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS8T3D_fix3
	.byte	21                      # Abbrev [21] 0x7c7:0x14 DW_TAG_variable
	.long	.Linfo_string140        # DW_AT_name
	.long	2011                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS9Tinty_f2D
	.byte	16                      # Abbrev [16] 0x7db:0xc DW_TAG_array_type
	.long	996                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x7e0:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	11                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	21                      # Abbrev [21] 0x7e7:0x14 DW_TAG_variable
	.long	.Linfo_string141        # DW_AT_name
	.long	2043                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	9                       # DW_AT_location
	.byte	3
	.quad	_ZTS10Tintxy_f3D
	.byte	16                      # Abbrev [16] 0x7fb:0xc DW_TAG_array_type
	.long	996                     # DW_AT_type
	.byte	17                      # Abbrev [17] 0x800:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	13                      # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	26                      # Abbrev [26] 0x807:0xd3 DW_TAG_subprogram
	.quad	.Lfunc_begin0           # DW_AT_low_pc
	.quad	.Lfunc_end0             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string152        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	156                     # DW_AT_decl_line
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0x825:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc0            # DW_AT_location
	.long	.Linfo_string169        # DW_AT_name
	.long	4227                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0x832:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc1            # DW_AT_location
	.long	.Linfo_string170        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0x83f:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc2            # DW_AT_location
	.long	.Linfo_string171        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0x84c:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc3            # DW_AT_location
	.long	.Linfo_string172        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x859:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\220~"
	.long	.Linfo_string167        # DW_AT_name
	.long	4215                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x866:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\310}"
	.long	.Linfo_string168        # DW_AT_name
	.long	4215                    # DW_AT_type
	.byte	29                      # Abbrev [29] 0x873:0x66 DW_TAG_lexical_block
	.long	.Ldebug_ranges0         # DW_AT_ranges
	.byte	30                      # Abbrev [30] 0x878:0x9 DW_TAG_variable
	.long	.Linfo_string172        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x881:0x9 DW_TAG_variable
	.long	.Linfo_string169        # DW_AT_name
	.long	4227                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x88a:0xd DW_TAG_variable
	.long	.Ldebug_loc4            # DW_AT_location
	.long	.Linfo_string173        # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x897:0xd DW_TAG_variable
	.long	.Ldebug_loc5            # DW_AT_location
	.long	.Linfo_string171        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x8a4:0xd DW_TAG_variable
	.long	.Ldebug_loc6            # DW_AT_location
	.long	.Linfo_string170        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x8b1:0xd DW_TAG_variable
	.long	.Ldebug_loc7            # DW_AT_location
	.long	.Linfo_string174        # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x8be:0xd DW_TAG_variable
	.long	.Ldebug_loc8            # DW_AT_location
	.long	.Linfo_string175        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x8cb:0xd DW_TAG_variable
	.long	.Ldebug_loc9            # DW_AT_location
	.long	.Linfo_string176        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	26                      # Abbrev [26] 0x8da:0x113 DW_TAG_subprogram
	.quad	.Lfunc_begin1           # DW_AT_low_pc
	.quad	.Lfunc_end1             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string153        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	180                     # DW_AT_decl_line
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0x8f8:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc10           # DW_AT_location
	.long	.Linfo_string169        # DW_AT_name
	.long	4227                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0x905:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc11           # DW_AT_location
	.long	.Linfo_string170        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0x912:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc12           # DW_AT_location
	.long	.Linfo_string171        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0x91f:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc13           # DW_AT_location
	.long	.Linfo_string172        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x92c:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\310z"
	.long	.Linfo_string167        # DW_AT_name
	.long	4215                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0x939:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\220{"
	.long	.Linfo_string168        # DW_AT_name
	.long	4215                    # DW_AT_type
	.byte	32                      # Abbrev [32] 0x946:0xa6 DW_TAG_lexical_block
	.quad	.Ltmp51                 # DW_AT_low_pc
	.quad	.Ltmp158                # DW_AT_high_pc
	.byte	30                      # Abbrev [30] 0x957:0x9 DW_TAG_variable
	.long	.Linfo_string172        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	30                      # Abbrev [30] 0x960:0x9 DW_TAG_variable
	.long	.Linfo_string169        # DW_AT_name
	.long	4227                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x969:0xd DW_TAG_variable
	.long	.Ldebug_loc14           # DW_AT_location
	.long	.Linfo_string173        # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x976:0xd DW_TAG_variable
	.long	.Ldebug_loc15           # DW_AT_location
	.long	.Linfo_string171        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x983:0xd DW_TAG_variable
	.long	.Ldebug_loc16           # DW_AT_location
	.long	.Linfo_string170        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x990:0xd DW_TAG_variable
	.long	.Ldebug_loc17           # DW_AT_location
	.long	.Linfo_string174        # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x99d:0xd DW_TAG_variable
	.long	.Ldebug_loc18           # DW_AT_location
	.long	.Linfo_string177        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x9aa:0xd DW_TAG_variable
	.long	.Ldebug_loc19           # DW_AT_location
	.long	.Linfo_string178        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x9b7:0xd DW_TAG_variable
	.long	.Ldebug_loc20           # DW_AT_location
	.long	.Linfo_string176        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x9c4:0xd DW_TAG_variable
	.long	.Ldebug_loc21           # DW_AT_location
	.long	.Linfo_string179        # DW_AT_name
	.long	131                     # DW_AT_type
	.byte	31                      # Abbrev [31] 0x9d1:0xd DW_TAG_variable
	.long	.Ldebug_loc22           # DW_AT_location
	.long	.Linfo_string180        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0x9de:0xd DW_TAG_variable
	.long	.Ldebug_loc23           # DW_AT_location
	.long	.Linfo_string181        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	26                      # Abbrev [26] 0x9ed:0xcd DW_TAG_subprogram
	.quad	.Lfunc_begin2           # DW_AT_low_pc
	.quad	.Lfunc_end2             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string153        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	230                     # DW_AT_decl_line
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0xa0b:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc24           # DW_AT_location
	.long	.Linfo_string169        # DW_AT_name
	.long	4232                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xa18:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc25           # DW_AT_location
	.long	.Linfo_string170        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xa25:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc26           # DW_AT_location
	.long	.Linfo_string171        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xa32:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc27           # DW_AT_location
	.long	.Linfo_string182        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xa3f:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc28           # DW_AT_location
	.long	.Linfo_string183        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xa4c:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc29           # DW_AT_location
	.long	.Linfo_string172        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	32                      # Abbrev [32] 0xa59:0x60 DW_TAG_lexical_block
	.quad	.Ltmp159                # DW_AT_low_pc
	.quad	.Ltmp222                # DW_AT_high_pc
	.byte	31                      # Abbrev [31] 0xa6a:0xd DW_TAG_variable
	.long	.Ldebug_loc30           # DW_AT_location
	.long	.Linfo_string172        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xa77:0xd DW_TAG_variable
	.long	.Ldebug_loc31           # DW_AT_location
	.long	.Linfo_string183        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xa84:0xd DW_TAG_variable
	.long	.Ldebug_loc32           # DW_AT_location
	.long	.Linfo_string182        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xa91:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\200~"
	.long	.Linfo_string171        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xa9e:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\370~"
	.long	.Linfo_string170        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xaab:0xd DW_TAG_variable
	.long	.Ldebug_loc33           # DW_AT_location
	.long	.Linfo_string169        # DW_AT_name
	.long	4232                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	26                      # Abbrev [26] 0xaba:0x101 DW_TAG_subprogram
	.quad	.Lfunc_begin3           # DW_AT_low_pc
	.quad	.Lfunc_end3             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string153        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	248                     # DW_AT_decl_line
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	27                      # Abbrev [27] 0xad8:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc34           # DW_AT_location
	.long	.Linfo_string169        # DW_AT_name
	.long	4237                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xae5:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc35           # DW_AT_location
	.long	.Linfo_string170        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xaf2:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc36           # DW_AT_location
	.long	.Linfo_string171        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xaff:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc37           # DW_AT_location
	.long	.Linfo_string182        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xb0c:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc38           # DW_AT_location
	.long	.Linfo_string183        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xb19:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc39           # DW_AT_location
	.long	.Linfo_string184        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xb26:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc40           # DW_AT_location
	.long	.Linfo_string144        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xb33:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc41           # DW_AT_location
	.long	.Linfo_string172        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	32                      # Abbrev [32] 0xb40:0x7a DW_TAG_lexical_block
	.quad	.Ltmp223                # DW_AT_low_pc
	.quad	.Ltmp287                # DW_AT_high_pc
	.byte	31                      # Abbrev [31] 0xb51:0xd DW_TAG_variable
	.long	.Ldebug_loc42           # DW_AT_location
	.long	.Linfo_string172        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xb5e:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\360}"
	.long	.Linfo_string144        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	28                      # Abbrev [28] 0xb6b:0xd DW_TAG_variable
	.byte	3                       # DW_AT_location
	.byte	145
	.ascii	"\370~"
	.long	.Linfo_string184        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xb78:0xd DW_TAG_variable
	.long	.Ldebug_loc43           # DW_AT_location
	.long	.Linfo_string183        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xb85:0xd DW_TAG_variable
	.long	.Ldebug_loc44           # DW_AT_location
	.long	.Linfo_string182        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xb92:0xd DW_TAG_variable
	.long	.Ldebug_loc45           # DW_AT_location
	.long	.Linfo_string171        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xb9f:0xd DW_TAG_variable
	.long	.Ldebug_loc46           # DW_AT_location
	.long	.Linfo_string170        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xbac:0xd DW_TAG_variable
	.long	.Ldebug_loc47           # DW_AT_location
	.long	.Linfo_string169        # DW_AT_name
	.long	4237                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	10                      # Abbrev [10] 0xbbb:0x52 DW_TAG_structure_type
	.long	.Linfo_string142        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0xbc2:0xd DW_TAG_member
	.long	.Linfo_string5          # DW_AT_name
	.long	111                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	33                      # Abbrev [33] 0xbcf:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin4           # DW_AT_low_pc
	.quad	.Lfunc_end4             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string154        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	29                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	34                      # Abbrev [34] 0xbe9:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4263                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	32                      # Abbrev [32] 0xbf1:0x1a DW_TAG_lexical_block
	.quad	.Ltmp288                # DW_AT_low_pc
	.quad	.Ltmp289                # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0xc02:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4263                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	36                      # Abbrev [36] 0xc0d:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin5           # DW_AT_low_pc
	.quad	.Lfunc_end5             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string155        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	10                      # Abbrev [10] 0xc25:0x52 DW_TAG_structure_type
	.long	.Linfo_string143        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0xc2c:0xd DW_TAG_member
	.long	.Linfo_string5          # DW_AT_name
	.long	111                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	33                      # Abbrev [33] 0xc39:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin6           # DW_AT_low_pc
	.quad	.Lfunc_end6             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string156        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	30                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	34                      # Abbrev [34] 0xc53:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4268                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	32                      # Abbrev [32] 0xc5b:0x1a DW_TAG_lexical_block
	.quad	.Ltmp292                # DW_AT_low_pc
	.quad	.Ltmp293                # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0xc6c:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4268                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	36                      # Abbrev [36] 0xc77:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin7           # DW_AT_low_pc
	.quad	.Lfunc_end7             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string157        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	26                      # Abbrev [26] 0xc8f:0x5f DW_TAG_subprogram
	.quad	.Lfunc_begin8           # DW_AT_low_pc
	.quad	.Lfunc_end8             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string147        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	41                      # DW_AT_decl_line
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	37                      # Abbrev [37] 0xcad:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc48           # DW_AT_location
	.long	3489                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0xcb7:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc49           # DW_AT_location
	.long	.Linfo_string185        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	32                      # Abbrev [32] 0xcc4:0x29 DW_TAG_lexical_block
	.quad	.Ltmp297                # DW_AT_low_pc
	.quad	.Ltmp299                # DW_AT_high_pc
	.byte	38                      # Abbrev [38] 0xcd5:0xa DW_TAG_variable
	.long	.Ldebug_loc50           # DW_AT_location
	.long	3489                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	31                      # Abbrev [31] 0xcdf:0xd DW_TAG_variable
	.long	.Ldebug_loc51           # DW_AT_location
	.long	.Linfo_string185        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	39                      # Abbrev [39] 0xcee:0x83 DW_TAG_class_type
	.long	.Linfo_string148        # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0xcf8:0xf DW_TAG_inheritance
	.long	.Linfo_string142        # DW_AT_name
	.long	3441                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0xd07:0xf DW_TAG_member
	.long	.Linfo_string144        # DW_AT_name
	.long	3451                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0xd16:0xf DW_TAG_member
	.long	.Linfo_string145        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	42                      # Abbrev [42] 0xd25:0xe DW_TAG_member
	.long	.Linfo_string147        # DW_AT_name
	.long	3473                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	41                      # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	33                      # Abbrev [33] 0xd33:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin9           # DW_AT_low_pc
	.quad	.Lfunc_end9             # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string158        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	42                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	34                      # Abbrev [34] 0xd4d:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	3489                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	32                      # Abbrev [32] 0xd55:0x1a DW_TAG_lexical_block
	.quad	.Ltmp300                # DW_AT_low_pc
	.quad	.Ltmp301                # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0xd66:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	3489                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	43                      # Abbrev [43] 0xd71:0xa DW_TAG_class_type
	.long	.Linfo_string142        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	7                       # Abbrev [7] 0xd7b:0x5 DW_TAG_pointer_type
	.long	3456                    # DW_AT_type
	.byte	43                      # Abbrev [43] 0xd80:0xa DW_TAG_class_type
	.long	.Linfo_string143        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	9                       # Abbrev [9] 0xd8a:0x7 DW_TAG_base_type
	.long	.Linfo_string146        # DW_AT_name
	.byte	4                       # DW_AT_encoding
	.byte	8                       # DW_AT_byte_size
	.byte	44                      # Abbrev [44] 0xd91:0x10 DW_TAG_subroutine_type
	.long	3466                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0xd96:0x5 DW_TAG_formal_parameter
	.long	3489                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0xd9b:0x5 DW_TAG_formal_parameter
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0xda1:0x5 DW_TAG_pointer_type
	.long	3310                    # DW_AT_type
	.byte	36                      # Abbrev [36] 0xda6:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin10          # DW_AT_low_pc
	.quad	.Lfunc_end10            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string159        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	26                      # Abbrev [26] 0xdbe:0x79 DW_TAG_subprogram
	.quad	.Lfunc_begin11          # DW_AT_low_pc
	.quad	.Lfunc_end11            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string147        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	83                      # DW_AT_decl_line
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	37                      # Abbrev [37] 0xddc:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc52           # DW_AT_location
	.long	3806                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0xde6:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc53           # DW_AT_location
	.long	.Linfo_string145        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	27                      # Abbrev [27] 0xdf3:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc54           # DW_AT_location
	.long	.Linfo_string185        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	32                      # Abbrev [32] 0xe00:0x36 DW_TAG_lexical_block
	.quad	.Ltmp304                # DW_AT_low_pc
	.quad	.Ltmp306                # DW_AT_high_pc
	.byte	38                      # Abbrev [38] 0xe11:0xa DW_TAG_variable
	.long	.Ldebug_loc55           # DW_AT_location
	.long	3806                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	31                      # Abbrev [31] 0xe1b:0xd DW_TAG_variable
	.long	.Ldebug_loc56           # DW_AT_location
	.long	.Linfo_string145        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	31                      # Abbrev [31] 0xe28:0xd DW_TAG_variable
	.long	.Ldebug_loc57           # DW_AT_location
	.long	.Linfo_string185        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	39                      # Abbrev [39] 0xe37:0x83 DW_TAG_class_type
	.long	.Linfo_string151        # DW_AT_name
	.byte	24                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0xe41:0xf DW_TAG_inheritance
	.long	.Linfo_string143        # DW_AT_name
	.long	3456                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0xe50:0xf DW_TAG_member
	.long	.Linfo_string144        # DW_AT_name
	.long	3770                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0xe5f:0xf DW_TAG_member
	.long	.Linfo_string150        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	42                      # Abbrev [42] 0xe6e:0xe DW_TAG_member
	.long	.Linfo_string147        # DW_AT_name
	.long	3785                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	83                      # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	33                      # Abbrev [33] 0xe7c:0x3d DW_TAG_subprogram
	.quad	.Lfunc_begin12          # DW_AT_low_pc
	.quad	.Lfunc_end12            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string160        # DW_AT_name
	.byte	2                       # DW_AT_decl_file
	.byte	84                      # DW_AT_decl_line
	.byte	1                       # DW_AT_external
	.byte	34                      # Abbrev [34] 0xe96:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	3806                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	32                      # Abbrev [32] 0xe9e:0x1a DW_TAG_lexical_block
	.quad	.Ltmp307                # DW_AT_low_pc
	.quad	.Ltmp308                # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0xeaf:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	3806                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0xeba:0x5 DW_TAG_pointer_type
	.long	3775                    # DW_AT_type
	.byte	43                      # Abbrev [43] 0xebf:0xa DW_TAG_class_type
	.long	.Linfo_string149        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	44                      # Abbrev [44] 0xec9:0x15 DW_TAG_subroutine_type
	.long	3466                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0xece:0x5 DW_TAG_formal_parameter
	.long	3806                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0xed3:0x5 DW_TAG_formal_parameter
	.long	3466                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0xed8:0x5 DW_TAG_formal_parameter
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0xede:0x5 DW_TAG_pointer_type
	.long	3639                    # DW_AT_type
	.byte	36                      # Abbrev [36] 0xee3:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin13          # DW_AT_low_pc
	.quad	.Lfunc_end13            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string161        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	26                      # Abbrev [26] 0xefb:0x5f DW_TAG_subprogram
	.quad	.Lfunc_begin14          # DW_AT_low_pc
	.quad	.Lfunc_end14            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string147        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	226                     # DW_AT_decl_line
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	37                      # Abbrev [37] 0xf19:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc58           # DW_AT_location
	.long	4273                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0xf23:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc59           # DW_AT_location
	.long	.Linfo_string145        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	32                      # Abbrev [32] 0xf30:0x29 DW_TAG_lexical_block
	.quad	.Ltmp311                # DW_AT_low_pc
	.quad	.Ltmp315                # DW_AT_high_pc
	.byte	38                      # Abbrev [38] 0xf41:0xa DW_TAG_variable
	.long	.Ldebug_loc60           # DW_AT_location
	.long	4273                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	31                      # Abbrev [31] 0xf4b:0xd DW_TAG_variable
	.long	.Ldebug_loc61           # DW_AT_location
	.long	.Linfo_string145        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	45                      # Abbrev [45] 0xf5a:0x3b DW_TAG_subprogram
	.quad	.Lfunc_begin15          # DW_AT_low_pc
	.quad	.Lfunc_end15            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string162        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	34                      # Abbrev [34] 0xf72:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4273                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	32                      # Abbrev [32] 0xf7a:0x1a DW_TAG_lexical_block
	.quad	.Ltmp316                # DW_AT_low_pc
	.quad	.Ltmp317                # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0xf8b:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4273                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	36                      # Abbrev [36] 0xf95:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin16          # DW_AT_low_pc
	.quad	.Lfunc_end16            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string163        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	26                      # Abbrev [26] 0xfad:0x5f DW_TAG_subprogram
	.quad	.Lfunc_begin17          # DW_AT_low_pc
	.quad	.Lfunc_end17            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string147        # DW_AT_name
	.byte	1                       # DW_AT_decl_file
	.byte	244                     # DW_AT_decl_line
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_external
	.byte	37                      # Abbrev [37] 0xfcb:0xa DW_TAG_formal_parameter
	.long	.Ldebug_loc62           # DW_AT_location
	.long	4394                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	27                      # Abbrev [27] 0xfd5:0xd DW_TAG_formal_parameter
	.long	.Ldebug_loc63           # DW_AT_location
	.long	.Linfo_string150        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	32                      # Abbrev [32] 0xfe2:0x29 DW_TAG_lexical_block
	.quad	.Ltmp320                # DW_AT_low_pc
	.quad	.Ltmp323                # DW_AT_high_pc
	.byte	38                      # Abbrev [38] 0xff3:0xa DW_TAG_variable
	.long	.Ldebug_loc64           # DW_AT_location
	.long	4394                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	31                      # Abbrev [31] 0xffd:0xd DW_TAG_variable
	.long	.Ldebug_loc65           # DW_AT_location
	.long	.Linfo_string150        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	45                      # Abbrev [45] 0x100c:0x3b DW_TAG_subprogram
	.quad	.Lfunc_begin18          # DW_AT_low_pc
	.quad	.Lfunc_end18            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string164        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	34                      # Abbrev [34] 0x1024:0x8 DW_TAG_formal_parameter
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4394                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	32                      # Abbrev [32] 0x102c:0x1a DW_TAG_lexical_block
	.quad	.Ltmp324                # DW_AT_low_pc
	.quad	.Ltmp325                # DW_AT_high_pc
	.byte	35                      # Abbrev [35] 0x103d:0x8 DW_TAG_variable
	.byte	1                       # DW_AT_location
	.byte	85
	.long	4394                    # DW_AT_type
	.byte	1                       # DW_AT_artificial
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
	.byte	36                      # Abbrev [36] 0x1047:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin19          # DW_AT_low_pc
	.quad	.Lfunc_end19            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string165        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	36                      # Abbrev [36] 0x105f:0x18 DW_TAG_subprogram
	.quad	.Lfunc_begin20          # DW_AT_low_pc
	.quad	.Lfunc_end20            # DW_AT_high_pc
	.byte	1                       # DW_AT_frame_base
	.byte	86
	.long	.Linfo_string166        # DW_AT_name
	.byte	1                       # DW_AT_external
	.byte	16                      # Abbrev [16] 0x1077:0xc DW_TAG_array_type
	.long	3466                    # DW_AT_type
	.byte	17                      # Abbrev [17] 0x107c:0x6 DW_TAG_subrange_type
	.long	824                     # DW_AT_type
	.byte	9                       # DW_AT_count
	.byte	0                       # End Of Children Mark
	.byte	46                      # Abbrev [46] 0x1083:0x5 DW_TAG_reference_type
	.long	3003                    # DW_AT_type
	.byte	46                      # Abbrev [46] 0x1088:0x5 DW_TAG_reference_type
	.long	3109                    # DW_AT_type
	.byte	46                      # Abbrev [46] 0x108d:0x5 DW_TAG_reference_type
	.long	4242                    # DW_AT_type
	.byte	10                      # Abbrev [10] 0x1092:0x15 DW_TAG_structure_type
	.long	.Linfo_string149        # DW_AT_name
	.byte	8                       # DW_AT_byte_size
	.byte	8                       # DW_AT_alignment
	.byte	6                       # Abbrev [6] 0x1099:0xd DW_TAG_member
	.long	.Linfo_string5          # DW_AT_name
	.long	111                     # DW_AT_type
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x10a7:0x5 DW_TAG_pointer_type
	.long	3003                    # DW_AT_type
	.byte	7                       # Abbrev [7] 0x10ac:0x5 DW_TAG_pointer_type
	.long	3109                    # DW_AT_type
	.byte	7                       # Abbrev [7] 0x10b1:0x5 DW_TAG_pointer_type
	.long	4278                    # DW_AT_type
	.byte	39                      # Abbrev [39] 0x10b6:0x64 DW_TAG_class_type
	.long	.Linfo_string189        # DW_AT_name
	.byte	40                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0x10c0:0xf DW_TAG_inheritance
	.long	.Linfo_string142        # DW_AT_name
	.long	3441                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0x10cf:0xf DW_TAG_member
	.long	.Linfo_string144        # DW_AT_name
	.long	3451                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0x10de:0xf DW_TAG_member
	.long	.Linfo_string186        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	41                      # Abbrev [41] 0x10ed:0xf DW_TAG_member
	.long	.Linfo_string187        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	41                      # Abbrev [41] 0x10fc:0xf DW_TAG_member
	.long	.Linfo_string188        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	42                      # Abbrev [42] 0x110b:0xe DW_TAG_member
	.long	.Linfo_string147        # DW_AT_name
	.long	4378                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	226                     # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	44                      # Abbrev [44] 0x111a:0x10 DW_TAG_subroutine_type
	.long	3466                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x111f:0x5 DW_TAG_formal_parameter
	.long	4273                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x1124:0x5 DW_TAG_formal_parameter
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	7                       # Abbrev [7] 0x112a:0x5 DW_TAG_pointer_type
	.long	4399                    # DW_AT_type
	.byte	39                      # Abbrev [39] 0x112f:0x82 DW_TAG_class_type
	.long	.Linfo_string195        # DW_AT_name
	.byte	56                      # DW_AT_byte_size
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	4                       # DW_AT_APPLE_runtime_class
	.byte	8                       # DW_AT_alignment
	.byte	40                      # Abbrev [40] 0x1139:0xf DW_TAG_inheritance
	.long	.Linfo_string142        # DW_AT_name
	.long	3441                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	41                      # Abbrev [41] 0x1148:0xf DW_TAG_member
	.long	.Linfo_string144        # DW_AT_name
	.long	3770                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	8
	.byte	41                      # Abbrev [41] 0x1157:0xf DW_TAG_member
	.long	.Linfo_string190        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	16
	.byte	41                      # Abbrev [41] 0x1166:0xf DW_TAG_member
	.long	.Linfo_string191        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	24
	.byte	41                      # Abbrev [41] 0x1175:0xf DW_TAG_member
	.long	.Linfo_string192        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	32
	.byte	41                      # Abbrev [41] 0x1184:0xf DW_TAG_member
	.long	.Linfo_string193        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	40
	.byte	41                      # Abbrev [41] 0x1193:0xf DW_TAG_member
	.long	.Linfo_string194        # DW_AT_name
	.long	3466                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	1                       # DW_AT_decl_line
	.byte	8                       # DW_AT_alignment
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	48
	.byte	42                      # Abbrev [42] 0x11a2:0xe DW_TAG_member
	.long	.Linfo_string147        # DW_AT_name
	.long	4529                    # DW_AT_type
	.byte	1                       # DW_AT_decl_file
	.byte	244                     # DW_AT_decl_line
	.byte	2                       # DW_AT_data_member_location
	.byte	35
	.byte	0
	.byte	0                       # End Of Children Mark
	.byte	44                      # Abbrev [44] 0x11b1:0x10 DW_TAG_subroutine_type
	.long	3466                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x11b6:0x5 DW_TAG_formal_parameter
	.long	4394                    # DW_AT_type
	.byte	14                      # Abbrev [14] 0x11bb:0x5 DW_TAG_formal_parameter
	.long	3466                    # DW_AT_type
	.byte	0                       # End Of Children Mark
	.byte	0                       # End Of Children Mark
.Ldebug_info_end0:
	.section	.debug_ranges,"",@progbits
.Ldebug_ranges0:
	.quad	.Ltmp1-.Lfunc_begin0
	.quad	.Ltmp34-.Lfunc_begin0
	.quad	.Ltmp35-.Lfunc_begin0
	.quad	.Ltmp49-.Lfunc_begin0
	.quad	0
	.quad	0
	.section	.debug_pubnames,"",@progbits
	.long	.LpubNames_end0-.LpubNames_begin0 # Length of Public Names Info
.LpubNames_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	4546                    # Compilation Unit Length
	.long	503                     # DIE offset
	.asciz	"_ZSt12_S_showpoint"    # External Name
	.long	3989                    # DIE offset
	.asciz	"_ZN9Tinty_f2DD0Ev"     # External Name
	.long	1715                    # DIE offset
	.asciz	"_ZTI11T1DFunction"     # External Name
	.long	3085                    # DIE offset
	.asciz	"_ZN11T1DFunctionD0Ev"  # External Name
	.long	4108                    # DIE offset
	.asciz	"~Tintxy_f3D"           # External Name
	.long	47                      # DIE offset
	.asciz	"std"                   # External Name
	.long	517                     # DIE offset
	.asciz	"_ZSt9_S_skipws"        # External Name
	.long	1765                    # DIE offset
	.asciz	"_ZTI8T2D_fix1"         # External Name
	.long	463                     # DIE offset
	.asciz	"_ZSt11_S_internal"     # External Name
	.long	603                     # DIE offset
	.asciz	"_ZSt9_S_badbit"        # External Name
	.long	3494                    # DIE offset
	.asciz	"_ZN8T2D_fix1D0Ev"      # External Name
	.long	445                     # DIE offset
	.asciz	"_ZSt6_S_dec"           # External Name
	.long	4191                    # DIE offset
	.asciz	"__sti___25_backgroundfield_quadr_cpp_6e0a9365" # External Name
	.long	2055                    # DIE offset
	.asciz	"Romberg_simple"        # External Name
	.long	1929                    # DIE offset
	.asciz	"_ZTVN10__cxxabiv120__si_class_type_infoE" # External Name
	.long	451                     # DIE offset
	.asciz	"_ZSt8_S_fixed"         # External Name
	.long	1971                    # DIE offset
	.asciz	"_ZTS8T3D_fix3"         # External Name
	.long	1855                    # DIE offset
	.asciz	"_ZTVN10__cxxabiv117__class_type_infoE" # External Name
	.long	4167                    # DIE offset
	.asciz	"_ZN10Tintxy_f3DD0Ev"   # External Name
	.long	510                     # DIE offset
	.asciz	"_ZSt10_S_showpos"      # External Name
	.long	750                     # DIE offset
	.asciz	"_ZNSt8ios_base11imbue_eventE" # External Name
	.long	1745                    # DIE offset
	.asciz	"_ZTI11T2DFunction"     # External Name
	.long	3191                    # DIE offset
	.asciz	"_ZN11T2DFunctionD0Ev"  # External Name
	.long	457                     # DIE offset
	.asciz	"_ZSt6_S_hex"           # External Name
	.long	475                     # DIE offset
	.asciz	"_ZSt6_S_oct"           # External Name
	.long	496                     # DIE offset
	.asciz	"_ZSt11_S_showbase"     # External Name
	.long	1534                    # DIE offset
	.asciz	"_ZTV9Tinty_f2D"        # External Name
	.long	469                     # DIE offset
	.asciz	"_ZSt7_S_left"          # External Name
	.long	1502                    # DIE offset
	.asciz	"_ZTV11T1DFunction"     # External Name
	.long	1594                    # DIE offset
	.asciz	"_ZTV8T2D_fix1"         # External Name
	.long	579                     # DIE offset
	.asciz	"_ZSt19_S_ios_fmtflags_min" # External Name
	.long	569                     # DIE offset
	.asciz	"_ZSt19_S_ios_fmtflags_max" # External Name
	.long	639                     # DIE offset
	.asciz	"_ZSt18_S_ios_iostate_min" # External Name
	.long	629                     # DIE offset
	.asciz	"_ZSt18_S_ios_iostate_max" # External Name
	.long	609                     # DIE offset
	.asciz	"_ZSt9_S_eofbit"        # External Name
	.long	744                     # DIE offset
	.asciz	"_ZNSt8ios_base11erase_eventE" # External Name
	.long	1835                    # DIE offset
	.asciz	"_ZTI10Tintxy_f3D"      # External Name
	.long	482                     # DIE offset
	.asciz	"_ZSt8_S_right"         # External Name
	.long	2023                    # DIE offset
	.asciz	"_ZTS10Tintxy_f3D"      # External Name
	.long	1634                    # DIE offset
	.asciz	"__I___25_backgroundfield_quadr_cpp_6e0a9365" # External Name
	.long	554                     # DIE offset
	.asciz	"_ZSt13_S_floatfield"   # External Name
	.long	561                     # DIE offset
	.asciz	"_ZSt19_S_ios_fmtflags_end" # External Name
	.long	1554                    # DIE offset
	.asciz	"_ZTV10Tintxy_f3D"      # External Name
	.long	1795                    # DIE offset
	.asciz	"_ZTI8T3D_fix3"         # External Name
	.long	547                     # DIE offset
	.asciz	"_ZSt12_S_basefield"    # External Name
	.long	621                     # DIE offset
	.asciz	"_ZSt18_S_ios_iostate_end" # External Name
	.long	1574                    # DIE offset
	.asciz	"_ZTV11T2DFunction"     # External Name
	.long	489                     # DIE offset
	.asciz	"_ZSt13_S_scientific"   # External Name
	.long	1877                    # DIE offset
	.asciz	"_ZTS11T1DFunction"     # External Name
	.long	3129                    # DIE offset
	.asciz	"T2DFunction::~T2DFunction" # External Name
	.long	1705                    # DIE offset
	.asciz	"__dso_handle"          # External Name
	.long	1939                    # DIE offset
	.asciz	"_ZTS8T2D_fix1"         # External Name
	.long	3930                    # DIE offset
	.asciz	"~Tinty_f2D"            # External Name
	.long	3811                    # DIE offset
	.asciz	"_ZN8T3D_fix3D0Ev"      # External Name
	.long	3379                    # DIE offset
	.asciz	"T2D_fix1::~T2D_fix1"   # External Name
	.long	52                      # DIE offset
	.asciz	"std::_ZSt4cerr"        # External Name
	.long	1815                    # DIE offset
	.asciz	"_ZTI9Tinty_f2D"        # External Name
	.long	1657                    # DIE offset
	.asciz	"_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE" # External Name
	.long	439                     # DIE offset
	.asciz	"_ZSt12_S_boolalpha"    # External Name
	.long	540                     # DIE offset
	.asciz	"_ZSt14_S_adjustfield"  # External Name
	.long	1991                    # DIE offset
	.asciz	"_ZTS9Tinty_f2D"        # External Name
	.long	615                     # DIE offset
	.asciz	"_ZSt10_S_failbit"      # External Name
	.long	2746                    # DIE offset
	.asciz	"Romberg"               # External Name
	.long	524                     # DIE offset
	.asciz	"_ZSt10_S_unitbuf"      # External Name
	.long	1909                    # DIE offset
	.asciz	"_ZTS11T2DFunction"     # External Name
	.long	3023                    # DIE offset
	.asciz	"T1DFunction::~T1DFunction" # External Name
	.long	756                     # DIE offset
	.asciz	"_ZNSt8ios_base13copyfmt_eventE" # External Name
	.long	597                     # DIE offset
	.asciz	"_ZSt10_S_goodbit"      # External Name
	.long	3708                    # DIE offset
	.asciz	"T3D_fix3::~T3D_fix3"   # External Name
	.long	532                     # DIE offset
	.asciz	"_ZSt12_S_uppercase"    # External Name
	.long	4013                    # DIE offset
	.asciz	"call"                  # External Name
	.long	1614                    # DIE offset
	.asciz	"_ZTV8T3D_fix3"         # External Name
	.long	0                       # End Mark
.LpubNames_end0:
	.section	.debug_pubtypes,"",@progbits
	.long	.LpubTypes_end0-.LpubTypes_begin0 # Length of Public Types Info
.LpubTypes_begin0:
	.short	2                       # DWARF Version
	.long	.Lcu_begin0             # Offset of Compilation Unit Info
	.long	4546                    # Compilation Unit Length
	.long	1735                    # DIE offset
	.asciz	"__class_type_info"     # External Name
	.long	432                     # DIE offset
	.asciz	"_ZSt13_Ios_Fmtflags"   # External Name
	.long	590                     # DIE offset
	.asciz	"_ZSt12_Ios_Iostate"    # External Name
	.long	836                     # DIE offset
	.asciz	"_ZSt6locale"           # External Name
	.long	737                     # DIE offset
	.asciz	"_ZNSt8ios_base5eventE" # External Name
	.long	131                     # DIE offset
	.asciz	"int"                   # External Name
	.long	1130                    # DIE offset
	.asciz	"_ZSt5ctypeIcE"         # External Name
	.long	979                     # DIE offset
	.asciz	"unsigned long"         # External Name
	.long	3441                    # DIE offset
	.asciz	"T1DFunction"           # External Name
	.long	945                     # DIE offset
	.asciz	"_ZNSt6locale5facetE"   # External Name
	.long	807                     # DIE offset
	.asciz	"void"                  # External Name
	.long	4278                    # DIE offset
	.asciz	"Tinty_f2D"             # External Name
	.long	4242                    # DIE offset
	.asciz	"T3DFunction"           # External Name
	.long	1271                    # DIE offset
	.asciz	"__SO__NSt6locale5facetE" # External Name
	.long	996                     # DIE offset
	.asciz	"signed char"           # External Name
	.long	1676                    # DIE offset
	.asciz	"_ZNSt8ios_base4InitE"  # External Name
	.long	62                      # DIE offset
	.asciz	"std::ostream"          # External Name
	.long	3466                    # DIE offset
	.asciz	"double"                # External Name
	.long	655                     # DIE offset
	.asciz	"_ZNSt8ios_base14_Callback_listE" # External Name
	.long	258                     # DIE offset
	.asciz	"_ZSt8ios_base"         # External Name
	.long	862                     # DIE offset
	.asciz	"_ZNSt6locale5_ImplE"   # External Name
	.long	1455                    # DIE offset
	.asciz	"_ZSt7num_putIcSt19ostreambuf_iteratorIcSt11char_traitsIcEEE" # External Name
	.long	1481                    # DIE offset
	.asciz	"_ZSt7num_getIcSt19istreambuf_iteratorIcSt11char_traitsIcEEE" # External Name
	.long	138                     # DIE offset
	.asciz	"_ZSt9basic_iosIcSt11char_traitsIcEE" # External Name
	.long	3639                    # DIE offset
	.asciz	"T3D_fix3"              # External Name
	.long	1401                    # DIE offset
	.asciz	"__locale_data"         # External Name
	.long	4399                    # DIE offset
	.asciz	"Tintxy_f3D"            # External Name
	.long	3310                    # DIE offset
	.asciz	"T2D_fix1"              # External Name
	.long	1413                    # DIE offset
	.asciz	"unsigned short"        # External Name
	.long	1013                    # DIE offset
	.asciz	"_ZSt15basic_streambufIcSt11char_traitsIcEE" # External Name
	.long	425                     # DIE offset
	.asciz	"long"                  # External Name
	.long	3456                    # DIE offset
	.asciz	"T2DFunction"           # External Name
	.long	1310                    # DIE offset
	.asciz	"__locale_struct"       # External Name
	.long	768                     # DIE offset
	.asciz	"_ZNSt8ios_base6_WordsE" # External Name
	.long	76                      # DIE offset
	.asciz	"_ZSo"                  # External Name
	.long	1785                    # DIE offset
	.asciz	"__si_class_type_info"  # External Name
	.long	0                       # End Mark
.LpubTypes_end0:
	.section	".note.GNU-stack","",@progbits
	.section	.debug_line,"",@progbits
.Lline_table_start0:
