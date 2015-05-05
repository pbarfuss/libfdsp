        .fpu neon
        .syntax unified
.eabi_attribute 24, 1
.eabi_attribute 25, 1
.macro function name, export=0
.macro endfunc
.size \name, . - \name
        .endfunc
        .purgem endfunc
.endm
        .text
        .global \name
\name:
        .type \name, %function
        .func \name
\name:
.endm

        .text
function vector_fmul_neon, export=1
1:
        vld1.32 {q0,q1}, [r1]!
        vld1.32 {q8,q9}, [r2]!
        vmul.f32 q10, q0, q8
        vmul.f32 q11, q1, q9
        vst1.32 {q10,q11}, [r0]!
        subs r3, r3, #8
        beq 2f
        b 1b
2:
        bx lr
endfunc
function vector_fmac_scalar_neon, export=1
        vdup.32         q15, r2
        bics            r12, r3, #7
        mov             r2, r0
        beq             3f
1:
        vld1.32         {q0},     [r1]!
        vld1.32         {q8},     [r2]!
        vld1.32         {q1},     [r1]!
        vld1.32         {q9},     [r2]!
        vmla.f32        q8,  q0,  q15
        vmla.f32        q9,  q1,  q15
        vst1.32         {q8},     [r0]!
        vst1.32         {q9},     [r0]!

        subs            r12, r12, #8
        beq             2f
        b               1b
2:
        ands            r3, r3, #7
        it              eq
        bxeq            lr
3:
        vld1.32         {q0},     [r1]!
        vld1.32         {q8},     [r2]!
        vmla.f32        q8,  q0,  q15
        vst1.32         {q8},     [r0]!
        subs            r3, r3, #4
        bgt             3b
        bx              lr
endfunc
function vector_fmul_scalar_neon, export=1
        vdup.32  q8, r2
1:
        vld1.32 {q0},[r1]!
        vmul.f32 q0, q0, q8
        vst1.32 {q0},[r0]!
        subs r3, r3, #4
        bgt 1b
        bx lr
endfunc
function vector_fmul_window_neon, export=1
        push {r4,r5,lr}
        ldr lr, [sp, #12]
        sub r2, r2, #8
        sub r5, lr, #2
        add r2, r2, r5, lsl #2
        add r4, r3, r5, lsl #3
        add ip, r0, r5, lsl #3
        mov r5, #-16
1:
        vld1.32 {d0,d1}, [r1]!
        vld1.32 {d2,d3}, [r2], r5
        vld1.32 {d4,d5}, [r3]!
        vld1.32 {d6,d7}, [r4], r5
        vmul.f32 d22, d0, d4
        vrev64.32 q3, q3
        vmul.f32 d23, d1, d5
        vrev64.32 q1, q1
        vmul.f32 d20, d0, d7
        vmul.f32 d21, d1, d6
        vmla.f32 d22, d3, d7
        vmla.f32 d23, d2, d6
        vmls.f32 d20, d3, d4
        vmls.f32 d21, d2, d5
        vrev64.32 q11, q11
        vswp d22, d23
        vst1.32 {d20,d21},[r0]!
        vst1.32 {d22,d23},[ip], r5
        subs lr, lr, #4
        beq 2f
        b 1b
2:
        pop {r4,r5,pc}
endfunc
function vector_fmul_reverse_neon, export=1
        add r2, r2, r3, lsl #2
        sub r2, r2, #32
        mov r12, #-32
1:
        vld1.32 {q0-q1}, [r1]!
        vld1.32 {q2-q3}, [r2], r12
        pld [r1, #32]
        pld [r2, #-32]
        vrev64.32 q2, q2
        vrev64.32 q3, q3
        vmul.f32 d16, d0, d7
        vmul.f32 d17, d1, d6
        vmul.f32 d18, d2, d5
        vmul.f32 d19, d3, d4
        vst1.32 {q8-q9}, [r0]!
        subs r3, r3, #8
        beq 2f
        b 1b
2:
        bx lr
endfunc
function vector_fmul_copy_neon, export=1
1:
        vld1.32 {q0,q1},[r1]!
        vld1.32 {q2,q3},[r1]!
        vst1.32 {q0,q1},[r0]!
        vst1.32 {q2,q3},[r0]!
        subs r2, r2, #16
        beq 2f
        b 1b
2:
        bx lr
endfunc
function vector_fmul_add_neon, export=1
        ldr r12, [sp]
1:
        vld1.32 {q0,q1}, [r1]!
        vld1.32 {q8,q9}, [r2]!
        vld1.32 {q2,q3}, [r3]!
        vmul.f32 q10, q0, q8
        vmul.f32 q11, q1, q9
        vadd.f32 q12, q2, q10
        vadd.f32 q13, q3, q11
        vst1.32 {q12,q13}, [r0]!
        subs r12, r12, #8
        beq 2f
        b 1b
2:
        bx lr
endfunc
function butterflies_float_neon, export=1
1:
        vld1.32 {q0},[r0]
        vld1.32 {q1},[r1]
        vsub.f32 q2, q0, q1
        vadd.f32 q1, q0, q1
        vst1.32 {q2},[r1]!
        vst1.32 {q1},[r0]!
        subs r2, r2, #4
        bgt 1b
        bx lr
endfunc
function scalarproduct_float_neon, export=1
        vmov.f32 q2, #0.0
1:
        vld1.32 {q0},[r0]!
        vld1.32 {q1},[r1]!
        vmla.f32 q2, q0, q1
        subs r2, r2, #4
        bgt 1b
        vadd.f32 d0, d4, d5
        vpadd.f32 d0, d0, d0
        vmov.32 r0, d0[0]
        bx lr
endfunc
function scalarproduct_symmetric_fir_float_neon, export=1
        add r3, r1, r2, lsl #3
        mov r12, #-12
        vmov.f32 q2, #0.0
1:
        vld1.32 {q0},[r0]!
        vld1.32 {q1},[r3], r12
        vld1.32 {q2},[r1]!
        vmla.f32 q3, q0, q2
        vmla.f32 q3, q1, q2
        subs r2, r2, #4
        bgt 1b
        vmov.f32 q1, #0.0
        vld1.32 {d2[0]}, [r0]
        vld1.32 {d3[0]}, [r1]
        vadd.f32 d0, d6, d7
        vpadd.f32 d0, d0, d0
        vmla.f32 d0, d2, d3
        vmov.32 r0, d0[0]
        bx lr
endfunc
function vector_clipf_neon, export=1
        vdup.32         q0,  r2
        vdup.32         q1,  r3
        ldr             r2,  [sp]
1:      
        vld1.f32        {q2},[r1]!
        vmin.f32        q10, q2,  q1
        vmax.f32        q10, q10, q0
        vld1.f32        {q3},[r1]!
        vmin.f32        q11, q3,  q1
        vmax.f32        q11, q11, q0
        vst1.f32        {q10},[r0]!
        vst1.f32        {q11},[r0]!
        subs            r2,  r2,  #8
        bgt 1b
        bx              lr
endfunc

function sbr_sum64x5_neon, export=1
        push {lr}
        add r1, r0, # 64*4
        add r2, r0, #128*4
        add r3, r0, #192*4
        add lr, r0, #256*4
        mov r12, #64
1:
        vld1.32 {q0}, [r0]
        vld1.32 {q1}, [r1]!
        vld1.32 {q2}, [r2]!
        vld1.32 {q3}, [r3]!
        vld1.32 {q8}, [lr]!
        vadd.f32 q0, q0, q1
        vadd.f32 q0, q0, q2
        vadd.f32 q0, q0, q3
        vadd.f32 q0, q0, q8
        vst1.32 {q0}, [r0]!
        subs r12, #4
        bgt 1b
        pop {pc}
endfunc
function sbr_qmf_pre_shuffle_neon, export=1
        add r1, r0, #60*4
        add r2, r0, #64*4
        vld1.32 {d0}, [r0]!
        vst1.32 {d0}, [r2]!
        mov r3, #-16
        mov r12, #24
        vmov.i32 q8, #1<<31
        vld1.32 {q0}, [r1], r3
        vld1.32 {d2}, [r0]!
1:
        vld1.32 {d3,d4}, [r0]!
        vrev64.32 q0, q0
        vld1.32 {q9}, [r1], r3
        veor q0, q0, q8
        vld1.32 {d5,d6}, [r0]!
        vswp d0, d1
        vrev64.32 q9, q9
        vst2.32 {q0,q1}, [r2]!
        vmov q10, q2
        veor q9, q9, q8
        vmov d2, d6
        vswp d18, d19
        vld1.32 {q0}, [r1], r3
        vst2.32 {q9,q10}, [r2]!
        subs r12, r12, #8
        bgt 1b
        vld1.32 {d3,d4}, [r0]!
        vrev64.32 q0, q0
        vld1.32 {q9}, [r1], r3
        veor q0, q0, q8
        vld1.32 {d5}, [r0]!
        vswp d0, d1
        vrev64.32 q9, q9
        vst2.32 {q0,q1}, [r2]!
        vswp d4, d5
        veor q1, q9, q8
        vst2.32 {d3,d5}, [r2]!
        vst2.32 {d2[0],d4[0]}, [r2]!
        bx lr
endfunc
function sbr_qmf_post_shuffle_neon, export=1
        add r2, r1, #60*4
        mov r3, #-16
        mov r12, #32
        vmov.i32 q8, #1<<31
1:
        pld [r2, #-32]
        vld1.32 {q1}, [r2], r3
        vld1.32 {q0}, [r1]!
        vrev64.32 q1, q1
        vswp d0, d1
        veor q0, q0, q8
        vst2.32 {d1,d3}, [r0]!
        vst2.32 {d0,d2}, [r0]!
        subs r12, r12, #4
        bgt 1b
        bx lr
endfunc
function sbr_qmf_deint_bfly_neon, export=1
        push            {lr}
        add             r2,  r2,  #60*4
        add             r3,  r0,  #124*4
        mov             r12, #64
        mov             lr,  #-16
1:
        vld1.32         {q0},     [r1]!
        vld1.32         {q1},     [r2], lr
        vrev64.32       q2,  q0
        vrev64.32       q3,  q1
        vadd.f32        d3,  d4,  d3
        vadd.f32        d2,  d5,  d2
        vsub.f32        d0,  d7,  d0
        vsub.f32        d1,  d6,  d1
        vst1.32         {q1},     [r3], lr
        vst1.32         {q0},     [r0]!
        subs            r12, r12, #4
        bgt             1b
        pop             {pc}
endfunc
function sbr_hf_g_filt_neon, export=1
        ldr             r12, [sp]
        add             r1,  r1,  r12, lsl #3
        mov             r12, #40*2*4
        sub             r3,  r3,  #1
        vld2.32         {d2[],d3[]}, [r2]!
        vld1.32         {d0},        [r1], r12
1:
        vld1.32         {d1},        [r1], r12
        vmul.f32        q3,  q0,  q1
        vld2.32         {d2[],d3[]}, [r2]!
        vld1.32         {d0},        [r1], r12
        vst1.32         {q3},        [r0]!
        subs            r3,  r3,  #2
        bgt             1b
        it              lt
        bxlt            lr
        vmul.f32        d0,  d0,  d2
        vst1.32         {d0},        [r0]!
        bx              lr
endfunc
function sbrenc_sum128x5_neon, export=1
        push {lr}
        add r1, r0, #128*4
        add r2, r0, #256*4
        add r3, r0, #384*4
        add lr, r0, #512*4
        mov r12, #128
1:
        vld1.32 {q0}, [r0]
        vld1.32 {q1}, [r1]!
        vld1.32 {q2}, [r2]!
        vld1.32 {q3}, [r3]!
        vld1.32 {q8}, [lr]!
        vadd.f32 q0, q0, q1
        vadd.f32 q0, q0, q2
        vadd.f32 q0, q0, q3
        vadd.f32 q0, q0, q8
        vst1.32 {q0}, [r0]!
        subs r12, #4
        bgt 1b
        pop {pc}
endfunc
function sbrenc_qmf_deint_bfly_neon, export=1
        push            {lr}
        mov             r2,  r1
        add             r2,  r2,  #124*4
        add             r3,  r0,  #124*4
        mov             r12, #64
        mov             lr,  #-16
1:
        vld1.32         {q0},     [r1]!
        vld1.32         {q1},     [r2], lr
        vrev64.32       q2,  q0
        vrev64.32       q3,  q1
        vadd.f32        d3,  d4,  d3
        vadd.f32        d2,  d5,  d2
        vsub.f32        d0,  d0,  d7
        vsub.f32        d1,  d1,  d6
        vst1.32         {q1},     [r3], lr
        vst1.32         {q0},     [r0]!
        subs            r12, r12, #4
        bgt             1b
        pop             {pc}
endfunc
function sbr_hf_gen_neon, export=1
        vld1.32         {q0}, [r2]

        vtrn.32         d0,  d1
        vneg.f32        d18, d1
        vtrn.32         d18, d1
        mov             r2, r3
        ldr             r3, [sp]
        add             r0, r0, r2, lsl #3
        add             r1, r1, r2, lsl #3
        sub             r1, r1, #2*8
        sub             r3, r3, r2
        vld1.32         {q1},     [r1]!
1:
        vld1.32         {q3},     [r1]!
        vrev64.32       q2,  q1
        vmov            q8,  q3
        vrev64.32       d20, d3
        vrev64.32       d21, d6
        vmla.f32        q3,  q1,  d0[0]
        vmla.f32        d6,  d4,  d18
        vmla.f32        d6,  d3,  d0[1]
        vmla.f32        d7,  d20, d18
        vmla.f32        d7,  d16, d0[1]
        vmla.f32        d6,  d5,  d1
        vmla.f32        d7,  d21, d1
        vmov            q1,  q8
        vst1.32         {q3},     [r0]!
        subs            r3,  r3,  #2
        bgt             1b
        bx              lr
endfunc
function sbr_qmf_synthesis_window_neon, export=1
        push {r4,r5,lr}
        mov r5, #64*4
1:
        mov r4, r2
        vld1.32 {q0}, [r1]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q10, q0, q1
        add r12, r1, #192*4
        vld1.32 {q2}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q2, q2, q1
        vadd.f32 q10, q2, q10
        add r12, r1, #256*4
        vld1.32 {q8}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q8, q8, q1
        vadd.f32 q10, q8, q10
        add r12, r1, #448*4
        vld1.32 {q2}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q2, q2, q1
        vadd.f32 q10, q2, q10
        add r12, r1, #512*4
        vld1.32 {q8}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q8, q8, q1
        vadd.f32 q10, q8, q10
        add r12, r1, #704*4
        vld1.32 {q2}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q2, q2, q1
        vadd.f32 q10, q2, q10
        add r12, r1, #768*4
        vld1.32 {q8}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q8, q8, q1
        vadd.f32 q10, q8, q10
        add r12, r1, #960*4
        vld1.32 {q2}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q2, q2, q1
        vadd.f32 q10, q2, q10
        add r12, r1, #1024*4
        vld1.32 {q8}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q8, q8, q1
        vadd.f32 q10, q8, q10
        add r12, r1, #1216*4
        vld1.32 {q2}, [r12]
        vld1.32 {q1}, [r4], r5
        vmul.f32 q2, q2, q1
        vadd.f32 q10, q2, q10
        vst1.32 {q10}, [r0]!

        subs r3, r3, #4
        add r1, r1, #16
        add r2, r2, #16
        beq 2f
        b 1b
2:
        pop {r4,r5,pc}
endfunc

function aacenc_calc_expspec_neon, export=1
        vmov.f32 q3, #0.5
        vmov.i32 q10, #0x7fffffff
1:
        vld1.32 {q0},[r1]!
        vand.f32 q0, q10, q0

        vrsqrte.f32 q1, q0
        vmul.f32 q2, q1, q1
        vmul.f32 q2, q2, q1
        vmul.f32 q2, q2, q0
        vsub.f32 q2, q1, q2
        vmul.f32 q2, q2, q3
        vadd.f32 q1, q1, q2
        vmul.f32 q0, q0, q1

        vrsqrte.f32 q1, q0
        vmul.f32 q2, q1, q1
        vmul.f32 q2, q2, q1
        vmul.f32 q2, q2, q0
        vsub.f32 q2, q1, q2
        vmul.f32 q2, q2, q3
        vadd.f32 q1, q1, q2
        vmul.f32 q1, q0, q1

        vmul.f32 q1, q0, q1
        vst1.32 {q1},[r0]!
        subs r2, r2, #4
        bgt 1b
        bx lr
endfunc
function vorbis_inverse_coupling_neon, export=1
        vmov.i32        q10, #1<<31
        mov             r3,  r0
        mov             r12, r1
1:
        vld1.32         {d2-d3},  [r1]!
        vld1.32         {d0-d1},  [r0]!
        vcle.s32        q8,  q1,  #0
        vand            q9,  q0,  q10
        veor            q1,  q1,  q9
        vand            q2,  q1,  q8
        vbic            q3,  q1,  q8
        vadd.f32        q1,  q0,  q2
        vsub.f32        q0,  q0,  q3
        vst1.32         {d2-d3},  [r3]!
        vst1.32         {d0-d1},  [r12]!
        subs            r2,  r2,  #4
        ble             2f
        b               1b
2:      
        bx              lr
endfunc

function sbr_qmf_deint_neg_neon, export=1
        add             r1,  r1,  #60*4
        add             r2,  r0,  #62*4
        mov             r3,  #-16
        mov             r12, #32
        vmov.i32        d2,  #1<<31
1:
        vld2.32         {d0,d1},  [r1,:128], r3
        veor            d0,  d0,  d2
        vrev64.32       d1,  d1
        vst1.32         {d0},     [r2,:64]
        vst1.32         {d1},     [r0,:64]!
        sub             r2,  r2,  #8
        subs            r12, r12, #2
        bgt             1b
        bx              lr
endfunc
function sbr_autocorrelate_neon, export=1
        vld1.32         {q0},     [r0]!
        vmov.f32        d0,  #0.0
        vmov.f32        q1,  #0.0
        vmov.f32        q3,  #0.0
        vmov.f32        d20, #0.0
        vmul.f32        d21, d1,  d1
        vmov            q8,  q0
        vmov            q11, q0
        mov             r12, #36
1:
        vld1.32         {q2},     [r0]!
        vrev64.32       q12, q2
        vmla.f32        q10, q2,  q2
        vmla.f32        d2,  d1,  d4
        vmla.f32        d3,  d1,  d24
        vmla.f32        d6,  d0,  d4
        vmla.f32        d7,  d0,  d24
        vmla.f32        d2,  d4,  d5
        vmla.f32        d3,  d4,  d25
        vmla.f32        d6,  d1,  d5
        vmla.f32        d7,  d1,  d25
        vmov            q0,  q2
        subs            r12, r12, #2
        bgt             1b
        vld1.32         {q2},     [r0]!
        vrev64.32       q12, q2
        vmla.f32        d2,  d1,  d4
        vmla.f32        d3,  d1,  d24
        vmla.f32        d6,  d0,  d4
        vmla.f32        d7,  d0,  d24
        vadd.f32        d20, d20, d21
        vrev64.32       d18, d17
        vneg.f32        s15, s15
        vmov            d21, d20
        vpadd.f32       d0,  d0,  d2
        vpadd.f32       d7,  d6,  d7
        vtrn.32         d1,  d3
        vsub.f32        d6,  d1,  d3
        vtrn.32         d0,  d6
        vpadd.f32       d20, d20, d21
        vst1.32         {d6},     [r1]!
        vst1.32         {d7[0]},  [r1]!
        vst1.32         {d7[1]},  [r1]!
        vst1.32         {d20[0]}, [r1]
        bx              lr
endfunc

function conv_fltp_to_flt_2ch_neon, export=1
	ldr	r12, [r1, #0]
	ldr	r1, [r1, #4]
1:
    vld1.32 {q0,q1}, [r1]!
    vld1.32 {q2,q3}, [r12]!
    vzip.32 q0, q2
    vzip.32 q1, q3
    vrev64.32 q0, q0
    vrev64.32 q2, q2
    vrev64.32 q1, q1
    vrev64.32 q3, q3
    vswp.32 q1, q2
    vst1.32 {q0,q1}, [r0]!
    vst1.32 {q2,q3}, [r0]!
    subs r2, r2, #8
    beq 2f
    b 1b
2:
    bx lr
endfunc

function conv_flt_to_fltp_2ch_neon, export=1
	ldr	r12, [r0, #0]
	ldr	r0, [r0, #4]
1:
    vld1.32 {q0,q1}, [r1]!
    vtrn.32 d0, d2
    vtrn.32 d1, d3
    vtrn.32 d0, d1
    vtrn.32 d2, d3
    vst1.32 {q0}, [r12]!
    vst1.32 {q1}, [r0]!
    subs r2, r2, #4
    beq 2f
    b 1b
2:
    bx lr
endfunc

