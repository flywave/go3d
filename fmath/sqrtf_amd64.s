// +build !noasm,!gccgo,!safe

#include "textflag.h"

// func Sqrtf(x float32) float32
TEXT ·Sqrtf(SB),NOSPLIT,$0
	SQRTSS x+0(FP), X0
	MOVSS X0, ret+8(FP)
	RET
