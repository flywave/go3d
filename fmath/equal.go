package fmath

import "math"

func AlmostEqual32(a, b, tol float32) bool {
	return math.Abs(float64(a-b)) < float64(tol)
}
