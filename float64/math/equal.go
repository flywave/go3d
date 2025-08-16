package math

import "math"

func AlmostEqual64(a, b, tol float64) bool {
	return math.Abs((a - b)) < (tol)
}
