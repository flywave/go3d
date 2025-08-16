package vec3

import (
	"math"
	"testing"
)

const (
	epsilon   = 1e-6
	epsilon32 = 1e-6
)

func TestConstants(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want T
	}{
		{"Zero", Zero, T{0, 0, 0}},
		{"UnitX", UnitX, T{1, 0, 0}},
		{"UnitY", UnitY, T{0, 1, 0}},
		{"UnitZ", UnitZ, T{0, 0, 1}},
		{"UnitXYZ", UnitXYZ, T{1, 1, 1}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if tt.vec != tt.want {
				t.Errorf("%s = %v, want %v", tt.name, tt.vec, tt.want)
			}
		})
	}
}

func TestLength(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want float32
	}{
		{"Zero vector", T{0, 0, 0}, 0},
		{"Unit X", T{1, 0, 0}, 1},
		{"Unit Y", T{0, 1, 0}, 1},
		{"Unit Z", T{0, 0, 1}, 1},
		{"3-4-5 triangle", T{3, 4, 0}, 5},
		{"Arbitrary vector", T{1, 2, 3}, float32(math.Sqrt(14))},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.vec.Length()
			if math.Abs(float64(got-tt.want)) > epsilon {
				t.Errorf("Length() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestLengthSqr(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want float32
	}{
		{"Zero vector", T{0, 0, 0}, 0},
		{"Unit X", T{1, 0, 0}, 1},
		{"3-4-5 triangle", T{3, 4, 0}, 25},
		{"Arbitrary vector", T{1, 2, 3}, 14},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.vec.LengthSqr()
			if got != tt.want {
				t.Errorf("LengthSqr() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestNormalize(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want T
	}{
		{"Already normalized", T{1, 0, 0}, T{1, 0, 0}},
		{"Simple vector", T{2, 0, 0}, T{1, 0, 0}},
		{"Arbitrary vector", T{1, 1, 1}, T{1 / float32(math.Sqrt(3)), 1 / float32(math.Sqrt(3)), 1 / float32(math.Sqrt(3))}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			vec := tt.vec
			vec.Normalize()
			if !vec.AlmostEqual(tt.want, epsilon32) {
				t.Errorf("Normalize() = %v, want %v", vec, tt.want)
			}
		})
	}
}

func TestNormalized(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want T
	}{
		{"Already normalized", T{1, 0, 0}, T{1, 0, 0}},
		{"Simple vector", T{2, 0, 0}, T{1, 0, 0}},
		{"Arbitrary vector", T{1, 1, 1}, T{1 / float32(math.Sqrt(3)), 1 / float32(math.Sqrt(3)), 1 / float32(math.Sqrt(3))}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.vec.Normalized()
			if !got.AlmostEqual(tt.want, epsilon32) {
				t.Errorf("Normalized() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestAdd(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		want T
	}{
		{"Zero vectors", T{0, 0, 0}, T{0, 0, 0}, T{0, 0, 0}},
		{"Add to zero", T{0, 0, 0}, T{1, 2, 3}, T{1, 2, 3}},
		{"Simple addition", T{1, 2, 3}, T{4, 5, 6}, T{5, 7, 9}},
		{"Negative numbers", T{-1, -2, -3}, T{1, 2, 3}, T{0, 0, 0}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Add function
			got := Add(&tt.a, &tt.b)
			if got != tt.want {
				t.Errorf("Add() = %v, want %v", got, tt.want)
			}

			// Test Add method
			vec := tt.a
			vec.Add(&tt.b)
			if vec != tt.want {
				t.Errorf("Add method = %v, want %v", vec, tt.want)
			}
		})
	}
}

func TestSub(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		want T
	}{
		{"Zero vectors", T{0, 0, 0}, T{0, 0, 0}, T{0, 0, 0}},
		{"Subtract from zero", T{0, 0, 0}, T{1, 2, 3}, T{-1, -2, -3}},
		{"Simple subtraction", T{5, 7, 9}, T{1, 2, 3}, T{4, 5, 6}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Sub function
			got := Sub(&tt.a, &tt.b)
			if got != tt.want {
				t.Errorf("Sub() = %v, want %v", got, tt.want)
			}

			// Test Sub method
			vec := tt.a
			vec.Sub(&tt.b)
			if vec != tt.want {
				t.Errorf("Sub method = %v, want %v", vec, tt.want)
			}
		})
	}
}

func TestDot(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		want float32
	}{
		{"Orthogonal vectors", T{1, 0, 0}, T{0, 1, 0}, 0},
		{"Parallel vectors", T{1, 0, 0}, T{1, 0, 0}, 1},
		{"Opposite vectors", T{1, 0, 0}, T{-1, 0, 0}, -1},
		{"Simple dot product", T{1, 2, 3}, T{4, 5, 6}, 32},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := Dot(&tt.a, &tt.b)
			if got != tt.want {
				t.Errorf("Dot() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestCross(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		want T
	}{
		{"Standard basis", T{1, 0, 0}, T{0, 1, 0}, T{0, 0, 1}},
		{"Reverse order", T{0, 1, 0}, T{1, 0, 0}, T{0, 0, -1}},
		{"Same vector", T{1, 2, 3}, T{1, 2, 3}, T{0, 0, 0}},
		{"Arbitrary vectors", T{1, 2, 3}, T{4, 5, 6}, T{-3, 6, -3}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := Cross(&tt.a, &tt.b)
			if got != tt.want {
				t.Errorf("Cross() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestScale(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		f    float32
		want T
	}{
		{"Scale by 1", T{1, 2, 3}, 1, T{1, 2, 3}},
		{"Scale by 0", T{1, 2, 3}, 0, T{0, 0, 0}},
		{"Scale by 2", T{1, 2, 3}, 2, T{2, 4, 6}},
		{"Scale by -1", T{1, 2, 3}, -1, T{-1, -2, -3}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Scale method
			vec := tt.vec
			vec.Scale(tt.f)
			if vec != tt.want {
				t.Errorf("Scale() = %v, want %v", vec, tt.want)
			}

			// Test Scaled method
			got := tt.vec.Scaled(tt.f)
			if got != tt.want {
				t.Errorf("Scaled() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestInvert(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want T
	}{
		{"Zero vector", T{0, 0, 0}, T{0, 0, 0}},
		{"Positive vector", T{1, 2, 3}, T{-1, -2, -3}},
		{"Negative vector", T{-1, -2, -3}, T{1, 2, 3}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Invert method
			vec := tt.vec
			vec.Invert()
			if vec != tt.want {
				t.Errorf("Invert() = %v, want %v", vec, tt.want)
			}

			// Test Inverted method
			got := tt.vec.Inverted()
			if got != tt.vec {
				t.Errorf("Inverted() = %v, want %v", got, tt.vec)
			}
		})
	}
}

func TestAbs(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want T
	}{
		{"All positive", T{1, 2, 3}, T{1, 2, 3}},
		{"All negative", T{-1, -2, -3}, T{1, 2, 3}},
		{"Mixed signs", T{-1, 2, -3}, T{1, 2, 3}},
		{"Zero vector", T{0, 0, 0}, T{0, 0, 0}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Abs method
			vec := tt.vec
			vec.Abs()
			if vec != tt.want {
				t.Errorf("Abs() = %v, want %v", vec, tt.want)
			}

			// Test Absed method
			got := tt.vec.Absed()
			if got != tt.want {
				t.Errorf("Absed() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestDistance(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		want float32
	}{
		{"Same point", T{1, 2, 3}, T{1, 2, 3}, 0},
		{"Simple distance", T{0, 0, 0}, T{1, 0, 0}, 1},
		{"3D distance", T{0, 0, 0}, T{1, 2, 2}, 3},
		{"Negative coordinates", T{-1, -1, -1}, T{1, 1, 1}, float32(math.Sqrt(12))},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Distance function
			got := Distance(&tt.a, &tt.b)
			if math.Abs(float64(got-tt.want)) > epsilon {
				t.Errorf("Distance() = %v, want %v", got, tt.want)
			}

			// Test SquareDistance function
			wantSqr := tt.want * tt.want
			gotSqr := SquareDistance(&tt.a, &tt.b)
			if gotSqr != wantSqr {
				t.Errorf("SquareDistance() = %v, want %v", gotSqr, wantSqr)
			}
		})
	}
}

func TestInterpolate(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		t    float32
		want T
	}{
		{"t=0", T{0, 0, 0}, T{10, 10, 10}, 0, T{0, 0, 0}},
		{"t=1", T{0, 0, 0}, T{10, 10, 10}, 1, T{10, 10, 10}},
		{"t=0.5", T{0, 0, 0}, T{10, 10, 10}, 0.5, T{5, 5, 5}},
		{"t=0.25", T{0, 0, 0}, T{10, 10, 10}, 0.25, T{2.5, 2.5, 2.5}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := Interpolate(&tt.a, &tt.b, tt.t)
			if !got.AlmostEqual(tt.want, epsilon32) {
				t.Errorf("Interpolate() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestClamp(t *testing.T) {
	min := T{0, 0, 0}
	max := T{10, 10, 10}

	tests := []struct {
		name string
		vec  T
		want T
	}{
		{"Within bounds", T{5, 5, 5}, T{5, 5, 5}},
		{"Below min", T{-5, -5, -5}, T{0, 0, 0}},
		{"Above max", T{15, 15, 15}, T{10, 10, 10}},
		{"Mixed bounds", T{-5, 5, 15}, T{0, 5, 10}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Clamp method
			vec := tt.vec
			vec.Clamp(&min, &max)
			if vec != tt.want {
				t.Errorf("Clamp() = %v, want %v", vec, tt.want)
			}

			// Test Clamped method
			got := tt.vec.Clamped(&min, &max)
			if got != tt.want {
				t.Errorf("Clamped() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestClamp01(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want T
	}{
		{"Within bounds", T{0.5, 0.5, 0.5}, T{0.5, 0.5, 0.5}},
		{"Below zero", T{-1, -1, -1}, T{0, 0, 0}},
		{"Above one", T{2, 2, 2}, T{1, 1, 1}},
		{"Mixed bounds", T{-1, 0.5, 2}, T{0, 0.5, 1}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Clamp01 method
			vec := tt.vec
			vec.Clamp01()
			if vec != tt.want {
				t.Errorf("Clamp01() = %v, want %v", vec, tt.want)
			}

			// Test Clamped01 method
			got := tt.vec.Clamped01()
			if got != tt.want {
				t.Errorf("Clamped01() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestAngle(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		want float32
	}{
		{"Same direction", T{1, 0, 0}, T{1, 0, 0}, 0},
		{"Orthogonal", T{1, 0, 0}, T{0, 1, 0}, float32(math.Pi / 2)},
		{"Opposite", T{1, 0, 0}, T{-1, 0, 0}, float32(math.Pi)},
		{"45 degrees", T{1, 0, 0}, T{1, 1, 0}, float32(math.Pi / 4)},
		{"Arbitrary angle", T{1, 2, 3}, T{4, 5, 6}, float32(math.Acos(32 / (math.Sqrt(14) * math.Sqrt(77))))},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := Angle(&tt.a, &tt.b)
			if math.Abs(float64(got-tt.want)) > epsilon {
				t.Errorf("Angle() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestMinMax(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		min  T
		max  T
	}{
		{"Simple min/max", T{1, 5, 3}, T{4, 2, 6}, T{1, 2, 3}, T{4, 5, 6}},
		{"Negative numbers", T{-1, -5, 3}, T{-4, 2, -6}, T{-4, -5, -6}, T{-1, 2, 3}},
		{"Equal vectors", T{1, 2, 3}, T{1, 2, 3}, T{1, 2, 3}, T{1, 2, 3}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			min := Min(&tt.a, &tt.b)
			if min != tt.min {
				t.Errorf("Min() = %v, want %v", min, tt.min)
			}

			max := Max(&tt.a, &tt.b)
			if max != tt.max {
				t.Errorf("Max() = %v, want %v", max, tt.max)
			}
		})
	}
}

func TestMul(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		want T
	}{
		{"Simple multiplication", T{1, 2, 3}, T{4, 5, 6}, T{4, 10, 18}},
		{"With zero", T{1, 2, 3}, T{0, 1, 2}, T{0, 2, 6}},
		{"Negative numbers", T{-1, 2, -3}, T{2, -3, 4}, T{-2, -6, -12}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test Mul function
			got := Mul(&tt.a, &tt.b)
			if got != tt.want {
				t.Errorf("Mul() = %v, want %v", got, tt.want)
			}

			// Test Mul method
			vec := tt.a
			vec.Mul(&tt.b)
			if vec != tt.want {
				t.Errorf("Mul method = %v, want %v", vec, tt.want)
			}
		})
	}
}

func TestAlmostEqual(t *testing.T) {
	tests := []struct {
		name string
		a    T
		b    T
		tol  float32
		want bool
	}{
		{"Exactly equal", T{1, 2, 3}, T{1, 2, 3}, 1e-6, true},
		{"Within tolerance", T{1, 2, 3}, T{1.0001, 2.0001, 3.0001}, 1e-3, true},
		{"Outside tolerance", T{1, 2, 3}, T{1.1, 2.1, 3.1}, 1e-3, false},
		{"Different components", T{1, 2, 3}, T{4, 5, 6}, 1e-6, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.a.AlmostEqual(tt.b, tt.tol)
			if got != tt.want {
				t.Errorf("AlmostEqual() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestLerp(t *testing.T) {
	tests := []struct {
		name  string
		v     T
		other T
		t     float32
		want  T
	}{
		{"t=0", T{0, 0, 0}, T{10, 10, 10}, 0, T{0, 0, 0}},
		{"t=1", T{0, 0, 0}, T{10, 10, 10}, 1, T{10, 10, 10}},
		{"t=0.5", T{0, 0, 0}, T{10, 10, 10}, 0.5, T{5, 5, 5}},
		{"t=0.25", T{0, 0, 0}, T{10, 10, 10}, 0.25, T{2.5, 2.5, 2.5}},
		{"Negative t", T{0, 0, 0}, T{10, 10, 10}, -0.5, T{-5, -5, -5}},
		{"t>1", T{0, 0, 0}, T{10, 10, 10}, 1.5, T{15, 15, 15}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.v.Lerp(&tt.other, tt.t)
			if !got.AlmostEqual(tt.want, epsilon32) {
				t.Errorf("Lerp() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestRotated(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		axis T
		rad  float32
		want T
	}{
		{"Rotate around Z axis 90 degrees", T{1, 0, 0}, T{0, 0, 1}, float32(math.Pi / 2), T{0, 1, 0}},
		{"Rotate around X axis 90 degrees", T{0, 1, 0}, T{1, 0, 0}, float32(math.Pi / 2), T{0, 0, 1}},
		{"Rotate around Y axis 90 degrees", T{1, 0, 0}, T{0, 1, 0}, float32(math.Pi / 2), T{0, 0, -1}},
		{"Rotate 180 degrees", T{1, 0, 0}, T{0, 0, 1}, float32(math.Pi), T{-1, 0, 0}},
		{"Rotate 360 degrees", T{1, 0, 0}, T{0, 0, 1}, float32(2 * math.Pi), T{1, 0, 0}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := Rotated(&tt.vec, &tt.axis, tt.rad)
			if !got.AlmostEqual(tt.want, 1e-3) {
				t.Errorf("Rotated() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestNormal(t *testing.T) {
	tests := []struct {
		name string
		vec  T
		want T
	}{
		{"X axis", T{1, 0, 0}, T{0, 1, 0}},
		{"Y axis", T{0, 1, 0}, T{1, 0, 0}},
		{"Z axis", T{0, 0, 1}, T{1, 0, 0}},
		{"Arbitrary vector", T{1, 1, 1}, T{0.70710677, -0.70710677, 0}},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.vec.Normal()
			// 检查是否为正交向量（点积接近0）
			if math.Abs(float64(Dot(&got, &tt.vec))) > 1e-6 {
				t.Errorf("Normal() = %v, not orthogonal to %v", got, tt.vec)
			}
			// 检查是否为标准化向量
			if math.Abs(float64(got.Length()-1)) > 1e-6 {
				t.Errorf("Normal() = %v, not normalized (length=%v)", got, got.Length())
			}
		})
	}
}

func TestParse(t *testing.T) {
	tests := []struct {
		name    string
		input   string
		want    T
		wantErr bool
	}{
		{"Valid input", "1 2 3", T{1, 2, 3}, false},
		{"With decimals", "1.5 2.5 3.5", T{1.5, 2.5, 3.5}, false},
		{"Negative numbers", "-1 -2 -3", T{-1, -2, -3}, false},
		{"Invalid format", "1 2", T{0, 0, 0}, true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := Parse(tt.input)
			if (err != nil) != tt.wantErr {
				t.Errorf("Parse() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !tt.wantErr && got != tt.want {
				t.Errorf("Parse() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestString(t *testing.T) {
	vec := T{1, 2, 3}
	got := vec.String()
	want := "1 2 3"
	if got != want {
		t.Errorf("String() = %v, want %v", got, want)
	}
}

func TestBoxIntersection(t *testing.T) {
	bb1 := Box{T{0, 0, 0}, T{1, 1, 1}}
	bb2 := Box{T{1, 1, 1}, T{2, 2, 2}}
	if !bb1.Intersects(&bb2) {
		t.Fail()
	}

	bb3 := Box{T{1, 2, 1}, T{2, 3, 2}}
	if bb1.Intersects(&bb3) {
		t.Fail()
	}
}
