// Package vec2 contains a 2D float64 vector type T and functions.
package vec2

import (
	"fmt"
	"math"

	"github.com/flywave/go3d/float64/generic"
)

var (
	// Zero holds a zero vector.
	Zero = T{}

	// UnitX holds a vector with X set to one.
	UnitX = T{1, 0}
	// UnitY holds a vector with Y set to one.
	UnitY = T{0, 1}
	// UnitXY holds a vector with X and Y set to one.
	UnitXY = T{1, 1}

	// MinVal holds a vector with the smallest possible component values.
	MinVal = T{-math.MaxFloat64, -math.MaxFloat64}
	// MaxVal holds a vector with the highest possible component values.
	MaxVal = T{+math.MaxFloat64, +math.MaxFloat64}
)

// T represents a 2D vector.
type T [2]float64

// From copies a T from a generic.T implementation.
func From(other generic.T) T {
	return T{other.Get(0, 0), other.Get(0, 1)}
}

// Parse parses T from a string. See also String()
func Parse(s string) (r T, err error) {
	_, err = fmt.Sscan(s, &r[0], &r[1])
	return r, err
}

// String formats T as string. See also Parse().
func (vec *T) String() string {
	return fmt.Sprint(vec[0], vec[1])
}

// Rows returns the number of rows of the vector.
func (vec *T) Rows() int {
	return 2
}

// Cols returns the number of columns of the vector.
func (vec *T) Cols() int {
	return 1
}

// Size returns the number elements of the vector.
func (vec *T) Size() int {
	return 2
}

// Slice returns the elements of the vector as slice.
func (vec *T) Slice() []float64 {
	return vec[:]
}

// Get returns one element of the vector.
func (vec *T) Get(col, row int) float64 {
	return vec[row]
}

// IsZero checks if all elements of the vector are zero.
func (vec *T) IsZero() bool {
	return vec[0] == 0 && vec[1] == 0
}

// Length returns the length of the vector.
// See also LengthSqr and Normalize.
func (vec *T) Length() float64 {
	return math.Hypot(vec[0], vec[1])
}

// LengthSqr returns the squared length of the vector.
// See also Length and Normalize.
func (vec *T) LengthSqr() float64 {
	return vec[0]*vec[0] + vec[1]*vec[1]
}

// Scale multiplies all element of the vector by f and returns vec.
func (vec *T) Scale(f float64) *T {
	vec[0] *= f
	vec[1] *= f
	return vec
}

// Scaled returns a copy of vec with all elements multiplies by f.
func (vec *T) Scaled(f float64) T {
	return T{vec[0] * f, vec[1] * f}
}

// Invert inverts the vector.
func (vec *T) Invert() *T {
	vec[0] = -vec[0]
	vec[1] = -vec[1]
	return vec
}

// Inverted returns an inverted copy of the vector.
func (vec *T) Inverted() T {
	return T{-vec[0], -vec[1]}
}

// Normalize normalizes the vector to unit length.
func (vec *T) Normalize() *T {
	sl := vec.LengthSqr()
	if sl == 0 || sl == 1 {
		return vec
	}
	return vec.Scale(1 / math.Sqrt(sl))
}

// Normalized returns a unit length normalized copy of the vector.
func (vec *T) Normalized() T {
	v := *vec
	v.Normalize()
	return v
}

// Add adds another vector to vec.
func (vec *T) Add(v *T) *T {
	vec[0] += v[0]
	vec[1] += v[1]
	return vec
}

// Sub subtracts another vector from vec.
func (vec *T) Sub(v *T) *T {
	vec[0] -= v[0]
	vec[1] -= v[1]
	return vec
}

// Mul multiplies the components of the vector with the respective components of v.
func (vec *T) Mul(v *T) *T {
	vec[0] *= v[0]
	vec[1] *= v[1]
	return vec
}

// Rotate rotates the vector counter-clockwise by angle.
func (vec *T) Rotate(angle float64) *T {
	*vec = vec.Rotated(angle)
	return vec
}

// Rotated returns a counter-clockwise rotated copy of the vector.
func (vec *T) Rotated(angle float64) T {
	sinus := math.Sin(angle)
	cosinus := math.Cos(angle)
	return T{
		vec[0]*cosinus - vec[1]*sinus,
		vec[0]*sinus + vec[1]*cosinus,
	}
}

// RotateAroundPoint rotates the vector counter-clockwise around a point.
func (vec *T) RotateAroundPoint(point *T, angle float64) *T {
	return vec.Sub(point).Rotate(angle).Add(point)
}

// Rotate90DegLeft rotates the vector 90 degrees left (counter-clockwise).
func (vec *T) Rotate90DegLeft() *T {
	temp := vec[0]
	vec[0] = -vec[1]
	vec[1] = temp
	return vec
}

// Rotate90DegRight rotates the vector 90 degrees right (clockwise).
func (vec *T) Rotate90DegRight() *T {
	temp := vec[0]
	vec[0] = vec[1]
	vec[1] = -temp
	return vec
}

// Angle returns the counter-clockwise angle of the vector from the x axis.
func (vec *T) Angle() float64 {
	return math.Atan2(vec[1], vec[0])
}

func (v *T) Lerp(other *T, t float64) T {
	return T{
		v[0] + t*(other[0]-v[0]),
		v[1] + t*(other[1]-v[1]),
	}
}

// Squared Distance between two vectors
func SquareDistance(a, b *T) float64 {
	d := Sub(a, b)
	return d.LengthSqr()
}

// Distance between two vectors
func Distance(a, b *T) float64 {
	d := Sub(a, b)
	return d.Length()
}

// Add returns the sum of two vectors.
func Add(a, b *T) T {
	return T{a[0] + b[0], a[1] + b[1]}
}

// Sub returns the difference of two vectors.
func Sub(a, b *T) T {
	return T{a[0] - b[0], a[1] - b[1]}
}

// Mul returns the component wise product of two vectors.
func Mul(a, b *T) T {
	return T{a[0] * b[0], a[1] * b[1]}
}

// Dot returns the dot product of two vectors.
func Dot(a, b *T) float64 {
	return a[0]*b[0] + a[1]*b[1]
}

// Cross returns the cross product of two vectors.
func Cross(a, b *T) T {
	return T{
		a[1]*b[0] - a[0]*b[1],
		a[0]*b[1] - a[1]*b[0],
	}
}

// Angle returns the angle between two vectors.
func Angle(a, b *T) float64 {
	v := Dot(a, b) / (a.Length() * b.Length())
	// prevent NaN
	if v > 1. {
		v = v - 2
	} else if v < -1. {
		v = v + 2
	}
	return math.Acos(v)
}

// IsLeftWinding returns if the angle from a to b is left winding.
func IsLeftWinding(a, b *T) bool {
	ab := b.Rotated(-a.Angle())
	return ab.Angle() > 0
}

// IsRightWinding returns if the angle from a to b is right winding.
func IsRightWinding(a, b *T) bool {
	ab := b.Rotated(-a.Angle())
	return ab.Angle() < 0
}

// Min returns the component wise minimum of two vectors.
func Min(a, b *T) T {
	min := *a
	if b[0] < min[0] {
		min[0] = b[0]
	}
	if b[1] < min[1] {
		min[1] = b[1]
	}
	return min
}

// Max returns the component wise maximum of two vectors.
func Max(a, b *T) T {
	max := *a
	if b[0] > max[0] {
		max[0] = b[0]
	}
	if b[1] > max[1] {
		max[1] = b[1]
	}
	return max
}

// Interpolate interpolates between a and b at t (0,1).
func Interpolate(a, b *T, t float64) T {
	t1 := 1 - t
	return T{
		a[0]*t1 + b[0]*t,
		a[1]*t1 + b[1]*t,
	}
}

// Clamp clamps the vector's components to be in the range of min to max.
func (vec *T) Clamp(min, max *T) *T {
	for i := range vec {
		if vec[i] < min[i] {
			vec[i] = min[i]
		} else if vec[i] > max[i] {
			vec[i] = max[i]
		}
	}
	return vec
}

// Clamped returns a copy of the vector with the components clamped to be in the range of min to max.
func (vec *T) Clamped(min, max *T) T {
	result := *vec
	result.Clamp(min, max)
	return result
}

// Clamp01 clamps the vector's components to be in the range of 0 to 1.
func (vec *T) Clamp01() *T {
	return vec.Clamp(&Zero, &UnitXY)
}

// Clamped01 returns a copy of the vector with the components clamped to be in the range of 0 to 1.
func (vec *T) Clamped01() T {
	result := *vec
	result.Clamp01()
	return result
}

func (vec *T) SetMin(c T) {
	if c[0] < vec[0] {
		vec[0] = c[0]
	}
	if c[1] < vec[1] {
		vec[1] = c[1]
	}
}
func (vec *T) SetMax(c T) {
	if c[0] > vec[0] {
		vec[0] = c[0]
	}
	if c[1] > vec[1] {
		vec[1] = c[1]
	}
}

func PointSegmentDistance(p1 *T, x1 *T, x2 *T) float64 {
	v := Sub(x2, x1)
	w := Sub(p1, x1)
	if w.IsZero() {
		return 0
	}

	vn := v.Normalized()
	wn := w.Normalized()
	cn := Dot(&vn, &wn)
	if math.Abs(math.Abs(cn)-1) < 1e-8 {
		return 0
	}
	return w.Length() * math.Sqrt((1 - cn*cn))
}

func PointSegmentVerticalPoint(p1 *T, x1 *T, x2 *T) *T {
	v := Sub(x2, x1)
	w := Sub(p1, x1)
	if w.IsZero() {
		return &T{p1[0], p1[1]}
	}

	vn := v.Normalized()
	wn := w.Normalized()
	cn := Dot(&vn, &wn)
	if math.Abs(math.Abs(cn)-1) < 1e-8 {
		return &T{p1[0], p1[1]}
	}
	c1 := Dot(&w, &v)
	c2 := Dot(&v, &v)
	b := c1 / c2
	v1 := v.Scaled(b)
	res := Add(x1, &v1)
	return &res
}
