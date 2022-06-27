// Package mat4 contains a 4x4 float32 matrix type T and functions.
package mat4

import (
	"fmt"

	math "github.com/flywave/go3d/fmath"

	"github.com/flywave/go3d/generic"
	"github.com/flywave/go3d/mat2"
	"github.com/flywave/go3d/mat3"
	"github.com/flywave/go3d/quaternion"
	"github.com/flywave/go3d/vec3"
	"github.com/flywave/go3d/vec4"
)

var (
	// Zero holds a zero matrix.
	Zero = T{}

	// Ident holds an ident matrix.
	Ident = T{
		vec4.T{1, 0, 0, 0},
		vec4.T{0, 1, 0, 0},
		vec4.T{0, 0, 1, 0},
		vec4.T{0, 0, 0, 1},
	}
)

// T represents a 4x4 matrix as 4 column vectors.
type T [4]vec4.T

// From copies a T from a generic.T implementation.
func From(other generic.T) T {
	r := Ident
	cols := other.Cols()
	rows := other.Rows()
	if !((cols == 2 && rows == 2) || (cols == 3 && rows == 3) || (cols == 4 && rows == 4)) {
		panic("Unsupported type")
	}
	for col := 0; col < cols; col++ {
		for row := 0; row < rows; row++ {
			r[col][row] = other.Get(col, row)
		}
	}
	return r
}

func FromArray(mtw [16]float32) T {
	mt := T{
		vec4.T{mtw[0], mtw[1], mtw[2], mtw[3]},
		vec4.T{mtw[4], mtw[5], mtw[6], mtw[7]},
		vec4.T{mtw[8], mtw[9], mtw[10], mtw[11]},
		vec4.T{mtw[12], mtw[13], mtw[14], mtw[15]},
	}
	return mt
}

// Parse parses T from a string. See also String()
func Parse(s string) (r T, err error) {
	_, err = fmt.Sscan(s,
		&r[0][0], &r[0][1], &r[0][2], &r[0][3],
		&r[1][0], &r[1][1], &r[1][2], &r[1][3],
		&r[2][0], &r[2][1], &r[2][2], &r[2][3],
		&r[3][0], &r[3][1], &r[3][2], &r[3][3],
	)
	return r, err
}

// String formats T as string. See also Parse().
func (mat *T) String() string {
	return fmt.Sprintf("%s %s %s %s", mat[0].String(), mat[1].String(), mat[2].String(), mat[3].String())
}

// Rows returns the number of rows of the matrix.
func (mat *T) Rows() int {
	return 4
}

// Cols returns the number of columns of the matrix.
func (mat *T) Cols() int {
	return 4
}

// Size returns the number elements of the matrix.
func (mat *T) Size() int {
	return 16
}

// Slice returns the elements of the matrix as slice.
func (mat *T) Slice() []float32 {
	return mat.Array()[:]
}

// Get returns one element of the matrix.
func (mat *T) Get(col, row int) float32 {
	return mat[col][row]
}

// IsZero checks if all elements of the matrix are zero.
func (mat *T) IsZero() bool {
	return *mat == Zero
}

// Scale multiplies the diagonal scale elements by f returns mat.
func (mat *T) Scale(f float32) *T {
	mat[0][0] *= f
	mat[1][1] *= f
	mat[2][2] *= f
	return mat
}

// Scaled returns a copy of the matrix with the diagonal scale elements multiplied by f.
func (mat *T) Scaled(f float32) T {
	result := *mat
	result.Scale(f)
	return result
}

// Mul multiplies every element by f and returns mat.
func (mat *T) Mul(f float32) *T {
	for i, col := range mat {
		for j := range col {
			mat[i][j] *= f
		}
	}
	return mat
}

// Muled returns a copy of the matrix with every element multiplied by f.
func (mat *T) Muled(f float32) T {
	result := *mat
	result.Mul(f)
	return result
}

// Mult multiplies this matrix with the given matrix m and saves the result in this matrix.
func (mat *T) MultMatrix(m *T) *T {
	// iterate over the rows of mat
	for i := range mat {
		row := vec4.T{mat[0][i], mat[1][i], mat[2][i], mat[3][i]}
		mat[0][i] = vec4.Dot4(&row, &m[0])
		mat[1][i] = vec4.Dot4(&row, &m[1])
		mat[2][i] = vec4.Dot4(&row, &m[2])
		mat[3][i] = vec4.Dot4(&row, &m[3])
	}
	return mat
}

// Trace returns the trace value for the matrix.
func (mat *T) Trace() float32 {
	return mat[0][0] + mat[1][1] + mat[2][2] + mat[3][3]
}

// Trace3 returns the trace value for the 3x3 sub-matrix.
func (mat *T) Trace3() float32 {
	return mat[0][0] + mat[1][1] + mat[2][2]
}

// AssignMat2x2 assigns a 2x2 sub-matrix and sets the rest of the matrix to the ident value.
func (mat *T) AssignMat2x2(m *mat2.T) *T {
	*mat = T{
		vec4.T{m[0][0], m[1][0], 0, 0},
		vec4.T{m[0][1], m[1][1], 0, 0},
		vec4.T{0, 0, 1, 0},
		vec4.T{0, 0, 0, 1},
	}
	return mat
}

// AssignMat3x3 assigns a 3x3 sub-matrix and sets the rest of the matrix to the ident value.
func (mat *T) AssignMat3x3(m *mat3.T) *T {
	*mat = T{
		vec4.T{m[0][0], m[1][0], m[2][0], 0},
		vec4.T{m[0][1], m[1][1], m[2][1], 0},
		vec4.T{m[0][2], m[1][2], m[2][2], 0},
		vec4.T{0, 0, 0, 1},
	}
	return mat
}

// AssignMul multiplies a and b and assigns the result to mat.
func (mat *T) AssignMul(a, b *T) *T {
	mat[0] = a.MulVec4(&b[0])
	mat[1] = a.MulVec4(&b[1])
	mat[2] = a.MulVec4(&b[2])
	mat[3] = a.MulVec4(&b[3])
	return mat
}

// MulVec4 multiplies v with mat and returns a new vector v' = M * v.
func (mat *T) MulVec4(v *vec4.T) vec4.T {
	return vec4.T{
		mat[0][0]*v[0] + mat[1][0]*v[1] + mat[2][0]*v[2] + mat[3][0]*v[3],
		mat[0][1]*v[0] + mat[1][1]*v[1] + mat[2][1]*v[2] + mat[3][1]*v[3],
		mat[0][2]*v[0] + mat[1][2]*v[1] + mat[2][2]*v[2] + mat[3][2]*v[3],
		mat[0][3]*v[0] + mat[1][3]*v[1] + mat[2][3]*v[2] + mat[3][3]*v[3],
	}
}

// TransformVec4 multiplies v with mat and saves the result in v.
func (mat *T) TransformVec4(v *vec4.T) {
	// Use intermediate variables to not alter further computations.
	x := mat[0][0]*v[0] + mat[1][0]*v[1] + mat[2][0]*v[2] + mat[3][0]*v[3]
	y := mat[0][1]*v[0] + mat[1][1]*v[1] + mat[2][1]*v[2] + mat[3][1]*v[3]
	z := mat[0][2]*v[0] + mat[1][2]*v[1] + mat[2][2]*v[2] + mat[3][2]*v[3]
	v[3] = mat[0][3]*v[0] + mat[1][3]*v[1] + mat[2][3]*v[2] + mat[3][3]*v[3]
	v[0] = x
	v[1] = y
	v[2] = z
}

// MulVec3 multiplies v (converted to a vec4 as (v_1, v_2, v_3, 1))
// with mat and divides the result by w. Returns a new vec3.
func (mat *T) MulVec3(v *vec3.T) vec3.T {
	v4 := vec4.FromVec3(v)
	v4 = mat.MulVec4(&v4)
	return v4.Vec3DividedByW()
}

// TransformVec3 multiplies v (converted to a vec4 as (v_1, v_2, v_3, 1))
// with mat, divides the result by w and saves the result in v.
func (mat *T) TransformVec3(v *vec3.T) {
	x := mat[0][0]*v[0] + mat[1][0]*v[1] + mat[2][0]*v[2] + mat[3][0]
	y := mat[0][1]*v[0] + mat[1][1]*v[1] + mat[2][1]*v[2] + mat[3][1]
	z := mat[0][2]*v[0] + mat[1][2]*v[1] + mat[2][2]*v[2] + mat[3][2]
	w := mat[0][3]*v[0] + mat[1][3]*v[1] + mat[2][3]*v[2] + mat[3][3]
	oow := 1 / w
	v[0] = x * oow
	v[1] = y * oow
	v[2] = z * oow
}

// MulVec3W multiplies v with mat with w as fourth component of the vector.
// Useful to differentiate between vectors (w = 0) and points (w = 1)
// without transforming them to vec4.
func (mat *T) MulVec3W(v *vec3.T, w float32) vec3.T {
	result := *v
	mat.TransformVec3W(&result, w)
	return result
}

// TransformVec3W multiplies v with mat with w as fourth component of the vector and
// saves the result in v.
// Useful to differentiate between vectors (w = 0) and points (w = 1)
// without transforming them to vec4.
func (mat *T) TransformVec3W(v *vec3.T, w float32) {
	// use intermediate variables to not alter further computations
	x := mat[0][0]*v[0] + mat[1][0]*v[1] + mat[2][0]*v[2] + mat[3][0]*w
	y := mat[0][1]*v[0] + mat[1][1]*v[1] + mat[2][1]*v[2] + mat[3][1]*w
	v[2] = mat[0][2]*v[0] + mat[1][2]*v[1] + mat[2][2]*v[2] + mat[3][2]*w
	v[0] = x
	v[1] = y
}

// SetTranslation sets the translation elements of the matrix.
func (mat *T) SetTranslation(v *vec3.T) *T {
	mat[3][0] = v[0]
	mat[3][1] = v[1]
	mat[3][2] = v[2]
	return mat
}

// Translate adds v to the translation part of the matrix.
func (mat *T) Translate(v *vec3.T) *T {
	mat[3][0] += v[0]
	mat[3][1] += v[1]
	mat[3][2] += v[2]
	return mat
}

// TranslateX adds dx to the X-translation element of the matrix.
func (mat *T) TranslateX(dx float32) *T {
	mat[3][0] += dx
	return mat
}

// TranslateY adds dy to the Y-translation element of the matrix.
func (mat *T) TranslateY(dy float32) *T {
	mat[3][1] += dy
	return mat
}

// TranslateZ adds dz to the Z-translation element of the matrix.
func (mat *T) TranslateZ(dz float32) *T {
	mat[3][2] += dz
	return mat
}

// Scaling returns the scaling diagonal of the matrix.
func (mat *T) Scaling() vec4.T {
	return vec4.T{mat[0][0], mat[1][1], mat[2][2], mat[3][3]}
}

// SetScaling sets the scaling diagonal of the matrix.
func (mat *T) SetScaling(s *vec4.T) *T {
	mat[0][0] = s[0]
	mat[1][1] = s[1]
	mat[2][2] = s[2]
	mat[3][3] = s[3]
	return mat
}

// ScaleVec3 multiplies the scaling diagonal of the matrix by s.
func (mat *T) ScaleVec3(s *vec3.T) *T {
	mat[0][0] *= s[0]
	mat[1][1] *= s[1]
	mat[2][2] *= s[2]
	return mat
}

func (mat *T) Quaternion() quaternion.T {
	// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
	// assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)
	m11, m12, m13 := mat[0][0], mat[1][0], mat[2][0]
	m21, m22, m23 := mat[0][1], mat[1][1], mat[2][1]
	m31, m32, m33 := mat[0][2], mat[1][2], mat[2][2]

	trace := m11 + m22 + m33
	var s, _w, _x, _y, _z float32

	if trace > 0 {
		s = 0.5 / math.Sqrt(trace+1.0)

		_w = 0.25 / s
		_x = (m32 - m23) * s
		_y = (m13 - m31) * s
		_z = (m21 - m12) * s
	} else if m11 > m22 && m11 > m33 {
		s = 2.0 * math.Sqrt(1.0+m11-m22-m33)

		_w = (m32 - m23) / s
		_x = 0.25 * s
		_y = (m12 + m21) / s
		_z = (m13 + m31) / s
	} else if m22 > m33 {
		s = 2.0 * math.Sqrt(1.0+m22-m11-m33)

		_w = (m13 - m31) / s
		_x = (m12 + m21) / s
		_y = 0.25 * s
		_z = (m23 + m32) / s
	} else {
		s = 2.0 * math.Sqrt(1.0+m33-m11-m22)

		_w = (m21 - m12) / s
		_x = (m13 + m31) / s
		_y = (m23 + m32) / s
		_z = 0.25 * s
	}
	return quaternion.T{
		_x,
		_y,
		_z,
		_w,
	}
}

// AssignQuaternion assigns a quaternion to the rotations part of the matrix and sets the other elements to their ident value.
func (mat *T) AssignQuaternion(q *quaternion.T) *T {
	xx := q[0] * q[0] * 2
	yy := q[1] * q[1] * 2
	zz := q[2] * q[2] * 2
	xy := q[0] * q[1] * 2
	xz := q[0] * q[2] * 2
	yz := q[1] * q[2] * 2
	wx := q[3] * q[0] * 2
	wy := q[3] * q[1] * 2
	wz := q[3] * q[2] * 2

	mat[0][0] = 1 - (yy + zz)
	mat[1][0] = xy - wz
	mat[2][0] = xz + wy
	mat[3][0] = 0

	mat[0][1] = xy + wz
	mat[1][1] = 1 - (xx + zz)
	mat[2][1] = yz - wx
	mat[3][1] = 0

	mat[0][2] = xz - wy
	mat[1][2] = yz + wx
	mat[2][2] = 1 - (xx + yy)
	mat[3][2] = 0

	mat[0][3] = 0
	mat[1][3] = 0
	mat[2][3] = 0
	mat[3][3] = 1

	return mat
}

// AssignXRotation assigns a rotation around the x axis to the rotation part of the matrix and sets the remaining elements to their ident value.
func (mat *T) AssignXRotation(angle float32) *T {
	cosine := math.Cos(angle)
	sine := math.Sin(angle)

	mat[0][0] = 1
	mat[1][0] = 0
	mat[2][0] = 0
	mat[3][0] = 0

	mat[0][1] = 0
	mat[1][1] = cosine
	mat[2][1] = -sine
	mat[3][1] = 0

	mat[0][2] = 0
	mat[1][2] = sine
	mat[2][2] = cosine
	mat[3][2] = 0

	mat[0][3] = 0
	mat[1][3] = 0
	mat[2][3] = 0
	mat[3][3] = 1

	return mat
}

// AssignYRotation assigns a rotation around the y axis to the rotation part of the matrix and sets the remaining elements to their ident value.
func (mat *T) AssignYRotation(angle float32) *T {
	cosine := math.Cos(angle)
	sine := math.Sin(angle)

	mat[0][0] = cosine
	mat[1][0] = 0
	mat[2][0] = sine
	mat[3][0] = 0

	mat[0][1] = 0
	mat[1][1] = 1
	mat[2][1] = 0
	mat[3][1] = 0

	mat[0][2] = -sine
	mat[1][2] = 0
	mat[2][2] = cosine
	mat[3][2] = 0

	mat[0][3] = 0
	mat[1][3] = 0
	mat[2][3] = 0
	mat[3][3] = 1

	return mat
}

// AssignZRotation assigns a rotation around the z axis to the rotation part of the matrix and sets the remaining elements to their ident value.
func (mat *T) AssignZRotation(angle float32) *T {
	cosine := math.Cos(angle)
	sine := math.Sin(angle)

	mat[0][0] = cosine
	mat[1][0] = -sine
	mat[2][0] = 0
	mat[3][0] = 0

	mat[0][1] = sine
	mat[1][1] = cosine
	mat[2][1] = 0
	mat[3][1] = 0

	mat[0][2] = 0
	mat[1][2] = 0
	mat[2][2] = 1
	mat[3][2] = 0

	mat[0][3] = 0
	mat[1][3] = 0
	mat[2][3] = 0
	mat[3][3] = 1

	return mat
}

// AssignCoordinateSystem assigns the rotation of a orthogonal coordinates system to the rotation part of the matrix and sets the remaining elements to their ident value.
func (mat *T) AssignCoordinateSystem(x, y, z *vec3.T) *T {
	mat[0][0] = x[0]
	mat[1][0] = x[1]
	mat[2][0] = x[2]
	mat[3][0] = 0

	mat[0][1] = y[0]
	mat[1][1] = y[1]
	mat[2][1] = y[2]
	mat[3][1] = 0

	mat[0][2] = z[0]
	mat[1][2] = z[1]
	mat[2][2] = z[2]
	mat[3][2] = 0

	mat[0][3] = 0
	mat[1][3] = 0
	mat[2][3] = 0
	mat[3][3] = 1

	return mat
}

// AssignEulerRotation assigns Euler angle rotations to the rotation part of the matrix and sets the remaining elements to their ident value.
func (mat *T) AssignEulerRotation(yHead, xPitch, zRoll float32) *T {
	sinH := math.Sin(yHead)
	cosH := math.Cos(yHead)
	sinP := math.Sin(xPitch)
	cosP := math.Cos(xPitch)
	sinR := math.Sin(zRoll)
	cosR := math.Cos(zRoll)

	mat[0][0] = cosR*cosH - sinR*sinP*sinH
	mat[1][0] = -sinR * cosP
	mat[2][0] = cosR*sinH + sinR*sinP*cosH
	mat[3][0] = 0

	mat[0][1] = sinR*cosH + cosR*sinP*sinH
	mat[1][1] = cosR * cosP
	mat[2][1] = sinR*sinH - cosR*sinP*cosH
	mat[3][1] = 0

	mat[0][2] = -cosP * sinH
	mat[1][2] = sinP
	mat[2][2] = cosP * cosH
	mat[3][2] = 0

	mat[0][3] = 0
	mat[1][3] = 0
	mat[2][3] = 0
	mat[3][3] = 1

	return mat
}

// ExtractEulerAngles extracts the rotation part of the matrix as Euler angle rotation values.
func (mat *T) ExtractEulerAngles() (yHead, xPitch, zRoll float32) {
	xPitch = math.Asin(mat[1][2])
	f12 := math.Abs(mat[1][2])
	if f12 > (1.0-0.0001) && f12 < (1.0+0.0001) { // f12 == 1.0
		yHead = 0.0
		zRoll = math.Atan2(mat[0][1], mat[0][0])
	} else {
		yHead = math.Atan2(-mat[0][2], mat[2][2])
		zRoll = math.Atan2(-mat[1][0], mat[1][1])
	}
	return yHead, xPitch, zRoll
}

// AssignPerspectiveProjection assigns a perspective projection transformation.
func (mat *T) AssignPerspectiveProjection(left, right, bottom, top, znear, zfar float32) *T {
	near2 := znear + znear
	ooFarNear := 1 / (zfar - znear)

	mat[0][0] = near2 / (right - left)
	mat[1][0] = 0
	mat[2][0] = (right + left) / (right - left)
	mat[3][0] = 0

	mat[0][1] = 0
	mat[1][1] = near2 / (top - bottom)
	mat[2][1] = (top + bottom) / (top - bottom)
	mat[3][1] = 0

	mat[0][2] = 0
	mat[1][2] = 0
	mat[2][2] = -(zfar + znear) * ooFarNear
	mat[3][2] = -2 * zfar * znear * ooFarNear

	mat[0][3] = 0
	mat[1][3] = 0
	mat[2][3] = -1
	mat[3][3] = 0

	return mat
}

// AssignOrthogonalProjection assigns an orthogonal projection transformation.
func (mat *T) AssignOrthogonalProjection(left, right, bottom, top, znear, zfar float32) *T {
	ooRightLeft := 1 / (right - left)
	ooTopBottom := 1 / (top - bottom)
	ooFarNear := 1 / (zfar - znear)

	mat[0][0] = 2 * ooRightLeft
	mat[1][0] = 0
	mat[2][0] = 0
	mat[3][0] = -(right + left) * ooRightLeft

	mat[0][1] = 0
	mat[1][1] = 2 * ooTopBottom
	mat[2][1] = 0
	mat[3][1] = -(top + bottom) * ooTopBottom

	mat[0][2] = 0
	mat[1][2] = 0
	mat[2][2] = -2 * ooFarNear
	mat[3][2] = -(zfar + znear) * ooFarNear

	mat[0][3] = 0
	mat[1][3] = 0
	mat[2][3] = 0
	mat[3][3] = 1

	return mat
}

// Determinant3x3 returns the determinant of the 3x3 sub-matrix.
func (mat *T) Determinant3x3() float32 {
	return mat[0][0]*mat[1][1]*mat[2][2] +
		mat[1][0]*mat[2][1]*mat[0][2] +
		mat[2][0]*mat[0][1]*mat[1][2] -
		mat[2][0]*mat[1][1]*mat[0][2] -
		mat[1][0]*mat[0][1]*mat[2][2] -
		mat[0][0]*mat[2][1]*mat[1][2]
}

func (mat *T) Determinant() float32 {
	s1 := mat[0][0]
	det1 := mat[1][1]*mat[2][2]*mat[3][3] +
		mat[2][1]*mat[3][2]*mat[1][3] +
		mat[3][1]*mat[1][2]*mat[2][3] -
		mat[3][1]*mat[2][2]*mat[1][3] -
		mat[2][1]*mat[1][2]*mat[3][3] -
		mat[1][1]*mat[3][2]*mat[2][3]

	s2 := mat[0][1]
	det2 := mat[1][0]*mat[2][2]*mat[3][3] +
		mat[2][0]*mat[3][2]*mat[1][3] +
		mat[3][0]*mat[1][2]*mat[2][3] -
		mat[3][0]*mat[2][2]*mat[1][3] -
		mat[2][0]*mat[1][2]*mat[3][3] -
		mat[1][0]*mat[3][2]*mat[2][3]
	s3 := mat[0][2]
	det3 := mat[1][0]*mat[2][1]*mat[3][3] +
		mat[2][0]*mat[3][1]*mat[1][3] +
		mat[3][0]*mat[1][1]*mat[2][3] -
		mat[3][0]*mat[2][1]*mat[1][3] -
		mat[2][0]*mat[1][1]*mat[3][3] -
		mat[1][0]*mat[3][1]*mat[2][3]
	s4 := mat[0][3]
	det4 := mat[1][0]*mat[2][1]*mat[3][2] +
		mat[2][0]*mat[3][1]*mat[1][2] +
		mat[3][0]*mat[1][1]*mat[2][2] -
		mat[3][0]*mat[2][1]*mat[1][2] -
		mat[2][0]*mat[1][1]*mat[3][2] -
		mat[1][0]*mat[3][1]*mat[2][2]
	return s1*det1 - s2*det2 + s3*det3 - s4*det4
}

// IsReflective returns true if the matrix can be reflected by a plane.
func (mat *T) IsReflective() bool {
	return mat.Determinant3x3() < 0
}

func swap(a, b *float32) {
	*a, *b = *b, *a
}

// Transpose transposes the matrix.
func (mat *T) Transpose() *T {
	swap(&mat[3][0], &mat[0][3])
	swap(&mat[3][1], &mat[1][3])
	swap(&mat[3][2], &mat[2][3])
	return mat.Transpose3x3()
}

// Transposed returns a transposed copy of the matrix.
func (mat *T) Transposed() T {
	result := *mat
	result.Transpose()
	return result
}

// Transpose3x3 transposes the 3x3 sub-matrix.
func (mat *T) Transpose3x3() *T {
	swap(&mat[1][0], &mat[0][1])
	swap(&mat[2][0], &mat[0][2])
	swap(&mat[2][1], &mat[1][2])
	return mat
}

// Adjugate computes the adjugate of this matrix and returns mat
func (mat *T) Adjugate() *T {
	matOrig := *mat
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			// - 1 for odd i+j, 1 for even i+j
			sign := float32(((i+j)%2)*-2 + 1)
			mat[i][j] = matOrig.maskedBlock(i, j).Determinant() * sign
		}
	}
	return mat.Transpose()
}

// Adjugated returns an adjugated copy of the matrix.
func (mat *T) Adjugated() T {
	result := *mat
	result.Adjugate()
	return result
}

// returns a 3x3 matrix without the i-th column and j-th row
func (mat *T) maskedBlock(blockI, blockJ int) *mat3.T {
	var m mat3.T
	m_i := 0
	for i := 0; i < 4; i++ {
		if i == blockI {
			continue
		}
		m_j := 0
		for j := 0; j < 4; j++ {
			if j == blockJ {
				continue
			}
			m[m_i][m_j] = mat[i][j]
			m_j++
		}
		m_i++
	}
	return &m
}

// Inverts the given matrix.
// Does not check if matrix is singular and may lead to strange results!
func (mat *T) Invert() *T {
	initialDet := mat.Determinant()
	mat.Adjugate()
	mat.Mul(1 / initialDet)
	return mat
}

// Inverted returns an inverted copy of the matrix.
// Does not check if matrix is singular and may lead to strange results!
func (mat *T) Inverted() T {
	result := *mat
	result.Invert()
	return result
}

func v3Combine(a *vec3.T, b *vec3.T, result *vec3.T, ascl float32, bscl float32) {

	result[0] = (ascl * a[0]) + (bscl * b[0])
	result[1] = (ascl * a[1]) + (bscl * b[1])
	result[2] = (ascl * a[2]) + (bscl * b[2])
}

func Decompose(mat *T) (*vec3.T, *quaternion.T, *vec3.T) {
	sx := (&vec3.T{mat[0][0], mat[0][1], mat[0][2]}).Length()
	sy := (&vec3.T{mat[1][0], mat[1][1], mat[1][2]}).Length()
	sz := (&vec3.T{mat[2][0], mat[2][1], mat[2][2]}).Length()

	// if determine is negative, we need to invert one scale
	det := mat.Determinant()
	if det < 0 {
		sx = -sx
	}

	position := vec3.T{mat[3][0], mat[3][1], mat[3][2]}

	// scale the rotation part
	invSX, invSY, invSZ := 1.0/sx, 1.0/sy, 1.0/sz
	matrix := T{
		vec4.T{mat[0][0] * invSX, mat[0][1] * invSX, mat[0][2] * invSX, 0},
		vec4.T{mat[1][0] * invSY, mat[1][1] * invSY, mat[1][2] * invSY, 0},
		vec4.T{mat[2][0] * invSZ, mat[2][1] * invSZ, mat[2][2] * invSZ, 0},
		vec4.T{0, 0, 0, 1},
	}
	quat := matrix.Quaternion()

	scale := vec3.T{sx, sy, sz}

	return &position, &quat, &scale
}

func Compose(pos *vec3.T, quat *quaternion.T, scale *vec3.T) *T {
	posMat := Ident
	posMat.SetTranslation(pos)
	quatMat := Ident
	quatMat.AssignQuaternion(quat)
	scaleMat := Ident
	scaleMat.ScaleVec3(scale)

	qsMat := Ident
	qsMat.AssignMul(&quatMat, &scaleMat)

	result := Ident
	result.AssignMul(&posMat, &qsMat)

	return &result
}

func (mat *T) LookAt(eye, target, up vec3.T) *T {
	_z := vec3.Sub(&eye, &target)

	if _z.LengthSqr() == 0 {
		// eye and target are in the same position
		_z[2] = 1
	}

	_z.Normalize()
	_x := vec3.Cross(&up, &_z)

	if _x.LengthSqr() == 0 {
		// up and z are parallel
		if math.Abs(up[2]) == 1 {
			_z[0] += 0.0001
		} else {
			_z[2] += 0.0001
		}

		_z.Normalize()
		_x = vec3.Cross(&up, &_z)
	}

	_x.Normalize()
	_y := vec3.Cross(&_z, &_x)

	mat[0][0] = _x[0]
	mat[0][1] = _y[0]
	mat[0][2] = _z[0]
	mat[1][0] = _x[1]
	mat[1][1] = _y[1]
	mat[1][2] = _z[1]
	mat[2][0] = _x[2]
	mat[2][1] = _y[2]
	mat[2][2] = _z[2]
	mat.Transpose()
	return mat
}

func AssignMul(a, b *T) *T {
	mat := Ident
	mat[0] = a.MulVec4(&b[0])
	mat[1] = a.MulVec4(&b[1])
	mat[2] = a.MulVec4(&b[2])
	mat[3] = a.MulVec4(&b[3])
	return &mat
}
