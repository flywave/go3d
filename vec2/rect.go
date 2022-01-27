package vec2

import (
	"fmt"
	"math"
)

// Rect is a coordinate system aligned rectangle defined by a Min and Max vector.
type Rect struct {
	Min T
	Max T
}

// ParseRect parses a Rect from a string. See also String()
func ParseRect(s string) (r Rect, err error) {
	_, err = fmt.Sscan(s, &r.Min[0], &r.Min[1], &r.Max[0], &r.Max[1])
	return r, err
}

func (rect *Rect) Width() float32 {
	return rect.Min[0] - rect.Max[0]
}

func (rect *Rect) Height() float32 {
	return rect.Min[1] - rect.Max[1]
}

func (rect *Rect) Size() float32 {
	width := rect.Width()
	height := rect.Height()
	return float32(math.Max(float64(width), float64(height)))
}

// Slice returns the elements of the vector as slice.
func (rect *Rect) Slice() []float32 {
	return rect.Array()[:]
}

func (rect *Rect) Array() *[4]float32 {
	return &[...]float32{
		rect.Min[0], rect.Min[1],
		rect.Max[0], rect.Max[1],
	}
}

// String formats Rect as string. See also ParseRect().
func (rect *Rect) String() string {
	return rect.Min.String() + " " + rect.Max.String()
}

// ContainsPoint returns if a point is contained within the rectangle.
func (rect *Rect) ContainsPoint(p *T) bool {
	return p[0] >= rect.Min[0] && p[0] <= rect.Max[0] &&
		p[1] >= rect.Min[1] && p[1] <= rect.Max[1]
}

func (rect *Rect) Contains(other *Rect) bool {
	return other.Min[0] >= rect.Min[0] &&
		other.Max[0] <= rect.Max[0] &&
		other.Min[1] >= rect.Min[1] &&
		other.Max[1] <= rect.Max[1]
}

func (rect *Rect) Intersects(other *Rect) bool {
	return other.Max[0] >= rect.Min[0] &&
		other.Min[0] <= rect.Max[0] &&
		other.Max[1] >= rect.Min[1] &&
		other.Min[1] <= rect.Max[1]
}

func (rect *Rect) Join(other *Rect) {
	rect.Min = Min(&rect.Min, &other.Min)
	rect.Max = Max(&rect.Max, &other.Max)
}

func (rect *Rect) Extend(p *T) {
	rect.Min = Min(&rect.Min, p)
	rect.Max = Max(&rect.Max, p)
}

// Joined returns the minimal rectangle containing both a and b.
func Joined(a, b *Rect) (rect Rect) {
	rect.Min = Min(&a.Min, &b.Min)
	rect.Max = Max(&a.Max, &b.Max)
	return rect
}
