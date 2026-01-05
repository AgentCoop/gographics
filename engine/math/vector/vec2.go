// vector/vector.go
package vector

import (
	"encoding/json"
	"fmt"
	"math"
)

// Vec2 represents a 2D vector
type Vec2[T Number] struct {
	X, Y T
}

// NewVec2 creates a new Vec2
func NewVec2[T Number](x, y T) Vec2[T] {
	return Vec2[T]{x, y}
}

// NewVec2Scalar creates a Vec2 with both components set to s
func NewVec2Scalar[T Number](s T) Vec2[T] {
	return Vec2[T]{s, s}
}

// Zero2 returns a zero vector
func Zero2[T Number]() Vec2[T] {
	return Vec2[T]{0, 0}
}

// One2 returns a vector with all components set to 1
func One2[T Number]() Vec2[T] {
	return Vec2[T]{1, 1}
}

// UnitX2 returns the X unit vector
func UnitX2[T Number]() Vec2[T] {
	return Vec2[T]{1, 0}
}

// UnitY2 returns the Y unit vector
func UnitY2[T Number]() Vec2[T] {
	return Vec2[T]{0, 1}
}

// Add adds two vectors
func (v Vec2[T]) Add(other Vec2[T]) Vec2[T] {
	return Vec2[T]{v.X + other.X, v.Y + other.Y}
}

// Sub subtracts two vectors
func (v Vec2[T]) Sub(other Vec2[T]) Vec2[T] {
	return Vec2[T]{v.X - other.X, v.Y - other.Y}
}

// Mul multiplies two vectors component-wise
func (v Vec2[T]) Mul(other Vec2[T]) Vec2[T] {
	return Vec2[T]{v.X * other.X, v.Y * other.Y}
}

// Div divides two vectors component-wise
func (v Vec2[T]) Div(other Vec2[T]) Vec2[T] {
	return Vec2[T]{v.X / other.X, v.Y / other.Y}
}

// AddScalar adds a scalar to each component
func (v Vec2[T]) AddScalar(s T) Vec2[T] {
	return Vec2[T]{v.X + s, v.Y + s}
}

// SubScalar subtracts a scalar from each component
func (v Vec2[T]) SubScalar(s T) Vec2[T] {
	return Vec2[T]{v.X - s, v.Y - s}
}

// MulScalar multiplies each component by a scalar
func (v Vec2[T]) MulScalar(s T) Vec2[T] {
	return Vec2[T]{v.X * s, v.Y * s}
}

// DivScalar divides each component by a scalar
func (v Vec2[T]) DivScalar(s T) Vec2[T] {
	return Vec2[T]{v.X / s, v.Y / s}
}

// Dot returns the dot product of two vectors
func (v Vec2[T]) Dot(other Vec2[T]) T {
	return v.X*other.X + v.Y*other.Y
}

// Cross returns the cross product (scalar in 2D)
func (v Vec2[T]) Cross(other Vec2[T]) T {
	return v.X*other.Y - v.Y*other.X
}

// Length returns the magnitude of the vector
func (v Vec2[T]) Length() T {
	return T(math.Sqrt(float64(v.X*v.X + v.Y*v.Y)))
}

// LengthSquared returns the squared magnitude
func (v Vec2[T]) LengthSquared() T {
	return v.X*v.X + v.Y*v.Y
}

// Normalize returns a normalized copy of the vector
func (v Vec2[T]) Normalize() Vec2[T] {
	length := v.Length()
	if float64(length) < EPSILON {
		return Zero2[T]()
	}
	invLength := 1 / length
	return Vec2[T]{v.X * invLength, v.Y * invLength}
}

// Distance returns the distance between two vectors
func (v Vec2[T]) Distance(other Vec2[T]) T {
	dx := v.X - other.X
	dy := v.Y - other.Y
	return T(math.Sqrt(float64(dx*dx + dy*dy)))
}

// DistanceSquared returns the squared distance
func (v Vec2[T]) DistanceSquared(other Vec2[T]) T {
	dx := v.X - other.X
	dy := v.Y - other.Y
	return dx*dx + dy*dy
}

// Lerp linearly interpolates between two vectors
func (v Vec2[T]) Lerp(other Vec2[T], t T) Vec2[T] {
	return Vec2[T]{
		v.X + (other.X-v.X)*t,
		v.Y + (other.Y-v.Y)*t,
	}
}

// Slerp spherically interpolates between two vectors
func (v Vec2[T]) Slerp(other Vec2[T], t T) Vec2[T] {
	dot := v.Dot(other)
	// Clamp dot to [-1, 1]
	if dot > 1 {
		dot = 1
	} else if dot < -1 {
		dot = -1
	}

	theta := T(math.Acos(float64(dot))) * t
	relative := other.Sub(v.MulScalar(dot)).Normalize()

	return v.MulScalar(T(math.Cos(float64(theta)))).Add(relative.MulScalar(T(math.Sin(float64(theta)))))
}

// Rotate rotates the vector by angle in radians
func (v Vec2[T]) Rotate(angle T) Vec2[T] {
	cos := T(math.Cos(float64(angle)))
	sin := T(math.Sin(float64(angle)))
	return Vec2[T]{
		v.X*cos - v.Y*sin,
		v.X*sin + v.Y*cos,
	}
}

// Angle returns the angle of the vector in radians
func (v Vec2[T]) Angle() T {
	return T(math.Atan2(float64(v.Y), float64(v.X)))
}

// AngleTo returns the angle between two vectors
func (v Vec2[T]) AngleTo(other Vec2[T]) T {
	dot := v.Dot(other)
	det := v.Cross(other)
	return T(math.Atan2(float64(det), float64(dot)))
}

// Perpendicular returns a perpendicular vector
func (v Vec2[T]) Perpendicular() Vec2[T] {
	return Vec2[T]{-v.Y, v.X}
}

// Reflect reflects the vector across a normal
func (v Vec2[T]) Reflect(normal Vec2[T]) Vec2[T] {
	dot := v.Dot(normal)
	return Vec2[T]{
		v.X - 2*dot*normal.X,
		v.Y - 2*dot*normal.Y,
	}
}

// Project projects the vector onto another
func (v Vec2[T]) Project(onto Vec2[T]) Vec2[T] {
	lengthSq := onto.LengthSquared()
	if float64(lengthSq) < EPSILON {
		return Zero2[T]()
	}
	return onto.MulScalar(v.Dot(onto) / lengthSq)
}

// Floor returns the vector with components rounded down
func (v Vec2[T]) Floor() Vec2[T] {
	return Vec2[T]{
		T(math.Floor(float64(v.X))),
		T(math.Floor(float64(v.Y))),
	}
}

// Ceil returns the vector with components rounded up
func (v Vec2[T]) Ceil() Vec2[T] {
	return Vec2[T]{
		T(math.Ceil(float64(v.X))),
		T(math.Ceil(float64(v.Y))),
	}
}

// Round returns the vector with components rounded
func (v Vec2[T]) Round() Vec2[T] {
	return Vec2[T]{
		T(math.Round(float64(v.X))),
		T(math.Round(float64(v.Y))),
	}
}

// Abs returns the vector with absolute values
func (v Vec2[T]) Abs() Vec2[T] {
	return Vec2[T]{
		T(math.Abs(float64(v.X))),
		T(math.Abs(float64(v.Y))),
	}
}

// Min returns component-wise minimum
func (v Vec2[T]) Min(other Vec2[T]) Vec2[T] {
	return Vec2[T]{
		T(math.Min(float64(v.X), float64(other.X))),
		T(math.Min(float64(v.Y), float64(other.Y))),
	}
}

// Max returns component-wise maximum
func (v Vec2[T]) Max(other Vec2[T]) Vec2[T] {
	return Vec2[T]{
		T(math.Max(float64(v.X), float64(other.X))),
		T(math.Max(float64(v.Y), float64(other.Y))),
	}
}

// Clamp clamps the vector components between min and max
func (v Vec2[T]) Clamp(min, max Vec2[T]) Vec2[T] {
	return Vec2[T]{
		T(math.Max(float64(min.X), math.Min(float64(max.X), float64(v.X)))),
		T(math.Max(float64(min.Y), math.Min(float64(max.Y), float64(v.Y)))),
	}
}

// ClampScalar clamps the vector components between min and max scalars
func (v Vec2[T]) ClampScalar(min, max T) Vec2[T] {
	return Vec2[T]{
		T(math.Max(float64(min), math.Min(float64(max), float64(v.X)))),
		T(math.Max(float64(min), math.Min(float64(max), float64(v.Y)))),
	}
}

// ClampLength clamps the vector length
func (v Vec2[T]) ClampLength(min, max T) Vec2[T] {
	length := v.Length()
	if float64(length) < EPSILON {
		return Zero2[T]()
	}
	if length < min {
		return v.Normalize().MulScalar(min)
	}
	if length > max {
		return v.Normalize().MulScalar(max)
	}
	return v
}

// Equals checks if two vectors are approximately equal
func (v Vec2[T]) Equals(other Vec2[T], epsilon ...float64) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(v.X-other.X)) <= eps && math.Abs(float64(v.Y-other.Y)) <= eps
}

// IsZero checks if the vector is approximately zero
func (v Vec2[T]) IsZero(epsilon ...float64) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(v.X)) <= eps && math.Abs(float64(v.Y)) <= eps
}

// Zero sets the vector to zero
func (v *Vec2[T]) Zero() {
	var zero T
	v.X, v.Y = zero, zero
}

// Set sets the vector components
func (v *Vec2[T]) Set(x, y T) {
	v.X, v.Y = x, y
}

// Copy copies another vector
func (v *Vec2[T]) Copy(other Vec2[T]) {
	v.X, v.Y = other.X, other.Y
}

// Clone returns a copy of the vector
func (v Vec2[T]) Clone() Vec2[T] {
	return Vec2[T]{v.X, v.Y}
}

// ToArray returns the vector as a slice
func (v Vec2[T]) ToArray() []T {
	return []T{v.X, v.Y}
}

// ToVec3 converts to Vec3 with z = 0
func (v Vec2[T]) ToVec3() Vec3[T] {
	var zero T
	return Vec3[T]{v.X, v.Y, zero}
}

// ToVec4 converts to Vec4 with z = 0, w = 1
func (v Vec2[T]) ToVec4() Vec4[T] {
	var zero T
	one := T(1)
	return Vec4[T]{v.X, v.Y, zero, one}
}

// String returns a string representation
func (v Vec2[T]) String() string {
	return fmt.Sprintf("Vec2[%T](%.6f, %.6f)", v.X, v.X, v.Y)
}

// MarshalJSON implements json.Marshaler
func (v Vec2[T]) MarshalJSON() ([]byte, error) {
	return []byte(fmt.Sprintf("[%f,%f]", v.X, v.Y)), nil
}

// UnmarshalJSON implements json.Unmarshaler
func (v *Vec2[T]) UnmarshalJSON(data []byte) error {
	var arr [2]T
	if err := json.Unmarshal(data, &arr); err != nil {
		return err
	}
	v.X, v.Y = arr[0], arr[1]
	return nil
}

// Transform transforms the vector by a 2x2 matrix
func (v Vec2[T]) Transform(m [4]T) Vec2[T] {
	return Vec2[T]{
		v.X*m[0] + v.Y*m[2],
		v.X*m[1] + v.Y*m[3],
	}
}
