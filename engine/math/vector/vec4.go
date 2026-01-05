package vector

import (
	"encoding/json"
	"fmt"
	"math"
)

// Vec4 represents a 4D vector or homogeneous coordinates
type Vec4[T Number] struct {
	X, Y, Z, W T
}

// NewVec4 creates a new Vec4
func NewVec4[T Number](x, y, z, w T) Vec4[T] {
	return Vec4[T]{x, y, z, w}
}

// NewVec4Scalar creates a Vec4 with all components set to s
func NewVec4Scalar[T Number](s T) Vec4[T] {
	return Vec4[T]{s, s, s, s}
}

// Zero4 returns a zero vector
func Zero4[T Number]() Vec4[T] {
	return Vec4[T]{0, 0, 0, 0}
}

// One4 returns a vector with all components set to 1
func One4[T Number]() Vec4[T] {
	return Vec4[T]{1, 1, 1, 1}
}

// NewPosition creates a position vector (w = 1)
func NewPosition[T Number](x, y, z T) Vec4[T] {
	one := T(1)
	return Vec4[T]{x, y, z, one}
}

// NewDirection creates a direction vector (w = 0)
func NewDirection[T Number](x, y, z T) Vec4[T] {
	var zero T
	return Vec4[T]{x, y, z, zero}
}

// FromVec2To4 creates a Vec4 from Vec2 with z and w components
func FromVec2To4[T Number](v Vec2[T], z, w T) Vec4[T] {
	return Vec4[T]{v.X, v.Y, z, w}
}

// FromVec3 creates a Vec4 from Vec3 with w component
func FromVec3[T Number](v Vec3[T], w T) Vec4[T] {
	return Vec4[T]{v.X, v.Y, v.Z, w}
}

// XYZ returns the XYZ components as Vec3
func (v Vec4[T]) XYZ() Vec3[T] {
	return Vec3[T]{v.X, v.Y, v.Z}
}

// XY returns the XY components as Vec2
func (v Vec4[T]) XY() Vec2[T] {
	return Vec2[T]{v.X, v.Y}
}

// Add adds two vectors
func (v Vec4[T]) Add(other Vec4[T]) Vec4[T] {
	return Vec4[T]{
		v.X + other.X,
		v.Y + other.Y,
		v.Z + other.Z,
		v.W + other.W,
	}
}

// Sub subtracts two vectors
func (v Vec4[T]) Sub(other Vec4[T]) Vec4[T] {
	return Vec4[T]{
		v.X - other.X,
		v.Y - other.Y,
		v.Z - other.Z,
		v.W - other.W,
	}
}

// Mul multiplies two vectors component-wise
func (v Vec4[T]) Mul(other Vec4[T]) Vec4[T] {
	return Vec4[T]{
		v.X * other.X,
		v.Y * other.Y,
		v.Z * other.Z,
		v.W * other.W,
	}
}

// Div divides two vectors component-wise
func (v Vec4[T]) Div(other Vec4[T]) Vec4[T] {
	return Vec4[T]{
		v.X / other.X,
		v.Y / other.Y,
		v.Z / other.Z,
		v.W / other.W,
	}
}

// AddScalar adds a scalar to each component
func (v Vec4[T]) AddScalar(s T) Vec4[T] {
	return Vec4[T]{v.X + s, v.Y + s, v.Z + s, v.W + s}
}

// SubScalar subtracts a scalar from each component
func (v Vec4[T]) SubScalar(s T) Vec4[T] {
	return Vec4[T]{v.X - s, v.Y - s, v.Z - s, v.W - s}
}

// MulScalar multiplies each component by a scalar
func (v Vec4[T]) MulScalar(s T) Vec4[T] {
	return Vec4[T]{v.X * s, v.Y * s, v.Z * s, v.W * s}
}

// DivScalar divides each component by a scalar
func (v Vec4[T]) DivScalar(s T) Vec4[T] {
	return Vec4[T]{v.X / s, v.Y / s, v.Z / s, v.W / s}
}

// Dot returns the dot product of two vectors
func (v Vec4[T]) Dot(other Vec4[T]) T {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z + v.W*other.W
}

// Length returns the magnitude of the vector
func (v Vec4[T]) Length() T {
	return T(math.Sqrt(float64(v.X*v.X + v.Y*v.Y + v.Z*v.Z + v.W*v.W)))
}

// LengthSquared returns the squared magnitude
func (v Vec4[T]) LengthSquared() T {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z + v.W*v.W
}

// Normalize returns a normalized copy of the vector
func (v Vec4[T]) Normalize() Vec4[T] {
	length := v.Length()
	if float64(length) < EPSILON {
		return Zero4[T]()
	}
	invLength := 1 / length
	return Vec4[T]{v.X * invLength, v.Y * invLength, v.Z * invLength, v.W * invLength}
}

// Distance returns the distance between two vectors
func (v Vec4[T]) Distance(other Vec4[T]) T {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	dw := v.W - other.W
	return T(math.Sqrt(float64(dx*dx + dy*dy + dz*dz + dw*dw)))
}

// DistanceSquared returns the squared distance
func (v Vec4[T]) DistanceSquared(other Vec4[T]) T {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	dw := v.W - other.W
	return dx*dx + dy*dy + dz*dz + dw*dw
}

// Lerp linearly interpolates between two vectors
func (v Vec4[T]) Lerp(other Vec4[T], t T) Vec4[T] {
	return Vec4[T]{
		v.X + (other.X-v.X)*t,
		v.Y + (other.Y-v.Y)*t,
		v.Z + (other.Z-v.Z)*t,
		v.W + (other.W-v.W)*t,
	}
}

// Homogenize divides by W component for perspective division
func (v Vec4[T]) Homogenize() Vec4[T] {
	if float32(math.Abs(float64(v.W))) < EPSILON {
		return Zero4[T]()
	}
	invW := 1.0 / v.W
	one := T(1)
	return Vec4[T]{v.X * invW, v.Y * invW, v.Z * invW, one}
}

// IsPoint returns true if W == 1 (or approximately 1)
func (v Vec4[T]) IsPoint(epsilon ...float64) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	one := T(1)
	return math.Abs(float64(v.W-one)) <= eps
}

// IsDirection returns true if W == 0 (or approximately 0)
func (v Vec4[T]) IsDirection(epsilon ...float64) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(v.W)) <= eps
}

// ToPosition converts to position (W = 1)
func (v Vec4[T]) ToPosition() Vec4[T] {
	if v.IsDirection() {
		one := T(1)
		return Vec4[T]{v.X, v.Y, v.Z, one}
	}
	return v.Homogenize()
}

// ToDirection converts to direction (W = 0)
func (v Vec4[T]) ToDirection() Vec4[T] {
	if v.IsDirection() {
		return v
	}
	var zero T
	return Vec4[T]{v.X, v.Y, v.Z, zero}
}

// Floor returns the vector with components rounded down
func (v Vec4[T]) Floor() Vec4[T] {
	return Vec4[T]{
		T(math.Floor(float64(v.X))),
		T(math.Floor(float64(v.Y))),
		T(math.Floor(float64(v.Z))),
		T(math.Floor(float64(v.W))),
	}
}

// Ceil returns the vector with components rounded up
func (v Vec4[T]) Ceil() Vec4[T] {
	return Vec4[T]{
		T(math.Ceil(float64(v.X))),
		T(math.Ceil(float64(v.Y))),
		T(math.Ceil(float64(v.Z))),
		T(math.Ceil(float64(v.W))),
	}
}

// Round returns the vector with components rounded
func (v Vec4[T]) Round() Vec4[T] {
	return Vec4[T]{
		T(math.Round(float64(v.X))),
		T(math.Round(float64(v.Y))),
		T(math.Round(float64(v.Z))),
		T(math.Round(float64(v.W))),
	}
}

// Abs returns the vector with absolute values
func (v Vec4[T]) Abs() Vec4[T] {
	return Vec4[T]{
		T(math.Abs(float64(v.X))),
		T(math.Abs(float64(v.Y))),
		T(math.Abs(float64(v.Z))),
		T(math.Abs(float64(v.W))),
	}
}

// Min returns component-wise minimum
func (v Vec4[T]) Min(other Vec4[T]) Vec4[T] {
	return Vec4[T]{
		T(math.Min(float64(v.X), float64(other.X))),
		T(math.Min(float64(v.Y), float64(other.Y))),
		T(math.Min(float64(v.Z), float64(other.Z))),
		T(math.Min(float64(v.W), float64(other.W))),
	}
}

// Max returns component-wise maximum
func (v Vec4[T]) Max(other Vec4[T]) Vec4[T] {
	return Vec4[T]{
		T(math.Max(float64(v.X), float64(other.X))),
		T(math.Max(float64(v.Y), float64(other.Y))),
		T(math.Max(float64(v.Z), float64(other.Z))),
		T(math.Max(float64(v.W), float64(other.W))),
	}
}

// Clamp clamps the vector components between min and max
func (v Vec4[T]) Clamp(min, max Vec4[T]) Vec4[T] {
	return Vec4[T]{
		T(math.Max(float64(min.X), math.Min(float64(max.X), float64(v.X)))),
		T(math.Max(float64(min.Y), math.Min(float64(max.Y), float64(v.Y)))),
		T(math.Max(float64(min.Z), math.Min(float64(max.Z), float64(v.Z)))),
		T(math.Max(float64(min.W), math.Min(float64(max.W), float64(v.W)))),
	}
}

// ClampScalar clamps the vector components between min and max scalars
func (v Vec4[T]) ClampScalar(min, max T) Vec4[T] {
	return Vec4[T]{
		T(math.Max(float64(min), math.Min(float64(max), float64(v.X)))),
		T(math.Max(float64(min), math.Min(float64(max), float64(v.Y)))),
		T(math.Max(float64(min), math.Min(float64(max), float64(v.Z)))),
		T(math.Max(float64(min), math.Min(float64(max), float64(v.W)))),
	}
}

// ClampLength clamps the vector length
func (v Vec4[T]) ClampLength(min, max T) Vec4[T] {
	length := v.Length()
	if float64(length) < EPSILON {
		return Zero4[T]()
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
func (v Vec4[T]) Equals(other Vec4[T], epsilon ...float64) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(v.X-other.X)) <= eps &&
		math.Abs(float64(v.Y-other.Y)) <= eps &&
		math.Abs(float64(v.Z-other.Z)) <= eps &&
		math.Abs(float64(v.W-other.W)) <= eps
}

// IsZero checks if the vector is approximately zero
func (v Vec4[T]) IsZero(epsilon ...float64) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(v.X)) <= eps &&
		math.Abs(float64(v.Y)) <= eps &&
		math.Abs(float64(v.Z)) <= eps &&
		math.Abs(float64(v.W)) <= eps
}

// Zero sets the vector to zero
func (v *Vec4[T]) Zero() {
	var zero T
	v.X, v.Y, v.Z, v.W = zero, zero, zero, zero
}

// Set sets the vector components
func (v *Vec4[T]) Set(x, y, z, w T) {
	v.X, v.Y, v.Z, v.W = x, y, z, w
}

// Copy copies another vector
func (v *Vec4[T]) Copy(other Vec4[T]) {
	v.X, v.Y, v.Z, v.W = other.X, other.Y, other.Z, other.W
}

// Clone returns a copy of the vector
func (v Vec4[T]) Clone() Vec4[T] {
	return Vec4[T]{v.X, v.Y, v.Z, v.W}
}

// ToArray returns the vector as a slice
func (v Vec4[T]) ToArray() []T {
	return []T{v.X, v.Y, v.Z, v.W}
}

// ToVec2 converts to Vec2, dropping Z and W
func (v Vec4[T]) ToVec2() Vec2[T] {
	return Vec2[T]{v.X, v.Y}
}

// ToVec3 converts to Vec3, dropping W
func (v Vec4[T]) ToVec3() Vec3[T] {
	if v.IsDirection() {
		return Vec3[T]{v.X, v.Y, v.Z}
	}
	// Homogenize first for points
	homogenized := v.Homogenize()
	return Vec3[T]{homogenized.X, homogenized.Y, homogenized.Z}
}

// String returns a string representation
func (v Vec4[T]) String() string {
	return fmt.Sprintf("Vec4[%T](%.6f, %.6f, %.6f, %.6f)", v.X, v.X, v.Y, v.Z, v.W)
}

// MarshalJSON implements json.Marshaler
func (v Vec4[T]) MarshalJSON() ([]byte, error) {
	return []byte(fmt.Sprintf("[%f,%f,%f,%f]", v.X, v.Y, v.Z, v.W)), nil
}

// UnmarshalJSON implements json.Unmarshaler
func (v *Vec4[T]) UnmarshalJSON(data []byte) error {
	var arr [4]T
	if err := json.Unmarshal(data, &arr); err != nil {
		return err
	}
	v.X, v.Y, v.Z, v.W = arr[0], arr[1], arr[2], arr[3]
	return nil
}

// Transform transforms the vector by a 4x4 matrix
func (v Vec4[T]) Transform(m [16]T) Vec4[T] {
	return Vec4[T]{
		v.X*m[0] + v.Y*m[4] + v.Z*m[8] + v.W*m[12],
		v.X*m[1] + v.Y*m[5] + v.Z*m[9] + v.W*m[13],
		v.X*m[2] + v.Y*m[6] + v.Z*m[10] + v.W*m[14],
		v.X*m[3] + v.Y*m[7] + v.Z*m[11] + v.W*m[15],
	}
}

// TransformAffine transforms by a 4x4 matrix assuming affine transformation
func (v Vec4[T]) TransformAffine(m [16]T) Vec4[T] {
	return Vec4[T]{
		v.X*m[0] + v.Y*m[4] + v.Z*m[8] + v.W*m[12],
		v.X*m[1] + v.Y*m[5] + v.Z*m[9] + v.W*m[13],
		v.X*m[2] + v.Y*m[6] + v.Z*m[10] + v.W*m[14],
		v.W, // W unchanged for affine transformations
	}
}
