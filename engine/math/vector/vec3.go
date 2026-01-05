package vector

import (
	"encoding/json"
	"fmt"
	"math"
)

// Vec3 represents a 3D vector
type Vec3[T Number] struct {
	X, Y, Z T
}

// NewVec3 creates a new Vec3
func NewVec3[T Number](x, y, z T) Vec3[T] {
	return Vec3[T]{x, y, z}
}

// NewVec3Scalar creates a Vec3 with all components set to s
func NewVec3Scalar[T Number](s T) Vec3[T] {
	return Vec3[T]{s, s, s}
}

// Zero3 returns a zero vector
func Zero3[T Number]() Vec3[T] {
	return Vec3[T]{0, 0, 0}
}

// One3 returns a vector with all components set to 1
func One3[T Number]() Vec3[T] {
	return Vec3[T]{1, 1, 1}
}

// UnitX3 returns the X unit vector
func UnitX3[T Number]() Vec3[T] {
	return Vec3[T]{1, 0, 0}
}

// UnitY3 returns the Y unit vector
func UnitY3[T Number]() Vec3[T] {
	return Vec3[T]{0, 1, 0}
}

// UnitZ3 returns the Z unit vector
func UnitZ3[T Number]() Vec3[T] {
	return Vec3[T]{0, 0, 1}
}

// Up returns the up vector (0, 1, 0)
func Up[T Number]() Vec3[T] {
	return Vec3[T]{0, 1, 0}
}

// Down returns the down vector (0, -1, 0)
func Down[T Number]() Vec3[T] {
	return Vec3[T]{0, -1, 0}
}

// Left returns the left vector (-1, 0, 0)
func Left[T Number]() Vec3[T] {
	return Vec3[T]{-1, 0, 0}
}

// Right returns the right vector (1, 0, 0)
func Right[T Number]() Vec3[T] {
	return Vec3[T]{1, 0, 0}
}

// Forward returns the forward vector (0, 0, -1) for OpenGL
func Forward[T Number]() Vec3[T] {
	return Vec3[T]{0, 0, -1}
}

// Back returns the back vector (0, 0, 1)
func Back[T Number]() Vec3[T] {
	return Vec3[T]{0, 0, 1}
}

// FromVec2To3 creates a Vec3 from Vec2 with z component
func FromVec2To3[T Number](v Vec2[T], z T) Vec3[T] {
	return Vec3[T]{v.X, v.Y, z}
}

// XY returns the XY components as Vec2
func (v Vec3[T]) XY() Vec2[T] {
	return Vec2[T]{v.X, v.Y}
}

// Add adds two vectors
func (v Vec3[T]) Add(other Vec3[T]) Vec3[T] {
	return Vec3[T]{v.X + other.X, v.Y + other.Y, v.Z + other.Z}
}

// Sub subtracts two vectors
func (v Vec3[T]) Sub(other Vec3[T]) Vec3[T] {
	return Vec3[T]{v.X - other.X, v.Y - other.Y, v.Z - other.Z}
}

// Mul multiplies two vectors component-wise
func (v Vec3[T]) Mul(other Vec3[T]) Vec3[T] {
	return Vec3[T]{v.X * other.X, v.Y * other.Y, v.Z * other.Z}
}

// Div divides two vectors component-wise
func (v Vec3[T]) Div(other Vec3[T]) Vec3[T] {
	return Vec3[T]{v.X / other.X, v.Y / other.Y, v.Z / other.Z}
}

// AddScalar adds a scalar to each component
func (v Vec3[T]) AddScalar(s T) Vec3[T] {
	return Vec3[T]{v.X + s, v.Y + s, v.Z + s}
}

// SubScalar subtracts a scalar from each component
func (v Vec3[T]) SubScalar(s T) Vec3[T] {
	return Vec3[T]{v.X - s, v.Y - s, v.Z - s}
}

// MulScalar multiplies each component by a scalar
func (v Vec3[T]) MulScalar(s T) Vec3[T] {
	return Vec3[T]{v.X * s, v.Y * s, v.Z * s}
}

// DivScalar divides each component by a scalar
func (v Vec3[T]) DivScalar(s T) Vec3[T] {
	return Vec3[T]{v.X / s, v.Y / s, v.Z / s}
}

// Dot returns the dot product of two vectors
func (v Vec3[T]) Dot(other Vec3[T]) T {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}

// Cross returns the cross product of two vectors
func (v Vec3[T]) Cross(other Vec3[T]) Vec3[T] {
	return Vec3[T]{
		v.Y*other.Z - v.Z*other.Y,
		v.Z*other.X - v.X*other.Z,
		v.X*other.Y - v.Y*other.X,
	}
}

// Length returns the magnitude of the vector
func (v Vec3[T]) Length() T {
	return T(math.Sqrt(float64(v.X*v.X + v.Y*v.Y + v.Z*v.Z)))
}

// LengthSquared returns the squared magnitude
func (v Vec3[T]) LengthSquared() T {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

// Normalize returns a normalized copy of the vector
func (v Vec3[T]) Normalize() Vec3[T] {
	length := v.Length()
	if float64(length) < EPSILON {
		return Zero3[T]()
	}
	invLength := 1 / length
	return Vec3[T]{v.X * invLength, v.Y * invLength, v.Z * invLength}
}

// Distance returns the distance between two vectors
func (v Vec3[T]) Distance(other Vec3[T]) T {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	return T(math.Sqrt(float64(dx*dx + dy*dy + dz*dz)))
}

// DistanceSquared returns the squared distance
func (v Vec3[T]) DistanceSquared(other Vec3[T]) T {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	return dx*dx + dy*dy + dz*dz
}

// Lerp linearly interpolates between two vectors
func (v Vec3[T]) Lerp(other Vec3[T], t T) Vec3[T] {
	return Vec3[T]{
		v.X + (other.X-v.X)*t,
		v.Y + (other.Y-v.Y)*t,
		v.Z + (other.Z-v.Z)*t,
	}
}

// Slerp spherically interpolates between two vectors
func (v Vec3[T]) Slerp(other Vec3[T], t T) Vec3[T] {
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

// Angle returns the angle between two vectors in radians
func (v Vec3[T]) Angle(other Vec3[T]) T {
	dot := v.Dot(other)
	lengthProduct := v.Length() * other.Length()
	if float64(lengthProduct) < EPSILON {
		return 0
	}
	return T(math.Acos(float64(dot / lengthProduct)))
}

// Reflect reflects the vector across a normal
func (v Vec3[T]) Reflect(normal Vec3[T]) Vec3[T] {
	dot := v.Dot(normal)
	return Vec3[T]{
		v.X - 2*dot*normal.X,
		v.Y - 2*dot*normal.Y,
		v.Z - 2*dot*normal.Z,
	}
}

// Refract refracts the vector
func (v Vec3[T]) Refract(normal Vec3[T], eta T) Vec3[T] {
	dot := v.Dot(normal)
	k := 1.0 - eta*eta*(1.0-dot*dot)
	if k < 0 {
		return Zero3[T]()
	}
	return v.MulScalar(eta).Sub(normal.MulScalar(eta*dot + T(math.Sqrt(float64(k)))))
}

// Project projects the vector onto another
func (v Vec3[T]) Project(onto Vec3[T]) Vec3[T] {
	lengthSq := onto.LengthSquared()
	if float64(lengthSq) < EPSILON {
		return Zero3[T]()
	}
	return onto.MulScalar(v.Dot(onto) / lengthSq)
}

// Reject returns the rejection of the vector from another
func (v Vec3[T]) Reject(onto Vec3[T]) Vec3[T] {
	return v.Sub(v.Project(onto))
}

// Floor returns the vector with components rounded down
func (v Vec3[T]) Floor() Vec3[T] {
	return Vec3[T]{
		T(math.Floor(float64(v.X))),
		T(math.Floor(float64(v.Y))),
		T(math.Floor(float64(v.Z))),
	}
}

// Ceil returns the vector with components rounded up
func (v Vec3[T]) Ceil() Vec3[T] {
	return Vec3[T]{
		T(math.Ceil(float64(v.X))),
		T(math.Ceil(float64(v.Y))),
		T(math.Ceil(float64(v.Z))),
	}
}

// Round returns the vector with components rounded
func (v Vec3[T]) Round() Vec3[T] {
	return Vec3[T]{
		T(math.Round(float64(v.X))),
		T(math.Round(float64(v.Y))),
		T(math.Round(float64(v.Z))),
	}
}

// Abs returns the vector with absolute values
func (v Vec3[T]) Abs() Vec3[T] {
	return Vec3[T]{
		T(math.Abs(float64(v.X))),
		T(math.Abs(float64(v.Y))),
		T(math.Abs(float64(v.Z))),
	}
}

// Min returns component-wise minimum
func (v Vec3[T]) Min(other Vec3[T]) Vec3[T] {
	return Vec3[T]{
		T(math.Min(float64(v.X), float64(other.X))),
		T(math.Min(float64(v.Y), float64(other.Y))),
		T(math.Min(float64(v.Z), float64(other.Z))),
	}
}

// Max returns component-wise maximum
func (v Vec3[T]) Max(other Vec3[T]) Vec3[T] {
	return Vec3[T]{
		T(math.Max(float64(v.X), float64(other.X))),
		T(math.Max(float64(v.Y), float64(other.Y))),
		T(math.Max(float64(v.Z), float64(other.Z))),
	}
}

// Clamp clamps the vector components between min and max
func (v Vec3[T]) Clamp(min, max Vec3[T]) Vec3[T] {
	return Vec3[T]{
		T(math.Max(float64(min.X), math.Min(float64(max.X), float64(v.X)))),
		T(math.Max(float64(min.Y), math.Min(float64(max.Y), float64(v.Y)))),
		T(math.Max(float64(min.Z), math.Min(float64(max.Z), float64(v.Z)))),
	}
}

// ClampScalar clamps the vector components between min and max scalars
func (v Vec3[T]) ClampScalar(min, max T) Vec3[T] {
	return Vec3[T]{
		T(math.Max(float64(min), math.Min(float64(max), float64(v.X)))),
		T(math.Max(float64(min), math.Min(float64(max), float64(v.Y)))),
		T(math.Max(float64(min), math.Min(float64(max), float64(v.Z)))),
	}
}

// ClampLength clamps the vector length
func (v Vec3[T]) ClampLength(min, max T) Vec3[T] {
	length := v.Length()
	if float64(length) < EPSILON {
		return Zero3[T]()
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
func (v Vec3[T]) Equals(other Vec3[T], epsilon ...float64) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(v.X-other.X)) <= eps &&
		math.Abs(float64(v.Y-other.Y)) <= eps &&
		math.Abs(float64(v.Z-other.Z)) <= eps
}

// IsZero checks if the vector is approximately zero
func (v Vec3[T]) IsZero(epsilon ...float64) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(v.X)) <= eps &&
		math.Abs(float64(v.Y)) <= eps &&
		math.Abs(float64(v.Z)) <= eps
}

// Zero sets the vector to zero
func (v *Vec3[T]) Zero() {
	var zero T
	v.X, v.Y, v.Z = zero, zero, zero
}

// Set sets the vector components
func (v *Vec3[T]) Set(x, y, z T) {
	v.X, v.Y, v.Z = x, y, z
}

// Copy copies another vector
func (v *Vec3[T]) Copy(other Vec3[T]) {
	v.X, v.Y, v.Z = other.X, other.Y, other.Z
}

// Clone returns a copy of the vector
func (v Vec3[T]) Clone() Vec3[T] {
	return Vec3[T]{v.X, v.Y, v.Z}
}

// ToArray returns the vector as a slice
func (v Vec3[T]) ToArray() []T {
	return []T{v.X, v.Y, v.Z}
}

// ToVec2 converts to Vec2, dropping Z
func (v Vec3[T]) ToVec2() Vec2[T] {
	return Vec2[T]{v.X, v.Y}
}

// ToVec4 converts to Vec4 with w = 1
func (v Vec3[T]) ToVec4() Vec4[T] {
	one := T(1)
	return Vec4[T]{v.X, v.Y, v.Z, one}
}

// String returns a string representation
func (v Vec3[T]) String() string {
	return fmt.Sprintf("Vec3[%T](%.6f, %.6f, %.6f)", v.X, v.X, v.Y, v.Z)
}

// MarshalJSON implements json.Marshaler
func (v Vec3[T]) MarshalJSON() ([]byte, error) {
	return []byte(fmt.Sprintf("[%f,%f,%f]", v.X, v.Y, v.Z)), nil
}

// UnmarshalJSON implements json.Unmarshaler
func (v *Vec3[T]) UnmarshalJSON(data []byte) error {
	var arr [3]T
	if err := json.Unmarshal(data, &arr); err != nil {
		return err
	}
	v.X, v.Y, v.Z = arr[0], arr[1], arr[2]
	return nil
}

// Transform transforms the vector by a 3x3 matrix
func (v Vec3[T]) Transform(m [9]T) Vec3[T] {
	return Vec3[T]{
		v.X*m[0] + v.Y*m[3] + v.Z*m[6],
		v.X*m[1] + v.Y*m[4] + v.Z*m[7],
		v.X*m[2] + v.Y*m[5] + v.Z*m[8],
	}
}
