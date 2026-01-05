package vector

import (
	"encoding/json"
	"fmt"
	"math"
)

// Vec3 represents a 3D vector
type Vec3 struct {
	X, Y, Z float32
}

// NewVec3 creates a new Vec3
func NewVec3(x, y, z float32) Vec3 {
	return Vec3{x, y, z}
}

// NewVec3Scalar creates a Vec3 with all components set to s
func NewVec3Scalar(s float32) Vec3 {
	return Vec3{s, s, s}
}

// Zero3 returns a zero vector
func Zero3() Vec3 {
	return Vec3{0, 0, 0}
}

// One3 returns a vector with all components set to 1
func One3() Vec3 {
	return Vec3{1, 1, 1}
}

// UnitX3 returns the X unit vector
func UnitX3() Vec3 {
	return Vec3{1, 0, 0}
}

// UnitY3 returns the Y unit vector
func UnitY3() Vec3 {
	return Vec3{0, 1, 0}
}

// UnitZ3 returns the Z unit vector
func UnitZ3() Vec3 {
	return Vec3{0, 0, 1}
}

// Up returns the up vector (0, 1, 0)
func Up() Vec3 {
	return Vec3{0, 1, 0}
}

// Down returns the down vector (0, -1, 0)
func Down() Vec3 {
	return Vec3{0, -1, 0}
}

// Left returns the left vector (-1, 0, 0)
func Left() Vec3 {
	return Vec3{-1, 0, 0}
}

// Right returns the right vector (1, 0, 0)
func Right() Vec3 {
	return Vec3{1, 0, 0}
}

// Forward returns the forward vector (0, 0, -1) for OpenGL
func Forward() Vec3 {
	return Vec3{0, 0, -1}
}

// Back returns the back vector (0, 0, 1)
func Back() Vec3 {
	return Vec3{0, 0, 1}
}

// Add adds two vectors
func (v Vec3) Add(other Vec3) Vec3 {
	return Vec3{v.X + other.X, v.Y + other.Y, v.Z + other.Z}
}

// Sub subtracts two vectors
func (v Vec3) Sub(other Vec3) Vec3 {
	return Vec3{v.X - other.X, v.Y - other.Y, v.Z - other.Z}
}

// Mul multiplies two vectors component-wise
func (v Vec3) Mul(other Vec3) Vec3 {
	return Vec3{v.X * other.X, v.Y * other.Y, v.Z * other.Z}
}

// Div divides two vectors component-wise
func (v Vec3) Div(other Vec3) Vec3 {
	return Vec3{v.X / other.X, v.Y / other.Y, v.Z / other.Z}
}

// AddScalar adds a scalar to each component
func (v Vec3) AddScalar(s float32) Vec3 {
	return Vec3{v.X + s, v.Y + s, v.Z + s}
}

// SubScalar subtracts a scalar from each component
func (v Vec3) SubScalar(s float32) Vec3 {
	return Vec3{v.X - s, v.Y - s, v.Z - s}
}

// MulScalar multiplies each component by a scalar
func (v Vec3) MulScalar(s float32) Vec3 {
	return Vec3{v.X * s, v.Y * s, v.Z * s}
}

// DivScalar divides each component by a scalar
func (v Vec3) DivScalar(s float32) Vec3 {
	return Vec3{v.X / s, v.Y / s, v.Z / s}
}

// Dot returns the dot product of two vectors
func (v Vec3) Dot(other Vec3) float32 {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}

// Cross returns the cross product
func (v Vec3) Cross(other Vec3) Vec3 {
	return Vec3{
		v.Y*other.Z - v.Z*other.Y,
		v.Z*other.X - v.X*other.Z,
		v.X*other.Y - v.Y*other.X,
	}
}

// Length returns the magnitude of the vector
func (v Vec3) Length() float32 {
	return float32(math.Sqrt(float64(v.X*v.X + v.Y*v.Y + v.Z*v.Z)))
}

// LengthSquared returns the squared magnitude
func (v Vec3) LengthSquared() float32 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z
}

// Normalize returns a normalized copy of the vector
func (v Vec3) Normalize() Vec3 {
	length := v.Length()
	if length < EPSILON {
		return Zero3()
	}
	invLength := 1.0 / length
	return Vec3{v.X * invLength, v.Y * invLength, v.Z * invLength}
}

// Distance returns the distance between two vectors
func (v Vec3) Distance(other Vec3) float32 {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	return float32(math.Sqrt(float64(dx*dx + dy*dy + dz*dz)))
}

// DistanceSquared returns the squared distance
func (v Vec3) DistanceSquared(other Vec3) float32 {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	return dx*dx + dy*dy + dz*dz
}

// Lerp linearly interpolates between two vectors
func (v Vec3) Lerp(other Vec3, t float32) Vec3 {
	return Vec3{
		v.X + (other.X-v.X)*t,
		v.Y + (other.Y-v.Y)*t,
		v.Z + (other.Z-v.Z)*t,
	}
}

// Slerp spherically interpolates between two vectors
func (v Vec3) Slerp(other Vec3, t float32) Vec3 {
	dot := v.Dot(other)
	// Clamp dot to [-1, 1]
	if dot > 1 {
		dot = 1
	} else if dot < -1 {
		dot = -1
	}

	theta := float32(math.Acos(float64(dot))) * t
	relative := other.Sub(v.MulScalar(dot)).Normalize()

	return v.MulScalar(float32(math.Cos(float64(theta)))).Add(relative.MulScalar(float32(math.Sin(float64(theta)))))
}

// Angle returns the angle between two vectors in radians
func (v Vec3) Angle(other Vec3) float32 {
	dot := v.Dot(other)
	lengthProduct := v.Length() * other.Length()
	if lengthProduct < EPSILON {
		return 0
	}
	return float32(math.Acos(float64(dot / lengthProduct)))
}

// Reflect reflects the vector across a normal
func (v Vec3) Reflect(normal Vec3) Vec3 {
	dot := v.Dot(normal)
	return Vec3{
		v.X - 2*dot*normal.X,
		v.Y - 2*dot*normal.Y,
		v.Z - 2*dot*normal.Z,
	}
}

// Refract refracts the vector
func (v Vec3) Refract(normal Vec3, eta float32) Vec3 {
	dot := v.Dot(normal)
	k := 1.0 - eta*eta*(1.0-dot*dot)
	if k < 0 {
		return Zero3()
	}
	return v.MulScalar(eta).Sub(normal.MulScalar(eta*dot + float32(math.Sqrt(float64(k)))))
}

// Project projects the vector onto another
func (v Vec3) Project(onto Vec3) Vec3 {
	lengthSq := onto.LengthSquared()
	if lengthSq < EPSILON {
		return Zero3()
	}
	return onto.MulScalar(v.Dot(onto) / lengthSq)
}

// Reject returns the rejection of the vector from another
func (v Vec3) Reject(onto Vec3) Vec3 {
	return v.Sub(v.Project(onto))
}

// Floor returns the vector with components rounded down
func (v Vec3) Floor() Vec3 {
	return Vec3{
		float32(math.Floor(float64(v.X))),
		float32(math.Floor(float64(v.Y))),
		float32(math.Floor(float64(v.Z))),
	}
}

// Ceil returns the vector with components rounded up
func (v Vec3) Ceil() Vec3 {
	return Vec3{
		float32(math.Ceil(float64(v.X))),
		float32(math.Ceil(float64(v.Y))),
		float32(math.Ceil(float64(v.Z))),
	}
}

// Round returns the vector with components rounded
func (v Vec3) Round() Vec3 {
	return Vec3{
		float32(math.Round(float64(v.X))),
		float32(math.Round(float64(v.Y))),
		float32(math.Round(float64(v.Z))),
	}
}

// Abs returns the vector with absolute values
func (v Vec3) Abs() Vec3 {
	return Vec3{
		float32(math.Abs(float64(v.X))),
		float32(math.Abs(float64(v.Y))),
		float32(math.Abs(float64(v.Z))),
	}
}

// Min returns component-wise minimum
func (v Vec3) Min(other Vec3) Vec3 {
	return Vec3{
		float32(math.Min(float64(v.X), float64(other.X))),
		float32(math.Min(float64(v.Y), float64(other.Y))),
		float32(math.Min(float64(v.Z), float64(other.Z))),
	}
}

// Max returns component-wise maximum
func (v Vec3) Max(other Vec3) Vec3 {
	return Vec3{
		float32(math.Max(float64(v.X), float64(other.X))),
		float32(math.Max(float64(v.Y), float64(other.Y))),
		float32(math.Max(float64(v.Z), float64(other.Z))),
	}
}

// Clamp clamps the vector components between min and max
func (v Vec3) Clamp(min, max Vec3) Vec3 {
	return Vec3{
		float32(math.Max(float64(min.X), math.Min(float64(max.X), float64(v.X)))),
		float32(math.Max(float64(min.Y), math.Min(float64(max.Y), float64(v.Y)))),
		float32(math.Max(float64(min.Z), math.Min(float64(max.Z), float64(v.Z)))),
	}
}

// ClampScalar clamps the vector components between min and max scalars
func (v Vec3) ClampScalar(min, max float32) Vec3 {
	return Vec3{
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.X)))),
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.Y)))),
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.Z)))),
	}
}

// ClampLength clamps the vector length
func (v Vec3) ClampLength(min, max float32) Vec3 {
	length := v.Length()
	if length < EPSILON {
		return Zero3()
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
func (v Vec3) Equals(other Vec3, epsilon ...float32) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return float32(math.Abs(float64(v.X-other.X))) <= eps &&
		float32(math.Abs(float64(v.Y-other.Y))) <= eps &&
		float32(math.Abs(float64(v.Z-other.Z))) <= eps
}

// IsZero checks if the vector is approximately zero
func (v Vec3) IsZero(epsilon ...float32) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return float32(math.Abs(float64(v.X))) <= eps &&
		float32(math.Abs(float64(v.Y))) <= eps &&
		float32(math.Abs(float64(v.Z))) <= eps
}

// Zero sets the vector to zero
func (v *Vec3) Zero() {
	v.X, v.Y, v.Z = 0, 0, 0
}

// Set sets the vector components
func (v *Vec3) Set(x, y, z float32) {
	v.X, v.Y, v.Z = x, y, z
}

// Copy copies another vector
func (v *Vec3) Copy(other Vec3) {
	v.X, v.Y, v.Z = other.X, other.Y, other.Z
}

// Clone returns a copy of the vector
func (v Vec3) Clone() Vec3 {
	return Vec3{v.X, v.Y, v.Z}
}

// ToArray returns the vector as a slice
func (v Vec3) ToArray() []float32 {
	return []float32{v.X, v.Y, v.Z}
}

// ToVec2 converts to Vec2, dropping Z
func (v Vec3) ToVec2() Vec2 {
	return Vec2{v.X, v.Y}
}

// ToVec4 converts to Vec4 with w = 1
func (v Vec3) ToVec4() Vec4 {
	return Vec4{v.X, v.Y, v.Z, 1}
}

// String returns a string representation
func (v Vec3) String() string {
	return fmt.Sprintf("Vec3(%.6f, %.6f, %.6f)", v.X, v.Y, v.Z)
}

// MarshalJSON implements json.Marshaler
func (v Vec3) MarshalJSON() ([]byte, error) {
	return []byte(fmt.Sprintf("[%f,%f,%f]", v.X, v.Y, v.Z)), nil
}

// UnmarshalJSON implements json.Unmarshaler
func (v *Vec3) UnmarshalJSON(data []byte) error {
	var arr [3]float32
	if err := json.Unmarshal(data, &arr); err != nil {
		return err
	}
	v.X, v.Y, v.Z = arr[0], arr[1], arr[2]
	return nil
}

// Transform transforms the vector by a 3x3 matrix
func (v Vec3) Transform(m [9]float32) Vec3 {
	return Vec3{
		v.X*m[0] + v.Y*m[3] + v.Z*m[6],
		v.X*m[1] + v.Y*m[4] + v.Z*m[7],
		v.X*m[2] + v.Y*m[5] + v.Z*m[8],
	}
}
