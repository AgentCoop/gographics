package vector

import (
	"encoding/json"
	"fmt"
	"math"
)

// Vec4 represents a 4D vector or homogeneous coordinates
type Vec4 struct {
	X, Y, Z, W float32
}

// NewVec4 creates a new Vec4
func NewVec4(x, y, z, w float32) Vec4 {
	return Vec4{x, y, z, w}
}

// NewVec4Scalar creates a Vec4 with all components set to s
func NewVec4Scalar(s float32) Vec4 {
	return Vec4{s, s, s, s}
}

// Zero4 returns a zero vector
func Zero4() Vec4 {
	return Vec4{0, 0, 0, 0}
}

// One4 returns a vector with all components set to 1
func One4() Vec4 {
	return Vec4{1, 1, 1, 1}
}

// NewPosition creates a position vector (w = 1)
func NewPosition(x, y, z float32) Vec4 {
	return Vec4{x, y, z, 1}
}

// NewDirection creates a direction vector (w = 0)
func NewDirection(x, y, z float32) Vec4 {
	return Vec4{x, y, z, 0}
}

// Add adds two vectors
func (v Vec4) Add(other Vec4) Vec4 {
	return Vec4{
		v.X + other.X,
		v.Y + other.Y,
		v.Z + other.Z,
		v.W + other.W,
	}
}

// Sub subtracts two vectors
func (v Vec4) Sub(other Vec4) Vec4 {
	return Vec4{
		v.X - other.X,
		v.Y - other.Y,
		v.Z - other.Z,
		v.W - other.W,
	}
}

// Mul multiplies two vectors component-wise
func (v Vec4) Mul(other Vec4) Vec4 {
	return Vec4{
		v.X * other.X,
		v.Y * other.Y,
		v.Z * other.Z,
		v.W * other.W,
	}
}

// Div divides two vectors component-wise
func (v Vec4) Div(other Vec4) Vec4 {
	return Vec4{
		v.X / other.X,
		v.Y / other.Y,
		v.Z / other.Z,
		v.W / other.W,
	}
}

// AddScalar adds a scalar to each component
func (v Vec4) AddScalar(s float32) Vec4 {
	return Vec4{v.X + s, v.Y + s, v.Z + s, v.W + s}
}

// SubScalar subtracts a scalar from each component
func (v Vec4) SubScalar(s float32) Vec4 {
	return Vec4{v.X - s, v.Y - s, v.Z - s, v.W - s}
}

// MulScalar multiplies each component by a scalar
func (v Vec4) MulScalar(s float32) Vec4 {
	return Vec4{v.X * s, v.Y * s, v.Z * s, v.W * s}
}

// DivScalar divides each component by a scalar
func (v Vec4) DivScalar(s float32) Vec4 {
	return Vec4{v.X / s, v.Y / s, v.Z / s, v.W / s}
}

// Dot returns the dot product of two vectors
func (v Vec4) Dot(other Vec4) float32 {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z + v.W*other.W
}

// Length returns the magnitude of the vector
func (v Vec4) Length() float32 {
	return float32(math.Sqrt(float64(v.X*v.X + v.Y*v.Y + v.Z*v.Z + v.W*v.W)))
}

// LengthSquared returns the squared magnitude
func (v Vec4) LengthSquared() float32 {
	return v.X*v.X + v.Y*v.Y + v.Z*v.Z + v.W*v.W
}

// Normalize returns a normalized copy of the vector
func (v Vec4) Normalize() Vec4 {
	length := v.Length()
	if length < EPSILON {
		return Zero4()
	}
	invLength := 1.0 / length
	return Vec4{v.X * invLength, v.Y * invLength, v.Z * invLength, v.W * invLength}
}

// Distance returns the distance between two vectors
func (v Vec4) Distance(other Vec4) float32 {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	dw := v.W - other.W
	return float32(math.Sqrt(float64(dx*dx + dy*dy + dz*dz + dw*dw)))
}

// DistanceSquared returns the squared distance
func (v Vec4) DistanceSquared(other Vec4) float32 {
	dx := v.X - other.X
	dy := v.Y - other.Y
	dz := v.Z - other.Z
	dw := v.W - other.W
	return dx*dx + dy*dy + dz*dz + dw*dw
}

// Lerp linearly interpolates between two vectors
func (v Vec4) Lerp(other Vec4, t float32) Vec4 {
	return Vec4{
		v.X + (other.X-v.X)*t,
		v.Y + (other.Y-v.Y)*t,
		v.Z + (other.Z-v.Z)*t,
		v.W + (other.W-v.W)*t,
	}
}

// Homogenize divides by W component for perspective division
func (v Vec4) Homogenize() Vec4 {
	if float32(math.Abs(float64(v.W))) < EPSILON {
		return Zero4()
	}
	invW := 1.0 / v.W
	return Vec4{v.X * invW, v.Y * invW, v.Z * invW, 1}
}

// IsPoint returns true if W == 1 (or approximately 1)
func (v Vec4) IsPoint(epsilon ...float32) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return float32(math.Abs(float64(v.W-1))) <= eps
}

// IsDirection returns true if W == 0 (or approximately 0)
func (v Vec4) IsDirection(epsilon ...float32) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return float32(math.Abs(float64(v.W))) <= eps
}

// ToPosition converts to position (W = 1)
func (v Vec4) ToPosition() Vec4 {
	if v.IsDirection() {
		return Vec4{v.X, v.Y, v.Z, 1}
	}
	return v.Homogenize()
}

// ToDirection converts to direction (W = 0)
func (v Vec4) ToDirection() Vec4 {
	if v.IsDirection() {
		return v
	}
	return Vec4{v.X, v.Y, v.Z, 0}
}

// Floor returns the vector with components rounded down
func (v Vec4) Floor() Vec4 {
	return Vec4{
		float32(math.Floor(float64(v.X))),
		float32(math.Floor(float64(v.Y))),
		float32(math.Floor(float64(v.Z))),
		float32(math.Floor(float64(v.W))),
	}
}

// Ceil returns the vector with components rounded up
func (v Vec4) Ceil() Vec4 {
	return Vec4{
		float32(math.Ceil(float64(v.X))),
		float32(math.Ceil(float64(v.Y))),
		float32(math.Ceil(float64(v.Z))),
		float32(math.Ceil(float64(v.W))),
	}
}

// Round returns the vector with components rounded
func (v Vec4) Round() Vec4 {
	return Vec4{
		float32(math.Round(float64(v.X))),
		float32(math.Round(float64(v.Y))),
		float32(math.Round(float64(v.Z))),
		float32(math.Round(float64(v.W))),
	}
}

// Abs returns the vector with absolute values
func (v Vec4) Abs() Vec4 {
	return Vec4{
		float32(math.Abs(float64(v.X))),
		float32(math.Abs(float64(v.Y))),
		float32(math.Abs(float64(v.Z))),
		float32(math.Abs(float64(v.W))),
	}
}

// Min returns component-wise minimum
func (v Vec4) Min(other Vec4) Vec4 {
	return Vec4{
		float32(math.Min(float64(v.X), float64(other.X))),
		float32(math.Min(float64(v.Y), float64(other.Y))),
		float32(math.Min(float64(v.Z), float64(other.Z))),
		float32(math.Min(float64(v.W), float64(other.W))),
	}
}

// Max returns component-wise maximum
func (v Vec4) Max(other Vec4) Vec4 {
	return Vec4{
		float32(math.Max(float64(v.X), float64(other.X))),
		float32(math.Max(float64(v.Y), float64(other.Y))),
		float32(math.Max(float64(v.Z), float64(other.Z))),
		float32(math.Max(float64(v.W), float64(other.W))),
	}
}

// Clamp clamps the vector components between min and max
func (v Vec4) Clamp(min, max Vec4) Vec4 {
	return Vec4{
		float32(math.Max(float64(min.X), math.Min(float64(max.X), float64(v.X)))),
		float32(math.Max(float64(min.Y), math.Min(float64(max.Y), float64(v.Y)))),
		float32(math.Max(float64(min.Z), math.Min(float64(max.Z), float64(v.Z)))),
		float32(math.Max(float64(min.W), math.Min(float64(max.W), float64(v.W)))),
	}
}

// ClampScalar clamps the vector components between min and max scalars
func (v Vec4) ClampScalar(min, max float32) Vec4 {
	return Vec4{
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.X)))),
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.Y)))),
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.Z)))),
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.W)))),
	}
}

// ClampLength clamps the vector length
func (v Vec4) ClampLength(min, max float32) Vec4 {
	length := v.Length()
	if length < EPSILON {
		return Zero4()
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
func (v Vec4) Equals(other Vec4, epsilon ...float32) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return float32(math.Abs(float64(v.X-other.X))) <= eps &&
		float32(math.Abs(float64(v.Y-other.Y))) <= eps &&
		float32(math.Abs(float64(v.Z-other.Z))) <= eps &&
		float32(math.Abs(float64(v.W-other.W))) <= eps
}

// IsZero checks if the vector is approximately zero
func (v Vec4) IsZero(epsilon ...float32) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return float32(math.Abs(float64(v.X))) <= eps &&
		float32(math.Abs(float64(v.Y))) <= eps &&
		float32(math.Abs(float64(v.Z))) <= eps &&
		float32(math.Abs(float64(v.W))) <= eps
}

// Zero sets the vector to zero
func (v *Vec4) Zero() {
	v.X, v.Y, v.Z, v.W = 0, 0, 0, 0
}

// Set sets the vector components
func (v *Vec4) Set(x, y, z, w float32) {
	v.X, v.Y, v.Z, v.W = x, y, z, w
}

// Copy copies another vector
func (v *Vec4) Copy(other Vec4) {
	v.X, v.Y, v.Z, v.W = other.X, other.Y, other.Z, other.W
}

// Clone returns a copy of the vector
func (v Vec4) Clone() Vec4 {
	return Vec4{v.X, v.Y, v.Z, v.W}
}

// ToArray returns the vector as a slice
func (v Vec4) ToArray() []float32 {
	return []float32{v.X, v.Y, v.Z, v.W}
}

// ToVec2 converts to Vec2, dropping Z and W
func (v Vec4) ToVec2() Vec2 {
	return Vec2{v.X, v.Y}
}

// ToVec3 converts to Vec3, dropping W
func (v Vec4) ToVec3() Vec3 {
	if v.IsDirection() {
		return Vec3{v.X, v.Y, v.Z}
	}
	// Homogenize first for points
	homogenized := v.Homogenize()
	return Vec3{homogenized.X, homogenized.Y, homogenized.Z}
}

// String returns a string representation
func (v Vec4) String() string {
	return fmt.Sprintf("Vec4(%.6f, %.6f, %.6f, %.6f)", v.X, v.Y, v.Z, v.W)
}

// MarshalJSON implements json.Marshaler
func (v Vec4) MarshalJSON() ([]byte, error) {
	return []byte(fmt.Sprintf("[%f,%f,%f,%f]", v.X, v.Y, v.Z, v.W)), nil
}

// UnmarshalJSON implements json.Unmarshaler
func (v *Vec4) UnmarshalJSON(data []byte) error {
	var arr [4]float32
	if err := json.Unmarshal(data, &arr); err != nil {
		return err
	}
	v.X, v.Y, v.Z, v.W = arr[0], arr[1], arr[2], arr[3]
	return nil
}

// Transform transforms the vector by a 4x4 matrix
func (v Vec4) Transform(m [16]float32) Vec4 {
	return Vec4{
		v.X*m[0] + v.Y*m[4] + v.Z*m[8] + v.W*m[12],
		v.X*m[1] + v.Y*m[5] + v.Z*m[9] + v.W*m[13],
		v.X*m[2] + v.Y*m[6] + v.Z*m[10] + v.W*m[14],
		v.X*m[3] + v.Y*m[7] + v.Z*m[11] + v.W*m[15],
	}
}

// TransformAffine transforms by a 4x4 matrix assuming affine transformation
func (v Vec4) TransformAffine(m [16]float32) Vec4 {
	return Vec4{
		v.X*m[0] + v.Y*m[4] + v.Z*m[8] + v.W*m[12],
		v.X*m[1] + v.Y*m[5] + v.Z*m[9] + v.W*m[13],
		v.X*m[2] + v.Y*m[6] + v.Z*m[10] + v.W*m[14],
		v.W, // W unchanged for affine transformations
	}
}
