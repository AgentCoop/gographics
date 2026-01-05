package vector

import (
	"encoding/json"
	"fmt"
	"math"
)

// Vec2 represents a 2D vector
type Vec2 struct {
	X, Y float32
}

// NewVec2 creates a new Vec2
func NewVec2(x, y float32) Vec2 {
	return Vec2{x, y}
}

// NewVec2Scalar creates a Vec2 with both components set to s
func NewVec2Scalar(s float32) Vec2 {
	return Vec2{s, s}
}

// Zero2 returns a zero vector
func Zero2() Vec2 {
	return Vec2{0, 0}
}

// One2 returns a vector with all components set to 1
func One2() Vec2 {
	return Vec2{1, 1}
}

// UnitX2 returns the X unit vector
func UnitX2() Vec2 {
	return Vec2{1, 0}
}

// UnitY2 returns the Y unit vector
func UnitY2() Vec2 {
	return Vec2{0, 1}
}

// Add adds two vectors
func (v Vec2) Add(other Vec2) Vec2 {
	return Vec2{v.X + other.X, v.Y + other.Y}
}

// Sub subtracts two vectors
func (v Vec2) Sub(other Vec2) Vec2 {
	return Vec2{v.X - other.X, v.Y - other.Y}
}

// Mul multiplies two vectors component-wise
func (v Vec2) Mul(other Vec2) Vec2 {
	return Vec2{v.X * other.X, v.Y * other.Y}
}

// Div divides two vectors component-wise
func (v Vec2) Div(other Vec2) Vec2 {
	return Vec2{v.X / other.X, v.Y / other.Y}
}

// AddScalar adds a scalar to each component
func (v Vec2) AddScalar(s float32) Vec2 {
	return Vec2{v.X + s, v.Y + s}
}

// SubScalar subtracts a scalar from each component
func (v Vec2) SubScalar(s float32) Vec2 {
	return Vec2{v.X - s, v.Y - s}
}

// MulScalar multiplies each component by a scalar
func (v Vec2) MulScalar(s float32) Vec2 {
	return Vec2{v.X * s, v.Y * s}
}

// DivScalar divides each component by a scalar
func (v Vec2) DivScalar(s float32) Vec2 {
	return Vec2{v.X / s, v.Y / s}
}

// Dot returns the dot product of two vectors
func (v Vec2) Dot(other Vec2) float32 {
	return v.X*other.X + v.Y*other.Y
}

// Cross returns the cross product (scalar in 2D)
func (v Vec2) Cross(other Vec2) float32 {
	return v.X*other.Y - v.Y*other.X
}

// Length returns the magnitude of the vector
func (v Vec2) Length() float32 {
	return float32(math.Sqrt(float64(v.X*v.X + v.Y*v.Y)))
}

// LengthSquared returns the squared magnitude
func (v Vec2) LengthSquared() float32 {
	return v.X*v.X + v.Y*v.Y
}

// Normalize returns a normalized copy of the vector
func (v Vec2) Normalize() Vec2 {
	length := v.Length()
	if length < EPSILON {
		return Zero2()
	}
	invLength := 1.0 / length
	return Vec2{v.X * invLength, v.Y * invLength}
}

// Distance returns the distance between two vectors
func (v Vec2) Distance(other Vec2) float32 {
	return v.Sub(other).Length()
}

// DistanceSquared returns the squared distance
func (v Vec2) DistanceSquared(other Vec2) float32 {
	dx := v.X - other.X
	dy := v.Y - other.Y
	return dx*dx + dy*dy
}

// Lerp linearly interpolates between two vectors
func (v Vec2) Lerp(other Vec2, t float32) Vec2 {
	return Vec2{
		v.X + (other.X-v.X)*t,
		v.Y + (other.Y-v.Y)*t,
	}
}

// Slerp spherically interpolates between two vectors
func (v Vec2) Slerp(other Vec2, t float32) Vec2 {
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

// Rotate rotates the vector by angle in radians
func (v Vec2) Rotate(angle float32) Vec2 {
	cos := float32(math.Cos(float64(angle)))
	sin := float32(math.Sin(float64(angle)))
	return Vec2{
		v.X*cos - v.Y*sin,
		v.X*sin + v.Y*cos,
	}
}

// Angle returns the angle of the vector in radians
func (v Vec2) Angle() float32 {
	return float32(math.Atan2(float64(v.Y), float64(v.X)))
}

// AngleTo returns the angle between two vectors
func (v Vec2) AngleTo(other Vec2) float32 {
	dot := v.Dot(other)
	det := v.Cross(other)
	return float32(math.Atan2(float64(det), float64(dot)))
}

// Perpendicular returns a perpendicular vector
func (v Vec2) Perpendicular() Vec2 {
	return Vec2{-v.Y, v.X}
}

// Reflect reflects the vector across a normal
func (v Vec2) Reflect(normal Vec2) Vec2 {
	dot := v.Dot(normal)
	return Vec2{
		v.X - 2*dot*normal.X,
		v.Y - 2*dot*normal.Y,
	}
}

// Project projects the vector onto another
func (v Vec2) Project(onto Vec2) Vec2 {
	lengthSq := onto.LengthSquared()
	if lengthSq < EPSILON {
		return Zero2()
	}
	return onto.MulScalar(v.Dot(onto) / lengthSq)
}

// Floor returns the vector with components rounded down
func (v Vec2) Floor() Vec2 {
	return Vec2{float32(math.Floor(float64(v.X))), float32(math.Floor(float64(v.Y)))}
}

// Ceil returns the vector with components rounded up
func (v Vec2) Ceil() Vec2 {
	return Vec2{float32(math.Ceil(float64(v.X))), float32(math.Ceil(float64(v.Y)))}
}

// Round returns the vector with components rounded
func (v Vec2) Round() Vec2 {
	return Vec2{float32(math.Round(float64(v.X))), float32(math.Round(float64(v.Y)))}
}

// Abs returns the vector with absolute values
func (v Vec2) Abs() Vec2 {
	return Vec2{float32(math.Abs(float64(v.X))), float32(math.Abs(float64(v.Y)))}
}

// Min returns component-wise minimum
func (v Vec2) Min(other Vec2) Vec2 {
	return Vec2{float32(math.Min(float64(v.X), float64(other.X))), float32(math.Min(float64(v.Y), float64(other.Y)))}
}

// Max returns component-wise maximum
func (v Vec2) Max(other Vec2) Vec2 {
	return Vec2{float32(math.Max(float64(v.X), float64(other.X))), float32(math.Max(float64(v.Y), float64(other.Y)))}
}

// Clamp clamps the vector components between min and max
func (v Vec2) Clamp(min, max Vec2) Vec2 {
	return Vec2{
		float32(math.Max(float64(min.X), math.Min(float64(max.X), float64(v.X)))),
		float32(math.Max(float64(min.Y), math.Min(float64(max.Y), float64(v.Y)))),
	}
}

// ClampScalar clamps the vector components between min and max scalars
func (v Vec2) ClampScalar(min, max float32) Vec2 {
	return Vec2{
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.X)))),
		float32(math.Max(float64(min), math.Min(float64(max), float64(v.Y)))),
	}
}

// ClampLength clamps the vector length
func (v Vec2) ClampLength(min, max float32) Vec2 {
	length := v.Length()
	if length < EPSILON {
		return Zero2()
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
func (v Vec2) Equals(other Vec2, epsilon ...float32) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return float32(math.Abs(float64(v.X-other.X))) <= eps && float32(math.Abs(float64(v.Y-other.Y))) <= eps
}

// IsZero checks if the vector is approximately zero
func (v Vec2) IsZero(epsilon ...float32) bool {
	eps := EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return float32(math.Abs(float64(v.X))) <= eps && float32(math.Abs(float64(v.Y))) <= eps
}

// Zero sets the vector to zero
func (v *Vec2) Zero() {
	v.X, v.Y = 0, 0
}

// Set sets the vector components
func (v *Vec2) Set(x, y float32) {
	v.X, v.Y = x, y
}

// Copy copies another vector
func (v *Vec2) Copy(other Vec2) {
	v.X, v.Y = other.X, other.Y
}

// Clone returns a copy of the vector
func (v Vec2) Clone() Vec2 {
	return Vec2{v.X, v.Y}
}

// ToArray returns the vector as a slice
func (v Vec2) ToArray() []float32 {
	return []float32{v.X, v.Y}
}

// ToVec3 converts to Vec3 with z = 0
func (v Vec2) ToVec3() Vec3 {
	return Vec3{v.X, v.Y, 0}
}

// ToVec4 converts to Vec4 with z = 0, w = 1
func (v Vec2) ToVec4() Vec4 {
	return Vec4{v.X, v.Y, 0, 1}
}

// String returns a string representation
func (v Vec2) String() string {
	return fmt.Sprintf("Vec2(%.6f, %.6f)", v.X, v.Y)
}

// MarshalJSON implements json.Marshaler
func (v Vec2) MarshalJSON() ([]byte, error) {
	return []byte(fmt.Sprintf("[%f,%f]", v.X, v.Y)), nil
}

// UnmarshalJSON implements json.Unmarshaler
func (v *Vec2) UnmarshalJSON(data []byte) error {
	var arr [2]float32
	if err := json.Unmarshal(data, &arr); err != nil {
		return err
	}
	v.X, v.Y = arr[0], arr[1]
	return nil
}

// Transform transforms the vector by a 2x2 matrix
func (v Vec2) Transform(m [4]float32) Vec2 {
	return Vec2{
		v.X*m[0] + v.Y*m[2],
		v.X*m[1] + v.Y*m[3],
	}
}
