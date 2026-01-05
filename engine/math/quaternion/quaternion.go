package quaternion

import (
	"encoding/json"
	"fmt"
	"math"

	"github.com/AgentCooper/gographics/engine/math/vector"
)

// Quaternion represents a 3D rotation using a 4D hypercomplex number
// q = w + xi + yj + zk = (w, v) where v = (x, y, z)
type Quaternion[T vector.Number] struct {
	W, X, Y, Z T
}

// ============= Creation Functions =============

// NewQuaternion creates a new quaternion from components
func NewQuaternion[T vector.Number](w, x, y, z T) Quaternion[T] {
	return Quaternion[T]{w, x, y, z}
}

// IdentityQuaternion returns the identity quaternion (no rotation)
func IdentityQuaternion[T vector.Number]() Quaternion[T] {
	var one T = 1
	var zero T = 0
	return Quaternion[T]{one, zero, zero, zero}
}

// ZeroQuaternion returns a zero quaternion
func ZeroQuaternion[T vector.Number]() Quaternion[T] {
	var zero T = 0
	return Quaternion[T]{zero, zero, zero, zero}
}

// FromAxisAngle creates a quaternion from an axis-angle representation
// axis: unit vector (x, y, z)
// angle: rotation angle in radians
func FromAxisAngle[T vector.Number](axis vector.Vec3[T], angle T) Quaternion[T] {
	halfAngle := angle / 2
	sinHalf := T(math.Sin(float64(halfAngle)))
	cosHalf := T(math.Cos(float64(halfAngle)))

	return Quaternion[T]{
		W: cosHalf,
		X: axis.X * sinHalf,
		Y: axis.Y * sinHalf,
		Z: axis.Z * sinHalf,
	}
}

// FromEuler creates a quaternion from Euler angles (yaw, pitch, roll)
// yaw (y): rotation around Y axis (up/down)
// pitch (x): rotation around X axis (left/right)
// roll (z): rotation around Z axis (forward/backward)
// Uses YXZ convention (common in 3D graphics)
func FromEuler[T vector.Number](yaw, pitch, roll T) Quaternion[T] {
	cy := T(math.Cos(float64(yaw * 0.5)))
	sy := T(math.Sin(float64(yaw * 0.5)))
	cp := T(math.Cos(float64(pitch * 0.5)))
	sp := T(math.Sin(float64(pitch * 0.5)))
	cr := T(math.Cos(float64(roll * 0.5)))
	sr := T(math.Sin(float64(roll * 0.5)))

	return Quaternion[T]{
		W: cr*cp*cy + sr*sp*sy,
		X: sr*cp*cy - cr*sp*sy,
		Y: cr*sp*cy + sr*cp*sy,
		Z: cr*cp*sy - sr*sp*cy,
	}
}

// FromRotationMatrix creates a quaternion from a 3x3 rotation matrix
func FromRotationMatrix[T vector.Number](m [9]T) Quaternion[T] {
	trace := m[0] + m[4] + m[8]

	var q Quaternion[T]

	if trace > 0 {
		s := T(math.Sqrt(float64(trace+1))) * 2
		q.W = 0.25 * s
		q.X = (m[5] - m[7]) / s
		q.Y = (m[6] - m[2]) / s
		q.Z = (m[1] - m[3]) / s
	} else if m[0] > m[4] && m[0] > m[8] {
		s := T(math.Sqrt(float64(1+m[0]-m[4]-m[8]))) * 2
		q.W = (m[5] - m[7]) / s
		q.X = 0.25 * s
		q.Y = (m[1] + m[3]) / s
		q.Z = (m[6] + m[2]) / s
	} else if m[4] > m[8] {
		s := T(math.Sqrt(float64(1+m[4]-m[0]-m[8]))) * 2
		q.W = (m[6] - m[2]) / s
		q.X = (m[1] + m[3]) / s
		q.Y = 0.25 * s
		q.Z = (m[5] + m[7]) / s
	} else {
		s := T(math.Sqrt(float64(1+m[8]-m[0]-m[4]))) * 2
		q.W = (m[1] - m[3]) / s
		q.X = (m[6] + m[2]) / s
		q.Y = (m[5] + m[7]) / s
		q.Z = 0.25 * s
	}

	return q.Normalize()
}

// LookAt creates a quaternion that rotates from forward direction to target direction
func LookAt[T vector.Number](forward, target, up vector.Vec3[T]) Quaternion[T] {
	// Normalize inputs
	f := forward.Normalize()
	t := target.Normalize()
	u := up.Normalize()

	// Handle edge cases
	if f.Equals(t) {
		return IdentityQuaternion[T]()
	}

	if f.Equals(t.MulScalar(-1)) {
		// 180 degree rotation around up axis
		return FromAxisAngle(u, T(math.Pi))
	}

	// Calculate rotation axis
	axis := f.Cross(t)
	if axis.LengthSquared() < vector.EPSILON {
		// Vectors are colinear but not opposite
		axis = f.Cross(u)
	}
	axis = axis.Normalize()

	// Calculate rotation angle
	cosAngle := f.Dot(t)
	angle := T(math.Acos(float64(cosAngle)))

	return FromAxisAngle(axis, angle)
}

// ============= Basic Operations =============

// Add adds two quaternions
func (q Quaternion[T]) Add(other Quaternion[T]) Quaternion[T] {
	return Quaternion[T]{
		W: q.W + other.W,
		X: q.X + other.X,
		Y: q.Y + other.Y,
		Z: q.Z + other.Z,
	}
}

// Sub subtracts two quaternions
func (q Quaternion[T]) Sub(other Quaternion[T]) Quaternion[T] {
	return Quaternion[T]{
		W: q.W - other.W,
		X: q.X - other.X,
		Y: q.Y - other.Y,
		Z: q.Z - other.Z,
	}
}

// Mul multiplies two quaternions (Hamilton product)
// Note: Quaternion multiplication is NOT commutative
func (q Quaternion[T]) Mul(other Quaternion[T]) Quaternion[T] {
	return Quaternion[T]{
		W: q.W*other.W - q.X*other.X - q.Y*other.Y - q.Z*other.Z,
		X: q.W*other.X + q.X*other.W + q.Y*other.Z - q.Z*other.Y,
		Y: q.W*other.Y - q.X*other.Z + q.Y*other.W + q.Z*other.X,
		Z: q.W*other.Z + q.X*other.Y - q.Y*other.X + q.Z*other.W,
	}
}

// MulScalar multiplies quaternion by a scalar
func (q Quaternion[T]) MulScalar(s T) Quaternion[T] {
	return Quaternion[T]{
		W: q.W * s,
		X: q.X * s,
		Y: q.Y * s,
		Z: q.Z * s,
	}
}

// DivScalar divides quaternion by a scalar
func (q Quaternion[T]) DivScalar(s T) Quaternion[T] {
	return Quaternion[T]{
		W: q.W / s,
		X: q.X / s,
		Y: q.Y / s,
		Z: q.Z / s,
	}
}

// Conjugate returns the conjugate quaternion (w, -x, -y, -z)
func (q Quaternion[T]) Conjugate() Quaternion[T] {
	return Quaternion[T]{
		W: q.W,
		X: -q.X,
		Y: -q.Y,
		Z: -q.Z,
	}
}

// Inverse returns the inverse quaternion (conjugate / norm²)
func (q Quaternion[T]) Inverse() Quaternion[T] {
	normSq := q.NormSquared()
	if normSq == 0 {
		return ZeroQuaternion[T]()
	}
	return q.Conjugate().DivScalar(normSq)
}

// Norm returns the magnitude (length) of the quaternion
func (q Quaternion[T]) Norm() T {
	return T(math.Sqrt(float64(q.W*q.W + q.X*q.X + q.Y*q.Y + q.Z*q.Z)))
}

// NormSquared returns the squared magnitude
func (q Quaternion[T]) NormSquared() T {
	return q.W*q.W + q.X*q.X + q.Y*q.Y + q.Z*q.Z
}

// Normalize returns a unit quaternion (norm = 1)
func (q Quaternion[T]) Normalize() Quaternion[T] {
	norm := q.Norm()
	if norm < vector.EPSILON {
		return IdentityQuaternion[T]()
	}
	return q.DivScalar(norm)
}

// IsUnit returns true if the quaternion is a unit quaternion (within epsilon)
func (q Quaternion[T]) IsUnit(epsilon ...float64) bool {
	eps := vector.EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(q.NormSquared()-1)) <= eps
}

// IsIdentity returns true if this is the identity quaternion (within epsilon)
func (q Quaternion[T]) IsIdentity(epsilon ...float64) bool {
	eps := vector.EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	identity := IdentityQuaternion[T]()
	return math.Abs(float64(q.W-identity.W)) <= eps &&
		math.Abs(float64(q.X-identity.X)) <= eps &&
		math.Abs(float64(q.Y-identity.Y)) <= eps &&
		math.Abs(float64(q.Z-identity.Z)) <= eps
}

// ============= Rotation Operations =============

// Rotate rotates a 3D vector by this quaternion
// For a unit quaternion q, rotation is: v' = q * v * q⁻¹
func (q Quaternion[T]) Rotate(v vector.Vec3[T]) vector.Vec3[T] {
	// Ensure we have a unit quaternion
	unitQ := q.Normalize()

	// Convert vector to pure quaternion (w=0)
	p := Quaternion[T]{W: 0, X: v.X, Y: v.Y, Z: v.Z}

	// Apply rotation: v' = q * v * q⁻¹
	// Since q is unit, q⁻¹ = q*
	// Important: Multiply in correct order: q * v * q⁻¹
	rotated := unitQ.Mul(p).Mul(unitQ.Conjugate())

	return vector.Vec3[T]{X: rotated.X, Y: rotated.Y, Z: rotated.Z}
}

// Slerp spherically interpolates between two quaternions
func (q Quaternion[T]) Slerp(other Quaternion[T], t T) Quaternion[T] {
	// Ensure both are unit quaternions
	q1 := q.Normalize()
	q2 := other.Normalize()

	// Calculate cosine of the angle between quaternions
	cosOmega := q1.Dot(q2)

	// If cosine is negative, negate one quaternion to take shorter path
	if cosOmega < 0 {
		q2 = q2.MulScalar(-1)
		cosOmega = -cosOmega
	}

	// Check if quaternions are very close
	k0, k1 := T(1.0), T(0.0)
	if cosOmega < 0.9999 {
		// Spherical interpolation
		sinOmega := T(math.Sqrt(float64(1 - cosOmega*cosOmega)))
		omega := T(math.Atan2(float64(sinOmega), float64(cosOmega)))
		oneOverSinOmega := 1.0 / sinOmega

		k0 = T(math.Sin(float64((1-t)*omega))) * oneOverSinOmega
		k1 = T(math.Sin(float64(t*omega))) * oneOverSinOmega
	} else {
		// Linear interpolation for very close quaternions
		k0 = 1 - t
		k1 = t
	}

	// Interpolate and normalize
	result := q1.MulScalar(k0).Add(q2.MulScalar(k1))
	return result.Normalize()
}

// Nlerp linearly interpolates and normalizes (faster but not constant speed)
func (q Quaternion[T]) Nlerp(other Quaternion[T], t T) Quaternion[T] {
	// Ensure both are unit quaternions
	q1 := q.Normalize()
	q2 := other.Normalize()

	// If dot product is negative, take the shorter path
	if q1.Dot(q2) < 0 {
		q2 = q2.MulScalar(-1)
	}

	// Linear interpolation
	result := q1.MulScalar(1 - t).Add(q2.MulScalar(t))
	return result.Normalize()
}

// Dot returns the dot product of two quaternions
func (q Quaternion[T]) Dot(other Quaternion[T]) T {
	return q.W*other.W + q.X*other.X + q.Y*other.Y + q.Z*other.Z
}

// Angle returns the angle between two quaternions in radians
func (q Quaternion[T]) Angle(other Quaternion[T]) T {
	dot := math.Abs(float64(q.Normalize().Dot(other.Normalize())))
	// Clamp to [-1, 1] due to floating point errors
	if dot > 1 {
		dot = 1
	}
	return T(2 * math.Acos(dot))
}

// ============= Conversion Functions =============

// ToAxisAngle converts quaternion to axis-angle representation
// Returns: (axis Vec3, angle float64)
func (q Quaternion[T]) ToAxisAngle() (vector.Vec3[T], T) {
	q = q.Normalize()

	angle := 2 * T(math.Acos(float64(q.W)))

	// Handle the singularity at angle = 0
	if angle < vector.EPSILON {
		return vector.Vec3[T]{X: 1, Y: 0, Z: 0}, 0
	}

	// Calculate axis
	s := T(math.Sqrt(float64(1 - q.W*q.W)))
	if s < vector.EPSILON {
		return vector.Vec3[T]{X: 1, Y: 0, Z: 0}, angle
	}

	axis := vector.Vec3[T]{
		X: q.X / s,
		Y: q.Y / s,
		Z: q.Z / s,
	}

	return axis.Normalize(), angle
}

// ToEuler converts quaternion to Euler angles (yaw, pitch, roll)
func (q Quaternion[T]) ToEuler() (yaw, pitch, roll T) {
	q = q.Normalize()

	// Roll (x-axis rotation)
	sinrCosp := 2 * (q.W*q.X + q.Y*q.Z)
	cosrCosp := 1 - 2*(q.X*q.X+q.Y*q.Y)
	roll = T(math.Atan2(float64(sinrCosp), float64(cosrCosp)))

	// Pitch (y-axis rotation)
	sinp := 2 * (q.W*q.Y - q.Z*q.X)
	if math.Abs(float64(sinp)) >= 1 {
		// Use 90 degrees if out of range
		pitch = T(math.Copysign(math.Pi/2, float64(sinp)))
	} else {
		pitch = T(math.Asin(float64(sinp)))
	}

	// Yaw (z-axis rotation)
	sinyCosp := 2 * (q.W*q.Z + q.X*q.Y)
	cosyCosp := 1 - 2*(q.Y*q.Y+q.Z*q.Z)
	yaw = T(math.Atan2(float64(sinyCosp), float64(cosyCosp)))

	return yaw, pitch, roll
}

// ToRotationMatrix converts quaternion to a 3x3 rotation matrix
func (q Quaternion[T]) ToRotationMatrix() [9]T {
	q = q.Normalize()

	xx := q.X * q.X
	yy := q.Y * q.Y
	zz := q.Z * q.Z
	xy := q.X * q.Y
	xz := q.X * q.Z
	yz := q.Y * q.Z
	wx := q.W * q.X
	wy := q.W * q.Y
	wz := q.W * q.Z

	return [9]T{
		1 - 2*(yy+zz), 2 * (xy + wz), 2 * (xz - wy),
		2 * (xy - wz), 1 - 2*(xx+zz), 2 * (yz + wx),
		2 * (xz + wy), 2 * (yz - wx), 1 - 2*(xx+yy),
	}
}

// ToRotationMatrix4x4 converts quaternion to a 4x4 homogeneous rotation matrix
func (q Quaternion[T]) ToRotationMatrix4x4() [16]T {
	q = q.Normalize()

	xx := q.X * q.X
	yy := q.Y * q.Y
	zz := q.Z * q.Z
	xy := q.X * q.Y
	xz := q.X * q.Z
	yz := q.Y * q.Z
	wx := q.W * q.X
	wy := q.W * q.Y
	wz := q.W * q.Z

	return [16]T{
		1 - 2*(yy+zz), 2 * (xy - wz), 2 * (xz + wy), 0,
		2 * (xy + wz), 1 - 2*(xx+zz), 2 * (yz - wx), 0,
		2 * (xz - wy), 2 * (yz + wx), 1 - 2*(xx+yy), 0,
		0, 0, 0, 1,
	}
}

// ============= Utility Functions =============

// Equals checks if two quaternions are approximately equal
func (q Quaternion[T]) Equals(other Quaternion[T], epsilon ...float64) bool {
	eps := vector.EPSILON
	if len(epsilon) > 0 {
		eps = epsilon[0]
	}
	return math.Abs(float64(q.W-other.W)) <= eps &&
		math.Abs(float64(q.X-other.X)) <= eps &&
		math.Abs(float64(q.Y-other.Y)) <= eps &&
		math.Abs(float64(q.Z-other.Z)) <= eps
}

// Clone returns a copy of the quaternion
func (q Quaternion[T]) Clone() Quaternion[T] {
	return Quaternion[T]{q.W, q.X, q.Y, q.Z}
}

// Set sets the quaternion components
func (q *Quaternion[T]) Set(w, x, y, z T) {
	q.W, q.X, q.Y, q.Z = w, x, y, z
}

// Copy copies another quaternion
func (q *Quaternion[T]) Copy(other Quaternion[T]) {
	q.W, q.X, q.Y, q.Z = other.W, other.X, other.Y, other.Z
}

// Zero sets the quaternion to zero
func (q *Quaternion[T]) Zero() {
	var zero T = 0
	q.W, q.X, q.Y, q.Z = zero, zero, zero, zero
}

// Identity sets the quaternion to identity
func (q *Quaternion[T]) Identity() {
	var zero T = 0
	var one T = 1
	q.W, q.X, q.Y, q.Z = one, zero, zero, zero
}

// String returns a string representation
func (q Quaternion[T]) String() string {
	return fmt.Sprintf("Quaternion[%T](%.4f, %.4f, %.4f, %.4f)", q.W, q.W, q.X, q.Y, q.Z)
}

// ToArray returns the quaternion as a slice
func (q Quaternion[T]) ToArray() []T {
	return []T{q.W, q.X, q.Y, q.Z}
}

// ============= Factory Functions for Common Rotations =============

// RotationX creates a quaternion representing rotation around X axis
func RotationX[T vector.Number](angle T) Quaternion[T] {
	halfAngle := angle / 2
	return Quaternion[T]{
		W: T(math.Cos(float64(halfAngle))),
		X: T(math.Sin(float64(halfAngle))),
		Y: 0,
		Z: 0,
	}
}

// RotationY creates a quaternion representing rotation around Y axis
func RotationY[T vector.Number](angle T) Quaternion[T] {
	halfAngle := angle / 2
	return Quaternion[T]{
		W: T(math.Cos(float64(halfAngle))),
		X: 0,
		Y: T(math.Sin(float64(halfAngle))),
		Z: 0,
	}
}

// RotationZ creates a quaternion representing rotation around Z axis
func RotationZ[T vector.Number](angle T) Quaternion[T] {
	halfAngle := angle / 2
	return Quaternion[T]{
		W: T(math.Cos(float64(halfAngle))),
		X: 0,
		Y: 0,
		Z: T(math.Sin(float64(halfAngle))),
	}
}

// RotationBetween creates a quaternion that rotates from vector a to vector b
func RotationBetween[T vector.Number](a, b vector.Vec3[T]) Quaternion[T] {
	a = a.Normalize()
	b = b.Normalize()

	dot := a.Dot(b)

	// Check for colinear vectors
	if dot > 0.999999 {
		// Same direction
		return IdentityQuaternion[T]()
	}

	if dot < -0.999999 {
		// Opposite direction, rotate 180 degrees around any orthogonal axis
		axis := a.Cross(vector.Vec3[T]{X: 1, Y: 0, Z: 0})
		if axis.LengthSquared() < vector.EPSILON {
			axis = a.Cross(vector.Vec3[T]{X: 0, Y: 1, Z: 0})
		}
		axis = axis.Normalize()
		return FromAxisAngle(axis, T(math.Pi))
	}

	// Normal case
	axis := a.Cross(b)
	s := T(math.Sqrt(float64((1 + dot) * 2)))
	invS := 1.0 / s

	return Quaternion[T]{
		W: s * 0.5,
		X: axis.X * invS,
		Y: axis.Y * invS,
		Z: axis.Z * invS,
	}
}

// ============= JSON Support =============

// MarshalJSON implements json.Marshaler
func (q Quaternion[T]) MarshalJSON() ([]byte, error) {
	return []byte(fmt.Sprintf("[%f,%f,%f,%f]", q.W, q.X, q.Y, q.Z)), nil
}

// UnmarshalJSON implements json.Unmarshaler
func (q *Quaternion[T]) UnmarshalJSON(data []byte) error {
	var arr [4]T
	if err := json.Unmarshal(data, &arr); err != nil {
		return err
	}
	q.W, q.X, q.Y, q.Z = arr[0], arr[1], arr[2], arr[3]
	return nil
}
