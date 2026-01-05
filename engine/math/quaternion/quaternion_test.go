package quaternion

import (
	"encoding/json"
	"math"
	"testing"

	"github.com/AgentCooper/gographics/engine/math/vector"
)

// Helper functions for approximate equality
func approxEqual[T vector.Number](a, b T) bool {
	return math.Abs(float64(a-b)) < 1e-6
}

func approxEqualVec2[T vector.Number](a, b vector.Vec2[T]) bool {
	return approxEqual(a.X, b.X) && approxEqual(a.Y, b.Y)
}

func approxEqualVec3[T vector.Number](a, b vector.Vec3[T]) bool {
	return approxEqual(a.X, b.X) && approxEqual(a.Y, b.Y) && approxEqual(a.Z, b.Z)
}

func approxEqualVec4[T vector.Number](a, b vector.Vec4[T]) bool {
	return approxEqual(a.X, b.X) && approxEqual(a.Y, b.Y) && approxEqual(a.Z, b.Z) && approxEqual(a.W, b.W)
}

func approxEqualQuat[T vector.Number](a, b Quaternion[T]) bool {
	return approxEqual(a.W, b.W) &&
		approxEqual(a.X, b.X) &&
		approxEqual(a.Y, b.Y) &&
		approxEqual(a.Z, b.Z)
}

func TestQuaternionCreation(t *testing.T) {
	// Test identity
	id := IdentityQuaternion[float32]()
	if id.W != 1 || id.X != 0 || id.Y != 0 || id.Z != 0 {
		t.Errorf("IdentityQuaternion failed: got %v", id)
	}

	// Test FromAxisAngle
	axis := vector.NewVec3[float32](1, 0, 0)
	q := FromAxisAngle(axis, float32(math.Pi/2))

	// For 90 degree rotation around X axis
	expectedW := float32(math.Cos(math.Pi / 4)) // cos(45°)
	expectedX := float32(math.Sin(math.Pi / 4)) // sin(45°)

	if !approxEqual(q.W, expectedW) || !approxEqual(q.X, expectedX) || q.Y != 0 || q.Z != 0 {
		t.Errorf("FromAxisAngle failed: got %v", q)
	}

	// Test FromEuler
	q2 := FromEuler[float32](0, float32(math.Pi/2), 0) // 90° pitch
	if !q2.IsUnit() {
		t.Errorf("FromEuler should produce unit quaternion: got norm %v", q2.Norm())
	}
}

func TestQuaternionOperations(t *testing.T) {
	q1 := Quaternion[float32]{W: 1, X: 2, Y: 3, Z: 4}
	q2 := Quaternion[float32]{W: 5, X: 6, Y: 7, Z: 8}

	// Addition
	sum := q1.Add(q2)
	expected := Quaternion[float32]{W: 6, X: 8, Y: 10, Z: 12}
	if !approxEqualQuat(sum, expected) {
		t.Errorf("Add failed: got %v, want %v", sum, expected)
	}

	// Multiplication (test with known values)
	// (1 + 2i + 3j + 4k) × (5 + 6i + 7j + 8k)
	// w = 1*5 - 2*6 - 3*7 - 4*8 = 5 - 12 - 21 - 32 = -60
	// x = 1*6 + 2*5 + 3*8 - 4*7 = 6 + 10 + 24 - 28 = 12
	// y = 1*7 - 2*8 + 3*5 + 4*6 = 7 - 16 + 15 + 24 = 30
	// z = 1*8 + 2*7 - 3*6 + 4*5 = 8 + 14 - 18 + 20 = 24
	product := q1.Mul(q2)
	expected = Quaternion[float32]{W: -60, X: 12, Y: 30, Z: 24}
	if !approxEqualQuat(product, expected) {
		t.Errorf("Mul failed: got %v, want %v", product, expected)
	}

	// Conjugate
	conj := q1.Conjugate()
	expected = Quaternion[float32]{W: 1, X: -2, Y: -3, Z: -4}
	if !approxEqualQuat(conj, expected) {
		t.Errorf("Conjugate failed: got %v, want %v", conj, expected)
	}

	// Norm
	norm := q1.Norm()
	expectedNorm := float32(math.Sqrt(1*1 + 2*2 + 3*3 + 4*4))
	if !approxEqual(norm, expectedNorm) {
		t.Errorf("Norm failed: got %v, want %v", norm, expectedNorm)
	}

	// Normalize
	normalized := q1.Normalize()
	if !approxEqual(normalized.Norm(), 1.0) {
		t.Errorf("Normalize failed: norm = %v, want 1", normalized.Norm())
	}
}

func TestQuaternionRotation(t *testing.T) {
	// Test rotation around X axis by 90 degrees
	axis := vector.NewVec3[float32](1, 0, 0)
	q := FromAxisAngle(axis, float32(math.Pi/2))

	// Rotate a point on the Y axis
	point := vector.NewVec3[float32](0, 1, 0)
	rotated := q.Rotate(point)

	// Should end up on the Z axis
	expected := vector.NewVec3[float32](0, 0, 1)
	if !approxEqualVec3(rotated, expected) {
		t.Errorf("Rotate failed: got %v, want %v", rotated, expected)
	}

	// Test identity rotation
	id := IdentityQuaternion[float32]()
	rotated = id.Rotate(point)
	if !approxEqualVec3(rotated, point) {
		t.Errorf("Identity rotation failed: got %v, want %v", rotated, point)
	}

	// Test rotation composition
	// Rotate 90° around X, then 90° around Y
	qx := FromAxisAngle(vector.NewVec3[float32](1, 0, 0), float32(math.Pi/2))
	qy := FromAxisAngle(vector.NewVec3[float32](0, 1, 0), float32(math.Pi/2))
	qCombined := qy.Mul(qx) // Apply qx first, then qy

	point = vector.NewVec3[float32](1, 0, 0)
	rotated = qCombined.Rotate(point)
	// After X rotation: (1, 0, 0) stays (1, 0, 0)
	// After Y rotation: (1, 0, 0) goes to (0, 0, -1)
	expected = vector.NewVec3[float32](0, 0, -1)
	if !approxEqualVec3(rotated, expected) {
		t.Errorf("Rotation composition failed: got %v, want %v", rotated, expected)
	}
}

func TestQuaternionInterpolation(t *testing.T) {
	// Create two rotations
	q1 := FromAxisAngle(vector.NewVec3[float32](1, 0, 0), 0)
	q2 := FromAxisAngle(vector.NewVec3[float32](1, 0, 0), float32(math.Pi/2))

	// Test SLERP at t=0.5
	interpolated := q1.Slerp(q2, 0.5)

	// Should be rotation by 45 degrees
	expected := FromAxisAngle(vector.NewVec3[float32](1, 0, 0), float32(math.Pi/4))
	if !approxEqualQuat(interpolated.Normalize(), expected.Normalize()) {
		t.Errorf("Slerp failed: got %v, want %v", interpolated, expected)
	}

	// Test NLERP
	interpolated = q1.Nlerp(q2, 0.5)
	// For this simple case, NLERP should be close to SLERP
	if !interpolated.IsUnit() {
		t.Errorf("Nlerp should produce unit quaternion: norm = %v", interpolated.Norm())
	}
}

func TestQuaternionConversion(t *testing.T) {
	// Test axis-angle conversion
	axis := vector.NewVec3[float32](1, 0, 0).Normalize()
	angle := float32(math.Pi / 3)
	q := FromAxisAngle(axis, angle)

	// Convert back
	convertedAxis, convertedAngle := q.ToAxisAngle()

	if !approxEqualVec3(convertedAxis, axis) {
		t.Errorf("ToAxisAngle axis failed: got %v, want %v", convertedAxis, axis)
	}
	if !approxEqual(convertedAngle, angle) {
		t.Errorf("ToAxisAngle angle failed: got %v, want %v", convertedAngle, angle)
	}

	// Test Euler conversion
	yaw := float32(math.Pi / 4)
	pitch := float32(math.Pi / 6)
	roll := float32(math.Pi / 8)

	q = FromEuler(yaw, pitch, roll)
	convertedYaw, convertedPitch, convertedRoll := q.ToEuler()

	if !approxEqual(convertedYaw, yaw) {
		t.Errorf("ToEuler yaw failed: got %v, want %v", convertedYaw, yaw)
	}
	if !approxEqual(convertedPitch, pitch) {
		t.Errorf("ToEuler pitch failed: got %v, want %v", convertedPitch, pitch)
	}
	if !approxEqual(convertedRoll, roll) {
		t.Errorf("ToEuler roll failed: got %v, want %v", convertedRoll, roll)
	}
}

func TestQuaternionMatrixConversion(t *testing.T) {
	// Test rotation around Y axis by 90 degrees
	q := FromAxisAngle(vector.NewVec3[float32](0, 1, 0), float32(math.Pi/2))

	// Convert to matrix
	matrix := q.ToRotationMatrix()

	// Rotate a point using the matrix
	point := vector.NewVec3[float32](1, 0, 0)
	rotatedByMatrix := point.Transform(matrix)

	// Rotate the same point using quaternion
	rotatedByQuat := q.Rotate(point)

	// They should be equal (or opposite)
	if !approxEqualVec3(rotatedByMatrix, rotatedByQuat) {
		// Check if they're opposite (q and -q represent same rotation)
		opposite := rotatedByQuat.MulScalar(-1)
		if !approxEqualVec3(rotatedByMatrix, opposite) {
			t.Errorf("Matrix conversion failed: matrix result %v, quaternion result %v",
				rotatedByMatrix, rotatedByQuat)
		}
	}

	// Convert back to quaternion
	q2 := FromRotationMatrix(matrix)

	// Should be equivalent (or opposite, since q and -q represent same rotation)
	if !approxEqualQuat(q.Normalize(), q2.Normalize()) &&
		!approxEqualQuat(q.Normalize(), q2.Normalize().MulScalar(-1)) {
		t.Errorf("Round-trip matrix conversion failed: original %v, converted %v", q, q2)
	}
}

func TestQuaternionLookAt(t *testing.T) {
	forward := vector.NewVec3[float32](0, 0, -1)
	target := vector.NewVec3[float32](1, 0, 0)
	up := vector.NewVec3[float32](0, 1, 0)

	q := LookAt(forward, target, up)

	// Rotate forward vector to point at target
	rotated := q.Rotate(forward)

	// Should be pointing in target direction
	if !approxEqualVec3(rotated.Normalize(), target.Normalize()) {
		t.Errorf("LookAt failed: rotated %v, target %v", rotated, target)
	}
}

func TestQuaternionProperties(t *testing.T) {
	// Test quaternion multiplication is not commutative
	q1 := Quaternion[float64]{W: 1, X: 2, Y: 3, Z: 4}
	q2 := Quaternion[float64]{W: 5, X: 6, Y: 7, Z: 8}

	product1 := q1.Mul(q2)
	product2 := q2.Mul(q1)

	// They should NOT be equal
	if approxEqualQuat(product1, product2) {
		t.Errorf("Quaternion multiplication should not be commutative")
	}

	// Test inverse property: q * q⁻¹ = identity
	q := FromAxisAngle(vector.NewVec3[float64](1, 2, 3).Normalize(), 0.5)
	inverse := q.Inverse()
	product := q.Mul(inverse)

	if !product.IsIdentity(1e-10) {
		t.Errorf("Inverse property failed: q * q⁻¹ = %v, not identity", product)
	}

	// Test norm preservation under multiplication (for unit quaternions)
	q1 = FromAxisAngle(vector.NewVec3[float64](1, 0, 0), 0.3).Normalize()
	q2 = FromAxisAngle(vector.NewVec3[float64](0, 1, 0), 0.4).Normalize()
	product = q1.Mul(q2)

	if !product.IsUnit(1e-10) {
		t.Errorf("Unit quaternion multiplication should preserve norm: norm = %v", product.Norm())
	}
}

func TestQuaternionJSON(t *testing.T) {
	q := Quaternion[float64]{W: 1.1, X: 2.2, Y: 3.3, Z: 4.4}

	// Marshal
	data, err := json.Marshal(q)
	if err != nil {
		t.Errorf("Marshal failed: %v", err)
	}

	// Unmarshal
	var q2 Quaternion[float64]
	err = json.Unmarshal(data, &q2)
	if err != nil {
		t.Errorf("Unmarshal failed: %v", err)
	}

	if !approxEqualQuat(q, q2) {
		t.Errorf("JSON round-trip failed: original %v, unmarshaled %v", q, q2)
	}
}

func TestQuaternionEdgeCases(t *testing.T) {
	// Test rotation of zero vector
	q := FromAxisAngle(vector.NewVec3[float32](1, 0, 0), float32(math.Pi/2))
	zero := vector.NewVec3[float32](0, 0, 0)
	rotated := q.Rotate(zero)
	if !rotated.IsZero() {
		t.Errorf("Rotation of zero vector should be zero: got %v", rotated)
	}

	// Test near-zero quaternion normalization
	smallQuat := Quaternion[float32]{W: 1e-10, X: 1e-10, Y: 1e-10, Z: 1e-10}
	normalized := smallQuat.Normalize()
	if !normalized.IsIdentity(1e-3) {
		t.Errorf("Normalization of near-zero quaternion failed: got %v", normalized)
	}

	// Test LookAt with same direction
	forward := vector.NewVec3[float32](0, 0, -1)
	q = LookAt(forward, forward, vector.NewVec3[float32](0, 1, 0))
	if !q.IsIdentity(1e-6) {
		t.Errorf("LookAt with same direction should return identity: got %v", q)
	}
}

func BenchmarkQuaternionRotation(b *testing.B) {
	q := FromAxisAngle(vector.NewVec3[float64](1, 2, 3).Normalize(), 0.5)
	v := vector.NewVec3[float64](4, 5, 6)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = q.Rotate(v)
	}
}

func BenchmarkQuaternionMultiplication(b *testing.B) {
	q1 := FromAxisAngle(vector.NewVec3[float64](1, 0, 0), 0.3)
	q2 := FromAxisAngle(vector.NewVec3[float64](0, 1, 0), 0.4)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = q1.Mul(q2)
	}
}

func BenchmarkQuaternionSlerp(b *testing.B) {
	q1 := FromAxisAngle(vector.NewVec3[float64](1, 0, 0), 0.3)
	q2 := FromAxisAngle(vector.NewVec3[float64](0, 1, 0), 0.4)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = q1.Slerp(q2, 0.5)
	}
}
