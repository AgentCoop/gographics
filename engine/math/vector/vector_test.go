// vector/vector_test.go
package vector

import (
	"encoding/json"
	"math"
	"testing"
)

// Helper functions for approximate equality
func approxEqual[T Number](a, b T) bool {
	return math.Abs(float64(a-b)) < 1e-6
}

func approxEqualVec2[T Number](a, b Vec2[T]) bool {
	return approxEqual(a.X, b.X) && approxEqual(a.Y, b.Y)
}

func approxEqualVec3[T Number](a, b Vec3[T]) bool {
	return approxEqual(a.X, b.X) && approxEqual(a.Y, b.Y) && approxEqual(a.Z, b.Z)
}

func approxEqualVec4[T Number](a, b Vec4[T]) bool {
	return approxEqual(a.X, b.X) && approxEqual(a.Y, b.Y) && approxEqual(a.Z, b.Z) && approxEqual(a.W, b.W)
}

// ============= Vector2 Tests =============

func TestVec2Creation(t *testing.T) {
	// Test float32
	v32 := NewVec2[float32](1.5, 2.5)
	if v32.X != 1.5 || v32.Y != 2.5 {
		t.Errorf("NewVec2[float32] failed: got (%v, %v), want (1.5, 2.5)", v32.X, v32.Y)
	}

	// Test float64
	v64 := NewVec2[float64](1.5, 2.5)
	if v64.X != 1.5 || v64.Y != 2.5 {
		t.Errorf("NewVec2[float64] failed: got (%v, %v), want (1.5, 2.5)", v64.X, v64.Y)
	}

	// Test scalar creation
	vScalar32 := NewVec2Scalar[float32](3.0)
	if vScalar32.X != 3.0 || vScalar32.Y != 3.0 {
		t.Errorf("NewVec2Scalar[float32] failed: got (%v, %v), want (3.0, 3.0)", vScalar32.X, vScalar32.Y)
	}

	// Test zero and one
	zero32 := Zero2[float32]()
	if zero32.X != 0 || zero32.Y != 0 {
		t.Errorf("Zero2[float32] failed")
	}

	one32 := One2[float32]()
	if one32.X != 1 || one32.Y != 1 {
		t.Errorf("One2[float32] failed")
	}

	// Test unit vectors
	unitX := UnitX2[float32]()
	if unitX.X != 1 || unitX.Y != 0 {
		t.Errorf("UnitX2 failed")
	}

	unitY := UnitY2[float32]()
	if unitY.X != 0 || unitY.Y != 1 {
		t.Errorf("UnitY2 failed")
	}
}

func TestVec2Operations(t *testing.T) {
	// Test with float32
	v1 := NewVec2[float32](1, 2)
	v2 := NewVec2[float32](3, 4)

	// Addition
	sum := v1.Add(v2)
	expected := NewVec2[float32](4, 6)
	if !approxEqualVec2(sum, expected) {
		t.Errorf("Add failed: got %v, want %v", sum, expected)
	}

	// Subtraction
	diff := v1.Sub(v2)
	expected = NewVec2[float32](-2, -2)
	if !approxEqualVec2(diff, expected) {
		t.Errorf("Sub failed: got %v, want %v", diff, expected)
	}

	// Multiplication (component-wise)
	mul := v1.Mul(v2)
	expected = NewVec2[float32](3, 8)
	if !approxEqualVec2(mul, expected) {
		t.Errorf("Mul failed: got %v, want %v", mul, expected)
	}

	// Division (component-wise)
	v3 := NewVec2[float32](6, 8)
	v4 := NewVec2[float32](2, 4)
	div := v3.Div(v4)
	expected = NewVec2[float32](3, 2)
	if !approxEqualVec2(div, expected) {
		t.Errorf("Div failed: got %v, want %v", div, expected)
	}

	// Scalar operations
	v5 := NewVec2[float32](2, 3)

	// Add scalar
	addScalar := v5.AddScalar(1)
	expected = NewVec2[float32](3, 4)
	if !approxEqualVec2(addScalar, expected) {
		t.Errorf("AddScalar failed: got %v, want %v", addScalar, expected)
	}

	// Multiply scalar
	mulScalar := v5.MulScalar(2)
	expected = NewVec2[float32](4, 6)
	if !approxEqualVec2(mulScalar, expected) {
		t.Errorf("MulScalar failed: got %v, want %v", mulScalar, expected)
	}
}

func TestVec2DotProduct(t *testing.T) {
	v1 := NewVec2[float32](1, 2)
	v2 := NewVec2[float32](3, 4)

	dot := v1.Dot(v2)
	expected := float32(1*3 + 2*4) // 11
	if !approxEqual(dot, expected) {
		t.Errorf("Dot failed: got %v, want %v", dot, expected)
	}

	// Test with float64 for precision
	v1d := NewVec2[float64](1.1, 2.2)
	v2d := NewVec2[float64](3.3, 4.4)
	dot64 := v1d.Dot(v2d)
	expected64 := 1.1*3.3 + 2.2*4.4
	if !approxEqual(dot64, expected64) {
		t.Errorf("Dot[float64] failed: got %v, want %v", dot64, expected64)
	}
}

func TestVec2CrossProduct(t *testing.T) {
	v1 := NewVec2[float32](1, 2)
	v2 := NewVec2[float32](3, 4)

	cross := v1.Cross(v2)
	expected := float32(1*4 - 2*3) // -2
	if !approxEqual(cross, expected) {
		t.Errorf("Cross failed: got %v, want %v", cross, expected)
	}
}

func TestVec2Length(t *testing.T) {
	v := NewVec2[float32](3, 4)

	length := v.Length()
	expected := float32(5.0) // sqrt(3² + 4²) = 5
	if !approxEqual(length, expected) {
		t.Errorf("Length failed: got %v, want %v", length, expected)
	}

	lengthSq := v.LengthSquared()
	expectedSq := float32(25.0)
	if !approxEqual(lengthSq, expectedSq) {
		t.Errorf("LengthSquared failed: got %v, want %v", lengthSq, expectedSq)
	}
}

func TestVec2Normalize(t *testing.T) {
	v := NewVec2[float32](3, 4)
	normalized := v.Normalize()

	// Check length is 1
	length := normalized.Length()
	if !approxEqual(length, 1.0) {
		t.Errorf("Normalize length failed: got %v, want 1", length)
	}

	// Check direction is preserved
	expected := NewVec2[float32](0.6, 0.8) // (3/5, 4/5)
	if !approxEqualVec2(normalized, expected) {
		t.Errorf("Normalize direction failed: got %v, want %v", normalized, expected)
	}

	// Test zero vector
	zero := Zero2[float32]()
	normalizedZero := zero.Normalize()
	if !normalizedZero.IsZero() {
		t.Errorf("Normalize zero vector failed: got %v, want zero", normalizedZero)
	}
}

func TestVec2Distance(t *testing.T) {
	v1 := NewVec2[float32](1, 2)
	v2 := NewVec2[float32](4, 6)

	distance := v1.Distance(v2)
	expected := float32(5.0) // sqrt((4-1)² + (6-2)²) = sqrt(9 + 16) = 5
	if !approxEqual(distance, expected) {
		t.Errorf("Distance failed: got %v, want %v", distance, expected)
	}

	distanceSq := v1.DistanceSquared(v2)
	expectedSq := float32(25.0)
	if !approxEqual(distanceSq, expectedSq) {
		t.Errorf("DistanceSquared failed: got %v, want %v", distanceSq, expectedSq)
	}
}

func TestVec2Lerp(t *testing.T) {
	v1 := NewVec2[float32](0, 0)
	v2 := NewVec2[float32](10, 10)

	// t = 0.5
	result := v1.Lerp(v2, 0.5)
	expected := NewVec2[float32](5, 5)
	if !approxEqualVec2(result, expected) {
		t.Errorf("Lerp(0.5) failed: got %v, want %v", result, expected)
	}

	// t = 0 (should return v1)
	result = v1.Lerp(v2, 0)
	if !approxEqualVec2(result, v1) {
		t.Errorf("Lerp(0) failed: got %v, want %v", result, v1)
	}

	// t = 1 (should return v2)
	result = v1.Lerp(v2, 1)
	if !approxEqualVec2(result, v2) {
		t.Errorf("Lerp(1) failed: got %v, want %v", result, v2)
	}
}

func TestVec2Rotate(t *testing.T) {
	v := NewVec2[float32](1, 0)

	// Rotate 90 degrees (π/2 radians)
	rotated := v.Rotate(float32(math.Pi / 2))
	expected := NewVec2[float32](0, 1)
	if !approxEqualVec2(rotated, expected) {
		t.Errorf("Rotate(90°) failed: got %v, want %v", rotated, expected)
	}

	// Rotate 180 degrees
	rotated = v.Rotate(float32(math.Pi))
	expected = NewVec2[float32](-1, 0)
	if !approxEqualVec2(rotated, expected) {
		t.Errorf("Rotate(180°) failed: got %v, want %v", rotated, expected)
	}

	// Rotate 360 degrees (should return original)
	rotated = v.Rotate(float32(2 * math.Pi))
	if !approxEqualVec2(rotated, v) {
		t.Errorf("Rotate(360°) failed: got %v, want %v", rotated, v)
	}
}

func TestVec2Angle(t *testing.T) {
	// Test angle of vector (1, 0)
	v1 := NewVec2[float32](1, 0)
	angle := v1.Angle()
	if !approxEqual(angle, 0.0) {
		t.Errorf("Angle of (1,0) failed: got %v, want 0", angle)
	}

	// Test angle of vector (0, 1)
	v2 := NewVec2[float32](0, 1)
	angle = v2.Angle()
	expected := float32(math.Pi / 2)
	if !approxEqual(angle, expected) {
		t.Errorf("Angle of (0,1) failed: got %v, want %v", angle, expected)
	}

	// Test angle between vectors
	v3 := NewVec2[float32](1, 0)
	v4 := NewVec2[float32](0, 1)
	angleBetween := v3.AngleTo(v4)
	expected = float32(math.Pi / 2)
	if !approxEqual(angleBetween, expected) {
		t.Errorf("AngleTo failed: got %v, want %v", angleBetween, expected)
	}
}

func TestVec2Reflect(t *testing.T) {
	// Reflecting across the x-axis
	v := NewVec2[float32](1, -1)
	normal := NewVec2[float32](0, 1) // Pointing up
	reflected := v.Reflect(normal)
	expected := NewVec2[float32](1, 1)
	if !approxEqualVec2(reflected, expected) {
		t.Errorf("Reflect failed: got %v, want %v", reflected, expected)
	}
}

func TestVec2Project(t *testing.T) {
	v := NewVec2[float32](3, 4)
	onto := NewVec2[float32](1, 0) // x-axis

	projected := v.Project(onto)
	expected := NewVec2[float32](3, 0) // Projection onto x-axis
	if !approxEqualVec2(projected, expected) {
		t.Errorf("Project failed: got %v, want %v", projected, expected)
	}

	// Project onto zero vector
	zero := Zero2[float32]()
	projectedZero := v.Project(zero)
	if !projectedZero.IsZero() {
		t.Errorf("Project onto zero failed: got %v, want zero", projectedZero)
	}
}

func TestVec2Clamp(t *testing.T) {
	v := NewVec2[float32](5, -5)
	min := NewVec2[float32](-2, -2)
	max := NewVec2[float32](2, 2)

	clamped := v.Clamp(min, max)
	expected := NewVec2[float32](2, -2)
	if !approxEqualVec2(clamped, expected) {
		t.Errorf("Clamp failed: got %v, want %v", clamped, expected)
	}

	// Scalar clamp
	clampedScalar := v.ClampScalar(-1, 1)
	expected = NewVec2[float32](1, -1)
	if !approxEqualVec2(clampedScalar, expected) {
		t.Errorf("ClampScalar failed: got %v, want %v", clampedScalar, expected)
	}
}

func TestVec2JSON(t *testing.T) {
	v := NewVec2[float64](1.5, 2.5)

	// Marshal
	data, err := json.Marshal(v)
	if err != nil {
		t.Errorf("Marshal failed: %v", err)
	}

	expectedJSON := `[1.500000,2.500000]`
	if string(data) != expectedJSON {
		t.Errorf("Marshal failed: got %s, want %s", data, expectedJSON)
	}

	// Unmarshal
	var v2 Vec2[float64]
	err = json.Unmarshal(data, &v2)
	if err != nil {
		t.Errorf("Unmarshal failed: %v", err)
	}

	if !approxEqualVec2(v2, v) {
		t.Errorf("Unmarshal failed: got %v, want %v", v2, v)
	}
}

// ============= Vector3 Tests =============

func TestVec3Creation(t *testing.T) {
	v32 := NewVec3[float32](1, 2, 3)
	if v32.X != 1 || v32.Y != 2 || v32.Z != 3 {
		t.Errorf("NewVec3[float32] failed")
	}

	v64 := NewVec3[float64](1, 2, 3)
	if v64.X != 1 || v64.Y != 2 || v64.Z != 3 {
		t.Errorf("NewVec3[float64] failed")
	}

	// Test unit vectors
	unitX := UnitX3[float32]()
	if unitX.X != 1 || unitX.Y != 0 || unitX.Z != 0 {
		t.Errorf("UnitX3 failed")
	}

	// Test direction helpers
	up := Up[float32]()
	if up.X != 0 || up.Y != 1 || up.Z != 0 {
		t.Errorf("Up failed")
	}

	forward := Forward[float32]()
	if forward.X != 0 || forward.Y != 0 || forward.Z != -1 {
		t.Errorf("Forward failed")
	}
}

func TestVec3Operations(t *testing.T) {
	v1 := NewVec3[float32](1, 2, 3)
	v2 := NewVec3[float32](4, 5, 6)

	// Cross product
	cross := v1.Cross(v2)
	expected := NewVec3[float32](
		2*6-3*5, // -3
		3*4-1*6, // 6
		1*5-2*4, // -3
	)
	if !approxEqualVec3(cross, expected) {
		t.Errorf("Cross failed: got %v, want %v", cross, expected)
	}

	// Test cross product properties
	// v1 × v2 = -(v2 × v1)
	cross2 := v2.Cross(v1)
	if !approxEqualVec3(cross, cross2.MulScalar(-1)) {
		t.Errorf("Cross anti-commutative property failed")
	}

	// v1 × v1 = 0
	crossSelf := v1.Cross(v1)
	if !crossSelf.IsZero() {
		t.Errorf("Cross self should be zero: got %v", crossSelf)
	}
}

func TestVec3Refract(t *testing.T) {
	// Simple refraction test (air to water, normal incidence)
	incident := NewVec3[float32](0, 0, -1) // Coming straight down
	normal := NewVec3[float32](0, 0, 1)    // Surface normal pointing up
	eta := float32(1.0 / 1.33)             // air to water

	refracted := incident.Refract(normal, eta)

	// For normal incidence, refracted ray should continue straight
	expected := NewVec3[float32](0, 0, -1)
	if !approxEqualVec3(refracted, expected) {
		t.Errorf("Refract failed: got %v, want %v", refracted, expected)
	}
}

func TestVec3Angle(t *testing.T) {
	v1 := NewVec3[float32](1, 0, 0)
	v2 := NewVec3[float32](0, 1, 0)

	angle := v1.Angle(v2)
	expected := float32(math.Pi / 2)
	if !approxEqual(angle, expected) {
		t.Errorf("Vec3 Angle failed: got %v, want %v", angle, expected)
	}

	// Same vector should have 0 angle
	angle = v1.Angle(v1)
	if !approxEqual(angle, 0) {
		t.Errorf("Vec3 Angle same vector failed: got %v, want 0", angle)
	}
}

// ============= Vector4 Tests =============

func TestVec4Creation(t *testing.T) {
	// Test position creation
	pos := NewPosition[float32](1, 2, 3)
	if pos.X != 1 || pos.Y != 2 || pos.Z != 3 || pos.W != 1 {
		t.Errorf("NewPosition failed: got %v", pos)
	}

	// Test direction creation
	dir := NewDirection[float32](1, 2, 3)
	if dir.X != 1 || dir.Y != 2 || dir.Z != 3 || dir.W != 0 {
		t.Errorf("NewDirection failed: got %v", dir)
	}

	// Test from Vec2
	v2 := NewVec2[float32](1, 2)
	v4 := FromVec2To4(v2, 3, 4)
	if v4.X != 1 || v4.Y != 2 || v4.Z != 3 || v4.W != 4 {
		t.Errorf("FromVec2 failed: got %v", v4)
	}
}

func TestVec4Homogenize(t *testing.T) {
	// Test position homogenization
	pos := NewVec4[float32](2, 4, 6, 2) // Represents (1, 2, 3, 1)
	homogenized := pos.Homogenize()
	expected := NewVec4[float32](1, 2, 3, 1)
	if !approxEqualVec4(homogenized, expected) {
		t.Errorf("Homogenize failed: got %v, want %v", homogenized, expected)
	}

	// Test direction (should remain unchanged or become zero)
	dir := NewDirection[float32](1, 2, 3)
	homogenizedDir := dir.Homogenize()
	if !homogenizedDir.IsZero() {
		t.Errorf("Homogenize direction failed: got %v, want zero", homogenizedDir)
	}
}

func TestVec4PointDirection(t *testing.T) {
	// Test IsPoint
	pos := NewPosition[float32](1, 2, 3)
	if !pos.IsPoint() {
		t.Errorf("IsPoint failed for position")
	}

	dir := NewDirection[float32](1, 2, 3)
	if !dir.IsDirection() {
		t.Errorf("IsDirection failed for direction")
	}

	// Test conversion
	pos2 := dir.ToPosition()
	if !pos2.IsPoint() {
		t.Errorf("ToPosition failed")
	}

	dir2 := pos.ToDirection()
	if !dir2.IsDirection() {
		t.Errorf("ToDirection failed")
	}
}

func TestVec4Transform(t *testing.T) {
	// Identity matrix
	identity := [16]float32{
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1,
	}

	v := NewVec4[float32](1, 2, 3, 1)
	transformed := v.Transform(identity)

	if !approxEqualVec4(transformed, v) {
		t.Errorf("Transform with identity failed: got %v, want %v", transformed, v)
	}

	// Translation matrix
	translate := [16]float32{
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		5, 6, 7, 1,
	}

	transformed = v.Transform(translate)
	expected := NewVec4[float32](1+5, 2+6, 3+7, 1)
	if !approxEqualVec4(transformed, expected) {
		t.Errorf("Transform with translation failed: got %v, want %v", transformed, expected)
	}
}

func TestVec4AffineTransform(t *testing.T) {
	// Translation matrix (affine)
	translate := [16]float32{
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		5, 6, 7, 1,
	}

	// Test with position
	pos := NewPosition[float32](1, 2, 3)
	transformed := pos.TransformAffine(translate)
	expected := NewVec4[float32](1+5, 2+6, 3+7, 1)
	if !approxEqualVec4(transformed, expected) {
		t.Errorf("TransformAffine position failed: got %v, want %v", transformed, expected)
	}

	// Test with direction (W should remain 0)
	dir := NewDirection[float32](1, 2, 3)
	transformedDir := dir.TransformAffine(translate)
	expectedDir := NewVec4[float32](1, 2, 3, 0) // Directions don't translate
	if !approxEqualVec4(transformedDir, expectedDir) {
		t.Errorf("TransformAffine direction failed: got %v, want %v", transformedDir, expectedDir)
	}
}

// ============= Benchmark Tests =============

func BenchmarkVec2Add(b *testing.B) {
	v1 := NewVec2[float64](1, 2)
	v2 := NewVec2[float64](3, 4)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = v1.Add(v2)
	}
}

func BenchmarkVec2Dot(b *testing.B) {
	v1 := NewVec2[float64](1, 2)
	v2 := NewVec2[float64](3, 4)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = v1.Dot(v2)
	}
}

func BenchmarkVec2Length(b *testing.B) {
	v := NewVec2[float64](3, 4)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = v.Length()
	}
}

func BenchmarkVec3Cross(b *testing.B) {
	v1 := NewVec3[float64](1, 2, 3)
	v2 := NewVec3[float64](4, 5, 6)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = v1.Cross(v2)
	}
}

func BenchmarkVec4Transform(b *testing.B) {
	v := NewVec4[float64](1, 2, 3, 1)
	m := [16]float64{
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		5, 6, 7, 1,
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = v.Transform(m)
	}
}

// ============= Edge Case Tests =============

func TestVec2EdgeCases(t *testing.T) {
	// Test with very small values
	v := NewVec2[float64](1e-10, 2e-10)
	length := v.Length()
	if length <= 0 {
		t.Errorf("Length failed for very small values: got %v", length)
	}

	// Test with very large values
	v = NewVec2[float64](1e10, 2e10)
	length = v.Length()
	if math.IsInf(float64(length), 0) || math.IsNaN(float64(length)) {
		t.Errorf("Length failed for very large values: got %v", length)
	}

	// Test NaN and Inf
	v = NewVec2[float64](math.NaN(), 1)
	length = v.Length()
	if !math.IsNaN(float64(length)) {
		t.Errorf("Length with NaN should return NaN: got %v", length)
	}
}

func TestVec3EdgeCases(t *testing.T) {
	// Test zero vector normalization
	zero := Zero3[float32]()
	normalized := zero.Normalize()
	if !normalized.IsZero() {
		t.Errorf("Normalize zero vector failed: got %v", normalized)
	}

	// Test near-zero vector
	epsilonVec := NewVec3[float32](1e-8, 1e-8, 1e-8)
	normalized = epsilonVec.Normalize()
	// Should either be zero or normalized
	if !normalized.IsZero(1e-3) && math.Abs(float64(normalized.Length()-1)) > 1e-3 {
		t.Errorf("Normalize near-zero vector failed: length = %v", normalized.Length())
	}
}

func TestVec4EdgeCases(t *testing.T) {
	// Test homogenization with near-zero W
	v := NewVec4[float64](1, 2, 3, 1e-10)
	homogenized := v.Homogenize()
	if !homogenized.IsZero(1e-3) {
		t.Errorf("Homogenize with near-zero W failed: got %v", homogenized)
	}

	// Test negative W
	v = NewVec4[float64](2, 4, 6, -2)
	homogenized = v.Homogenize()
	expected := NewVec4[float64](-1, -2, -3, 1)
	if !approxEqualVec4(homogenized, expected) {
		t.Errorf("Homogenize with negative W failed: got %v, want %v", homogenized, expected)
	}
}

// ============= Property-Based Tests =============

func TestVec2Properties(t *testing.T) {
	// Test commutative property of addition
	v1 := NewVec2[float64](math.Pi, math.E)
	v2 := NewVec2[float64](math.Sqrt2, math.Sqrt(3))

	sum1 := v1.Add(v2)
	sum2 := v2.Add(v1)
	if !approxEqualVec2(sum1, sum2) {
		t.Errorf("Addition commutative property failed")
	}

	// Test associative property of addition
	v3 := NewVec2[float64](1.1, 2.2)
	sumA := v1.Add(v2).Add(v3)
	sumB := v1.Add(v2.Add(v3))
	if !approxEqualVec2(sumA, sumB) {
		t.Errorf("Addition associative property failed")
	}

	// Test distributive property
	scalar := 2.5
	dist1 := v1.MulScalar(scalar).Add(v2.MulScalar(scalar))
	dist2 := v1.Add(v2).MulScalar(scalar)
	if !approxEqualVec2(dist1, dist2) {
		t.Errorf("Distributive property failed")
	}

	// Test dot product symmetry
	dot1 := v1.Dot(v2)
	dot2 := v2.Dot(v1)
	if !approxEqual(dot1, dot2) {
		t.Errorf("Dot product symmetry failed")
	}
}

func TestVec3Properties(t *testing.T) {
	// Test cross product anti-commutative property
	v1 := NewVec3[float64](1, 2, 3)
	v2 := NewVec3[float64](4, 5, 6)

	cross1 := v1.Cross(v2)
	cross2 := v2.Cross(v1)

	// v1 × v2 = -(v2 × v1)
	if !approxEqualVec3(cross1, cross2.MulScalar(-1)) {
		t.Errorf("Cross product anti-commutative property failed")
	}

	// Test Jacobi identity: a × (b × c) + b × (c × a) + c × (a × b) = 0
	v3 := NewVec3[float64](7, 8, 9)

	term1 := v1.Cross(v2.Cross(v3))
	term2 := v2.Cross(v3.Cross(v1))
	term3 := v3.Cross(v1.Cross(v2))

	sum := term1.Add(term2).Add(term3)
	if !sum.IsZero(1e-10) {
		t.Errorf("Jacobi identity failed: got %v", sum)
	}
}

// ============= Type Conversion Tests =============

func TestTypeConversions(t *testing.T) {
	// Vec2 to Vec3
	v2 := NewVec2[float32](1, 2)
	v3 := v2.ToVec3()
	expected3 := NewVec3[float32](1, 2, 0)
	if !approxEqualVec3(v3, expected3) {
		t.Errorf("Vec2.ToVec3 failed: got %v, want %v", v3, expected3)
	}

	// Vec2 to Vec4
	v4 := v2.ToVec4()
	expected4 := NewVec4[float32](1, 2, 0, 1)
	if !approxEqualVec4(v4, expected4) {
		t.Errorf("Vec2.ToVec4 failed: got %v, want %v", v4, expected4)
	}

	// Vec3 to Vec2
	v3b := NewVec3[float32](1, 2, 3)
	v2b := v3b.ToVec2()
	expected2 := NewVec2[float32](1, 2)
	if !approxEqualVec2(v2b, expected2) {
		t.Errorf("Vec3.ToVec2 failed: got %v, want %v", v2b, expected2)
	}

	// Vec3 to Vec4
	v4b := v3b.ToVec4()
	expected4b := NewVec4[float32](1, 2, 3, 1)
	if !approxEqualVec4(v4b, expected4b) {
		t.Errorf("Vec3.ToVec4 failed: got %v, want %v", v4b, expected4b)
	}

	// Vec4 to Vec2
	v4c := NewVec4[float32](1, 2, 3, 1)
	v2c := v4c.ToVec2()
	expected2c := NewVec2[float32](1, 2)
	if !approxEqualVec2(v2c, expected2c) {
		t.Errorf("Vec4.ToVec2 failed: got %v, want %v", v2c, expected2c)
	}

	// Vec4 to Vec3 (position)
	v4d := NewPosition[float32](2, 4, 6)
	v3c := v4d.ToVec3()
	expected3c := NewVec3[float32](2, 4, 6)
	if !approxEqualVec3(v3c, expected3c) {
		t.Errorf("Vec4.ToVec3 position failed: got %v, want %v", v3c, expected3c)
	}

	// Vec4 to Vec3 (direction)
	v4e := NewDirection[float32](1, 2, 3)
	v3d := v4e.ToVec3()
	expected3d := NewVec3[float32](1, 2, 3)
	if !approxEqualVec3(v3d, expected3d) {
		t.Errorf("Vec4.ToVec3 direction failed: got %v, want %v", v3d, expected3d)
	}
}

// ============= Utility Function Tests =============

func TestVec2String(t *testing.T) {
	v := NewVec2[float64](1.234567, 2.345678)
	str := v.String()
	// Check it contains the values
	if len(str) == 0 {
		t.Errorf("String() returned empty string")
	}
}

func TestVec2CopyAndClone(t *testing.T) {
	original := NewVec2[float32](1, 2)

	// Test Clone
	cloned := original.Clone()
	if !approxEqualVec2(cloned, original) {
		t.Errorf("Clone failed: got %v, want %v", cloned, original)
	}

	// Test that they're independent
	cloned.X = 99
	if original.X == 99 {
		t.Errorf("Clone not independent: modifying clone affected original")
	}

	// Test Copy method
	var copyVec Vec2[float32]
	copyVec.Copy(original)
	if !approxEqualVec2(copyVec, original) {
		t.Errorf("Copy failed: got %v, want %v", copyVec, original)
	}
}

func TestVec2SetAndZero(t *testing.T) {
	var v Vec2[float32]

	// Test Set
	v.Set(5, 6)
	if v.X != 5 || v.Y != 6 {
		t.Errorf("Set failed: got (%v, %v)", v.X, v.Y)
	}

	// Test Zero
	v.Zero()
	if !v.IsZero() {
		t.Errorf("Zero failed: got (%v, %v)", v.X, v.Y)
	}
}
