// vector/vector.go
package vector

// Number is a constraint for numeric types
type Number interface {
	~float32 | ~float64
}

// EPSILON is the default epsilon value for floating-point comparisons
const EPSILON = 1e-6
