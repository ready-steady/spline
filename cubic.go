package spline

import (
	"github.com/ready-steady/linear/system"
)

// Cubic is a cubic-spline interpolant.
type Cubic struct {
	nodes   []float64
	weights []float64
}

// NewCubic constructs a cubic-spline interpolant for a function y = f(x) given
// as a series of points (x, y). The x coordinates should be a strictly
// increasing sequence with at least two elements. The corresponding y
// coordinates can be multidimensional.
func NewCubic(x, y []float64) *Cubic {
	nn := len(x)
	nd := len(y) / nn

	dx := make([]float64, nn-1)
	dydx := make([]float64, (nn-1)*nd)
	for i := 0; i < (nn - 1); i++ {
		dx[i] = x[i+1] - x[i]
		for j := 0; j < nd; j++ {
			dydx[i*nd+j] = (y[(i+1)*nd+j] - y[i*nd+j]) / dx[i]
		}
	}

	s := &Cubic{}

	switch nn {
	case 2:
		s.nodes = []float64{x[0], x[1]}
		s.weights = make([]float64, nd*4)
		for j := 0; j < nd; j++ {
			s.weights[j*4+2] = dydx[j]
			s.weights[j*4+3] = y[j]
		}
	case 3:
		s.nodes = []float64{x[0], x[2]}
		s.weights = make([]float64, nd*4)
		for j := 0; j < nd; j++ {
			c1 := (dydx[nd+j] - dydx[j]) / (x[2] - x[0])
			s.weights[j*4+1] = c1
			s.weights[j*4+2] = dydx[j] - c1*dx[0]
			s.weights[j*4+3] = y[j]
		}
	default:
		xb := x[2] - x[0]
		xe := x[nn-1] - x[nn-3]

		a := make([]float64, nn)
		for i := 0; i < (nn - 2); i++ {
			a[i] = dx[i+1]
		}
		a[nn-2] = xe

		b := make([]float64, nn)
		b[0] = dx[1]
		for i := 1; i < (nn - 1); i++ {
			b[i] = 2 * (dx[i] + dx[i-1])
		}
		b[nn-1] = dx[nn-3]

		c := make([]float64, nn)
		c[1] = xb
		for i := 2; i < nn; i++ {
			c[i] = dx[i-2]
		}

		d := make([]float64, nd*nn)
		for j := 0; j < nd; j++ {
			d[j*nn] = ((dx[0]+2*xb)*dx[1]*dydx[j] + dx[0]*dx[0]*dydx[nd+j]) / xb
			for i := 1; i < (nn - 1); i++ {
				d[j*nn+i] = 3 * (dx[i]*dydx[(i-1)*nd+j] + dx[i-1]*dydx[i*nd+j])
			}
			d[j*nn+nn-1] = (dx[nn-2]*dx[nn-2]*dydx[(nn-3)*nd+j] +
				(2*xe+dx[nn-2])*dx[nn-3]*dydx[(nn-2)*nd+j]) / xe
		}

		slopes := system.ComputeTridiagonal(a, b, c, d)

		s.nodes = make([]float64, nn)
		copy(s.nodes, x)

		s.weights = make([]float64, (nn-1)*nd*4)
		for i, k := 0, 0; i < (nn - 1); i++ {
			for j := 0; j < nd; j++ {
				α := (dydx[i*nd+j] - slopes[j*nn+i]) / dx[i]
				β := (slopes[j*nn+i+1] - dydx[i*nd+j]) / dx[i]
				s.weights[k] = (β - α) / dx[i]
				s.weights[k+1] = 2*α - β
				s.weights[k+2] = slopes[j*nn+i]
				s.weights[k+3] = y[i*nd+j]
				k += 4
			}
		}
	}

	return s
}

// Evaluate interpolates the function values y = f(x). The x coordinates should
// be an increasing sequence whose elements belong to the range of the points
// that the interpolant has been constructed with.
func (s *Cubic) Evaluate(x []float64) []float64 {
	nn, np := len(s.nodes), len(x)
	nd := len(s.weights) / (4 * (nn - 1))

	y := make([]float64, np*nd)

	for i, l := 0, 0; i < np; i++ {
		for x[i] > s.nodes[l+1] && l < (nn-2) {
			l++
		}

		z := x[i] - s.nodes[l]
		for j, k := 0, l*nd*4; j < nd; j++ {
			y[i*nd+j] = z*(z*(z*s.weights[k]+s.weights[k+1])+s.weights[k+2]) + s.weights[k+3]
			k += 4
		}
	}

	return y
}
