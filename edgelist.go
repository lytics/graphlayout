package graphlayout

import (
	"errors"

	"gonum.org/v1/gonum/mat"
)

var (
	ErrMatrixAymmetric = errors.New("matrix not symmetric")
	ErrMatrixEmpty     = errors.New("matrix not symmetric")
	ErrNoEdges         = errors.New("graph contains no edges")
)

type edgelist struct {
	n    int
	m    int
	vals []float64
}

func isWeak(x, y float64) bool {
	return x > 0. || y > 0.
}

func isStrong(x, y float64) bool {
	return x > 0. && y > 0.
}

// vectorized + symmetrized
func matrixToEdgelist(mat *mat.Dense) (*edgelist, error) {
	r, c := mat.Dims()
	if r == 0 || c == 0 {
		return nil, ErrMatrixEmpty
	}
	if r != c {
		return nil, ErrMatrixAymmetric
	}

	var fi, fj float64

	// TODO: optimize capacity??
	cap := int(float64(r) / 2.)
	src := make([]float64, 0, cap)
	tgt := make([]float64, 0, cap)
	val := make([]float64, 0, cap)
	for i := 0; i < r; i++ {
		for j := 0; j < i; j++ {
			// since we're symmetrizing, look at either (i,j) or (j,i)
			val1 := mat.At(i, j)
			val2 := mat.At(j, i)
			// no edge, skipping
			if !isWeak(val1, val2) {
				continue
			}

			// take the midpoint of the edge weights
			v := (val1 + val2) / 2.
			fi = float64(i)
			fj = float64(j)

			// append (i, j)
			src = append(src, fi)
			tgt = append(tgt, fj)
			val = append(val, v)

			// append (j, i)
			src = append(src, fj)
			tgt = append(tgt, fi)
			val = append(val, v)
		}
	}

	if len(src) == 0 {
		return nil, ErrNoEdges
	}

	return &edgelist{
		n:    r,
		m:    len(src),
		vals: append(append(src, tgt...), val...),
	}, nil
}
