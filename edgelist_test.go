package graphlayout

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"gonum.org/v1/gonum/mat"
)

func TestEdgelist(t *testing.T) {
	rawedges := []float64{0., 1., 1., 1., 0., 0., 0., 1., 1., 0., 0., 1., 0., 0., 1., 1., 1., 0., 0., 0., 1., 0., 1., 1., 0.}
	symedges := []float64{1, 0, 2, 0, 2, 1, 3, 0, 3, 1, 4, 0, 4, 2, 4, 3, 0, 1, 0, 2, 1, 2, 0, 3, 1, 3, 0, 4, 2, 4, 3, 4, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 1, 1, 0.5, 0.5}
	m := mat.NewDense(5, 5, rawedges)
	edges, err := matrixToEdgelist(m)
	assert.Equal(t, err, nil)
	assert.Equal(t, symedges, edges.vals)
}
