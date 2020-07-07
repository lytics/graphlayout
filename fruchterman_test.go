package graphlayout

import (
	"context"
	"math"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
	"golang.org/x/sync/errgroup"
	"gonum.org/v1/gonum/mat"
)

func TestFruchtermanReingold(t *testing.T) {
	m := mat.NewDense(5, 5, nil)
	stop := make(chan struct{})

	// test empty graph
	{
		_, err := FruchtermanReingoldLayout(m, stop, nil)
		require.Equal(t, ErrNoEdges, err)
	}

	// test asymmetric graph
	{
		_, err := FruchtermanReingoldLayout(mat.NewDense(5, 6, nil), stop, nil)
		require.Equal(t, ErrMatrixAymmetric, err)
	}

	// test FR layout
	{
		//  Graph representation of `m`:
		//
		//  B --- D
		//  | \ / |
		//  |  A  |
		//  | / \ |
		//  C --- E
		//

		m = mat.NewDense(5, 5, []float64{0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0})

		layout, err := FruchtermanReingoldLayout(m, stop, nil)
		require.Equal(t, nil, err)

		xy := map[string]struct {
			x float64
			y float64
		}{
			"a": {layout.At(0, 0), layout.At(0, 1)},
			"b": {layout.At(1, 0), layout.At(1, 1)},
			"c": {layout.At(2, 0), layout.At(2, 1)},
			"d": {layout.At(3, 0), layout.At(3, 1)},
			"e": {layout.At(4, 0), layout.At(4, 1)},
		}
		distance := func(a, b string) float64 {
			return euclidean(xy[a].x-xy[b].x, xy[a].y-xy[b].y)
		}
		approxequal := func(x, y, margin float64) bool {
			if x == y {
				return true
			}
			return math.Abs(1-x/y) < margin
		}

		// because of the nature of the simulation, not *every* condition
		// below will be true *every* time, so allow some conditions to fail
		// sometimes
		var passing int
		if approxequal(distance("a", "b")+distance("a", "e"), distance("b", "e"), 0.6) {
			passing++
		}
		if approxequal(distance("c", "a")+distance("a", "d"), distance("c", "d"), 0.6) {
			passing++
		}
		if distance("b", "a") < distance("b", "e")*1.2 {
			passing++
		}
		if distance("c", "a") < distance("c", "d")*1.2 {
			passing++
		}

		assert.Truef(t, passing >= 2, "at least two distance checks must pass")
	}

	// test stopping
	{
		// run the layout in the background
		group, _ := errgroup.WithContext(context.Background())
		group.Go(func() error {
			// try performing a trillion iterations
			_, err := FruchtermanReingoldLayout(m, stop, &FruchtermanReingoldConfig{
				Niter: 1e9,
			})
			return err
		})
		// let run for 10ms, then stop it
		time.Sleep(time.Millisecond * 10)
		stop <- struct{}{}
		err := group.Wait()
		assert.Equal(t, ErrSimulationStopped, err)
	}
}
