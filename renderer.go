package graphlayout

import (
  "errors"
  "gonum.org/v1/gonum/mat"
)

type GraphRenderer func (*mat.Dense, interface{}, <-chan struct{}) (*mat.Dense, error) // adjacency to xy coords

var (
  ErrSimulationStopped = errors.New("layout simulation stopped")
)
