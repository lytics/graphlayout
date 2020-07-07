package graphlayout

import (
	"math"
	"math/rand"
	"time"

	"gonum.org/v1/gonum/mat"
)

type FruchtermanReingoldConfig struct {
	// number of simulation iterations to run
	Niter int

	// maximum change in position for any iteration
	MaxDelta float64

	// size of space to explore
	Area float64

	// cooling exponent for annealing
	CoolExp float64

	// radius at which vertex/vertex repulsion cancels out attraction
	// of adjacent vertices
	RepulseRad float64

	// to speed calculations on large graphs, the plot region is divided at each
	// iteration into `ncell`x`ncell` "cells", which are used to define
	// neighborhoods and force calculation. moderate numbers of cells result in
	// fastest calculations. too few can look odd, and too many can take too long
	Ncells float64

	// cell jitter: factor (in units of cell width) used in assigning vertices
	// to cells.  small values may generate "grid-like" anomalies for graphs
	// with many isolates
	CellJitter float64

	// cell pointpointrad: square "radius" (in units of cells) such that exact
	// point interaction calculations are used for all vertices belonging to any
	// two cells less than or equal to this distance apart.  higher values
	// approximate the true F/R solution, but increases computational cost.
	CPPR float64

	// cell pointcellrad: squared "radius" (in units of cells) such that approximate
	// point/cell interaction calculations are used for all vertices belonging to
	// two cells less than or equal to this distance apart (and not within the
	// point/point radius). higher values provide somewhat better approximations
	// to the true F/R solution at slightly increased computational cost.
	CPCR float64

	// cell cellcellrad: squared "radius" (in units of cells) such that the approximate
	// cell/cell interactino calculations are used for all vertices belonging to
	// any two cells less than or equal to this distance apart (and not within the
	// point/point or point/cell radii). higher values provide somewhat better
	// approximations to the true F/R solution at slightly increased computational
	// cost.  note that cells beyond this radius (if any) do not interact, save
	// through edge attraction.
	CCCR float64
}

func FruchtermanReingoldLayout(matrix *mat.Dense, stop <-chan struct{}, conf *FruchtermanReingoldConfig) (*mat.Dense, error) {
	// list for the stop chan
	shouldStop := func() bool {
		select {
		case <-stop:
			return true
		default:
			return false
		}
	}

	d, err := matrixToEdgelist(matrix)
	if err != nil {
		return nil, err
	}

	n := float64(d.n)
	nint := int(n)

	rnd := rand.New(rand.NewSource(time.Now().UnixNano()))
	rnorm := func(stddev float64) float64 {
		return rnd.NormFloat64() * stddev
	}

	// default simulation parameters
	if conf == nil {
		conf = &FruchtermanReingoldConfig{}
	}

	if conf.Niter == 0 {
		conf.Niter = 500
	}

	if conf.MaxDelta == 0. {
		conf.MaxDelta = float64(n)
	}

	if conf.Area == 0. {
		conf.Area = n * n
	}

	if conf.CoolExp == 0. {
		conf.CoolExp = 3.
	}

	if conf.RepulseRad == 0. {
		conf.RepulseRad = conf.Area * math.Log(n)
	}

	if conf.Ncells == 0. {
		conf.Ncells = math.Ceil(math.Pow(n, 0.4))
	}

	if conf.CellJitter == 0. {
		conf.CellJitter = 0.5
	}

	// conf.CPPR defaults to zero

	if conf.CPCR == 0. {
		conf.CPCR = 18.
	}

	if conf.CCCR == 0. {
		conf.CCCR = conf.Ncells * conf.Ncells
	}

	// temperature
	var t float64

	// dyadic euclidean distance
	var ded float64

	// current x/y distances
	var xd, yd float64

	// repuslive/attractive forces
	var rf, af float64

	// current x/y widths
	var xwid, ywid float64

	// distance between cells
	var celldis float64

	// iterators
	var i, j, k, l int

	// indexes
	var ix, iy, jx, jy float64

	// cells
	var vcells *vcell
	var p *vcell
	var p2 *vcell
	var vlp *vlist
	var vlp2 *vlist

	frk2 := conf.Area / float64(n)
	frk := math.Sqrt(frk2)

	// calculate repulsive force based on dyadic euclidean distance
	repulse := func(ded float64) float64 {
		return frk2 * (1/ded - ded*ded/conf.RepulseRad)
	}

	// calculate attractive force based on dyadic euclidean distance
	attract := func(weight, ded float64) float64 {
		return weight * ded * ded / frk
	}

	// cool the annealing temperature
	cool := func(i int) float64 {
		return conf.MaxDelta * math.Pow(float64(conf.Niter-i)/float64(conf.Niter), conf.CoolExp)
	}

	// bounds
	var xmin, ymin, xmax, ymax float64

	x := make([]float64, nint)
	y := make([]float64, nint)
	dx := make([]float64, nint)
	dy := make([]float64, nint)
	cellid := make([]int, nint)

	small := 1e-6

	// seed some xy coords
	twopi := 2 * math.Pi
	for j = 0; j < nint; j++ {
		tempa := rnd.Float64()
		x[j] = n / twopi * math.Sin(twopi*tempa)
		y[j] = n / twopi * math.Cos(twopi*tempa)
	}

	for i = 0; i < conf.Niter; i++ {
		if shouldStop() {
			return nil, ErrSimulationStopped
		}

		xmin = math.Inf(1)
		ymin = math.Inf(1)
		xmax = math.Inf(-1)
		ymax = math.Inf(-1)

		// get extrama to form cells
		for j = 0; j < nint; j++ {
			xmin = math.Min(xmin, x[j])
			ymin = math.Min(ymin, y[j])
			xmax = math.Max(xmax, x[j])
			ymax = math.Max(ymax, y[j])
		}
		xmin -= 0.0001 * (xmax - xmin)
		ymin -= 0.0001 * (ymax - ymin)
		xmax += 0.0001 * (xmax - xmin)
		ymax += 0.0001 * (ymax - ymin)
		xwid = (xmax - xmin) / conf.Ncells
		ywid = (ymax - ymin) / conf.Ncells

		vcells = &vcell{}

		// calculate attractive forces
		for j = 0; j < nint; j++ {
			if shouldStop() {
				return nil, ErrSimulationStopped
			}
			xjitter := small * xwid
			yjitter := small * ywid
			jx = math.Max(math.Min(x[j]+rnorm(xwid*conf.CellJitter), xmax-xjitter), xmin+xjitter)
			jy = math.Max(math.Min(y[j]+rnorm(ywid*conf.CellJitter), ymax-yjitter), ymin+xjitter)
			cellid[j] = int(math.Floor((jx-xmin)/xwid) + conf.Ncells*math.Floor((jy-ymin)/ywid))

			// find j's cell
			for p = vcells; p != nil && p.next != nil && p.id != cellid[j]; p = p.next {
			}

			if p == nil {
				p = &vcell{
					id: cellid[j],
				}
				vcells = p
			} else if p.id != cellid[j] {
				// got to the end, insert new element
				p.next = &vcell{}
				p = p.next
				p.id = cellid[j]
			}

			// add j to the membership stack for this cell
			p.count++
			vlp = &vlist{
				v:    j,
				next: p.memb,
			}
			p.memb = vlp
			p.xm = (p.xm*(p.count-1) + x[j]) / p.count
			p.ym = (p.ym*(p.count-1) + y[j]) / p.count
		}

		// limit the maximum displacement to the temperature t
		t = cool(i)

		// clear the deltas
		for j = 0; j < nint; j++ {
			dx[j] = 0.
			dy[j] = 0.
		}

		// increement deltas for general force effects, using cells
		for p = vcells; p != nil; p = p.next {
			if shouldStop() {
				return nil, ErrSimulationStopped
			}
			for p2 = p; p2 != nil; p2 = p2.next {
				ix = float64(p.id % int(conf.Ncells))
				jx = float64(p2.id % int(conf.Ncells))
				iy = math.Floor(float64(p.id) / conf.Ncells)
				jy = math.Floor(float64(p2.id) / conf.Ncells)
				celldis = float64((ix-jx)*(ix-jx) + (iy-jy)*(iy-jy))
				if celldis <= conf.CPPR+0.001 {
					// use point/point calculations
					for vlp = p.memb; vlp != nil; vlp = vlp.next {
						if p == p2 {
							vlp2 = vlp.next
						} else {
							vlp2 = p2.memb
						}
						for ; vlp2 != nil; vlp2 = vlp2.next {
							// create L2 differences
							xd = x[vlp.v] - x[vlp2.v]
							yd = y[vlp.v] - y[vlp2.v]
							ded = euclidean(xd, yd)
							if ded == 0. {
								ded = 1
							}
							xd /= ded
							yd /= ded
							rf = repulse(ded)
							dx[vlp.v] += xd * rf
							dx[vlp2.v] -= xd * rf
							dy[vlp.v] += yd * rf
							dy[vlp2.v] -= yd * rf
						}
					}
				} else if celldis <= conf.CPCR+0.001 {
					// use point/cell calculations
					for vlp = p.memb; vlp != nil; vlp = vlp.next {
						// L2 norm!!
						xd = x[vlp.v] - p2.xm
						yd = y[vlp.v] - p2.ym
						ded = euclidean(xd, yd)
						if ded == 0. {
							ded = 1
						}
						xd /= ded
						yd /= ded
						rf = repulse(ded)
						dx[vlp.v] += xd * rf * p2.count
						dy[vlp.v] += yd * rf * p2.count
					}
				} else if celldis <= conf.CCCR+0.001 {
					// use cell/cell calculations
					xd = p.xm - p2.xm
					yd = p.ym - p2.ym
					ded = euclidean(xd, yd)
					if ded == 0. {
						ded = 1
					}
					xd /= ded
					yd /= ded
					rf = repulse(ded)
					// increment force to each member of p and p2
					for vlp = p.memb; vlp != nil; vlp = vlp.next {
						dx[vlp.v] += xd * rf * p2.count
						dy[vlp.v] += yd * rf * p2.count
					}
					for vlp = p2.memb; vlp != nil; vlp = vlp.next {
						dx[vlp.v] -= xd * rf * p.count
						dy[vlp.v] -= yd * rf * p.count
					}
				}
			}
		}

		// calculate attraction along edges
		for j = 0; j < d.m; j++ {
			k = int(d.vals[j])
			l = int(d.vals[j+d.m])
			xd = x[k] - x[l]
			yd = y[k] - y[l]
			ded = euclidean(xd, yd)
			af = attract(d.vals[j+2*d.m], ded)
			dx[k] -= xd * af
			dx[l] += xd * af
			dy[k] -= yd * af
			dy[l] += yd * af
		}

		// dampen motion, if needed, and move the points
		for j = 0; j < nint; j++ {
			ded = euclidean(dx[j], dy[j])
			if ded > t {
				ded = t / ded
				dx[j] *= ded
				dy[j] *= ded
			}
			x[j] += dx[j]
			y[j] += dy[j]
		}
	}

	output := mat.NewDense(nint, 2, nil)
	for i = 0; i < nint; i++ {
		output.Set(i, 0, x[i])
		output.Set(i, 1, y[i])
	}
	return output, nil

}

func euclidean(x, y float64) float64 {
	return math.Sqrt(x*x + y*y)
}

type vlist struct {
	v    int
	next *vlist
}

type vcell struct {
	id    int
	count float64
	xm    float64
	ym    float64
	memb  *vlist
	next  *vcell
}
