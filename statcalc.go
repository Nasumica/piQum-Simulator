package main

// Author: Srbislav D. Nešić, srbislav.nesic@fincore.com

import "math"

// # Statistical calculator
type StatCalc struct {
	Cnt int     `json:"cnt"`           // category sample length
	Sum float64 `json:"sum"`           // sum of items
	Sqr float64 `json:"sqr,omitempty"` // sum of squares
	Min float64 `json:"min"`           // min (low)
	Max float64 `json:"max"`           // max (high)
	Avg float64 `json:"avg"`           // average (mean)
	Dev float64 `json:"dev,omitempty"` // standard deviation
	Nul int     `json:"nul,omitempty"` // zero items count
	Val float64 `json:"val"`           // last value
	Cat string  `json:"cat,omitempty"` // category name
}

// Reset statistical calculator.
func (s *StatCalc) Reset() {
	s.Cnt, s.Sum, s.Sqr, s.Min, s.Max = 0, 0, 0, 0, 0
	s.Avg, s.Dev, s.Nul, s.Val = 0, 0, 0, 0
}

// Add values to statistical calculator.
func (sc *StatCalc) Add(values ...float64) {
	for _, x := range values {
		if sc.Cnt == 0 {
			sc.Min = x
			sc.Max = x
		} else if sc.Min > x {
			sc.Min = x
		} else if sc.Max < x {
			sc.Max = x
		}
		sc.Cnt++
		if x == 0 {
			sc.Nul++
		} else {
			sc.Sum += x
			sc.Sqr += x * x
		}
		if sc.Min == sc.Max {
			sc.Avg = sc.Min
		} else {
			n := float64(sc.Cnt)
			sc.Avg = sc.Sum / n
			sc.Dev = math.Sqrt(math.Abs(n*sc.Sqr-sc.Sum*sc.Sum)) / n
		}
		sc.Val = x
	}
}

// Add integers to statistical calculator.
func (sc *StatCalc) Int(values ...int) {
	for _, i := range values {
		sc.Add(float64(i))
	}
}
