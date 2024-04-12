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
	s.Avg, s.Dev, s.Nul, s.Val, s.Cat = 0, 0, 0, 0, ""
}

// Add values to statistical calculator.
func (s *StatCalc) Add(values ...float64) {
	for _, x := range values {
		if s.Cnt == 0 {
			s.Min = x
			s.Max = x
		} else if s.Min > x {
			s.Min = x
		} else if s.Max < x {
			s.Max = x
		}
		s.Cnt++
		if x == 0 {
			s.Nul++
		} else {
			s.Sum += x
			s.Sqr += x * x
		}
		if s.Min == s.Max {
			s.Avg = s.Min
		} else {
			n := float64(s.Cnt)
			s.Avg = s.Sum / n
			s.Dev = math.Sqrt(math.Abs(n*s.Sqr-s.Sum*s.Sum)) / n
		}
		s.Val = x
	}
}
