package main

// Author: Srbislav D. Nešić, srbislav.nesic@fincore.com

import (
	"crypto/rand"
	"math"
	"math/big"
	"sync"
	"time"
)

// # Linear congruential pseudo-random numbers generator
type LCPRNG struct {
	seed uint64     // generator seed
	dog  sync.Mutex // watchdog Šarko
}

// Current seed.
func (rnd *LCPRNG) Seed() uint64 {
	rnd.dog.Lock()
	defer rnd.dog.Unlock()
	return rnd.seed
}

// Inititalize generator with seeds or system clock.
func (rnd *LCPRNG) Randomize(seeds ...uint64) (seed uint64) {
	xorshift := func(x uint64) uint64 { // by George Marsaglia
		x ^= x << 13
		x ^= x >> 7
		x ^= x << 17
		return x
	}
	if len(seeds) == 0 {
		seed = xorshift(uint64(time.Now().UnixNano()) / 100)
		if g, e := rand.Int(rand.Reader, new(big.Int).SetBit(new(big.Int), 64, 1)); e == nil {
			seed ^= g.Uint64() // seed from computer crypto entropy generator
		}
	} else {
		for _, s := range seeds {
			seed = xorshift(seed) ^ s
		}
	}
	rnd.dog.Lock()         // Meni je nekako logičnije
	defer rnd.dog.Unlock() // da bude obrnuto, ... :)
	rnd.seed = seed
	return
}

// Next random value from generator.
//
// (Cycle after 584554 years and 18 days for 1 million randoms per second.)
/*
	y = x * a + c (mod 2⁶⁴)
where constants
	a = 0x5851f42d4c957f2d
	c = 0x14057b7ef767814f
given by Knuth in MMIX RISC processor.
*/
func (rnd *LCPRNG) Next() uint64 {
	rnd.dog.Lock()         // ... da pustim kuče dok ja radim,
	defer rnd.dog.Unlock() // a da ga vežem kad rade drugi. :)
	const (
		a uint64 = 0x5851f42d4c957f2d // multiplier
		c uint64 = 0x14057b7ef767814f // incrementer
	)
	rnd.seed *= a // constants
	rnd.seed += c // by Knuth
	return rnd.seed
}

// Previous random value from generator.
/*
If
	y = x * a + c
then
	x = (y - c) / a
	  = y/a - c/a
	  = 1/a * y - 1/a * c
	  = b * y + d
where
	b = 1/a = MulInv64(a)
	d = -b * c
For given constants a and c
	b = 0xc097ef87329e28a5
	d = 0x9995b5b621535015
*/
func (rnd *LCPRNG) Prev() uint64 {
	rnd.dog.Lock()
	defer rnd.dog.Unlock()
	const (
		b uint64 = 0xc097ef87329e28a5 // multiplier
		d uint64 = 0x9995b5b621535015 // incrementer
	)
	rnd.seed *= b
	rnd.seed += d
	return rnd.seed
}

// Next random value from generator limited to range [0, n].
func (rnd *LCPRNG) Limited(n uint64) uint64 {
	if n != 0 {
		if n++; n == 0 { // n = 2⁶⁴
			n = rnd.Next()
		} else { // acceptance-rejection: p(reject) = 2⁶⁴ % n / 2⁶⁴
			for f, u := -n/n+1, n; n >= u; { // f = 2⁶⁴ / n
				n = rnd.Next() / f
			}
		}
	}
	return n
}

// Random integer in range [m, n].
func (rnd *LCPRNG) Int(m, n int) int {
	if m > n {
		m, n = n, m
	}
	return m + int(rnd.Limited(uint64(n-m)))
}

// Random integer in range [0, n) for positive n, -1 for n = 0, else in range [n, -1].
//
//	μ = (n - 1) / 2
//	σ² = (n² - 1) / 12
func (rnd *LCPRNG) Choice(n int) int {
	if n > 1 {
		return int(rnd.Limited(uint64(n - 1)))
	} else if n >= 0 {
		return n - 1
	} else {
		return -int(rnd.Limited(uint64(-n-1))) - 1
	}
}

// True with probability k/n.
func (rnd *LCPRNG) Choose(n, k int) bool {
	return (n > 0) && (k > 0) && (n <= k || rnd.Choice(n) < k)
}

// Knuth shuffle (Fisher-Yates).
func (rnd *LCPRNG) Shuffle(a *[]int) {
	for j, i := 0, len(*a); i > 1; (*a)[i], (*a)[j] = (*a)[j], (*a)[i] {
		j, i = rnd.Choice(i), i-1
	}
}

// Array of n integers in range [m, m + n) in random order.
func (rnd *LCPRNG) Fill(m, n int) (a []int) {
	if n > 0 {
		a = make([]int, n)
		for i := range a {
			j := rnd.Int(0, i)
			a[i], a[j] = a[j], m+i
		}
	}
	return
}

// Random permutation.
func (rnd *LCPRNG) Permutation(n int) []int {
	return rnd.Fill(0, n)
}

// Random combination k of n elements (quickpick).
func (rnd *LCPRNG) Combination(n, k int) (c []int) {
	for i := 0; n > 0 && k > 0; i, n = i+1, n-1 {
		if rnd.Choose(n, k) {
			c = append(c, i)
			k--
		}
	}
	return
}

// Random sample of k elements.
func (rnd *LCPRNG) Sample(k int, a *[]int) []int {
	s := rnd.Combination(len(*a), k)
	for i, j := range s {
		s[i] = (*a)[j]
	}
	return s
}

// Hypergeometric distribution random variable.
//
//	p(hits) = HypGeomDist(hits, draw, succ, size)
func (rnd *LCPRNG) HyperGeometric(draw, succ, size int) (hits int) {
	if size >= draw && size >= succ {
		for ; draw > 0 && succ > 0; draw-- {
			if size == succ {
				if draw < succ {
					hits += draw
				} else {
					hits += succ
				}
				break
			}
			if rnd.Choose(size, succ) {
				hits++
				succ--
			}
			size--
		}
	}
	return
}

// Random array index for non-empty array else -1.
func (rnd *LCPRNG) Index(a *[]int) int {
	return rnd.Choice(len(*a))
}

// Random item of array.
func (rnd *LCPRNG) Item(a *[]int) int {
	if i := rnd.Index(a); i < 0 {
		return i
	} else {
		return (*a)[i]
	}
}

// Random value from non-empty array else NaN.
func (rnd *LCPRNG) Value(values *[]float64) float64 {
	if n := rnd.Choice(len(*values)); n < 0 {
		return math.NaN()
	} else {
		return (*values)[n]
	}
}

// From precalculated cumulative mass function table c,
// calculate loaded uniform (weighted) random integer in range [0, len(c)).
//
// The sequence must be non-negative, non-decreasing.
// The last element of the sequence must be greater than 0.
func (rnd *LCPRNG) Loaded(c *[]int) int {
	r := len(*c) - 1
	if r > 0 { // data present and not single
		n := rnd.Choice((*c)[r]) // last c is "probabilityDown"
		for l := 0; l < r; {     // binary search
			m := (l + r) / 2
			if n < (*c)[m] { // c[m] = ∑ "probabilityUp" to m
				r = m
			} else {
				l = m + 1
			}
		} // at this point l = r
	}
	return r
}

// Weighted-uniform random variable.
//
//	ProbabilityUp[i] = w[i]
//	ProbabolityDown  = ∑ w
func (rnd *LCPRNG) Weighted(w *[]int) int {
	t := 0       // total mass (probabilityDown)
	c := []int{} // cumulative mass table
	for _, m := range *w {
		if m < 0 {
			return -1 // no negative mass
		}
		t += m
		c = append(c, t)
	}
	if t == 0 {
		return rnd.Index(&c) // random photon
	} else {
		return rnd.Loaded(&c)
	}
}

// Uniform random number in range (0, 1).
func (rnd *LCPRNG) Random() float64 {
	var n uint64
	for n == 0 {
		n = rnd.Next() >> 11 // trim to 53 bits mantissa
	}
	const ε float64 = 0x1p-53 // ε = 2⁻⁵³
	return ε * float64(n)
}

// Random angle (0, 2π).
func (rnd *LCPRNG) Angle() float64 {
	const τ float64 = 2 * math.Pi // τ = 2π = 0x3243F6A8885A3p-47
	return τ * rnd.Random()
}

// Uniform random number in range (0, x).
func (rnd *LCPRNG) Uniform(x float64) float64 {
	if x != 0 {
		x *= rnd.Random()
	}
	return x
}

// Uniform random number in range (a, b).
func (rnd *LCPRNG) Range(a, b float64) float64 {
	return a + rnd.Uniform(b-a)
}

// True with probability p.
func (rnd *LCPRNG) Bernoulli(p float64) bool {
	return (p >= 1) || (p > 0 && p > rnd.Random())
}

// Binomial distribution random variable.
func (rnd *LCPRNG) Binomial(n int, p float64) (b int) {
	if p <= 0 || n <= 0 {
		return 0
	} else if p >= 1 {
		return n
	}
	const limit = 50
	x, q := float64(n), 1-p
	if (n > limit) && (x*p > 9*q) && (x*q > 9*p) { // Central Limit Theorem
		x *= p           // mean
		q *= x           // variance
		q = math.Sqrt(q) // standard deviation
		for b = -1; (b < 0) || (b > n); {
			b = rnd.Discrete(x, q)
		}
	} else {
		for ; n > 0; n-- {
			if rnd.Bernoulli(p) {
				b++
			}
		}
	}
	return
}

// Exponential distribution random variable.
func (rnd *LCPRNG) Exponential(ƛ ...float64) float64 {
	e := -math.Log1p(-rnd.Random()) // domain = [0, 36.7368]
	if len(ƛ) > 0 {
		e /= ƛ[0]
	}
	return e
}

// Pascal (negative binomial) distribution random variable.
func (rnd *LCPRNG) Pascal(r int, p float64) (n float64) {
	if r <= 0 || p <= 0 {
		n = math.Inf(1)
	} else if p < 1 {
		for p = -math.Log1p(-p); r > 0; r-- {
			n += math.Floor(rnd.Exponential(p))
		}
	}
	return
}

// Geometric distribution random variable.
//
//	p(n) = p · (1 - p)ⁿ
func (rnd *LCPRNG) Geometric(p float64) float64 {
	return rnd.Pascal(1, p)
}

// Rayleigh distribution random variable.
func (rnd *LCPRNG) Rayleigh(σ float64) float64 {
	if σ != 0 {
		σ *= math.Sqrt(2 * rnd.Exponential()) // Box-Muller transform
	}
	return σ
}

// Arcus distribution random variable (-1, 1).
//
//	μ = 0
//	σ² = 1/2
func (rnd *LCPRNG) Arcus() float64 {
	return math.Sin(rnd.Angle())
}

// ArcSine distribution random variable (0, 1).
//
//	μ = 1/2
//	σ² = 1/8
func (rnd *LCPRNG) ArcSine() float64 {
	return (rnd.Arcus() + 1) / 2
}

// Gauss distribution random variable.
//
//	μ = 0
//	σ = 1
func (rnd *LCPRNG) Gauss() float64 {
	return rnd.Rayleigh(rnd.Arcus()) // domain = [±8.57167]
}

// Normal distribution random variable.
func (rnd *LCPRNG) Normal(μ, σ float64) float64 {
	if σ != 0 {
		σ *= rnd.Gauss()
	}
	return μ + σ
}

// Discrete normal distribution random variable.
func (rnd *LCPRNG) Discrete(μ, σ float64) int {
	return int(math.Round(rnd.Normal(μ, σ)))
}

// Skew-normal distribution random variable.
func (rnd *LCPRNG) SkewNormal(ξ, ω, ɑ float64) (s float64) {
	const limit = 1024
	if ω != 0 {
		switch { // some special cases
		case ɑ == 0:
			s = rnd.Gauss()
		case ɑ < -limit:
			s = -math.Abs(rnd.Gauss())
		case ɑ > +limit:
			s = +math.Abs(rnd.Gauss())
		default:
			u, v := rnd.Target(1)
			if u > v {
				u, v = v, u
			}
			switch {
			case ɑ == -1:
				s = u
			case ɑ == +1:
				s = v
			default:
				s = (u*(1-ɑ) + v*(1+ɑ)) / math.Sqrt(2*(1+ɑ*ɑ))
			}
		}
		s *= ω
	}
	s += ξ
	return
}

// Log-normal distribution random variable.
func (rnd *LCPRNG) LogNormal(μ, σ float64) float64 {
	return math.Exp(rnd.Normal(μ, σ))
}

// Exponentially modified normal distribution random variable.
func (rnd *LCPRNG) ExpNormal(μ, σ, ƛ float64) float64 {
	return rnd.Normal(μ, σ) + rnd.Exponential(ƛ)
}

// Laplace  distribution random variable.
func (rnd *LCPRNG) Laplace(μ, b float64) (l float64) {
	if b > 0 {
		l = b * rnd.Exponential()
		if rnd.Bernoulli(0.5) {
			l = -l
		}
	}
	l += μ
	return
}

// Cauchy  distribution random variable.
func (rnd *LCPRNG) Cauchy(x0, ɣ float64) (c float64) {
	if ɣ > 0 {
		c = ɣ * math.Tan(rnd.Angle()) // due to inexact π: tan(π/2) = 16331239353195392
	}
	c += x0
	return
}

// Tukey distribution random variable.
func (rnd *LCPRNG) Tukey(ƛ float64) float64 {
	p := rnd.Random()
	q := 1 - p
	if ƛ == 0 {
		return math.Log(p / q)
	} else {
		return (math.Pow(p, ƛ) - math.Pow(q, ƛ)) / ƛ
	}
}

// Logistic distribution random variable.
func (rnd *LCPRNG) Logistic(μ, s float64) (l float64) {
	if s > 0 {
		l = s * rnd.Tukey(0)
	}
	l += μ
	return
}

// Poisson distribution random variable.
//
//	p(n) = exp(-ƛ) · ƛⁿ / n!
func (rnd *LCPRNG) Poisson(ƛ float64) (n int) {
	const limit = 256
	if ƛ > 0 {
		if ƛ < limit { // Knuth method
			l := math.Exp(-ƛ)
			for p := rnd.Random(); p > l; n++ {
				p *= rnd.Random()
			}
		} else { // Variance stabilizing
			const adj = 0.25 // adjustment (empirical = variance)
			x := rnd.Normal(math.Sqrt(ƛ-adj), 0.5)
			n = int(math.Round(x * x))
		}
	}
	return
}

// Skellam distribution random variable.
func (rnd *LCPRNG) Skellam(μ1, μ2 float64) (n int) {
	if μ1 >= 0 && μ2 >= 0 {
		n = rnd.Poisson(μ1) - rnd.Poisson(μ2)
	}
	return
}

// Hermite distribution random variable.
func (rnd *LCPRNG) Hermite(ɑ1, ɑ2 float64) (n int) {
	if ɑ1 >= 0 && ɑ2 >= 0 {
		n = rnd.Poisson(ɑ1) + 2*rnd.Poisson(ɑ2)
	}
	return
}

// χ² distribution random variable with k degrees of freedom.
//
// Sum of k squared Gauss randoms.
func (rnd *LCPRNG) ChiSquared(k int) (x float64) {
	const limit = 256
	if k < limit {
		if k > 1 {
			for x = 1; k > 1; k -= 2 {
				x -= rnd.Uniform(x)
			}
			x = -2 * math.Log(x)
		}
		if k > 0 {
			x += rnd.Exponential() * (rnd.Arcus() + 1)
		}
	} else { // Central Limit Theorem
		x = float64(k)
		x = rnd.Normal(x, math.Sqrt(2*x))
	}
	return
}

// χ distribution random variable with k degrees of freedom.
//
// Length of k-dimensional vector with Gauss random coordinates.
func (rnd *LCPRNG) Chi(k int) (x float64) {
	switch { // speed up by special cases
	case k == 1:
		x = math.Abs(rnd.Gauss()) // Half-Normal distribution
	case k == 2:
		x = rnd.Rayleigh(1) // Rayleigh distribution
	case k >= 3:
		x = math.Sqrt(rnd.ChiSquared(k))
	}
	return
}

// Erlang distribution random variable.
func (rnd *LCPRNG) Erlang(k int, ƛ float64) float64 {
	return rnd.ChiSquared(2*k) / (2 * ƛ)
}

// Γ distribution random variable.
//
// Sum of α Exponential(β) randoms.
//
//	μ = ɑ / β
//	σ² = ɑ / β²
func (rnd *LCPRNG) Gamma(ɑ float64, β ...float64) (g float64) {
	if ɑ > 0 {
		t, a := math.Modf(2 * ɑ) // trunc & frac
		if a > 0 {
			a /= 2
			// Ahrens-Dieter acceptance-rejection method
			for u, r, l, f := a+math.E, 1/a, a-1, false; !f; {
				if rnd.Uniform(u) < math.E {
					g = math.Pow(rnd.Random(), r)
					f = rnd.Bernoulli(math.Exp(-g))
				} else {
					g = 1 + rnd.Exponential()
					f = rnd.Bernoulli(math.Pow(g, l))
				}
			}
		}
		g += rnd.ChiSquared(int(t)) / 2 // add integer part
		if len(β) > 0 {
			g /= β[0]
		}
	}
	return
}

// Β distribution random variable.
func (rnd *LCPRNG) Beta(ɑ, β float64) (b float64) {
	if ɑ > 0 && β > 0 {
		switch { // some special cases
		case ɑ == 1 && β == 1:
			b = rnd.Random() // uniform
		case ɑ == 1:
			b = 1 - math.Pow(1-rnd.Random(), 1/β) // maximum of β uniform randoms
		case β == 1:
			b = math.Pow(rnd.Random(), 1/ɑ) // minimum of α uniform randoms
		case ɑ == 0.5 && β == 0.5:
			b = rnd.ArcSine()
		default:
			b = rnd.Gamma(ɑ)
			if b != 0 {
				b /= b + rnd.Gamma(β)
			}
		}
	}
	return
}

// Β' distribution random variable.
//
//	Γ(α) / Γ(β) = Γ(α, Γ(β))
func (rnd *LCPRNG) BetaPrime(ɑ, β float64) (b float64) {
	b = rnd.Beta(ɑ, β)
	if b != 0 && b != 1 {
		b /= 1 - b
	}
	return
}

// Student's t-distribution random variable with ν degrees of freedom.
//
// Normal distribution with
//
//	μ = 0
//	σ² = ν / χ²(ν)
//
// For ν -> ∞, σ -> 1
//
//	StudentsT(∞) = Normal(0, 1) = Gauss()
func (rnd *LCPRNG) StudentsT(ν float64) (t float64) {
	if ν > 0 {
		t = rnd.Gauss()
		if !math.IsInf(ν, 1) {
			ν /= 2
			t *= math.Sqrt(ν / rnd.Gamma(ν))
		}
	}
	return
}

// Snedecor's F-distribution random variable.
//
//	(χ²(d₁) / d₁) / (χ²(d₂) / d₂)
func (rnd *LCPRNG) SnedecorsF(d1, d2 float64) (f float64) {
	if d1 > 0 && d2 > 0 {
		f = rnd.BetaPrime(d1/2, d2/2) * d2 / d1
	}
	return
}

// Dirichlet distribution random array which sum is equal to 1.
//
// Weighted random cuts.
//
//	s = Σ ɑ
//	μᵢ = ɑᵢ / s
//	σᵢ = sqrt(ɑᵢ ✶ (s - ɑᵢ) / (s + 1)) / s
func (rnd *LCPRNG) Dirichlet(ɑ ...float64) (d []float64) {
	if n := len(ɑ); n > 0 {
		d = make([]float64, n)
		if n == 1 {
			d[0] = 1
		} else {
			var s float64
			for i, a := range ɑ {
				d[i] = rnd.Gamma(a)
				s += d[i]
			}
			if s == 0 {
				for i := range d {
					d[i] = rnd.Exponential()
					s += d[i]
				}
			}
			for i := range d {
				d[i] /= s
			}
		}
	}
	return
}

// Maxwell–Boltzmann distribution random variable (3 degrees of freedom).
//
// Random speed of particle (m/s) in 3D space where
//
//	M = molar mass (g/mol)
//	T = temperature (°C)
//
// For example:
//
//	M (oxygen molecule O₂) = 16 · 2 = 32 g/mol
//	T (room temperature) = 25°C
//
//	Maxwellian(32, 25)
func (rnd *LCPRNG) Maxwellian(M, T float64) (v float64) {
	const (
		O = -273.15       // Absolute zero (°C)
		c = 299792458     // Speed of light (m/s)
		N = 6.02214076e23 // Avogadro constant (mol⁻¹)
		k = 1.380649e-23  // Boltzmann constant (J/K)
		R = N * k * 1000  // Ideal gas constant (scaled kg/g)
	)
	T -= O // temperature in Kelvin
	switch {
	case M < 0 || T < 0: // unknown (yet)
		v = math.NaN()
	case M == 0: // photon
		v = c
	case T == 0: // absolute zero
		v = 0
	default: // Brownian motion
		v = math.Sqrt(rnd.ChiSquared(3) * R * T / M)
		v = math.Min(v, c)
	}
	return
}

// Inverse Gausian distribution random variable.
func (rnd *LCPRNG) Wald(μ, ƛ float64) (w float64) {
	if μ > 0 && ƛ > 0 {
		ƛ *= 2
		w = μ * rnd.ChiSquared(1)
		w = μ * (w - math.Sqrt(w*(2*ƛ+w)))
		w = μ + w/ƛ
		if rnd.Uniform(μ+w) > μ {
			w = μ * μ / w
		}
	}
	return
}

// Pareto  distribution random variable.
func (rnd *LCPRNG) Pareto(xm, ɑ float64) (p float64) {
	if xm > 0 && ɑ > 0 {
		p = xm * math.Exp(rnd.Exponential(ɑ))
	}
	return
}

// Lomax  distribution random variable.
func (rnd *LCPRNG) Lomax(ɑ, ƛ float64) float64 {
	return rnd.Pareto(ƛ, ɑ) - ƛ
}

// Weibull  distribution random variable.
func (rnd *LCPRNG) Weibull(ƛ, k float64) (w float64) {
	if ƛ > 0 && k > 0 {
		w = ƛ * math.Pow(rnd.Exponential(), 1/k)
	}
	return
}

// Logarithmic-uniform random variable.
func (rnd *LCPRNG) Logarithmic(a, b float64) (l float64) {
	if a > b {
		a, b = b, a
	}
	if a > 0 {
		if a == b {
			l = a
		} else {
			l = math.Exp(rnd.Range(math.Log(a), math.Log(b)))
			if l < a {
				l = a
			} else if l > b {
				l = b
			}
		}
	}
	return
}

// Benford law random integer in range [m, n] for positive arguments.
func (rnd *LCPRNG) Benford(m, n int) (b int) {
	if m > n {
		m, n = n, m
	}
	if m > 0 {
		if m == n {
			b = m
		} else {
			b = int(math.Exp(rnd.Range(math.Log(float64(m)), math.Log(float64(n+1)))))
			if b < m {
				b = m
			} else if b > n {
				b = n
			}
		}
	}
	return
}

// Irwin-Hall distribution random variable.
func (rnd *LCPRNG) IrwinHall(x float64) (h float64) {
	if x > 0 {
		const limit = 64
		if x < limit {
			t, f := math.Modf(x)
			for n := int(t); n > 0; n-- {
				h += rnd.Random()
			}
			h += rnd.Uniform(f)
		} else { // Central Limit Theorem
			h = rnd.Normal(x/2, math.Sqrt(x/12))
		}
	}
	return
}

// Bates distribution random variable.
func (rnd *LCPRNG) Bates(x, a, b float64) (c float64) {
	if x > 0 {
		c = (b - a) * rnd.IrwinHall(x) / x
	}
	c += a
	return
}

// Triangulat distribution random variable.
func (rnd *LCPRNG) Triangular(a, b, mode float64) (t float64) {
	if a > b {
		a, b = b, a
	}
	if a <= mode && mode <= b {
		w, d := b-a, mode-a
		if x := rnd.Uniform(w); x < d {
			t = a + math.Sqrt(x*d)
		} else {
			t = b - math.Sqrt((w-x)*(b-mode))
		}
	}
	return
}

// Sort integers with random pivot.
func (rnd *LCPRNG) Sort(x *[]int) {
	const treshold = 16 // algorithm selection treshold

	type part struct {
		l, r int
	}
	q := part{0, len(*x) - 1}
	queue := []part{q} // partition queue

	var l, r, p int // left, right and pivot

	qsort := func() { // Quick Sort by Sir C. A. R. Hoare (1960)
		l, r = q.l, q.r
		// p = (*x)[(l+r)/2] // middle pivot
		p = (*x)[rnd.Int(l, r)] // random pivot
		for l <= r {
			for (*x)[l] < p {
				l++
			}
			for p < (*x)[r] {
				r--
			}
			if l <= r {
				(*x)[l], (*x)[r] = (*x)[r], (*x)[l]
				l++
				r--
			}
		}
		if q.l < r {
			queue = append(queue, part{q.l, r})
		}
		if l < q.r {
			queue = append(queue, part{l, q.r})
		}
	}

	isort := func() { // Insertion Sort (better for short arrays)
		for r = q.l + 1; r <= q.r; r++ {
			p = (*x)[r]
			for l = r; l > q.l && (*x)[l-1] > p; l-- {
				(*x)[l] = (*x)[l-1]
			}
			(*x)[l] = p
		}
	}

	for len(queue) > 0 { // main loop
		q, queue = queue[0], queue[1:]
		if (q.r - q.l) > treshold {
			qsort()
		} else {
			isort()
		}
	}
}

// Weighted-uniform random variation k of n elements
// where weights = tuning, n = len(tuning) and k = podium,
// calculated by race simulation standing list.
func (rnd *LCPRNG) Podium(podium int, tuning *[]int) (stand []int) { // not optimised, tested
	cars := len(*tuning) // number of cars

	// censor podium
	if podium < 0 {
		podium = 0
	} else if podium > cars {
		podium = cars
	}

	stand = make([]int, podium) // standing list
	stop := make([]bool, cars)  // is car stop
	total := 0                  // total distance remaining
	place := 0                  // battle for place
	count := 0                  // number of cars who completed the race
	finish := 0                 // finish line
	car := 0                    // current car
	speed := 0                  // and speed
	race := count < podium      // Gentlemen, start your engines!

	if race { // convoy head (favorites)
		for _, speed = range *tuning {
			if speed > 0 {
				total += speed
			}
		}
		for (total != 0) && race {
			finish = rnd.Choice(total)
			for car, speed = range *tuning {
				if !stop[car] && speed > 0 {
					finish -= speed
					if finish < 0 { // chequered flag
						stop[car] = true // stop the car
						total -= speed
						stand[place] = car
						place++
						count++
						race = count < podium
						break
					}
				}
			}
		}
	}

	if race { // convoy body (uniform)
		for _, speed = range *tuning {
			if speed == 0 {
				total++
			}
		}
		for (total != 0) && race {
			finish = rnd.Choice(total)
			for car, speed = range *tuning {
				if !stop[car] && speed == 0 {
					finish--
					if finish < 0 { // chequered flag
						stop[car] = true // stop the car
						total--
						stand[place] = car
						place++
						count++
						race = count < podium
						break
					}
				}
			}
		}
	}

	if race { // convoy tail
		for _, speed = range *tuning {
			if speed < 0 {
				total -= speed
			}
		}
		place = cars // backwards
		for (total != 0) && race {
			finish = rnd.Choice(total)
			for car, speed = range *tuning {
				if !stop[car] {
					finish += speed
					if finish < 0 { // chequered flag
						stop[car] = true // stop the car
						total += speed
						place--
						if place < podium {
							stand[place] = car
							count++
							race = count < podium
						}
						break
					}
				}
			}
		}
	}

	return
}

// Weighted random permutation.
func (rnd *LCPRNG) Race(tuning *[]int) []int {
	return rnd.Podium(len(*tuning), tuning)
}

// Random forest.
//
// Random string of nested parenthesis.
//
// TAOCP Vol 4a, p 453, Algorithm W.
func (rnd *LCPRNG) Forest(n int) (f string) {
	for p, q := n, n; q > 0; {
		if rnd.Choose((q+p)*(q-p+1), (q-p)*(q+1)) {
			f += ")"
			q--
		} else {
			f += "("
			p--
		}
	}
	return
}

// Cut deck of cards.
func (rnd *LCPRNG) CutDeck(deck *[]int) (l, r []int) {
	c := rnd.Binomial(len(*deck), 0.5)
	l, r = append(l, (*deck)[:c]...), append(r, (*deck)[c:]...)
	return
}

// Interleave cards from left and right hand.
//
// Gilbert-Shannon-Reeds model.
func (rnd *LCPRNG) DoveTail(l, r *[]int) (d []int) {
	i, j := len(*l), len(*r)
	n := i + j
	d = make([]int, n)
	for n > 0 {
		if rnd.Choose(n, i) {
			n--
			i--
			d[n] = (*l)[i]
		} else {
			n--
			j--
			d[n] = (*r)[j]
		}
	}
	return
}

// Riffle shuffle deck of cards.
func (rnd *LCPRNG) RiffleShuffle(deck *[]int) {
	if n := len(*deck); n > 1 {
		// by Bayer & Diaconis (n = 8 for standard deck)
		for n = int(math.Log2(float64(n)) * 1.5); n > 0; n-- {
			l, r := rnd.CutDeck(deck)
			(*deck) = rnd.DoveTail(&l, &r)
		}
	}
}

// Array of n random integers which sum is equal to s.
//
//	μ = s / n
//	σ = sqrt(s ✶ (n - 1)) / n
//
// Deli špil od s karata na n približno jednakih delova.
//
// Metoda određuje "koliko dinara će sakupiti svako dete",
// kada kum na "Kume, izgoreti kesa!" baci s dinara, a ispred crkve se nalazi n dece.
func (rnd *LCPRNG) Scatter(s, n int) (d []int) {
	if n > 0 {
		d = make([]int, n)
		if s != 0 {
			const l = 2 * 53 * math.Ln2
			if n > 1 && math.Abs(float64(s)) < float64(n-1)*l { // Bernoulli method
				for ; s > 0; s-- {
					d[rnd.Choice(n)]++
				}
				for ; s < 0; s++ {
					d[rnd.Choice(n)]--
				}
			} else { // Central Limit Theorem
				for n > 1 {
					t, c := float64(s), float64(n) // total and count
					n--
					d[n] = rnd.Discrete(t/c, math.Sqrt(math.Abs(t)*(c-1))/c)
					s -= d[n]
				}
				d[0] = s
			}
		}
	}
	return
}

// Random point on circle of radius r.
func (rnd *LCPRNG) Circle(r float64) (x, y float64) {
	if r != 0 {
		y, x = math.Sincos(rnd.Angle())
		x, y = r*x, r*y
	}
	return
}

// Uniform random point in disc of radius r.
func (rnd *LCPRNG) Disc(r float64) (x, y float64) {
	return rnd.Circle(r * math.Sqrt(rnd.Random()))
}

// Shooting on target bullet dispersion.
func (rnd *LCPRNG) Target(dispersion float64) (x, y float64) {
	return rnd.Circle(rnd.Rayleigh(dispersion))
}

// Color pixel random dither.
func (rnd *LCPRNG) Dither(r, g, b byte) bool {
	const ( // CIELAB white-point
		x = 0.212671232040624
		y = 0.715159645674898
		z = 1 - (x + y)
	)
	gray := float64(r)*x + float64(g)*y + float64(b)*z
	return rnd.Bernoulli(gray / 255)
}

// House edge random variable for given return-to-player.
//
//	μ = 1 - rtp
//	σ = rtp
func (rnd *LCPRNG) Edge(rtp float64) float64 {
	if rtp > 0 {
		return 1 - rtp*rnd.Exponential()
	} else {
		return 1
	}
}

// 3 dice roll.
//
// Returns random dice roll (111-666), virtue (1-56) and frequency (1, 3, 6).
//
// Ludus Clericalis, TAOCP 4b, pp 493-494.
func (rnd *LCPRNG) SicBo() (dice []int, virtue, freq int) {
	d := dice
	for roll := rnd.Choice(216); len(d) < 3; roll /= 6 {
		d = append(d, roll%6+1)
	}
	dice = append(dice, d...) // random variation
	rnd.Sort(&d)              // random combination
	virtue = 56 - ((6-d[0])*(7-d[0])*(8-d[0])/6 + (6-d[1])*(7-d[1])/2 + (6-d[2])/1)
	switch {
	case d[0] == d[2]: // triplet
		freq = 1
	case d[0] != d[1] && d[1] != d[2]: // singles
		freq = 6
	default: // double + single
		freq = 3
	}
	return
}

// Slot reels stop positions and grid.
func (rnd *LCPRNG) Slot(reels *[][]int, height ...int) (stop []int, grid [][]int) {
	l := len(height)
	for i, r := range *reels {
		s := rnd.Index(&r)
		stop = append(stop, s)
		if s >= 0 {
			r = append(r[s:], r[:s]...)
			if l > 0 {
				if h := height[i%l]; h < len(r) {
					r = r[:h]
				}
			}
		}
		grid = append(grid, r)
	}
	return
}

// 2-adic multiplicative inverse for odd o else 0.
//
//	o · r = 1 (mod 2⁶⁴)
//
// practically
//
//	r = 1 / o
func MulInv64(o uint64) (r uint64) {
	if o != 0 {
		o /= -o & o // trim right zeroes
		for m, b := uint64(0), uint64(1); b != 0; b <<= 1 {
			if m |= b; o*r&m != 1 {
				r |= b
			}
		}
	}
	return
}

// n!
func Factorial(n int) float64 {
	return math.Gamma(float64(n + 1))
}

// n! / (n - k)!
func FallFact(n, k int) (f float64) {
	if n < 0 || k <= n {
		for f = 1; k > 0; n, k = n-1, k-1 {
			f *= float64(n)
		}
	}
	return
}

// Binomial coefficient.
//
//	n! / (n - k)! / k!
func Binomial(n, k int) (b float64) {
	b = 1
	if n < 0 { // Newton extension
		if k <= n {
			k = n - k
		}
		if 0 <= k {
			n = k - n - 1
		}
		if k&1 != 0 {
			b = -b
		}
	}
	if 0 <= k && k <= n { // Pascal triangle
		if l := n - k; k > l {
			k = l
		}
		for i := 1; i <= k; i, n = i+1, n-1 {
			b = b * float64(n) / float64(i) // do not change
		}
	} else {
		b = 0
	}
	return
}

// Multinomial coefficient.
//
//	(k₀ + k₁ + k₂ + ··· )! / (k₀! ✶ k₁! ✶ k₂! ✶ ··· )
func Multinomial(k ...int) (m float64) {
	m = 1
	n := 0
	for _, j := range k {
		n += j
		if m *= Binomial(n, j); m == 0 {
			break
		}
	}
	return
}

// Hyper-geometric distribution probability.
//
// Equivalent to Excel
//
//	HYPGEOMDIST(hits, draw, succ, size).
func HypGeomDist(hits, draw, succ, size int) (prob float64) {
	if prob = Binomial(succ, hits); prob != 0 {
		if prob *= Binomial(size-succ, draw-hits); prob != 0 {
			prob /= Binomial(size, draw)
		}
	}
	return
}

// Calculate combination index and probability.
func Ludus(sides int, dice ...int) (total, index int, prob float64) {
	if sides >= 0 {
		{
			var r LCPRNG
			r.Randomize()
			r.Sort(&dice)
		}
		k := len(dice)
		n := sides + k - 1
		total = int(Binomial(n, k))
		index, prob = total, 1
		l, c := 0, 0
		for i, d := range dice {
			if d < 1 || d > sides {
				return 0, 0, 0 // inconsistent
			}
			index -= int(Binomial(n-i-d, k-i))
			if d != l {
				l, c = d, 0
			}
			c += sides
			prob = prob * float64(i+1) / float64(c)
		}
	}
	return
}

// # Cheap rng for non-rgs stuff
//
// Whole Sort Of General Mish-Mash (H2G2).
var RND LCPRNG

func init() {
	RND.Randomize()
}
