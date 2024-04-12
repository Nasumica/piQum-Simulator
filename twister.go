package main

import (
	"sync"
	"time"
)

// # Mersenne Twister 19937 - 64 bit random number generator
//
// https://itechlabs.com/certificates/Srbislav/RNG_Certificate_Srbislav_27022015.pdf
type MT_19937_64 struct {
	index int        // state index
	state []uint64   // state vector
	dog   sync.Mutex // Šarko
}

const (
	mt_size  int    = (19937 + 64 - 1) / 64 // state size in octabytes
	mt_lmask uint64 = 1<<31 - 1             // lower bits mask
	mt_umask uint64 = ^mt_lmask             // upper bits mask
	mt_mat_a uint64 = 0xB5026F5AA96619E9    // constant by Makoto Matsumoto 松本 眞
	mt_mult  uint64 = 6364136223846793005   // constant by Donald Knuth
)

func (mt *MT_19937_64) New(seed ...uint64) *MT_19937_64 {
	if len(mt.state) != mt_size { // init
		n := uint64(time.Now().UnixNano()) / 100 // UnixNano always ends with '00' (Windows)
		if len(seed) > 0 {
			n = seed[0]
		}
		mt.state = make([]uint64, mt_size)
		for i := range mt.state {
			n += uint64(i)
			mt.state[i] = n
			n = (n ^ (n >> 62)) * mt_mult
		}
		mt.index = mt_size
	}
	return mt
}

// Next random value from Mersenne Twister.
func (mt *MT_19937_64) Next() (n uint64) {
	mt.dog.Lock()         // Šarko mora biti vezan ...
	defer mt.dog.Unlock() // ... dok se izvršava funkcija.

	if len(mt.state) != mt_size { // init
		mt.New()
	}

	if mt.index++; mt.index >= mt_size { // twist
		for i, j, m := 0, 1, mt_size/2; i < mt_size; i, j, m = i+1, (j+1)%mt_size, (m+1)%mt_size {
			n = (mt.state[i] & mt_umask) | (mt.state[j] & mt_lmask)
			mt.state[i] = mt.state[m] ^ (n >> 1) ^ ((n & 1) * mt_mat_a)
		}
		mt.index = 0
	}

	n = mt.state[mt.index]

	// xorshift64 improvement by George Marsaglia
	n ^= (n >> 29) & 0x5555555555555555
	n ^= (n << 17) & 0x71D67FFFEDA60000
	n ^= (n << 37) & 0xFFF7EEE000000000
	n ^= (n >> 43)

	return
}

// Next random value from Mersenne Twister in range [0, n]
func (mt *MT_19937_64) Limited(m, n uint64) uint64 {
	if m > n {
		m, n = n, m // swap limits
	}
	if n -= m; n != 0 {
		if n++; n == 0 { // n = 2⁶⁴
			n = mt.Next()
		} else { // acceptance-rejection
			f, u := -n/n+1, n                // f = 2⁶⁴ / n, u = n\
			for n = mt.Next() / f; n >= u; { // p(reject) = 2⁶⁴ % n / 2⁶⁴
				n = mt.Next() / f
			}
		}
	}
	return m + n
}

// Uniform random integer in range [m, n].
func (mt *MT_19937_64) Int(m, n int) int {
	if m > n {
		m, n = n, m
	}
	return m + int(mt.Limited(0, uint64(n-m)))
}

// Uniform random number in range [0, 1).
func (mt *MT_19937_64) Random() float64 {
	const eps float64 = 0x1p-53 // ε = 1 · 2⁻⁵³
	n := mt.Next() >> 11        // trim to 53 bit mantissa
	return eps * float64(n)
}

// Uniform random integer in range [0, n).
func (mt *MT_19937_64) Choice(n int) int {
	if n > 1 {
		c := uint64(n) // c = n
		c = mt.Limited(0, (-c/c+1)*c-1) % c
		return int(c)
	} else if n >= 0 {
		return n - 1
	} else {
		return -mt.Choice(-n)
	}
}

func (mt *MT_19937_64) Uint(n int) int {
	if n > 1 {
		c := uint64(n)
		for f, n := -c/c+1, c; c >= n; { // f = 2⁶⁴ / n
			c = mt.Next() / f // p(c >= n) = 2⁶⁴ % n / 2⁶⁴
		}
		return int(c)
	} else if n >= 0 {
		return n - 1 // -1 for n = 0 (special case)
	} else {
		panic("Invalid argument.")
	}
}

// True with probability k/n.
func (mt *MT_19937_64) Choose(n, k int) bool {
	if n < 0 {
		n, k = -n, -k // normalize denominator
	}
	return k >= 0 && (n <= k || mt.Choice(n) < k)
}

// Random array index for non-empty array else negative.
func (mt *MT_19937_64) Index(a *[]int) int {
	return mt.Choice(len(*a))
}

// Random item of array.
func (mt *MT_19937_64) Item(a *[]int) int {
	if i := mt.Index(a); i < 0 {
		return i
	} else {
		return (*a)[i]
	}
}
