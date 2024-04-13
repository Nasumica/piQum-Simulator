package main

// Author: Srbislav D. Nešić, srbislav.nesic@fincore.com

import (
	"errors"
	"fmt"
	"time"
)

type spil = uint64

const standard_pack spil = (1 << 52) - 1 // 52 ones

/*
	AAAA KKKK QQQQ JJJJ TTTT 9999 8888 7777 6666 5555 4444 3333 2222
	♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠ ♣♥♦♠
	0100 0100 0100 0100 0100 0000 0000 0000 0000 0000 0000 0000 0000 = 0x4444400000000
	poker.Hand(0x4444400000000) = 9EDCBA 9 royal flush [A♥ K♥ Q♥ J♥ 10♥]
*/

// # BitPoker Game Lore
type BitPoker struct {
	pack  spil       // pack of cards
	order [10]int    // hands order
	power [10]int    // hands strength
	Wheel bool       // Is 5432A = straight?
	Kinds [13]string // kinds
	Suits [4]string  // suits
	Names [10]string // hands names
}

// Init poker game.
func (bp *BitPoker) Init(pack spil, wheel bool, order []int, maps ...[]string) (err error) {
	bp.pack = pack & standard_pack
	if l := bp.Length(bp.pack); l < 5 {
		err = errors.New("unexpectedError")
	}

	bp.Wheel = wheel

	copy(bp.order[:], order)
	const order_mask spil = 1<<10 - 1 // 10 ones
	if len(order) != 10 || bp.Deflate(order) != order_mask {
		err = errors.New("unexpectedError")
	} else {
		for p, o := range bp.order {
			bp.power[o] = p
		}
	}

	bp.Kinds = [13]string{"2", "3", "4", "5", "6", "7", "8", "9", "T", "J", "Q", "K", "A"}
	bp.Suits = [4]string{"♠", "♦", "♥", "♣"} // preferans order
	// bp.Suits = [4]string{"s", "d", "h", "c"} // preferans order
	bp.Names = [10]string{"high card", "1 pair", "2 pair", "3 of kind",
		"straight", "flush", "full house", "4 of kind", "straight flush", "royal flush"}
	if l := len(maps); l > 0 {
		copy(bp.Kinds[:], maps[0])
		if l > 1 {
			copy(bp.Suits[:], maps[1])
			if l > 2 {
				copy(bp.Names[:], maps[2])
			}
		}
	}

	return
}

// Pack of cards used in poker game.
func (bp *BitPoker) Pack() spil {
	return bp.pack
}

// Rest of cards.
func (bp *BitPoker) Rest(hold spil) spil {
	return bp.pack & ^hold
}

// Engine pack of cards used in poker game.
func (bp *BitPoker) PackOfCards() []int {
	return bp.Expand(bp.pack)
}

// Classic game of poker with 52 cards.
/*
	pack = 0xFFFFFFFFFFFFF
	9        4  royal flush
	8       36  straight flush
	7      624  4 of kind
	6     3744  full house
	5     5108  flush
	4    10200  straight
	3    54912  3 of kind
	2   123552  2 pair
	1  1098240  1 pair
	0  1302540  high card
	Σ  2598960  total
*/
func (bp *BitPoker) Classic() {
	bp.Init(standard_pack, true, []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) // standard hands order
}

// Game of poker with 36 cards (6+).
/*
	pack = 0xFFFFFFFFF0000
	9        4  royal flush
	8       16  straight flush
	7      288  4 of kind
	5      484  flush
	6     1728  full house
	4     5100  straight
	3    16128  3 of kind
	2    36288  2 pair
	1   193536  1 pair
	0   123420  high card
	Σ   376992  total
*/
func (bp *BitPoker) SixUp() {
	// remove first 16 cards from standard pack
	const pack = standard_pack >> 16 << 16
	bp.Init(pack, false, []int{0, 1, 2, 3, 4, 6, 5, 7, 8, 9}) // advance 5 (flush) 1 position
}

// Game of poker with 32 cards (7+).
/*
	pack = 0xFFFFFFFF00000
	9        4  royal flush
	8       12  straight flush
	5      208  flush
	7      224  4 of kind
	6     1344  full house
	4     4080  straight
	3    10752  3 of kind
	2    24192  2 pair
	1   107520  1 pair
	0    53040  high card
	Σ   201376  total
*/
func (bp *BitPoker) SevenUp() {
	// remove first 20 cards from standard pack
	const pack = standard_pack >> 20 << 20
	bp.Init(pack, false, []int{0, 1, 2, 3, 4, 6, 7, 5, 8, 9}) // advance 5 (flush) 2 positions
}

// Ensures non empty pack of cards.
//
// If not specified, Classic Poker is assumed.
func (bp *BitPoker) ensure() {
	if bp.pack&standard_pack == 0 {
		bp.Classic()
	}
}

/*
Returns poker hand as string in format

	"SXXXXX H hand name"

where S is one digit hand strength, H is the hand code (same as S for standard 52 cards  deck)

	9 = "royal flush"
	8 = "straight flush"
	7 = "4 of kind"
	6 = "full house"
	5 = "flush"
	4 = "straight"
	3 = "3 of kind"
	2 = "2 pair"
	1 = "1 pair"
	0 = "high card"

and XXXXX is 5 hecadecimal digits cards codes:
2 to 9, 10=A, J=B, Q=C, K=D, A=E or A=1 in bicycle,
sorted in order of significance.

	fmt.Println(poker.Hand(0x1111100000000)) // [A♠ K♠ Q♠ J♠ T♠]
	9EDCBA 9 royal flush

	fmt.Println(poker.Hand(0x1000000008421)) // [A♠ 5♣ 4♥ 3♦ 2♠]
	454321 4 straight

Hands can be compared:

	var poker game.BitPoker
	if a, b := poker.Hand(x), poker.Hand(y); a > b {
		fmt.Println("x wins")
	} else if a < b {
		fmt.Println("y wins")
	} else {
		fmt.Println("split")
	}

There are 7462 different hands:

	075432 high card (weakest)
	076432 high card
	...
	8DCBA9 straight flush
	9EDCBA royal flush (strongest)
*/
func (bp *BitPoker) Hand(hold spil) string {
	const (
		kenta int = 0b11111          // 5 consecutive kinds
		A5432 int = 0b10000000011110 // A5432 mask
		royal int = 0b11111000000000 // AKQJT mask
	)

	bp.ensure()

	result, cards := "", bp.Inflate(hold&bp.pack)

	if len(cards) == 5 {
		kc, km, sm := map[int]int{}, 0, 0 // cards, kinds & suits
		for _, c := range cards {
			k, s := c>>2+1, c&3
			kc[k]++      // kinds count
			km |= 1 << k // kinds mask
			sm |= 1 << s // suits mask
		}

		bicycle := bp.Wheel && km == A5432          // five-high straight
		straight := bicycle || km/(km&-km) == kenta // 5 consecutive kinds or A5432
		flush := sm/(sm&-sm) == 1                   // only one suit

		type rank struct {
			k int // kind
			c int // count
		}
		var list []rank
		for k, c := range kc { // sort list order by c desc, k desc
			if bicycle && k == 13 {
				k = 0 // in bicycle A is the lowest card
			}
			r, i := rank{k, c}, len(list)
			list = append(list, r)
			for j := i - 1; i > 0 && (list[j].c < c || (list[j].c == c && list[j].k < k)); j-- {
				list[i], i = list[j], j
			}
			list[i] = r
		}

		hand := 0 // poker hand

		if flush {
			if straight {
				if km == royal {
					hand = 9 // "royal flush"
				} else {
					hand = 8 // "straight flush"
				}
			} else {
				hand = 5 // "flush"
			}
		} else if straight {
			hand = 4 // "straight"
		} else {
			value := 0
			for _, o := range list {
				value = 10*value + o.c
			}
			switch value {
			case 41:
				hand = 7 // "4 of kind"
			case 32:
				hand = 6 // "full house"
			case 311:
				hand = 3 // "3 of kind"
			case 221:
				hand = 2 // "2 pair"
			case 2111:
				hand = 1 // "1 pair"
			default:
				hand = 0 // "high card"
			}
		}

		strength := bp.power[hand] // mapped strength

		code := []byte{digit(strength)}
		for _, o := range list {
			for d := digit(o.k + 1); o.c > 0; o.c-- {
				code = append(code, d)
			}
		}
		code = append(code, ' ', digit(hand))

		result = string(code) + " " + bp.Names[hand]
	}

	return result
}

// Returns poker hand strength.
func (bp *BitPoker) HandStrength(code string) int {
	if code == "" {
		return -1
	}
	s := 0
	for i := 0; i < 6; i++ {
		s = 16*s + number(code[i])
	}
	s = 10*s + number(code[7])
	return s
}

// Returns poker hand code only.
func (bp *BitPoker) HandCode(hand string) int {
	if hand == "" {
		return -1
	}
	return number(hand[7])
}

// Returns poker hand code only, from 0 to 9.
func (bp *BitPoker) Code(hold spil) int {
	return bp.HandCode(bp.Hand(hold))
}

// Best poker hand with keep number of cards from hold and rest from desk.
//
// Returns hand code and selected cards from hold and desk.
func (bp *BitPoker) Holdem(hold, desk spil, keep int) (mc string, mh, md spil) {
	bp.ensure()                             // ensure pack
	hold, desk = hold&bp.pack, desk&bp.pack // remove non-standard cards (if any)
	desk ^= hold & desk                     // remove duplicates from desk (if any)

	var sh, sd sampler // hold and desk samples

	if sh.Init(hold, keep) && sd.Init(desk, 5-keep) {
		mb := mh
		for !sh.Eof() {
			h := sh.Next() // keep cards from hold

			sd.Reset()
			for !sd.Eof() {
				d := sd.Next() // 5-keep cards from desk

				b := h | d // test for the best
				if c := bp.Hand(b); c > mc || (c == mc && b > mb) {
					mc, mh, md, mb = c, h, d, b // current best hand
				}
			}
		}
	}

	return // code, hold, desk
}

// Best poker hand from all hold cards and rest from desk.
//
// Returns hand code and selected cards from desk.
func (bp *BitPoker) ElGordo(hold, desk spil) (string, spil) {
	bp.ensure()
	code, _, draw := bp.Holdem(hold, desk, bp.Length(hold&bp.pack))
	return code, draw
}

// Best poker hand from all hold cards and rest of pack.
//
// Returns hand code and joker replacements.
func (bp *BitPoker) JokerWild(hold spil) (string, spil) {
	bp.ensure()
	return bp.ElGordo(hold, hold^bp.pack)
}

// Best poker hand from pile of 5 or more cards.
//
// Returns hand code and selected cards from pile.
func (bp *BitPoker) BestOf(pile spil) (string, spil) {
	return bp.ElGordo(0, pile)
}

// What could this be?
func (bp *BitPoker) WhatIs(pile spil) (string, spil, spil) {
	bp.ensure()
	pile &= bp.pack
	if l := bp.Length(pile); l > 5 {
		c, best := bp.BestOf(pile)
		return c, best, best ^ pile
	} else if l < 5 {
		c, wild := bp.JokerWild(pile)
		return c, pile | wild, wild
	} else {
		c := bp.Hand(pile)
		return c, pile, 0
	}
}

// Texas Holdem best poker hand from hold and desk.
//
// Returns hand code and selected cards from hold and desk.
func (bp *BitPoker) TexasHoldem(hold, desk spil) (string, spil, spil) {
	desk ^= hold & desk                   // no duplicates
	code, best := bp.BestOf(hold | desk)  // best hand
	return code, hold & best, desk & best // split hand
}

// Omaha Holdem best poker hand with 2 cards from hold and 3 from desk.
//
// Returns hand code and selected cards from hold and desk.
func (bp *BitPoker) OmahaHoldem(hold, desk spil) (string, spil, spil) {
	return bp.Holdem(hold, desk, 2)
}

// Create cards indices from pile of cards.
func (bp *BitPoker) Inflate(pile spil) []int {
	r := []int{}
	/*
		// Naive way
		for i := 0; pile != 0; i++ {
			if pile&1 != 0 {
				r = append(r, i)
			}
			pile >>= 1
		}
		return r
	*/
	for pile != 0 {
		r = append(r, bp.Rightmost(pile))
		pile &= pile - 1
	}
	return r
}

// Create pile of cards from cards indices.
func (bp *BitPoker) Deflate(r []int) spil {
	var n spil
	for _, b := range r {
		if 0 <= b && b < 64 {
			n |= 1 << b
		}
	}
	return n
}

// Arange pile as cards in order of significance.
func (bp *BitPoker) Arange(pile spil) []int {
	bp.ensure()
	pile &= bp.pack
	d := bp.Inflate(pile)
	if l := len(d); l == 5 {
		p := bp.Hand(pile)
		b := []byte(p[1:5])
		for i, c := range b {
			k := number(c) - 1 // get kind
			if k == 0 {
				k = 13 // bicycle
			}
			k--
			for j := 4; j > i; j-- { // find card
				if d[j]>>2 == k {
					k = d[j]
					for l := j; l > i; l-- { // make room
						d[l] = d[l-1]
					}
					d[i] = k // move card
					break
				}
			}
		}
	} else {
		for i, j := 0, l-1; i < j; i, j = i+1, j-1 { // reflect
			d[i], d[j] = d[j], d[i]
		}
	}
	return d
}

// Debug view of cards list.
func (bp *BitPoker) Stringify(d []int, maps ...[]string) []string {
	k, s := bp.Kinds, bp.Suits
	if len(maps) > 0 {
		copy(k[:], maps[0])
		if len(maps) > 1 {
			copy(s[:], maps[1])
		}
	}
	h := make([]string, len(d))
	for i, c := range d {
		h[i] = k[c>>2] + s[c&3]
	}
	return h
}

// Human view of arranged cards list.
func (bp *BitPoker) Humanize(pile spil, maps ...[]string) []string {
	return bp.Stringify(bp.Arange(pile), maps...)
}

// Create pile of cards from engine cards.
func (bp *BitPoker) Squeeze(cards ...int) spil {
	bp.ensure()
	var pile spil
	for _, c := range cards {
		if 1 <= c && c <= 52 {
			pile |= 1 << (c - 1)
		}
	}
	return pile & bp.pack
}

// Expand pile of cards as engine cards.
func (bp *BitPoker) Expand(pile spil) []int {
	a := bp.Arange(pile)
	for i := range a {
		a[i]++
	}
	return a
}

// Determines players standing list according to given rule.
func (bp *BitPoker) PlayPokerHand(players []spil, desk spil, rule func(hold, desk spil) (string, spil, spil)) (codes []string, stand [][]int) {
	// calculate hands and sort list
	type hand struct {
		p    int    // player
		c    string // code
		h, d spil   // hold, desk (future use)
	}
	hands := make([]hand, len(players))
	for p, hold := range players {
		c, h, d := rule(hold, desk)
		x := hand{p, c, h, d}
		for q := p - 1; p > 0 && hands[q].c < c; q-- { // insertion sort
			p, hands[p] = q, hands[q]
		}
		hands[p] = x
	}

	// calculate standongs
	c, p := "", -1
	for i, x := range hands {
		if c != x.c || i == 0 {
			c, p, codes, stand = x.c, p+1, append(codes, x.c), append(stand, []int{})
		}
		stand[p] = append(stand[p], x.p)
	}

	return
}

// Determines players standing list according to Texas Holdem rule.
func (bp *BitPoker) PlayTexasHand(players []spil, desk spil) ([]string, [][]int) {
	return bp.PlayPokerHand(players, desk, bp.TexasHoldem)
}

// Determines players standing list according to Omaha Holdem rule.
func (bp *BitPoker) PlayOmahaHand(players []spil, desk spil) ([]string, [][]int) {
	return bp.PlayPokerHand(players, desk, bp.OmahaHoldem)
}

// Determines players standing list according to Classic Poker rule.
func (bp *BitPoker) PlayClassicHand(players []spil) ([]string, [][]int) {
	return bp.PlayTexasHand(players, 0)
}

// Convert number to digit.
func digit(b int) byte {
	if b < 10 {
		return byte(b + '0')
	} else {
		return byte(b + 'A' - 10)
	}
}

// Convert digit to number.
func number(d byte) int {
	if d < 'A' {
		return int(d) - '0'
	} else {
		return int(d) - 'A' + 10
	}
}

/*
Bit-trik koji prožima većinu metoda je
	n &= n -1
Ako je n oblika
	n     = x...x 1 0...0
gde su desno od 1 svi ostali bitovi 0, odnosno, 1 je poslednja jedinica u n, onda je
	n - 1 = x...x 0 1...1
odnosno, poslednji bit 1 postaje 0, a nule desno od njega postaju 1
Drugim rečima, n - 1 je negacija poslednje jedinice i bitova desno od nje, dok levi bitovi ostaju nepromenjeni.
Kada se nad ove dve vrednosti primeni bitwise-and (&) dobiće se vrednost
	x...x 1 0...0 = n
	x...x 0 1...1 = n - 1
&	-------------
	x...x 0 0...0 = n & (n - 1)
odnosno, biće "obrisan" poslednji bit 1 u n.

Sledeći bit-trik je
	n & -n
Opet, ako je n oblika
	n = x...x 1 0...0
Vrednost -n je predstavljena u potpunom komplementu (2-adic) kao -n = ~n + 1
	y...y 0 1...1 = ~n
				1 = 1
+	-------------
	y...y 1 0...0 = -n
gde su vrednosti bitova y negacije (komoplement) vrednsoti bitova x, odakle sledi x...x & y...y = 0.
Praktično, -n vrši negaciju svih bitova levo od poslednje jedinice, desni bitovi ostaju neproomenjeni.
Onda je
	x...x 1 0...0 = n
	y...y 1 0...0 = -n
&	-------------
	0...0 1 0...0 = n & -n
odnosno, iz n će biti "obrisani" svi bitovi sem poslednje jedinice.
Ovaj broj predstavlja najveći mogući stepen broja 2 kojim n može biti podeljen bez ostatka.
Izraz
	n / (n & -n) = x...x 1 0...0 / 0...0 1 0...0 = x...x 1
"precrtava" sve desne nule iz n, odnosno, šiftuje n udesno najviše moguće, poput
	for n & 1 == 0 {
		n >>= 1
	}
ali bez korrišćenja petlje (loop-less).

GO, umesto simbola ~ kao unarnog operatora bitwise-not, koristi simbol ^ (xor), pa je -n = ^n + 1.

Još jedan bit-trik koji autor često koristi, ali koji se ovde ne pojavljuje,
za n tipa uint64 i n > 1
	f = -n / n + 1
Iako je, naizgled, vrednost izraza 0, treba uzeti u obzir da se prvo izračunava unarni minus.
Kako je u 2-adic -n = 2⁶⁴ - n, onda je
	f = (2⁶⁴ - n) / n + 1
	  = 2⁶⁴ / n - n / n + 1
	  = 2⁶⁴ / n - 1 + 1
	  = 2⁶⁴ / n
Naravno, 1 - n / n = 0, jer je, u ovom slučaju, minus binarna operacija.
*/

// Count bit 1.
func (bp *BitPoker) Length(n spil) int {
	// Keringhan bitcount (faster for sparse n)
	o := 0
	for ; n != 0; o++ {
		n &= n - 1
	}
	return o
	/*
		// TAOCP 4A, 143
		const (
			a  uint64 = 0x0101010101010101
			m0        = 0x55 * a
			m1        = 0x33 * a
			m2        = 0x0f * a
		)
		y := uint64(n)
		y = y - ((y >> 1) & m0)
		y = (y & m1) + ((y >> 2) & m1)
		y = (y + (y >> 4)) & m2
		y = (a * y) >> 56
		return int(y)
	*/
}

/*
// Decompose variable as array of different items with only one bit set.
func Explode(n int) []int {
	e := []int{}
	for n != 0 {
		e = append(e, n&-n)
		n &= n - 1
	}
	return e
}

// Inversion of Explode.
func Implode(e []int) int {
	i := 0
	for _, b := range e {
		i |= b
	}
	return i
}
*/

// Leftmost bit 1 position.
func (bp *BitPoker) Leftmost(n spil) int { // = floor(lg(n))
	if n == 0 {
		return -1
	}
	/*
		// Naive way
		l := 63
		for ; n > 0; l-- {
			n <<= 1
		}
		return l
	*/
	// Knuth, TAOCP 4A, 142
	l := 0
	for b := 32; b > 0; b >>= 1 {
		if m := n >> b; m != 0 {
			l += b
			n = m
		}
	}
	return l
}

// Rightmost bit 1 position.
func (bp *BitPoker) Rightmost(n spil) int { // = Length(n - 1) + 1 - Length(n)
	if n == 0 {
		return -1
	}
	/*
		// Naive way
		r := 0
		for ; n&1 == 0; r++ {
			n >>= 1
		}
		return r
	*/
	// Knuth, TAOCP 4A, 142.
	n = (-n & n) * 0x03f79d71b4ca8b09 >> 58
	return []int{ // de Brujin cycle
		0, 1, 56, 2, 57, 49, 28, 3, 61, 58, 42, 50, 38, 29, 17, 4,
		62, 47, 59, 36, 45, 43, 51, 22, 53, 39, 33, 30, 24, 18, 12, 5,
		63, 55, 48, 27, 60, 41, 37, 16, 46, 35, 44, 21, 52, 32, 23, 11,
		54, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19, 9, 13, 8, 7, 6,
	}[n]
}

// # Cards combinations generator
//
// Algorithm L, TAOCP 4A, 358.
type sampler struct {
	n int    // total cards
	k int    // sample size
	e bool   // eof?
	b []spil // bits
	c []int  // counters
}

// Initialise generator with pile of cards and sample size and return success.
func (s *sampler) Init(pile spil, size int) bool {
	s.b, s.c = []spil{}, []int{}
	for pile != 0 {
		s.b = append(s.b, -pile&pile)
		pile &= pile - 1
	}
	s.n, s.k = len(s.b), size
	if s.e = s.k < 0 || s.n < s.k; !s.e {
		s.c = make([]int, s.k+2)
		s.c[s.k], s.c[s.k+1] = s.n, 0
		s.Reset()
	}
	return !s.e
}

// Reset generator to 1st combination and return success.
func (s *sampler) Reset() bool {
	if s.e = s.k < 0 || s.n < s.k; !s.e {
		for j := 0; j < s.k; j++ {
			s.c[j] = j
		}
	}
	return !s.e
}

// Next combination.
func (s *sampler) Next() (c spil) {
	if !s.e {
		var j int
		for j = 0; j < s.k; j++ {
			c |= s.b[s.c[j]]
		}
		for j = 0; s.c[j]+1 == s.c[j+1]; j++ {
			s.c[j] = j
		}
		if s.e = j >= s.k; !s.e {
			s.c[j]++
		}
	}
	return c
}

// End of file?
func (s *sampler) Eof() bool {
	return s.e
}

// # Game of Poker Lore
//
// (wrapper around BitPoker)
type Poker struct{ bp BitPoker }

// Init poker game.
func (pok *Poker) Init(cards []int, wheel bool, order []int, maps ...[]string) error {
	return pok.bp.Init(pok.bp.Squeeze(cards...), wheel, order, maps...)
}

// Classic game of poker with 52 cards.
/*
	9        4  royal flush
	8       36  straight flush
	7      624  4 of kind
	6     3744  full house
	5     5108  flush
	4    10200  straight
	3    54912  3 of kind
	2   123552  2 pair
	1  1098240  1 pair
	0  1302540  high card
	Σ  2598960  total
*/
func (pok *Poker) Classic() {
	pok.bp.Classic()
}

// Game of poker with 36 cards (6+).
/*
	9        4  royal flush
	8       16  straight flush
	7      288  4 of kind
	5      484  flush
	6     1728  full house
	4     5100  straight
	3    16128  3 of kind
	2    36288  2 pair
	1   193536  1 pair
	0   123420  high card
	Σ   376992  total
*/
func (pok *Poker) SixUp() {
	pok.bp.SixUp()
}

// Game of poker with 32 cards (7+).
/*
	9        4  royal flush
	8       12  straight flush
	5      208  flush
	7      224  4 of kind
	6     1344  full house
	4     4080  straight
	3    10752  3 of kind
	2    24192  2 pair
	1   107520  1 pair
	0    53040  high card
	Σ   201376  total
*/
func (pok *Poker) SevenUp() {
	pok.bp.SevenUp()
}

// Engine pack of cards used in poker game.
func (pok *Poker) PackOfCards() []int {
	return pok.bp.PackOfCards()
}

func (pok *Poker) Hand(hold []int) string {
	return pok.bp.Hand(pok.bp.Squeeze(hold...))
}

// Returns poker hand code only.
func (pok *Poker) HandCode(hand string) int {
	if hand == "" {
		return -1
	} else {
		return number(hand[7])
	}
}

// Returns poker hand code only, from 0 to 9.
func (pok *Poker) Code(hold []int) int {
	return pok.HandCode(pok.Hand(hold))
}

func (pok *Poker) Holdem(hold, desk []int, keep int) (string, []int, []int) {
	c, h, d := pok.bp.Holdem(pok.bp.Squeeze(hold...), pok.bp.Squeeze(desk...), keep)
	return c, pok.bp.Expand(h), pok.bp.Expand(d)
}

// Best poker hand from all hold cards and rest from desk.
//
// Returns hand code and selected cards from desk.
func (pok *Poker) ElGordo(hold, desk []int) (string, []int) {
	h := pok.bp.Squeeze(hold...)
	l := pok.bp.Length(h)
	code, _, draw := pok.bp.Holdem(h, pok.bp.Squeeze(desk...), l)
	return code, pok.bp.Expand(draw)
}

// Best poker hand from all hold cards and rest of pack.
//
// Returns hand code and joker replacements.
func (pok *Poker) JokerWild(hold []int) (string, []int) {
	r := pok.bp.Squeeze(hold...) ^ pok.bp.pack
	return pok.ElGordo(hold, pok.bp.Expand(r))
}

// Best poker hand from pile of 5 or more cards.
//
// Returns hand code and selected cards from pile.
func (pok *Poker) BestOf(pile []int) (string, []int) {
	return pok.ElGordo([]int{}, pile)
}

// What could this be?
func (pok *Poker) WhatIs(pile []int) (string, []int, []int) {
	c, h, w := pok.bp.WhatIs(pok.bp.Squeeze(pile...))
	return c, pok.bp.Expand(h), pok.bp.Expand(w)
}

func (pok *Poker) Arange(pile []int) []int {
	hold := pok.bp.Arange(pok.bp.Squeeze(pile...))
	for i := range hold {
		hold[i]++
	}
	return hold
}

func (pok *Poker) Stringify(cards []int, maps ...[]string) []string {
	k, s := pok.bp.Kinds, pok.bp.Suits
	if len(maps) > 0 {
		copy(k[:], maps[0])
		if len(maps) > 1 {
			copy(s[:], maps[1])
		}
	}
	h := make([]string, len(cards))
	for i, c := range cards {
		if 1 <= c && c <= 52 {
			c--
			h[i] = k[c>>2] + s[c&3]
		} else {
			h[i] = "*"
		}
	}
	return h
}

func (pok *Poker) Humanize(pile []int, maps ...[]string) []string {
	return pok.Stringify(pok.Arange(pile), maps...)
}

// Texas Holdem best poker hand from hold and desk.
//
// Returns hand code and selected cards from hold and desk.
func (pok *Poker) TexasHoldem(hold, desk []int) (string, []int, []int) {
	c, h, d := pok.bp.TexasHoldem(pok.bp.Squeeze(hold...), pok.bp.Squeeze(desk...))
	return c, pok.bp.Expand(h), pok.bp.Expand(d)
}

// Omaha Holdem best poker hand with 2 cards from hold and 3 from desk.
//
// Returns hand code and selected cards from hold and desk.
func (pok *Poker) OmahaHoldem(hold, desk []int) (string, []int, []int) {
	c, h, d := pok.bp.OmahaHoldem(pok.bp.Squeeze(hold...), pok.bp.Squeeze(desk...))
	return c, pok.bp.Expand(h), pok.bp.Expand(d)
}

// Determines players standing list according to given rule.
func (pok *Poker) PlayPokerHand(players [][]int, desk []int, rule func(hold, desk []int) (string, []int, []int)) (codes []string, stand [][]int) {
	// calculate hands and sort list
	type hand struct {
		p    int    // player
		c    string // code
		h, d []int  // hold, desk (future use)
	}
	hands := make([]hand, len(players))
	for p, hold := range players {
		c, h, d := rule(hold, desk)
		x := hand{p, c, h, d}
		for q := p - 1; p > 0 && hands[q].c < c; q-- { // insertion sort
			p, hands[p] = q, hands[q]
		}
		hands[p] = x
	}

	// calculate standongs
	c, p := "", -1
	for i, x := range hands {
		if c != x.c || i == 0 {
			c, p, codes, stand = x.c, p+1, append(codes, x.c), append(stand, []int{})
		}
		stand[p] = append(stand[p], x.p)
	}

	return
}

// Determines players standing list according to Texas Holdem rule.
func (pok *Poker) PlayTexasHand(players [][]int, desk []int) ([]string, [][]int) {
	return pok.PlayPokerHand(players, desk, pok.TexasHoldem)
}

// Determines players standing list according to Omaha Holdem rule.
func (pok *Poker) PlayOmahaHand(players [][]int, desk []int) ([]string, [][]int) {
	return pok.PlayPokerHand(players, desk, pok.OmahaHoldem)
}

// Determines players standing list according to Classic Poker rule.
func (pok *Poker) PlayClassicHand(players [][]int) ([]string, [][]int) {
	return pok.PlayTexasHand(players, []int{})
}

func (pok *Poker) Length(hold []int) int {
	return pok.bp.Length(pok.bp.Squeeze(hold...))
}

/*
Card Charasters True Type Font mapping:
var suits = []string{"}", "[", "{", "]"}
var kinds = []string{"2", "3", "4", "5", "6", "7", "8", "9", "=", "J", "Q", "K", "A"}
*/

var WildCardComb = [6][10]int{
	{1302540, 1098240, 123552, 54912, 10200, 5108, 3744, 624, 36, 4},
	{0, 169848, 0, 82368, 10332, 2696, 2808, 2509, 144, 20},
	{0, 0, 0, 13320, 3840, 888, 0, 3796, 216, 40},
	{0, 0, 0, 0, 0, 0, 0, 1142, 144, 40},
	{0, 0, 0, 0, 0, 0, 0, 0, 32, 20},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}

func CountCombs(w int) (card, wild [10]int) {
	var poker BitPoker
	poker.Classic()

	var base = [10]int{1302540, 1098240, 123552, 54912, 10200, 5108, 3744, 624, 36, 4}

	var left, right sampler

	left.Init(poker.Pack(), 5-w)
	for !left.Eof() {
		l := left.Next()
		right.Init(poker.Rest(l), w)
		best := 0
		for !right.Eof() {
			r := right.Next()
			hand := l | r
			code := poker.Code(hand)
			card[code]++
			if best < code {
				best = code
			}
		}
		wild[best]++
	}

	for i, b := range base {
		card[i] /= b
	}

	return
}

func (bp *BitPoker) MinMax(hold spil, min int) (worst, best string) {
	bp.ensure() // ensure pack
	if min < 2 {
		min = 2
	}

	hold &= bp.pack // remove non-standard cards (if any)
	desk := bp.Rest(hold)

	var hand, wild sampler // hold and desk samples

	if l := bp.Length(hold); l >= min && l <= 5 {
		if hand.Init(hold, l) && wild.Init(desk, 5-l) {
			for !hand.Eof() {
				h := hand.Next() // keep cards from hold
				wild.Reset()
				for !wild.Eof() {
					d := wild.Next() // 5-keep cards from desk
					test := h | d    // test for the best
					c := bp.Hand(test)
					if best == "" {
						worst, best = c, c
					} else if worst > c {
						worst = c
					} else if best < c {
						best = c
					}
				}
			}
		}
	}

	return
}

func (bp *BitPoker) Likelihood(hold spil) (list []int, sum int) {
	bp.ensure() // ensure pack

	hold &= bp.pack // remove non-standard cards (if any)

	start := time.Now()
	if l := bp.Length(hold); l < 0 {
		list = WildCardComb[0][:]
	} else if l < 5 {
		var wild sampler
		if wild.Init(bp.Rest(hold), 5-l) {
			list = make([]int, 10)
			for !wild.Eof() {
				list[bp.Code(hold|wild.Next())]++
			}
		}
	} else {
		var pile sampler
		pile.Init(hold, 5)
		list = make([]int, 10)
		for !pile.Eof() {
			list[bp.Code(pile.Next())]++
		}
	}
	for _, l := range list {
		sum += l
	}
	elapsed := time.Since(start).Seconds()
	fmt.Println(elapsed)

	return
}

func (pok *Poker) Likelihood(hold []int) ([]int, int) {
	return pok.bp.Likelihood(pok.bp.Squeeze(hold...))
}

func GCD(u, v uint64) uint64 {
	for v != 0 {
		u, v = v, u%v
	}
	return u
}

func MaxWinProb(deck, m, n int) (uint64, uint64) {
	dn := uint64(deck)
	for i := 1; i <= n; i++ {
		dn *= uint64(deck)
	}
	g := GCD(dn, 20)
	return 20 / g, dn / g
}

func HandProb(wheel bool, h ...int) {
	var p Poker
	p.Classic()
	p.bp.Wheel = wheel

	l, s := p.Likelihood(h)
	fmt.Println()
	fmt.Print(p.Stringify(h))
	if wheel {
		fmt.Print("  -  wheel")
	}
	fmt.Println()
	for j := range l {
		i := 9 - j
		c := uint64(l[i])
		g := GCD(c, uint64(s))
		n := p.bp.Names[i]
		p := float64(c) / float64(s)
		fmt.Printf("%-20s", n)
		if p > 0 {
			fmt.Printf("  %12.9f%%  =  %d / %d", 100*p, c/g, uint64(s)/g)
		}
		fmt.Println()
	}
}

var TriFoilList = []uint16{15, 31, 47, 200, 206, 319, 575, 831}

type trifoil struct {
	hand []string
	like []int
}

func Trifoil() (list []int) {
	var p BitPoker
	p.Classic()
	p.Wheel = false

	tri := map[int]trifoil{}
	var t trifoil
	var s sampler
	s.Init(p.pack, 3)
	for !s.Eof() {
		h := s.Next()
		l, _ := p.Likelihood(h)
		c := 0
		for i := 9; i >= 0; i-- {
			c += c
			if l[i] > 0 {
				c++
			}
		}
		t.hand = p.Humanize(h)
		t.like = l
		tri[c] = t
	}
	for i := range tri {
		list = append(list, i)
	}
	var r LCPRNG
	r.Sort(&list)
	return
}
