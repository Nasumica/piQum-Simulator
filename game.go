package main

// Author: Srbislav D. Nešić, srbislav.nesic@fincore.com

import (
	"math"
)

var skill = 100 // in percent

type MetaVirtue struct {
	Free  bool      `json:"f,omitempty"`
	Mult  float64   `json:"m,omitempty"`
	Times []float64 `json:"x,omitempty"`
}

type Metamorphic struct {
	Level   int          `json:"level"`
	Virtue  []MetaVirtue `json:"virtue"`
	Poker   []int        `json:"poker,omitempty"`
	Walk    int          `json:"walk,omitempty"`
	Collect int          `json:"collect,omitempty"`
}

type GameData struct {
	Seed  int                 `json:"seed"`
	Deals int                 `json:"deals"`
	Games int                 `json:"games"`
	Bet   float64             `json:"bet"`
	Win   float64             `json:"win"`
	Cat   map[string]StatCalc `json:"category"`
	Lev   []StatCalc          `json:"level"`
	Pick  [][]int             `json:"pick"`
	Wild  []int               `json:"wild,omitempty"`
	Meta  Metamorphic         `json:"meta"`
}

var Data GameData

var poker Poker

var serva = false

var fp_weight = 2.5

var this struct {
	bet   float64
	games int
	deals int
	base  int
}

func addlev(lev int, val float64) {
	Data.Lev[lev].Add(val)
}

func addcat(cat string, val float64) {
	s := Data.Cat[cat]
	s.Cat = cat
	s.Add(val)
	Data.Cat[cat] = s
}

func Deal() (hand []CardType, jokers int) {
	this.deals++
	Data.Deals++
	hand = make([]CardType, JSON.Config.Deal[0].HandLength)
	jokers = 0
	d := &JSON.Config.Setup.Decks[Data.Meta.Level]
	l := len(*d) - 1
	for i := range hand {
		j := RND.Int(i, l)
		(*d)[i], (*d)[j] = (*d)[j], (*d)[i]
		n := (*d)[i]
		if n == 0 {
			jokers++
		}
		c := JSON.Config.Setup.Cards[n]
		c.Virtue.Val = n
		c.Virtue.Idx = i
		hand[i] = c
	}
	return
}

func QuickPick() (qp int) {
	qp = -1
	if Data.Meta.Level == baselevel && skill != 0 { // base level only
		list := []int{}
		if skill > 0 {
			if RND.Choose(100, skill) {
				best := math.SmallestNonzeroFloat64
				for i, v := range Data.Meta.Virtue {
					if v.Free {
						v.Mult = fp_weight // free-pick mult weight
					}
					if v.Mult >= best {
						if v.Mult > best {
							best, list = v.Mult, []int{} // better choice
						}
						list = append(list, i) // best choices list
					}
				}
			}
		} else {
			if RND.Choose(100, -skill) {
				worst := math.MaxFloat64
				for i, v := range Data.Meta.Virtue {
					if v.Free {
						v.Mult = 1 // free-pick mult weight
					}
					if v.Mult <= worst {
						if v.Mult < worst {
							worst, list = v.Mult, []int{} // worse choice
						}
						list = append(list, i) // worst choices list
						if worst == 0 {
							break
						}
					}
				}
			}
		}
		qp = RND.Item(&list) // skill pick
	}
	if qp < 0 {
		qp = RND.Choice(JSON.Config.Deal[0].HandLength) // random pick
	}
	return
}

func Play(bet float64) {
	level := Data.Meta.Level

	if level == 0 {
		this.base++
	}

	totalwin := 0.

	addwin := func(cat string, win float64) {
		addcat(cat, win)
		totalwin += win
	}

	cards, jokers := Deal()
	pick := QuickPick()
	card := cards[pick]
	virt := Data.Meta.Virtue[pick]
	Data.Pick[level][card.Virtue.Val]++

	if level == freelevel {
		addcat("jokers", 0)
	} else {
		this.games++
		Data.Games++
		if level == baselevel && virt.Free {
			// addcat("free", 0)
			this.bet += bet
			Data.Bet += bet
			addcat("bet", bet)
			addwin("free", bet)
		} else {
			this.bet += bet
			Data.Bet += bet
			addcat("bet", bet)
		}
		if level != baselevel {
			addcat("step", 0)
		}
	}

	if win := card.Virtue.Win; win > 0 { // royal win
		mult := math.Max(1, virt.Mult)
		addcat("mult", mult)
		win *= mult
		win *= bet
		addwin(card.Kind, win)
		addcat("royals", win)
		addcat("width", float64(len(virt.Times)))
		if virt.Mult > 0 {
			addcat("product", virt.Mult)
		}
	}

	if serva && len(JSON.Feature.Poker.Hands) == 0 && level == 0 {
		var hold spil
		for _, c := range cards {
			if v := c.Virtue.Val; 1 <= v && v <= 52 {
				hold |= 1 << (v - 1)
			}
		}
		if h := poker.bp.Hand(hold); h != "" {
			c := poker.HandCode(h)
			// w := [10]float64{0, 1, 2, 3, 5, 10, 20, 30, 50, 100}[c]
			w := [10]float64{0, 1, 2, 3, 10, 20, 30, 100, 200, 300}[c]
			if w > 0 {
				h = h[9:]
				w *= bet
				addwin(h, w)
				addcat("serva", w)
			}
		}
	}

	if win := card.Virtue.Bonus; win > 0 { // bonus win
		win *= bet
		addwin(card.Kind, win)
		addcat("bonus", win)
		addcat("walk", float64(Data.Meta.Walk))
		Data.Meta.Walk = 0
	}

	if level == 0 && JSON.Feature.Play.Collect { //poker hand
		if len(Data.Meta.Poker) >= 5 {
			Data.Meta.Poker = []int{} // reset poker hand
		}
		Data.Meta.Collect++
		draw := card.Virtue.Val
		if 0 <= draw && draw <= 52 { // joker and standard cards only
			hold := poker.bp.Squeeze(draw)
			pile := poker.bp.Squeeze(Data.Meta.Poker...)
			new := (draw > 0 && hold&pile == 0) || (draw == 0 && len(Data.Meta.Poker)+JSON.Feature.Play.Jokers >= 5)
			// new := (draw > 0 && hold&pile == 0)
			if new {
				Data.Meta.Poker = append(Data.Meta.Poker, draw)
				if len(Data.Meta.Poker) == 5 {
					addcat("collect", float64(Data.Meta.Collect))
					Data.Meta.Collect = 0
					best, wild := poker.JokerWild(Data.Meta.Poker) // get best hand and wild replacements list
					code := poker.HandCode(best)                   // extract hand code (0-9)
					hand := JSON.Feature.Poker.Hands[code]         // hand from feature
					win := hand.Win
					win *= bet
					addwin(hand.Name, win)
					if win > 0 {
						addcat("poker", win)
					}
					addcat("hands", win)
					if l := len(wild); l > 0 {
						if len(Data.Wild) == 0 {
							Data.Wild = make([]int, 53)
						}
						Data.Wild[0] += l
						for _, w := range wild {
							Data.Wild[w]++
						}
					}
				}
			}
		}
	}

	switch level { // calculate virtues
	case baselevel:
		for i, c := range cards {
			v := Data.Meta.Virtue[i]
			v.Free = c.Virtue.Free
			if c.Virtue.Mult == 0 {
				v.Times = []float64{}
				v.Mult = 0
			} else {
				l := len(v.Times)
				if l == 0 {
					v.Mult = 1
				}
				if JSON.Feature.Play.MaxMultLen == 0 || l < JSON.Feature.Play.MaxMultLen {
					v.Mult *= c.Virtue.Mult
				}
				if JSON.Feature.Play.MaxMult > 0 && v.Mult > JSON.Feature.Play.MaxMult {
					v.Mult = JSON.Feature.Play.MaxMult
				}
				v.Times = append(v.Times, c.Virtue.Mult)
			}
			Data.Meta.Virtue[i] = v
		}
	case freelevel:
		Data.Meta.Virtue = make([]MetaVirtue, JSON.Config.Deal[0].HandLength)
		v := Data.Meta.Virtue[pick]
		v.Free = card.Virtue.Free
		if card.Virtue.Mult > 0 {
			v.Times = append(v.Times, card.Virtue.Mult)
			v.Mult = card.Virtue.Mult
		}
		Data.Meta.Virtue[pick] = v
	}

	switch { // determine next level
	case jokers == JSON.Config.Deal[0].HandLength:
		Data.Meta.Level = freelevel
		// Data.Meta.Level = baselevel
	case card.Virtue.Adv == 0:
		Data.Meta.Level = baselevel
	default:
		Data.Meta.Level += card.Virtue.Adv
		Data.Meta.Walk++
	}

	if totalwin > 0 {
		addcat("total", totalwin)
	}
	addlev(level, totalwin)
	addcat("summary", totalwin)
	Data.Win += totalwin
}

func PokerCRate(deck, jokers int) float64 {
	d, j := float64(deck), float64(jokers)
	return d*257449/3248700 + d/(48+j)
}

func PokerWRate(deck, jokers int) float64 {
	d, j := float64(deck), float64(jokers)
	return (d * (15606252 + 257449*j)) / (12 * (990900 + 100877*j))
}

func PokerRtp(deck, jokers int) float64 {
	d, j := float64(deck), float64(jokers)
	return (24 * (3116500 + 638117*j)) / (d * (15606252 + 257449*j))
}

func CalcPokerRTP() (rtp, rate, μ, σ float64) {
	deck := JSON.Config.Setup.Decks[0]
	card := WildCardComb[0]
	wild := WildCardComb[1]

	jokers := 0
	for _, c := range deck {
		if c == 0 {
			jokers++
		}
	}

	var population, success, prize, squares float64
	for i, h := range JSON.Feature.Poker.Hands {
		c := float64(5*card[i] + jokers*wild[i])
		population += c
		if w := h.Win; w > 0 {
			success += c
			c *= w
			prize += c
			squares += c * w
		}
	}

	if success > 0 { // prob := success / population
		d, j := float64(len(deck)), float64(jokers)
		rate = (d/52 + d/51 + d/50 + d/49 + d/(48+j)) * population / success
		μ = prize / success
		σ = math.Sqrt(squares/success - μ*μ)
		rtp = μ / rate
	}

	return
}

func init() {
	poker.Classic()
}
