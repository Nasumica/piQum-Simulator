package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"strconv"
	"time"
)

// Author: Srbislav D. Nešić, srbislav.nesic@fincore.com

const (
	million         = 1000 * 1000        // million ^ 1
	billion         = million * million  // million ^ 2
	trillion        = million * billion  // million ^ 3
	quadrillion     = million * trillion // million ^ 4
	milliard        = 1000 * million
	billiard        = 1000 * billion
	trilliard       = 1000 * trillion
	quadrilliard    = 1000 * quadrillion
	usa_billion     = milliard
	usa_trillion    = billion
	usa_quadrillion = billiard
)

var result = "result"
var max_mult float64

func ShowDecks() {
	var rnd LCPRNG
	for n, deck := range JSON.Config.Setup.Decks {
		rnd.Sort(&deck)
		p := poker.bp.Squeeze(deck...) ^ (1<<52 - 1)
		s := poker.bp.Length(p)
		j := 0
		row := ""
		row += fmt.Sprintf("deck[%d] = ", n)
		row += "["
		d := ""
		var y StatCalc
		for _, n := range deck {
			c := JSON.Config.Setup.Cards[n]
			if max_mult < c.Virtue.Mult {
				max_mult = c.Virtue.Mult
			}
			y.Add(c.Virtue.Bonus)
			if n == 0 {
				j++
			} else if n > 52 {
				x := JSON.Config.Setup.Cards[n].Kind
				row += fmt.Sprintf("%s%s", d, x)
				d = ", "
			}
		}
		row += "]"
		if p == 0 && s == 0 {
			row += " + 52 standard"
		}
		if j > 0 {
			row += fmt.Sprintf(" + %d jokers", j)
		}
		row += fmt.Sprintf(" = %d", len(deck))
		row = fmt.Sprintf("%-80s", row)
		if y.Avg > 0 {
			row += fmt.Sprintf("%10.2f ± %.2f", y.Avg, y.Dev)
		}
		fmt.Println(row)
	}
}

func CalcProb(iter, hands int) {
	deck := JSON.Config.Setup.Decks[0]
	cards := len(deck)
	jokers := 0
	ups := 0
	mults := 0
	for _, n := range deck {
		c := JSON.Config.Setup.Cards[n]
		if n == 0 {
			jokers++
		}
		if c.Virtue.Adv == 1 {
			ups++
		}
		if c.Virtue.Mult > 0 {
			mults++
		}
	}
	// n := float64(iter) * float64(hands)
	p := HypGeomDist(5, 5, jokers, cards)
	q := 1 - p
	pjokers := 1 - math.Pow(q, float64(hands))
	pbonus := float64(ups) / float64(cards) * 5 / 3
	bonus := int(math.Round(float64(iter) * pbonus))
	base := iter
	freespin := int(math.Round(float64(base) * pjokers))
	fmt.Println(base, bonus, freespin)
}

func Cats(cat ...string) {
	for _, x := range cat {
		w := true
		if x != "" && x[:1] == ":" {
			w = false
			x = x[1:]
		}
		switch x {
		case "":
			fmt.Println()
		case "*":
			fmt.Println("category                 sum      rtp       count         rate     min     max           μ ± σ")
		case "-":
			fmt.Println("—————————————————————————————————————————————————————————————————————————————————————————————————————")
		case "=":
			fmt.Println("═════════════════════════════════════════════════════════════════════════════════════════════════════")
		default:
			var c StatCalc
			var e bool
			if i, err := strconv.Atoi(x); err == nil {
				c = Data.Lev[i]
				e = true
			} else {
				c, e = Data.Cat[x]
				e = e || w
			}
			if e {
				fmt.Printf("%-16s", x)
				if c.Sum > 0 {
					fmt.Printf("  %10.0f", c.Sum)
				} else {
					fmt.Printf("%12s", "")
				}
				if c.Sum > 0 && w {
					rtp := 100 * c.Sum / Data.Bet
					fmt.Printf("  %6.2f%%", rtp)
				} else {
					fmt.Printf("%9s", "")
				}
				// rate := float64(Data.Games) / float64(c.Cnt)
				rate := float64(Data.Lev[0].Cnt) / float64(c.Cnt)
				fmt.Printf("  %10d  %11.2f", c.Cnt, rate)
				if c.Max > 0 {
					fmt.Printf("  %6.0f  %6.0f", c.Min, c.Max)
					if c.Min != c.Max {
						fmt.Printf("  %10.2f ± %.2f", c.Avg, c.Dev)
					}
				}
				fmt.Println()
			}
		}
	}
}

func ShowResult() {
	exist := func(cat string) bool {
		_, e := Data.Cat[cat]
		return e
	}

	Cats("=", "*", "-")
	fmt.Printf("bet %24.0f %20d\n", Data.Bet, Data.Games)

	if exist("royals") {
		Cats("-", "J", "Q", "K", "A", "-", "royals")
	}

	if exist("bonus") {
		Cats("-", "mini", "minor", "major", "grand", "mega", "-", "bonus")
	}

	if exist("hands") {
		Cats("-")
		for _, h := range JSON.Feature.Poker.Hands {
			if h.Win > 0 {
				Cats(h.Name)
			}
		}
		// Cats("-", "poker", "hands")
		Cats("-", "poker")
	}

	if exist("serva") {
		f := false
		for _, h := range poker.bp.Names {
			if exist(h) {
				if !f {
					Cats("-")
					f = true
				}
				Cats(h)
			}
		}
		Cats("-", "serva")
	}

	Cats("-", "free", "-", "total", "=")

	for i := range JSON.Config.Setup.Decks {
		Cats(strconv.Itoa(i))
	}
	// Cats("-", "summary", "-")
	Cats("-")

	// Cats(":free", ":jokers", ":mult")
	Cats(":mult", ":product", ":width", ":walk", ":step", ":collect")

	Cats("=")
}

func LoadResult(game string) {
	// filename := "./result/" + game + "-result.json"
	filename := "./" + result + "/" + game + "-result.json"
	data, _ := os.ReadFile(filename)
	json.Unmarshal(data, &Data)
}

func SaveResult(game string) {
	// filename := "./result/" + game + "-result.json"
	filename := "./" + result + "/" + game + "-result.json"
	data, _ := json.MarshalIndent(Data, "", "\t")
	os.WriteFile(filename, data, 0644)
}

func Reset() {
	d := len(JSON.Config.Setup.Decks)
	c := len(JSON.Config.Setup.Cards)
	Data.Meta.Level = JSON.Config.Deal[0].DefaultDeck
	Data.Meta.Poker = []int{}
	Data.Meta.Virtue = make([]MetaVirtue, JSON.Config.Deal[0].HandLength)
	Data.Meta.Walk = 0
	Data.Cat = make(map[string]StatCalc)
	Data.Lev = make([]StatCalc, d)
	Data.Pick = make([][]int, d)
	for i := range Data.Pick {
		Data.Pick[i] = make([]int, c)
	}
	Data.Bet = 0
	Data.Win = 0
	Data.Deals = 0
	Data.Games = 0
	Data.Seed = int(RND.Seed())
}

func strategy() {
	type str struct {
		n string
		v float64
	}
	var s []str
	s = append(s, str{"none", 0})
	s = append(s, str{"fp", fp_weight})
	s = append(s, str{"2x", 2})
	s = append(s, str{"3x", 3})
	s = append(s, str{"4x", 4})
	s = append(s, str{"5x", 5})
	for i := 0; i < len(s)-1; i++ {
		for j := i + 1; j < len(s); j++ {
			if s[i].v > s[j].v {
				s[i], s[j] = s[j], s[i]
			}
		}
	}
	d := "strategy:  "
	for _, x := range s {
		fmt.Printf("%s%s", d, x.n)
		d = " < "
	}
	fmt.Println()
}

var bet float64 = 1

func CSV() (rows []string) {
	sum, _ := os.Create("xls_summary.csv")
	defer sum.Close()
	det, _ := os.Create("xls_details.csv")
	defer det.Close()

	b := Data.Cat["bet"]
	w := Data.Cat["total"]

	rtp := 100 * w.Sum / b.Sum
	pp := b.Avg
	base, bonus, free := 0, 0, 0
	for i, s := range Data.Lev {
		if i == 0 {
			base += s.Cnt
		} else if i+1 == len(Data.Lev) {
			free += s.Cnt
		} else {
			bonus += s.Cnt
		}
	}

	sum.WriteString("\"tickets\"\t\"bet\"\t\"win\"\t\"rtp\"\t\"base\"\t\"bonus\"\t\"free\"\t\"pp\"\n")
	sum.WriteString(fmt.Sprintf("\"%d\"\t", b.Cnt))
	sum.WriteString(fmt.Sprintf("\"%.0f\"\t", b.Sum))
	sum.WriteString(fmt.Sprintf("\"%.0f\"\t", w.Sum))
	sum.WriteString(fmt.Sprintf("\"%.2f\"\t", rtp))
	sum.WriteString(fmt.Sprintf("\"%d\"\t", base))
	sum.WriteString(fmt.Sprintf("\"%d\"\t", bonus))
	sum.WriteString(fmt.Sprintf("\"%d\"\t", free))
	sum.WriteString(fmt.Sprintf("\"%.2f\"\n", pp))

	rows = append(rows, "\"cat\"\t\"count\"\t\"sum\"\t\"min\"\t\"max\"\t\"sqr\"\t\"nul\"")

	one := func(h string, s StatCalc) {
		r := fmt.Sprintf("\"%s\"", h)
		r += fmt.Sprintf("\t\"%d\"", s.Cnt)
		r += fmt.Sprintf("\t\"%0.f\"", s.Sum)
		r += fmt.Sprintf("\t\"%0.f\"", s.Min)
		r += fmt.Sprintf("\t\"%0.f\"", s.Max)
		r += fmt.Sprintf("\t\"%0.f\"", s.Sqr)
		if s.Nul == 0 {
			r += "\t\"\""
		} else {
			r += fmt.Sprintf("\t\"%d\"", s.Nul)
		}
		rows = append(rows, r)
	}

	cat := func(head ...string) {
		for _, h := range head {
			s, e := Data.Cat[h]
			if e {
				one(h, s)
			}
		}
	}

	lev := func() {
		for i, s := range Data.Lev {
			h := fmt.Sprintf("%d", i)
			one(h, s)
		}
	}

	cat("mult")
	one("fp", Data.Cat["free"])
	cat("J", "Q", "K", "A")
	cat("mini", "minor", "major", "grand", "mega")

	for _, h := range JSON.Feature.Poker.Hands {
		if h.Win > 0 {
			cat(h.Name)
		}
	}

	lev()

	for _, r := range rows {
		det.WriteString(r + "\n")
	}

	return
}

func test(game string, iter int) {
	LoadConfig(game)
	Reset()

	LoadResult(game)

	RND.Randomize(uint64(Data.Seed))

	if iter < 0 {
		iter = -iter - Data.Lev[0].Cnt
	}

	start := time.Now()
	this.bet, this.games, this.deals, this.base = 0, 0, 0, 0
	// for this.games < iter {
	for this.base < iter {
		Play(bet)
	}
	for Data.Meta.Level != 0 {
		Play(bet)
	}
	elapsed := time.Since(start).Seconds()
	Data.Seed = int(RND.Seed())

	if iter > 0 {
		SaveResult(game)
	}

	rtp := 100 * Data.Cat["total"].Sum / Data.Bet
	fmt.Printf("%s  (%.2f%%)\n", game, rtp)
	strategy()
	fmt.Println()

	ShowDecks()

	fmt.Println()
	fmt.Println("price points:    ", JSON.Prces)

	fmt.Printf("max multipliers:  ")
	if JSON.Feature.Play.MaxMultLen > 0 {
		m := math.Pow(max_mult, float64(JSON.Feature.Play.MaxMultLen))
		if JSON.Feature.Play.MaxMult <= 0 || m < JSON.Feature.Play.MaxMult {
			JSON.Feature.Play.MaxMult = m
		}
		fmt.Printf("%d\n", JSON.Feature.Play.MaxMultLen)
	} else {
		fmt.Printf("no limit\n")
	}
	fmt.Printf("max product:      ")
	if JSON.Feature.Play.MaxMult > 0 {
		fmt.Printf("%.0f\n", JSON.Feature.Play.MaxMult)
	} else {
		fmt.Printf("no limit\n")
	}
	if JSON.Feature.Play.MaxMultLen > 0 {
		d, m, w := JSON.Config.Setup.Decks[0], 0, 0.
		for _, c := range d {
			if JSON.Config.Setup.Cards[c].Virtue.Mult > 0 {
				m++
			}
			if JSON.Config.Setup.Cards[c].Virtue.Win > w {
				w = JSON.Config.Setup.Cards[c].Virtue.Win
			}
		}
		w *= math.Pow(max_mult, float64(JSON.Feature.Play.MaxMultLen))
		pup, pdn := MaxWinProb(len(d), m, JSON.Feature.Play.MaxMultLen)
		p := float64(pup) / float64(pdn)
		// r := float64(pdn) / float64(pup)
		fmt.Printf("max odd:          %.0f  probability = %d / %d = %v\n", w, pup, pdn, p)
	}
	fmt.Println()

	if len(JSON.Feature.Poker.Hands) > 0 {
		rtp, rate, mu, sigma := CalcPokerRTP()
		fmt.Printf("expected poker rtp = %.2f%%   rate = %.2f   mean = %.2f ± %.2f\n", 100*rtp, rate, mu, sigma)
		fmt.Println()
	}

	ShowResult()

	fmt.Printf("(all rates are relative to the number of base games:  rate = %d / count)\n", Data.Lev[0].Cnt)

	if this.deals > 0 && elapsed > 0 {
		speed := float64(this.deals) / elapsed
		fmt.Printf("elapsed = %.3f\",  speed = %.0f deals / s\n", elapsed, speed)
	}

	// CSV()
}

func main() {
	games = "games/config"     // games folder
	result = "games/result"    // result folder
	skill = 100                // in percent (-100, 100)
	fp_weight = 2.5            // 1.5 or 2.5
	bet = 1                    // ticket bet
	deals := -1000 * million   // games to test
	game := "piqum-classic-95" // game config

	if len(os.Args) > 1 {
		game = os.Args[1]
		if len(os.Args) > 2 {
			if n, e := strconv.ParseInt(os.Args[2], 10, 64); e == nil {
				deals = int(n) * million
			}
		}
	}

	test(game, deals) // do it
}
