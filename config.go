package main

// Author: Srbislav D. Nešić, srbislav.nesic@fincore.com

import (
	"encoding/json"
	"fmt"
	"os"
)

type CardVirtue struct {
	Win   float64 `json:"win"`        // win amount
	Bonus float64 `json:"bonus"`      // win bonus
	Mult  float64 `json:"multiplier"` // win multiplier
	Free  bool    `json:"freepick"`   // free-pick
	Adv   int     `json:"advance"`    // level advance
	Level int     `json:"level"`      // level advance
	Val   int     `json:"-"`          // value
	Idx   int     `json:"-"`          // index
}

type CardType struct {
	Kind   string     `json:"kind"`             // card kind
	Suit   string     `json:"suit,omitempty"`   // card suit
	Virtue CardVirtue `json:"virtue,omitempty"` // card virtue
}

type CardsSetup struct {
	Decks [][]int    `json:"decks"`
	Cards []CardType `json:"cards,omitempty"`
}

type CardsDeal struct {
	HandLength       int      `json:"handLength"`               // hand length
	DefaultDeck      int      `json:"defaultDeck,omitempty"`    // default deck index
	MinPick          int      `json:"minCardsPicked,omitempty"` // min # of picked cards by player
	MaxPick          int      `json:"maxCardsPicked,omitempty"` // max # of picked cards by player
	AutPickoEnabled  bool     `json:"autoPickEnabled,omitempty"`
	FeaturesPreSpin  []string `json:"featuresPreSpin,omitempty"`
	FeaturesPostSpin []string `json:"featuresPostSpin,omitempty"`
}

type CardsConfiguration struct {
	Setup CardsSetup  `json:"setup"`
	Deal  []CardsDeal `json:"deals"`
}

type PlayFeature struct {
	Index      int     `json:"showDeckIndex"`   // show deck index
	History    int     `json:"handsHistory"`    // history length
	MaxMult    float64 `json:"maxProduct"`      // max mult length
	MaxMultLen int     `json:"maxMultiplier"`   // max mult length
	LevelMult  bool    `json:"levelMultiplier"` // rpyals level multiplier
	Collect    bool    `json:"pokerCollect"`    // poker collect
	Jokers     int     `json:"pokerJokers"`     // max number of wilds in poker
}

type PokerFeature struct {
	Kinds []string           `json:"kinds"` // cards kinds
	Suits []string           `json:"suits"` // cards suits
	Hands []PokerHandFeature `json:"hands"` // hands features
}

type PokerHandFeature struct {
	Name string  `json:"name"` // hand name
	Win  float64 `json:"win"`  // hand win
}

type Features struct {
	Play  PlayFeature  `json:"piqumplay"`
	Poker PokerFeature `json:"piqumpoker"`
}

type GameMeta struct {
	Domain string `json:"domain"`
}

type GameChecksum struct {
	Strategy string `json:"strategy"`
	Config   any    `json:"strategyConfiguratio"`
}

type GameConfiguration struct {
	GameId   string             `json:"gameId"`
	Config   CardsConfiguration `json:"cardsConfiguration"`
	Feature  Features           `json:"feature"`
	Meta     GameMeta           `json:"metamorphic"`
	Checksum GameChecksum       `json:"checksumConfiguration"`
	Prces    []float64          `json:"ticketPrice"`
}

var JSON GameConfiguration

var baselevel, freelevel int

var games = "games"

func LoadConfig(game string) {
	// filename := "./games/" + game + ".json"
	filename := "./" + games + "/" + game + ".json"
	data, err := os.ReadFile(filename)
	if err != nil {
		fmt.Println("No gane.", game)
		os.Exit(1)
	}
	json.Unmarshal(data, &JSON)
	baselevel = JSON.Config.Deal[0].DefaultDeck  // first deck
	freelevel = len(JSON.Config.Setup.Decks) - 1 // last deck
	if JSON.GameId != game {
		fmt.Println("Invalid game name.")
	}
	if JSON.Meta.Domain != game {
		fmt.Println("Invalid meta domain.")
	}
	if JSON.Checksum.Strategy == "" {
		fmt.Println("Invalid checksum.")
	}
}
