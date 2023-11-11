package main

import (
	"fmt"
	"testing"
)

func TestAdjustPos(t *testing.T) {
	fmt.Println(AdjustPos(80, "65M2D1I35M"))
}

func TestCountMatchSNPInMH(t *testing.T) {
	var pair = [][2]AlleleMH{
		{"A-T-A-G-T", "A-T-A-G-T"}, //Match == 5
		{"A-T-T-G-T", "A-T-A-G-T"}, // Never Match == -1
		{".-.-.-A-T", "A-T-A-G-T"}, // not overhang == -1
		{".-.-.-G-T", "A-T-A-G-T"}, // left overhang == 2
		{"A-T-A-.-.", "A-T-A-G-T"}, // right overhang == 3
	}
	for _, v := range pair {
		fmt.Printf("%s set to %d\n", v[0], CountMatchSNPInMH(v[0], v[1]))
	}
}
