package main

import (
	"image/color"
	"math/rand"
)

func Clip(n int, m int, M int) int {
	if n < m {
		n = m
	}
	if n > M {
		n = M
	}
	return n
}

func sqDiffUInt8(x, y uint8) uint64 {
	d := uint64(x) - uint64(y)
	return d * d
}

func RandomSampleGene(g []Gene, n int) []Gene {
	idx := rand.Perm(n)
	sample := make([]Gene, n)
	for i, id := range idx {
		sample[i] = g[id]
	}
	return sample
}
func RemoveGene(s []Gene, i int) []Gene {
	s[i] = s[len(s)-1]
	// We do not need to put s[i] at the end, as it will be discarded anyway
	return s[:len(s)-1]
}

func DeepCopy(g []Gene) []Gene {
	var tmpCenter *[2]int
	var tmpColor *color.RGBA
	genome := make([]Gene, len(g))
	for i, gene := range g {
		genome[i].Radius = gene.Radius
		tmpCenter = new([2]int)
		*tmpCenter = *&gene.Center
		tmpColor = new(color.RGBA)
		*tmpColor = *&gene.Color

		genome[i].Center = *tmpCenter
		genome[i].Color = *tmpColor
	}
	return genome
}
