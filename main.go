package main

import (
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"math/rand"
	"os"
	"sync"
	"time"
)

type Gene struct {
	Center [2]int
	Radius int
	Color  color.RGBA
}

type Result struct {
	fitness float64
	genome  []Gene
	out     image.Image
}

var (
	initialGenes int
	population   int
	mutationProb float64
	addProb      float64
	removeProb   float64

	minRadius int
	maxRadius int

	saveIter int

	width  int
	height int

	best Result

	img image.Image

	wg sync.WaitGroup

	genome []Gene
	out    image.Image
)

func (g Gene) init() Gene {
	g.Radius = rand.Intn(maxRadius-minRadius+1) + minRadius
	g.Center = [2]int{rand.Intn(width + 1), rand.Intn(height + 1)}
	g.Color.R, g.Color.G, g.Color.B, g.Color.A = uint8(rand.Intn(256)), uint8(rand.Intn(256)), uint8(rand.Intn(256)), 255
	return g
}

func (g Gene) MutateRadius(mutationSize float64) int {
	max := int(math.Floor(float64(g.Radius) * (1 + mutationSize)))
	min := int(math.Ceil(float64(g.Radius) * (1 - mutationSize)))
	if max-min == 0 {
		g.Radius = min
	} else {
		g.Radius = rand.Intn(max-min+1) + min
	}
	g.Radius = Clip(g.Radius, 1, 100)
	return g.Radius
}

func (g Gene) MutateCenter(mutationSize float64) [2]int {
	min0, max0 := int(math.Ceil(float64(g.Center[0])*(1-mutationSize))), int(math.Floor(float64(g.Center[0])*(1+mutationSize)))
	min1, max1 := int(math.Ceil(float64(g.Center[1])*(1-mutationSize))), int(math.Floor(float64(g.Center[1])*(1+mutationSize)))

	if max0-min0 == 0 {
		g.Center[0] = min0
	} else {
		g.Center[0] = rand.Intn(max0-min0+1) + min0
	}

	if max1-min1 == 0 {
		g.Center[1] = min1
	} else {
		g.Center[1] = rand.Intn(max1-min1+1) + min1
	}

	g.Center = [2]int{
		Clip(g.Center[0], 0, width),
		Clip(g.Center[1], 0, height)}

	return g.Center
}

func (g Gene) MutateColor(mutationSize float64) color.RGBA {
	cR := int(g.Color.R)
	cG := int(g.Color.G)
	cB := int(g.Color.B)
	minR, maxR := int(math.Ceil(float64(cR)*(1-mutationSize))), int(math.Floor(float64(cR)*(1+mutationSize)))
	minG, maxG := int(math.Ceil(float64(cG)*(1-mutationSize))), int(math.Floor(float64(cG)*(1+mutationSize)))
	minB, maxB := int(math.Ceil(float64(cB)*(1-mutationSize))), int(math.Floor(float64(cB)*(1+mutationSize)))

	if maxR-minR == 0 {
		g.Color.R = uint8(minR)
	} else {
		g.Color.R = uint8(rand.Intn(maxR-minR+1) + minR)
	}

	if maxG-minG == 0 {
		g.Color.G = uint8(minG)
	} else {
		g.Color.G = uint8(rand.Intn(maxG-minG+1) + minG)
	}

	if maxB-minB == 0 {
		g.Color.B = uint8(minB)
	} else {
		g.Color.B = uint8(rand.Intn(maxB-minB+1) + minB)
	}

	g.Color.R = uint8(Clip(int(g.Color.R), 0, 255))
	g.Color.G = uint8(Clip(int(g.Color.G), 0, 255))
	g.Color.B = uint8(Clip(int(g.Color.B), 0, 255))

	return g.Color
}

func (g Gene) Mutate() Gene {
	mutationSize := float64(math.Max(1, math.Round(rand.NormFloat64()*4+15))) / 100

	r := rand.Float64()

	if r < 0.33 {
		g.Radius = g.MutateRadius(mutationSize)
	} else if r < 0.66 {
		g.Center = g.MutateCenter(mutationSize)
	} else {
		g.Color = g.MutateColor(mutationSize)
	}
	return g
}

func ComputeFitness(genome []Gene) (float64, image.Image) {
	var fitness float64
	out := image.NewRGBA(image.Rect(0, 0, width, height))

	UnTransparent(out)

	for _, gene := range genome {
		DrawCircle(out, gene.Center[0], gene.Center[1], gene.Radius, gene.Color)
	}

	if img, ok := img.(*image.RGBA); ok {
		compVal, _ := CompareImage(img, out)
		fitness = 255 / float64(compVal)
	}
	return fitness, out
}

func ComputePopulation(gnm []Gene) Result {
	gene := Gene{}
	gene = gene.init()

	//genome := gnm

	genome := make([]Gene, len(gnm))

	for i, gn := range gnm {
		r := gn.Color.R
		g := gn.Color.G
		b := gn.Color.B

		cx := gn.Center[0]
		cy := gn.Center[1]
		rad := gn.Radius

		genome[i].Color = color.RGBA{r, g, b, 255}
		genome[i].Center = [2]int{cx, cy}
		genome[i].Radius = rad
	}

	if len(genome) < 200 {
		for _, g := range genome {
			if rand.Float64() < mutationProb {
				g = g.Mutate()
			}
		}
	} else {
		mut := RandomSampleGene(genome, int(float64(len(genome))*mutationProb))
		for _, g := range mut {
			g = g.Mutate()
		}
	}

	if rand.Float64() < addProb {
		genome = append(genome, gene)
	}

	if len(genome) > 0 && rand.Float64() < removeProb {
		genome = RemoveGene(genome, rand.Intn(len(genome)))
	}

	fitness, out := ComputeFitness(genome)

	result := Result{fitness, genome, out}
	return result
}

func worker(genomes chan []Gene, results chan<- Result) {
	defer wg.Done()
	for i := 0; i < 10; i++ {
		result := ComputePopulation(<-genomes)
		results <- result
		genomes <- result.genome
	}
}

func findBest(results chan Result) {
	max := (<-results).fitness
	for i := 0; i < population-1; i++ {
		result := <-results
		if result.fitness > max {
			max = result.fitness
			best = result
		}
	}
}

func main() {
	rand.Seed(time.Now().UTC().UnixNano())

	initialGenes = 50
	population = 50
	mutationProb = 0.05
	addProb = 0.3
	removeProb = 0.2

	minRadius = 5
	maxRadius = 15

	//imgPath := os.Args[1]
	imgPath := "test.png"
	imgF, _ := os.Open(imgPath)
	defer imgF.Close()
	img, _, _ = image.Decode(imgF)
	b := img.Bounds()

	width = b.Max.X
	height = b.Max.Y
	gen := 0
	workers := 5

	os.Mkdir("result", 777)

	genome = make([]Gene, initialGenes)

	for i, gene := range genome {
		genome[i] = gene.init()
	}

	_, out = ComputeFitness(genome)

	for gen < 1000 {

		genomes := make(chan []Gene, population)
		results := make(chan Result, population)

		for i := 0; i < population; i++ {
			j := genome
			genomes <- j
		}

		for w := 0; w < workers; w++ {
			wg.Add(1)
			go worker(genomes, results)
		}

		wg.Wait()

		findBest(results)

		genome = best.genome
		out = best.out

		close(results)
		close(genomes)

		gen++
		fmt.Printf("Currently on generation %d, fitness %.6f\n", gen, best.fitness)

		if gen%100 == 0 {
			genoutput, _ := os.Create(fmt.Sprintf("result/gen%d.png", gen))
			png.Encode(genoutput, out)
		}

	}
	output, _ := os.Create("testsave.png")
	err := png.Encode(output, out)

	if err != nil {
		fmt.Println(err)
	}
}
