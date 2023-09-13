package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"os/exec"
	"strconv"
	"strings"
)

type AlleleMH string

type MHMarker struct {
	ID      string
	Chr     string
	SNPs    []uint64
	Alleles map[AlleleMH]float64 // 单个个体每个基因型的统计深度

	Population map[string][2]AlleleMH //每个个体的基因型, 带"."的alleleMH都用单个"."表示
}

func (mh MHMarker) String() string {
	return ""
}
func (mh MHMarker) GetCHROM() string {
	return ""
}
func (mh MHMarker) GetPOS() uint64 {
	return 0
}

func (mh *MHMarker) IndividualGenotype(sample string) [2]AlleleMH {
	if diploid, ok := mh.Population[sample]; ok {
		return diploid
	} else {
		return [2]AlleleMH{".", "."}
	}
}

// 群体多态信息含量
func (mh *MHMarker) PIC() float64 {
	var (
		stat  = mh.allelePopulation()
		numP  = float64(len(mh.Population)) * 2
		total int
		freq  []float64
	)
	for _, count := range stat {
		freq = append(freq, float64(count)/numP)
		total += count
	}
	if total < 2*len(mh.Population) {
		freq = append(freq, numP-float64(total)/numP)
	}

	return pic(freq)
}
func pic(n []float64) float64 {
	// PIC = 1 - \sum_{i=1}^{n} p_i^2 - (\sum_{i=1}^{n} p_i^2)^2 + \sum_{i=1}^{n} p_i^4
	var sumOfSquares, sumOfBiquadratic float64
	for _, v := range n {
		sumOfSquares += math.Pow(v, 2)
		sumOfBiquadratic += math.Pow(v, 4)
	}
	return 1 - sumOfSquares - math.Pow(sumOfSquares, 2) + sumOfBiquadratic
}

// 群体等位基因数量统计
func (mh *MHMarker) allelePopulation() map[AlleleMH]int {
	var stat = make(map[AlleleMH]int)
	for _, diploid := range mh.Population {
		// 第一条染色体
		if _, ok := stat[diploid[0]]; !ok && diploid[0] != "." {
			stat[diploid[0]] = 1
		} else if ok && diploid[0] != "." {
			stat[diploid[0]] += 1
		}
		// 第二条染色体
		if _, ok := stat[diploid[1]]; !ok && diploid[1] != "." {
			stat[diploid[1]] = 1
		} else if ok && diploid[0] != "." {
			stat[diploid[1]] += 1
		}
	}

	return stat
}

// Ae refers to the effective number of alleles
func (mh *MHMarker) Ae() int {
	return len(mh.allelePopulation())
}

func (mh *MHMarker) appendSNP(n uint64) {
	mh.SNPs = append(mh.SNPs, n)
}

// Classify method decides which overlapping alleles belong to a complete allele.
// Each overlapping allele add 1/len(SNPs) to corresponding allele.
// Only Onw or Two alleles would be retained, other seldom alleles (frequency is lee than 10 %) would be abandoned.
func (mh *MHMarker) Classify() {
	// extract compete alleles
	var (
		alleles  []AlleleMH
		coverage float64
	)

	for k, _ := range mh.Alleles {
		if i := strings.Index(string(k), "."); i == -1 {
			alleles = append(alleles, AlleleMH(k))
		}
	}
	if len(alleles) == 0 {
		return
	}

	// classify
	len_SNP := len(alleles[0])/2 + 1
	for k, v := range mh.Alleles {
		coverage += v
		if i := strings.Index(string(k), "."); i != -1 {
			for _, allele := range alleles {
				if string(allele)[:i] == string(k)[:i] {
					mh.Alleles[allele] += v * (float64(i) / 2) * (1 / float64(len_SNP))
				}
			}
			delete(mh.Alleles, k)
		}
	}

	// abandon seldom alleles
	for k, v := range mh.Alleles {
		if v/coverage <= 0.1 {
			delete(mh.Alleles, k)
		}
	}
}

func (mh *MHMarker) SimpleString() string {
	var (
		s       string
		alleles []AlleleMH
	)

	for allele, _ := range mh.Alleles {
		alleles = append(alleles, allele)
	}

	switch len(alleles) {
	case 0:
		s = fmt.Sprintf("%s\t%s\t%s", mh.ID, ".", ".") // missing alignment
	case 1:
		s = fmt.Sprintf("%s\t%s\t%s", mh.ID, alleles[0], alleles[0])
	case 2:
		s = fmt.Sprintf("%s\t%s\t%s", mh.ID, alleles[0], alleles[1])
	default:
		s = fmt.Sprintf("#%s\t%v", mh.ID, mh.Alleles)
	}
	return s
}

func NewMHMarker(ID string, Chr string, SNPs []uint64) *MHMarker {
	return &MHMarker{ID: ID, Chr: Chr, SNPs: SNPs, Alleles: make(map[AlleleMH]float64)}
}

func NewMHMarkerCollection() map[string]MHMarker {
	var set = make(map[string]MHMarker, 0)
	return set
}

// input marker form file.
func NewMHMarkers(path string) []MHMarker {
	var (
		markers []MHMarker
		current *MHMarker
	)

	handle, err := os.Open(path)
	if err != nil {
		log.Fatalln(err)
	}
	scanner := bufio.NewScanner(handle)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) < 7 {
			continue
		}
		chr := fields[0]
		Pos, _ := strconv.ParseUint(fields[1], 10, 64)
		mhID := fields[6]
		if current == nil {
			current = NewMHMarker(mhID, chr, []uint64{Pos})
		} else if current.Chr != chr {
			markers = append(markers, *current)
			current = NewMHMarker(mhID, chr, []uint64{Pos})
		} else if current.ID == mhID {
			current.appendSNP(Pos)
		} else {
			markers = append(markers, *current)
			current = NewMHMarker(mhID, chr, []uint64{Pos})
		}
	}
	markers = append(markers, *current)
	return markers
}

var (
	SAM1 = "" // path
	MH1  = MHMarker{ID: "mhGP01", Chr: "Chr1", SNPs: []uint64{16574710, 16574718, 16574732}}
)

func ExtractMHAlleles(path string, marker *MHMarker) {
	handle, err := os.Open(path)
	if err != nil {
		log.Fatalln(err)
	}
	scanner := bufio.NewScanner(handle)
	for scanner.Scan() {
		SAMrecord := NewSAM(scanner.Text())
		if SAMrecord == nil {
			continue
		}
		allele := SAMrecord.TypingMH(*marker)
		if allele == "" {
			continue
		}
		// overlap
		if _, ok := marker.Alleles[allele]; ok {
			marker.Alleles[allele] += 1
		} else {
			marker.Alleles[allele] = 1
		}
	}
}

func CMDSamtools(bamfile string, mh *MHMarker) {
	cmd := exec.Command("samtools",
		"view",
		"-o", fmt.Sprintf("temp/%s.sam", mh.ID),
		bamfile,
		fmt.Sprintf("%s:%d-%d", mh.Chr, mh.SNPs[0]-100, mh.SNPs[len(mh.SNPs)-1]+100))
	//fmt.Println(cmd.Args, "okkkkk")
	cmd.Output()
}

func CMDRemoveTemp(mh *MHMarker) {
	cmd := exec.Command("rm",
		fmt.Sprintf("temp/%s.sam", mh.ID))
	cmd.Output()
}

/*
func main() {
	flag.Parse()
	if *SNP_path == "" || *BAM == "" {
		flag.Usage()
		os.Exit(1)
	}

	makers := NewMHMarkers(*SNP_path)
	for _, marker := range makers {
		CMDSamtools(*BAM, &marker)
		ExtractMHAlleles(fmt.Sprintf("temp/%s.sam", marker.ID), &marker)
		CMDRemoveTemp(&marker)
		marker.Classify()
		fmt.Println(marker.SimpleString())
	}
}
*/
