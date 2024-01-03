package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
)

type SNP struct {
	VCFFormat

	// The maximum number of each markers is four alleles.
	// The four elements of array represents the coverage of each alleles corresponding to ["A", "T", "C", "G"].
	Alleles [4]uint64
}

// BASE is a kind of nucleotide in DNA (A, T, G, C).
type BASE = string

var SortedBASE = [4]BASE{"A", "T", "C", "G"}

// Type assertion at compile time, to check SNP implements GeneticMarker interface.
var _ GeneticMarker = (*SNP)(nil)

func (snp SNP) GetPOS() uint64 {
	return snp.POS
}
func (snp SNP) GetCHROM() string {
	return snp.CHROM
}
func (snp SNP) GetAlleles() []uint64 {
	return snp.Alleles[:]
}

func (snp SNP) VerboseString() string {
	return fmt.Sprintf("%s\t%d\t%d\t%d\t%d",
		snp.ID, snp.Alleles[0], snp.Alleles[1], snp.Alleles[2], snp.Alleles[3])
}

func (snp SNP) String() string {
	var genotype = snp.DetermineGenotype()
	var genotypeSlice = genotype[:]
	// There are four compounds that make of DNA (Adenine, Guanine, Cytosine, Thymine).
	sort.Strings(genotypeSlice)
	return fmt.Sprintf("%s\t%s\t%s", snp.ID, genotypeSlice[0], genotypeSlice[1])
}

// DetermineGenotype removes less than three percent of BASE from four possibility.
func (snp SNP) DetermineGenotype() [2]BASE {
	var count = snp.Alleles[0] + snp.Alleles[1] + snp.Alleles[2] + snp.Alleles[3]
	var genotype []BASE
	for i, base := range SortedBASE {
		if float64(snp.Alleles[i])/float64(count) > *min_freq {
			genotype = append(genotype, base)
		}
	}

	// According to the number of retained alleles, homozygous, heterozygous or failure genotype would be determined.
	switch len(genotype) {
	case 1:
		return [2]BASE{genotype[0], genotype[0]}
	case 2:
		return [2]BASE{genotype[0], genotype[1]}
	default: // Void or three and more.
		return [2]BASE{}
	}
}

func ExtractSNP(file *os.File, marker *SNP) {
	//Debug #1: reset pointer offset. (2023-09-13)
	_, err := file.Seek(0, io.SeekStart)
	check(err)

	//fmt.Println(*marker)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		SAMrecord := NewSAM(scanner.Text())
		if SAMrecord == nil {
			continue
		}
		allele := SAMrecord.TypingSNP(*marker)
		if allele == "N" {
			continue
		}
		switch allele {
		case "A":
			marker.Alleles[0]++
		case "T":
			marker.Alleles[1]++
		case "C":
			marker.Alleles[2]++
		case "G":
			marker.Alleles[3]++
		}
	}
}
