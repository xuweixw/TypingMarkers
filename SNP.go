package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
)

type SNP struct {
	VCFFormat

	// The maximum number of each markers is four alleles.
	// The four elements of array represents the coverage of each alleles corresponding to ["A", "T", "C", "G"].
	Alleles [4]uint64
}

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

func (snp SNP) String() string {
	return fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%d",
		snp.ID, snp.Alleles[0], snp.Alleles[1], snp.Alleles[2], snp.Alleles[3], snp.POS)
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
