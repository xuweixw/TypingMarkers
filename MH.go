package main

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strings"
)

type AlleleMH = string

type MH struct {
	VCFFormat

	// A microhaplotype marker consists of two and more SNP sites.
	// The position of first SNP records in VCFFormat.POS, and others record in here.
	// the offset is relative to first SNP.
	OffSet []uint64

	Alleles     map[AlleleMH]float64 // 单个个体每个基因型的统计深度
	RareAlleles map[AlleleMH]float64
	Population  map[string][2]AlleleMH //每个个体的基因型, 带"."的alleleMH都用单个"."表示
}

func (mh MH) String() string {
	var genotype = mh.DetermineGenotype()
	var genotypeSlice = genotype[:]
	// There are four compounds that make of DNA (Adenine, Guanine, Cytosine, Thymine).
	sort.Strings(genotypeSlice)

	return fmt.Sprintf("%s\t%s\t%s", mh.ID, genotypeSlice[0], genotypeSlice[1])
}

// VerboseString print the details of each alleles of each markers, including coverage/count/depth, position.
func (mh MH) VerboseString() string {
	return fmt.Sprintf("%s\t%d\t%s\t%s\t%s", mh.CHROM, mh.POS, mh.ID, mapToString(mh.Alleles), mapToString(mh.RareAlleles))
}

// mapToString is generic function to print map type.
func mapToString[T int | float64 | float32](m map[string]T) string {
	var (
		s                = strings.Builder{}
		mapSortedByValue = make([][2]interface{}, 0, len(m))
	)
	for k, v := range m {
		mapSortedByValue = append(mapSortedByValue, [2]interface{}{k, v})
	}
	sort.Slice(mapSortedByValue, func(i, j int) bool {
		return mapSortedByValue[i][1].(T) > mapSortedByValue[j][1].(T)
	})
	for _, v := range mapSortedByValue {
		s.WriteString(fmt.Sprintf("%s:%0.f ", v[0].(string), v[1].(float64)))
	}
	return s.String()
}

func (mh MH) GetCHROM() string {
	return mh.CHROM
}
func (mh MH) GetPOS() uint64 {
	return mh.POS
}

// DetermineGenotype return genotype.
func (mh MH) DetermineGenotype() [2]AlleleMH {
	var (
		count    float64
		genotype []AlleleMH
	)
	// Calculate depth.
	for _, n := range mh.Alleles {
		count += n
	}
	// Omit alleles which of frequency are less than 3%.
	for allele, n := range mh.Alleles {
		if n/count > *min_freq {
			genotype = append(genotype, allele)
		}
	}

	// According to the number of retained alleles, homozygous, heterozygous or failure genotype would be determined.
	switch len(genotype) {
	case 1:
		return [2]AlleleMH{genotype[0], genotype[0]}
	case 2:
		return [2]AlleleMH{genotype[0], genotype[1]}
	default: // Void or three and more.
		return [2]AlleleMH{}
	}
}

func (mh *MH) IndividualGenotype(sample string) [2]AlleleMH {
	if diploid, ok := mh.Population[sample]; ok {
		return diploid
	} else {
		return [2]AlleleMH{".", "."}
	}
}

// 群体多态信息含量
func (mh *MH) PIC() float64 {
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
func (mh *MH) allelePopulation() map[AlleleMH]int {
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
func (mh *MH) Ae() int {
	return len(mh.allelePopulation())
}

//func (mh *MHMarker) appendSNP(n uint64) {
//	mh.SNPs = append(mh.SNPs, n)
//}

// Classify method decides which overlapping alleles belong to a complete allele.
// Each overlapping allele add 1/len(SNPs) to corresponding allele.
// Only Onw or Two alleles would be retained, other seldom alleles (frequency is lee than 10 %) would be abandoned.
func (mh *MH) Classify() {
	// extract compete alleles
	var (
		alleles  []AlleleMH
		coverage float64
	)

	for k := range mh.Alleles {
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

func (mh *MH) SimpleString() string {
	var (
		s       string
		alleles []AlleleMH
	)

	for allele := range mh.Alleles {
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

func NewMHCollection() map[string]MH {
	var set = make(map[string]MH, 0)
	return set
}

// input marker form file.
//func NewMHMarkers(path string) []MHMarker {
//	var (
//		markers []MHMarker
//		current *MHMarker
//	)
//
//	handle, err := os.Open(path)
//	if err != nil {
//		log.Fatalln(err)
//	}
//	scanner := bufio.NewScanner(handle)
//	for scanner.Scan() {
//		fields := strings.Split(scanner.Text(), "\t")
//		if len(fields) < 7 {
//			continue
//		}
//		chr := fields[0]
//		Pos, _ := strconv.ParseUint(fields[1], 10, 64)
//		mhID := fields[6]
//		if current == nil {
//			current = NewMHMarker(mhID, chr, []uint64{Pos})
//		} else if current.Chr != chr {
//			markers = append(markers, *current)
//			current = NewMHMarker(mhID, chr, []uint64{Pos})
//		} else if current.ID == mhID {
//			current.appendSNP(Pos)
//		} else {
//			markers = append(markers, *current)
//			current = NewMHMarker(mhID, chr, []uint64{Pos})
//		}
//	}
//	markers = append(markers, *current)
//	return markers
//}

//var (
//	SAM1 = "" // path
//	MH1  = MHMarker{ID: "mhGP01", Chr: "Chr1", SNPs: []uint64{16574710, 16574718, 16574732}}
//)

func ExtractMHAlleles(file *os.File, marker *MH) {
	//Debug #1: reset pointer offset. (2023-09-13)
	_, err := file.Seek(0, io.SeekStart)
	check(err)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		SAMrecord := NewSAM(scanner.Text())
		if SAMrecord == nil {
			continue
		}
		allele := SAMrecord.TypingMH(*marker)
		if allele == "" {
			continue
		}
		// MH.Alleles
		// MH.RareAlleles
		// allele

		// 如果当前等位基因完整，先判断是否已知，如已知，跳入下一个循环，如未知，加入罕见列表。
		// 如果当前等位基因残缺，先判断可否加入已知，已经累计到已知，跳入下一个循环，如未知，再判断可否加入未知。
		var commonAllele bool
		for existAllele := range marker.Alleles {
			if n := CountMatchSNPInMH(allele, existAllele); n > 0 {
				marker.Alleles[existAllele] += float64(n)
				commonAllele = true
			}
		}
		if !commonAllele && strings.Index(allele, ".") == -1 { // Un-match without Overhang
			if c, ok := marker.RareAlleles[allele]; ok {
				marker.RareAlleles[allele] = c + float64(len(allele)/2+1)
			} else {
				marker.RareAlleles[allele] = float64(len(allele)/2 + 1)
			}
		}
	}
}

// CountMatchSNPInMH returns complete and overhang match number of SNP in a MH allele or -1 for never match.
func CountMatchSNPInMH(new, exist AlleleMH) (n int) {
	if len(new) != len(exist) {
		return -1
	}
	for i := 0; i < len(new); i += 2 {
		switch {
		case new[i] == exist[i]:
			n++
		case new[i] == '.' && n == 0: // for left overhang
			n = 0
		case new[i] == '.' && n != 0: // for right overhang
			return
		default: // for never match
			return -1
		}
	}
	return
}
