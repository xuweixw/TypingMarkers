package main

import (
	"strconv"
	"strings"
	"unicode"
)

type SAM struct {
	seqID       string
	flag        uint64
	chr         string
	pos         uint64
	mapQ        uint64
	cigar       string
	refNext     string
	posNext     uint64
	templateLen int64
	seq         string
	qual        string
}

/*
NewSAM
A01045:788:HYM3WDSX2:4:1504:12536:25927 99
Chr1    16574508        42      150M    =       16574679        321
GCCATTTTCCCAGGAAATACTCCAATGAATCATTCCTAGGCACTTGATCTCAACAATTTCTATTTAT
CAATGTTTATCTCATTAGTGATTACAAAAAGCAGCTTTATCCTCCCCAACTCCCTCCCTTCATGTTTGTTTTTTTCTAATTAG
FFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF
FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFFFFFFF,,,FF:FFFFFF:FFFFF:FFFF,FFFFFF,FFFFF
AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:150        YS:i:-20        YT:Z:CP
*/
func NewSAM(record string) *SAM {
	var (
		align = new(SAM)
		err   error
	)

	fields := strings.Split(record, "\t")
	if len(fields) < 11 {
		return nil
	}
	align.seqID = fields[0]
	align.flag, err = strconv.ParseUint(fields[1], 10, 64)
	check(err)
	align.chr = fields[2]
	align.pos, err = strconv.ParseUint(fields[3], 10, 64)
	check(err)
	align.mapQ, err = strconv.ParseUint(fields[4], 10, 60)
	check(err)
	align.cigar = fields[5]
	align.refNext = fields[6]
	align.posNext, err = strconv.ParseUint(fields[7], 10, 64)
	check(err)
	align.templateLen, err = strconv.ParseInt(fields[8], 10, 64)
	check(err)
	align.seq = fields[9]
	align.qual = fields[10]

	return align
}

// TypingMH returns the allele. If the seq overlaps the microhaplotype, the missing SNPs represent to ".".
// If the seq doesn't overlap the microhaplotype, empty string was returned.
func (s *SAM) TypingMH(mh MHMarker) AlleleMH {
	var alleleSNP []string

	if s.chr != mh.Chr || s.cigar != "150M" || mh.SNPs[0]-s.pos > 150 || mh.SNPs[len(mh.SNPs)-1] < s.pos {
		return ""
	}
	for _, p := range mh.SNPs {
		if i := p - s.pos; i >= 0 && i < 150 {
			alleleSNP = append(alleleSNP, string(s.seq[i]))
		} else {
			alleleSNP = append(alleleSNP, ".")
		}
	}
	return AlleleMH(strings.Join(alleleSNP, "-"))
}

// TypingSNP returns the allele.
// vcf format start with 1st base having position 1
func (s *SAM) TypingSNP(marker SNP) string {
	if s.chr != marker.CHROM || marker.POS-s.pos > uint64(len(s.seq)) || marker.POS-1 < s.pos {
		return ""
	}
	// calculate offset
	pos := AdjustPos(marker.POS-s.pos, s.cigar)

	// 1st base having position 1
	return string(s.seq[pos-1])
}

func AdjustPos(pos uint64, s string) uint64 {
	var adjP, maximum uint64

	nucleotideGroups := strings.FieldsFunc(s, func(r rune) bool {
		return unicode.IsLetter(r)
	})
	totalGroups := len(nucleotideGroups)

	for i := 0; i >= 0; {
		i = strings.IndexFunc(s, func(r rune) bool {
			return unicode.IsLetter(r)
		})
		//  fmt.Println(totalGroups, nucleotideGroups, s, pos)
		currentGroup, _ := strconv.ParseUint(nucleotideGroups[len(nucleotideGroups)-totalGroups], 10, 64)
		// each "D" decreases 1
		// each "I" pluses 1
		switch s[i] {
		case 'M':
			maximum += currentGroup
			if maximum >= pos {
				return pos
			}
		case 'I':
			maximum += currentGroup
			pos += currentGroup
		case 'D':
			maximum -= currentGroup
			pos -= currentGroup
		}
		if adjP > maximum {
			adjP += maximum
		}
		totalGroups--
		s = s[:i] + "-" + s[i+1:]
	}
	return pos
}
