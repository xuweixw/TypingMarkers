package main

import (
	"fmt"
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
NewSAM parses a string to SAM struct.
seqID 99	CHR	POS	42	24M	=	16574679	321	sequence	quality	AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:150	YS:i:-20	YT:Z:CP
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
func (s *SAM) TypingMH(mh MH) AlleleMH {

	// the record don't overlap with MicroHaplotype marker.
	if s.chr != mh.CHROM || //s.cigar != "150M" ||
		mh.POS-s.pos > uint64(len(s.seq)) || mh.OffSet[len(mh.OffSet)-1]+mh.POS < s.pos {
		return ""
	}

	fmt.Println("okk", s.chr, mh)

	var alleleSNP = []string{string(s.seq[mh.POS-s.pos])} // The first SNP.
	for _, sub := range mh.OffSet {                       // The rest of SNPs.
		if i := (sub + mh.POS) - s.pos; i >= 0 && i < uint64(len(s.seq)) {
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
		return "N"
	}
	// calculate offset
	pos := AdjustPos(marker.POS-s.pos, s.cigar)
	// cigar sum is more than the length of sequences, return "N"
	if int(pos) > len(s.seq) {
		return "N"
	}

	// 1st base having position 1
	return string(s.seq[pos-1])
}

// AdjustPos offset sequence index according to cigar value 's'.
func AdjustPos(pos uint64, s string) uint64 {
	// extract cigar number into a slice.
	nucleotideGroups := strings.FieldsFunc(s, func(r rune) bool {
		return unicode.IsLetter(r)
	})
	restGroups := len(nucleotideGroups)

	for i := 0; i >= 0; {
		// extract cigar character one by one.
		i = strings.IndexFunc(s, func(r rune) bool {
			return unicode.IsLetter(r)
		})
		if i == -1 {
			return pos
		}

		currentGroup, _ := strconv.ParseUint(nucleotideGroups[len(nucleotideGroups)-restGroups], 10, 64)
		// each "D" decreases 1
		// each "I" pluses 1
		switch s[i] {
		//case 'M':
		case 'I':
			pos += currentGroup
		case 'D':
			pos -= currentGroup
		}
		restGroups--
		s = s[:i] + "-" + s[i+1:]
	}
	return pos
}
