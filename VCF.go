package main

import (
	"bufio"
	"fmt"

	"os"
	"strconv"
	"strings"
)

type INFOKey string
type INFOKeyCollection map[INFOKey]string

// InfoKeyCollection
var _ = INFOKeyCollection{
	"AF":     "allele frequency",
	"AN":     "the number of alternative allele",
	"OFFSET": "a particular field for microhaplotype",
}

type INFOValue []any

func (v INFOValue) String() string {
	var newV []string
	switch v[0].(type) {
	case float64, float32:
		for i := range v {
			newV = append(newV, fmt.Sprintf("%f", v[i]))
		}
	case uint64, int64, int32, uint:
		for i := range v {
			newV = append(newV, fmt.Sprintf("%d", v[i]))
		}
	}
	return strings.Join(newV, ",")
}

/*
VCFFormat struct hold a VCF record. The Position starts from one not zero.
The format follows as:

	#CHROM	POS	ID	REF	ALT	QUAL FILTER INFO	FORMAT	NA00001
	1	2827694	rs2376870	CGTGGATGCGGGGAC	C	.	PASS	SVTYPE=DEL;END=2827762;HOMLEN=1;HOMSEQ=G;SVLEN=-68	GT:GQ	1/1:13.9
	2	321682	.	T	<DEL>	6    PASS	SVTYPE=DEL;END=321887;SVLEN=-205;CIPOS=-56,20;CIEND=-10,62	GT:GQ	0/1:12
*/
type VCFFormat struct {
	CHROM  string
	POS    uint64
	ID     string
	REF    string
	ALT    string
	QUAL   string
	FILTER string
	INFO   map[INFOKey]INFOValue
}

func (v *VCFFormat) String() string {
	var properties []string

	for key, value := range v.INFO {
		properties = append(properties, string(key)+"="+value.String())
	}

	return fmt.Sprintf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
		v.CHROM, v.POS, v.ID,
		v.REF, v.ALT,
		v.QUAL, v.FILTER,
		strings.Join(properties, ";"))
}

// parseINFO parse INFO field from a string into a map.
func parseINFO(s string) map[INFOKey]INFOValue {
	var INFO = make(map[INFOKey]INFOValue)
	property := strings.Split(s, ";")
	if len(property) == 0 {
		return nil
	}
	for _, substring := range property {
		if index := strings.IndexAny(substring, "="); index != -1 {
			INFO[INFOKey(substring[0:index])] = parseINFOValue(substring[index+1:])
		}
	}
	return INFO
}

// parseINFOValue parse INFO value of each property into a slice, whatever one or more values.
func parseINFOValue(s string) INFOValue {
	var InfoValue []any
	if values := strings.Split(s, ","); len(values) == 0 {
		InfoValue = []any{s}
	} else {
		for _, v := range values {
			InfoValue = append(InfoValue, v)
		}
	}
	return INFOValue(InfoValue)
}

// NewVCFFormat import SNP site form VCF format file and return a array of SNP.
func NewVCFFormat(file *os.File) (records []GeneticMarker) {
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")

		if len(fields) < 8 || fields[0][0] == '#' {
			continue
		}
		Pos, _ := strconv.ParseUint(fields[1], 10, 64)
		record := VCFFormat{
			CHROM: fields[0], POS: Pos, ID: fields[2],
			REF: fields[3], ALT: fields[4],
			QUAL: fields[5], FILTER: fields[6],
			INFO: parseINFO(fields[7]),
		}
		if v, ok := record.INFO["OFFSET"]; ok {
			var offset []uint64
			for i := range v {
				// Cautious, there is a type assertion for converting 'any' to 'string'
				sub, err := strconv.ParseUint(v[i].(string), 10, 64)
				check(err)
				offset = append(offset, sub)
			}
			// The keys of MH.Alleles were initially set by VCF Ref and Alt fields.
			var alleles = map[AlleleMH]float64{record.REF: 0}
			for _, allele := range strings.Split(record.ALT, ",") {
				alleles[allele] = 0
			}
			records = append(records, MH{VCFFormat: record, OffSet: offset, Alleles: alleles, RareAlleles: make(map[AlleleMH]float64)})
		} else {
			records = append(records, SNP{VCFFormat: record})
		}
	}
	return
}

type GeneticMarker interface {
	GetCHROM() string
	GetPOS() uint64
	String() string
}
