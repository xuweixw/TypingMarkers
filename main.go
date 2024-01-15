package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"time"
)

var (
	OUT      = flag.String("OUT", "result", "specify the prefix of all output files")
	SAMPath  = flag.String("SAM", "", "specify SAM path")
	VCF      = flag.String("VCF", "", "specify SNP path")
	min_freq = flag.Float64("min_freq", 0.03, "specify minimum frequency of each "+
		"allele, ranging from 0 for high depth to 1 for low depth")
)

const (
	VERSION    = "v1.0.2"
	UpdateDate = "2024-01-10"
)

func main() {
	flag.Parse()
	if *VCF == "" || *SAMPath == "" {
		flag.Usage()
		fmt.Println(VERSION, UpdateDate)
		os.Exit(1)
	}
	//filterSAM()
	statAlleles()
	//plotStat()

	pointInTime := time.Now()
	fmt.Printf("Congratulations, the program has finished successfully! (now %d:%d)\n", pointInTime.Hour(), pointInTime.Minute())

}

//func filterSAM() {
//	cmd := exec.Command("samtools",
//		"view",
//		"-F", "4",
//		"-o", *OUT+"_valid.sam",
//		*SAMPath)
//	_, err := cmd.Output()
//	if err != nil {
//		log.Panic(err)
//	}
//}

func statAlleles() {
	outHandle, err := os.Create(*OUT + ".tab")
	outVerboseHandle, err := os.Create(*OUT + ".verbose.csv")
	defer func() {
		err = outHandle.Close()
		check(err)
		err := outVerboseHandle.Close()
		check(err)
	}()
	check(err)

	writer := bufio.NewWriter(outHandle)
	_, err = writer.WriteString(fmt.Sprintln("#Marker\tA\tT\tC\tG"))
	check(err)

	writerVerbose := bufio.NewWriter(outVerboseHandle)
	_, err = writerVerbose.WriteString(fmt.Sprintln("#Marker\tA\tT\tC\tG"))
	check(err)

	handleVCF, err := os.Open(*VCF)
	defer handleVCF.Close()
	check(err)

	handleSAM, err := os.Open(*SAMPath)
	defer handleSAM.Close()
	check(err)

	markers := NewVCFFormat(handleVCF)

	for _, marker := range markers {

		switch m := marker.(type) {
		case SNP:
			ExtractSNP(handleSAM, &m)
			_, err = writer.WriteString(m.String() + "\n")
			check(err)
			_, err = writerVerbose.WriteString(m.VerboseString() + "\n")
			check(err)
		case MH:
			ExtractMHAlleles(handleSAM, &m)
			_, err = writer.WriteString(m.String() + "\n")
			check(err)
			_, err = writerVerbose.WriteString(m.VerboseString() + "\n")
			check(err)
		}
	}
	err = writer.Flush()
	check(err)
	err = writerVerbose.Flush()
	check(err)
}

//func plotStat() {
//	cmd := exec.Command("Rscript",
//		"plotStat.R",
//		*OUT+".csv")
//	_, err := cmd.Output()
//	if err != nil {
//		log.Panic(err)
//	}
//}

// check function acts as a utility tool for entire error check.
func check(err error) {
	if err != nil {
		log.Panic(err)
	}
}
