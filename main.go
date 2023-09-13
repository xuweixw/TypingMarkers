package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"time"
)

var (
	//SNP_path = flag.String("SNP", "", "SNP path")
	OUT     = flag.String("OUT", "result", "specify the prefix of all output files")
	SAMPath = flag.String("SAM", "", "specify SAM path")
	VCF     = flag.String("VCF", "", "specify SNP path")
)

func main() {
	flag.Parse()
	if *VCF == "" || *SAMPath == "" {
		flag.Usage()
		os.Exit(1)
	}
	//filterSAM()
	statAlleles()
	//plotStat()

	pointInTime := time.Now()
	fmt.Printf("Congratulations, the program has finished successfully! (now %d:%d)\n", pointInTime.Hour(), pointInTime.Minute())

}

func filterSAM() {
	cmd := exec.Command("samtools",
		"view",
		"-F", "4",
		"-o", *OUT+"_valid.sam",
		*SAMPath)
	_, err := cmd.Output()
	if err != nil {
		log.Panic(err)
	}
}

func statAlleles() {
	outHandle, err := os.Create(*OUT + ".csv")
	defer func() {
		err = outHandle.Close()
		check(err)
	}()
	check(err)

	writer := bufio.NewWriter(outHandle)
	_, err = writer.WriteString(fmt.Sprintln("Marker\tA\tT\tC\tG"))
	check(err)

	handleVCF, err := os.Open(*VCF)
	check(err)

	handleSAM, err := os.Open(*SAMPath)
	check(err)

	markers := NewVCFFormat(handleVCF)

	for _, marker := range markers {

		switch m := marker.(type) {
		case SNP:
			ExtractSNP(handleSAM, &m)
			_, err = writer.WriteString(m.String() + "\n")
			check(err)
		case MHMarker:
			fmt.Println("ok")
		}
	}
	err = writer.Flush()
	check(err)
}

func plotStat() {
	cmd := exec.Command("Rscript",
		"plotStat.R",
		*OUT+".csv")
	_, err := cmd.Output()
	if err != nil {
		log.Panic(err)
	}
}

// check function acts as a utility tool for entire error check.
func check(err error) {
	if err != nil {
		log.Panic(err)
	}
}
