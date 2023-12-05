package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"
	"testing"
)

func TestAdjustPos(t *testing.T) {
	fmt.Println(AdjustPos(80, "65M2D1I35M"))
}

func TestCountMatchSNPInMH(t *testing.T) {
	var pair = [][2]AlleleMH{
		{"A-T-A-G-T", "A-T-A-G-T"}, //Match == 5
		{"A-T-T-G-T", "A-T-A-G-T"}, // Never Match == -1
		{".-.-.-A-T", "A-T-A-G-T"}, // not overhang == -1
		{".-.-.-G-T", "A-T-A-G-T"}, // left overhang == 2
		{"A-T-A-.-.", "A-T-A-G-T"}, // right overhang == 3
	}
	for _, v := range pair {
		fmt.Printf("%s set to %d\n", v[0], CountMatchSNPInMH(v[0], v[1]))
	}
}

var SAMPLES = [][2]string{
	{"1005", "1005"},
	{"1014", "1014"},
	{"1057", "1057"},
	{"1074", "1074"},
	{"1078", "1078"},
	{"1121", "1121"},
	{"1126", "1126"},
	{"1127", "1127"},
	{"1131", "1131"},
	{"1132", "1132"},
	{"1148", "1148"},
	{"1171", "1171"},
	{"1172", "1172"},
	{"1173", "1173"},
	{"1185", "1185"},
	{"1186", "1186"},
	{"1202", "1202"},
	{"1230", "1230"},
	{"1231", "1231"},
	{"1233", "1233"},
	{"1235", "1235"},
	{"1237", "1237"},
	{"1238", "1238"},
	{"1246", "1246"},
	{"1247", "1247"},
	{"1250", "1250"},
	{"1251", "1251"},
	{"1252", "1252"},
	{"386", "386"},
	{"386", "386"},
	{"386", "386"},
	{"425", "425"},
	{"491", "491"},
	{"494", "494"},
	{"515", "515"},
	{"520", "520"},
	{"530", "530"},
	{"573", "573"},
	{"593", "593"},
	{"614", "614"},
	{"627", "627"},
	{"630", "630"},
	{"670", "670"},
	{"685", "685"},
	{"711", "711"},
	{"732", "732"},
	{"738", "738"},
	{"763", "763"},
	{"792", "792"},
	{"801", "801"},
	{"811", "811"},
	{"813", "813"},
	{"814", "814"},
	{"824", "824"},
	{"839", "839"},
	{"855", "855"},
	{"882", "882"},
	{"917", "917"},
	{"920", "920"},
	{"926", "926"},
	{"947", "947"},
	{"965", "965"},
	{"966", "966"},
	{"985", "985"},
	{"987", "987"},
	{"988", "988"},
	{"991", "991"},
	{"997", "997"},
	{"998", "998"},
	{"AW", "AW"},
	{"BX", "BX"},
	{"FDD", "FDD"},
	{"FDSW202100238-1r_L3", "unnamed_jiaozi_progeny_2"},
	{"FDSW202100239-1r_L1", "unnamed_qiyuan_progeny_1"},
	{"FDSW202100239-1r_L3", "unnamed_qiyuan_progeny_1"},
	{"FDSW202100240-1r_L3", "unnamed_pingping"},
	{"FDSW202100241-1r_L1", "312-A"},
	{"FDSW202100241-1r_L3", "312-B"},
	{"FDSW202100242-1r_L3", "unnamed_cancan"},
	{"FDSW202100243-1r_L3", "unnamed_qinghe_progeny_1"},
	{"FDSW202100244-1r_L3", "unnamed_shishi_progeny_1"},
	{"FDSW202100245-1r_L3", "287"},
	{"FDSW202100246-1r_L3", "509"},
	{"FDSW202100247-1r_L3", "unnamed_aibang_progeny_1"},
	{"FDSW202100249-1r_L1", "unnamed_A4"},
	{"FDSW202100249-1r_L3", "unnamed_A4"},
	{"FDSW202100956-1r_L4", "870"},
	{"FDSW202100958-1r_L1", "858-A"},
	{"FDSW202100958-1r_L4", "858-B"},
	{"FDSW202100959-1r_L4", "857"},
	{"FDSW202100960-1r_L4", "725"},
	{"FDSW202100961-1r_L4", "678"},
	{"FDSW202100962-1r_L4", "871"},
	{"FDSW202100963-1r_L4", "681"},
	{"FDSW202143894-1r_L4", "1203"},
	{"FDSW202143895-1r_L1", "1224-A"},
	{"FDSW202143895-1r_L2", "1224-B"},
	{"FDSW202290695-1r_L1", "FDSW202290695-1r_L1"},
	{"FDSW202290695-1r_L2", "FDSW202290695-1r_L2"},
	{"FDSW202290696-1r_L1", "FDSW202290696-1r_L1"},
	{"FDSW202290696-1r_L2", "FDSW202290696-1r_L2"},
	{"FDSW210021237-1r_L1", "765"},
	{"FDSW210021238-1r_L1", "725"},
	{"FDSW210021239-1r_L1", "870-A"},
	{"FDSW210021239-1r_L2", "870-B"},
	{"FDSW210021240-1r_L1", "871"},
	{"FDSW210021241-1r_L1", "681"},
	{"FDSW210021242-1r_L1", "857"},
	{"FDSW210021243-1r_L1", "858"},
	{"FDSW210021244-1r_L1", "598-A"},
	{"FDSW210021244-1r_L2", "598-B"},
	{"FDSW210021700-1b_L2", "523"},
	{"FDSW210021701-1r_L1", "570"},
	{"FDSW210021702-1r_L1", "788"},
	{"FDSW210021703-1r_L1", "680"},
	{"FDSW210021704-1r_L1", "1032"},
	{"FDSW210021705-1r_L1", "1033"},
	{"FDSW210021706-1r_L1", "649"},
	{"FDSW210021707-1r_L1", "990"},
	{"FDSW210021708-1r_L1", "532"},
	{"FDSW210021709-1r_L1", "522"},
	{"FDSW210021710-1r_L1", "637"},
	{"FDSW210021711-1r_L2", "762"},
	{"FDSW210021712-1r_L2", "584"},
	{"FDSW210021713-1r_L2", "480"},
	{"FDSW210021714-1r_L2", "703"},
	{"FDSW210021715-1r_L2", "635"},
	{"FDSW210021716-1r_L2", "963"},
	{"FDSW210021717-1r_L2", "731"},
	{"FDSW210021718-1r_L2", "555"},
	{"FDSW210021719-1r_L2", "665"},
	{"FDSW210021720-1b_L2", "671"},
	{"FDSW210021721-1b_L2", "796"},
	{"FDSW210021722-1r_L2", "342"},
	{"FDSW210021723-1r_L3", "645"},
	{"FDSW210021724-1r_L3", "561"},
	{"FDSW210021725-1r_L3", "314"},
	{"FDSW210021726-1r_L3", "537"},
	{"FDSW210021727-1b_L2", "761"},
	{"FDSW210021728-1r_L3", "783"},
	{"FDSW210021729-1r_L3", "709"},
	{"FDSW210021730-1r_L3", "490"},
	{"FDSW210021731-1r_L3", "278"},
	{"FDSW210021732-1r_L3", "717"},
	{"FDSW210021733-2r_L3", "823"},
	{"FDSW210021734-1r_L4", "362"},
	{"FDSW210021735-1r_L4", "unnamed_jiumei"},
	{"FDSW210021736-1r_L4", "895"},
	{"FDSW210021737-1r_L4", "954"},
	{"FDSW210021738-1r_L4", "894"},
	{"FDSW210021739-1r_L4", "496"},
	{"FDSW210021740-1r_L4", "394"},
	{"FDSW210021741-1r_L4", "986"},
	{"FDSW210021742-1r_L4", "1082"},
	{"FDSW210021743-1r_L4", "1081"},
	{"FDSW210021744-1r_L2", "1077"},
	{"FDSW210021745-1r_L2", "948-A"},
	{"FDSW210021745-1r_L4", "948-B"},
	{"FDSW210021746-1r_L2", "994-A"},
	{"FDSW210021746-1r_L4", "994-B"},
	{"FDSW210021747-1r_L2", "1038"},
	{"FDSW210021748-1r_L2", "831"},
	{"FDSW210021749-1r_L2", "1016"},
	{"FDSW210021750-1r_L2", "1017-A"},
	{"FDSW210021750-1r_L4", "1017-B"},
	{"FDSW210021751-1r_L2", "1073"},
	{"FDSW210021752-1r_L2", "1056-A"},
	{"FDSW210021752-1r_L4", "1056-B"},
	{"FDSW210021753-1r_L2", "1122"},
	{"FDSW210021754-1r_L3", "1129"},
	{"FDSW210021755-1r_L3", "1013"},
	{"FDSW210021756-1r_L3", "540"},
	{"FDSW210021757-1b_L2", "401"},
	{"FDSW210021758-1r_L3", "1145"},
	{"FDSW210021759-1r_L3", "1124"},
	{"FDSW210021760-1r_L3", "1125"},
	{"FDSW210021761-1r_L3", "652"},
	{"FDSW210021762-1r_L3", "1119-A"},
	{"FDSW210021762-1r_L4", "1119-B"},
	{"FDSW210021763-1r_L3", "724-A"},
	{"FDSW210021763-1r_L4", "724-B"},
	{"FDSW210021764-1r_L3", "853"},
	{"FDSW210021765-1r_L3", "703"},
	{"FDSW210021766-1r_L2", "883"},
	{"FDSW210021767-1r_L3", "1008-A"},
	{"FDSW210021767-1r_L4", "1008-B"},
	{"FDSW210021768-1r_L3", "958-A"},
	{"FDSW210021768-1r_L4", "958-B"},
	{"FDSW210021769-1r_L3", "unnamed_lili_progny_1-A"},
	{"FDSW210021769-1r_L4", "unnamed_lili_progny_1-B"},
	{"FDSW210021770-1b_L2", "unnamed_jiaozi_progeny_1"},
	{"JX", "JX"},
	{"LH", "LH"},
	{"YS", "YS"},
}

var MARKERS = []string{
	"mh0XGP-001",
	"mh01GP-002",
	"mh01GP-003",
	"mh02GP-004",
	"mh02GP-005",
	"mh02GP-006",
	"mh02GP-007",
	"mh03GP-008",
	"mh03GP-009",
	"mh04GP-010",
	"mh05GP-011",
	"mh06GP-012",
	"mh06GP-013",
	"mh07GP-014",
	"mh07GP-015",
	"mh08GP-016",
	"mh08GP-017",
	"mh09GP-018",
	"mh10GP-019",
	"mh11GP-020",
	"mh12GP-021",
	"mh12GP-022",
	"mh13GP-023",
	"mh14GP-024",
	"mh15GP-025",
	"mh16GP-026",
	"mh16GP-027",
	"mh17GP-028",
	"mh17GP-029",
	"mh19GP-030",
	"mh19GP-031",
	"mh19GP-032",
	"mh20GP-034",
}

func TestRunCommand(t *testing.T) {
	for _, s := range SAMPLES {
		cmd := exec.Command("go", "run", "TypingMarkers",
			"-OUT", "out/"+s[1],
			"-SAM", "example/panda/"+s[0]+".sam",
			"-VCF", "example/microhaplotype-markers.vcf",
			"-min_freq", "0.2")
		n, _ := cmd.Output()

		fmt.Println("finished: ", s, string(n))
	}
}

type PanelMH struct {
	panel map[string][2]AlleleMH
}

func (p PanelMH) String() string {
	var sortedAllele []AlleleMH
	for _, m := range MARKERS {
		sortedAllele = append(sortedAllele, p.panel[m][0], p.panel[m][1])
	}
	return strings.Join(sortedAllele, "\t")
}

func Read(file *os.File) PanelMH {
	var panel = make(map[string][2]AlleleMH, len(MARKERS))
	reader := bufio.NewScanner(file)
	for reader.Scan() {
		if reader.Text()[0] == '#' {
			continue
		}
		if sub := strings.Split(reader.Text(), "\t"); len(sub) == 3 {
			panel[sub[0]] = [2]AlleleMH{sub[1], sub[2]}
		}
	}
	return PanelMH{panel: panel}
}

func TestReadAndMerge(t *testing.T) {
	var (
		header         = strings.Builder{}
		diploidMarkers []string
	)
	for i := 0; i < len(MARKERS); i++ {
		diploidMarkers = append(diploidMarkers, MARKERS[i], MARKERS[i])
	}
	header.WriteString("MARKER\t")
	header.WriteString(strings.Join(diploidMarkers, "\t"))
	fmt.Println(header.String())

	for _, s := range SAMPLES {
		file, err := os.Open("out/" + s[1] + ".tab")
		if err != nil {
			log.Panic(err)
		}
		panel := Read(file)
		fmt.Printf("%s\t%s\n", s[1], panel.String())
	}
}
