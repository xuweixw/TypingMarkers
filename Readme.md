## How to run this program

```bash
go run github.com/xuweixw/typing-markers -OUT demo -SAM demo.sam -VCF example/microhaplotype-markers.vcf -min_freq 0.2

go run TypingMarkers -OUT demo -SAM example/H28.sample.sam -VCF example/microhaplotype-markers.vcf -min_freq 0.2 
```

you can get a file with .csv suffix and specifying out prefix. In this example, it is demo.csv

debug #1. 当逐行遍历文件后，指针指向文件末尾，再一次读取时，需要把指针归向原处。pointer offset.

recall File.seek()

## Requirements

## Revise for low depth re-sequencing data

1. 等位基因必须是已知的，定义在vcf文件中；未定义的则为异常值抛出到BareAllele字典中。
2. 覆盖度为SNP组成数的倍数，不完全覆盖Read则累计至所有可能的已知等位基因中，但不加入到RareAllele。
3. overlap比对和插入缺失比对暂时不考虑。

? 检查mh20GP-034的OFFSET是20还是21。 21 checkbox

debug RareAllele统计频率不为SNP倍数的问题

goreleaser --snapshot --clean
