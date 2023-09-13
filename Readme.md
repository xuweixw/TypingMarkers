## How to run this program

```bash
go run github.com/xuweixw/typing-markers -OUT demo -SAM demo.sam -VCF demo.vcf
```

you can get a file with .csv suffix and specifying out prefix. In this example, it is demo.csv

debug #1. 当逐行遍历文件后，指针指向文件末尾，再一次读取时，需要把指针归向原处。pointer offset.

recall File.seek()


## Requirements
