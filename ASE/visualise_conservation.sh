# Download chain file from UCSC
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsHg38/mm39.hg38.all.chain.gz

gunzip mm39.hg38.all.chain.gz

# Download chromosome sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes

# Download UCSC scripts
rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/chainToBigChain ./
chmod +x chainToBigChain
rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/bedToBigBed ./
chmod +x bedToBigBed

./chainToBigChain ~/Downloads/mm39.hg38.all.chain ~/Downloads/mm39.hg38.all.bigChain ~/Downloads/mm39.hg38.all.bigLinkOut

awk 'BEGIN{OFS="\t"}{
  if ($12 ~ /e\+/) {
    $12 = sprintf("%.0f", $12)
  }
  print
}' /Users/graylachlan/LRZ Sync+Share/LGray/LINKS_study/IGV/mm39.hg38.all.bigChain \
> /Users/graylachlan/LRZ Sync+Share/LGray/LINKS_study/IGV/mm39.hg38.all.bigChain.fixed

## Create bigChain.as with the following content:
table bigChain
"pairwise alignment chains (chainToBigChain observed output)"
(
string chrom;        "Target chromosome"
uint chromStart;     "Target start"
uint chromEnd;       "Target end"
string name;         "Chain name"
uint score;          "Chain score"
char[1] strand;      "Target strand"
uint tSize;          "Target chromosome size"
string qName;        "Query chromosome"
uint qSize;          "Query chromosome size"
uint qStart;         "Query start"
uint qEnd;           "Query end"
uint chainId;        "Chain ID"
)

# Convert to bigChain format
./bedToBigBed \
  -type=bed6+6 \
  -as=/Users/graylachlan/LRZ Sync+Share/LGray/LINKS_study/IGV/bigChain.as \
  -tab \
  /Users/graylachlan/LRZ Sync+Share/LGray/LINKS_study/IGV/mm39.hg38.all.bigChain.fixed \
  /Users/graylachlan/LRZ Sync+Share/LGray/LINKS_study/IGV/mm39.chrom.sizes \
  /Users/graylachlan/LRZ Sync+Share/LGray/LINKS_study/IGV/mm39.hg38.all.bb