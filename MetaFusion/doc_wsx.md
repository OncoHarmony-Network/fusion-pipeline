## Run and test with hg38 references and data

> Recorded by WSX based on [metafusion使用代码.txt](metafusion使用代码.txt)

- Reference: <https://github.com/ccmbioinfo/MetaFusion/wiki/metafusion_sop>

Pull docker image

```sh
docker pull mapostolides/metafusion
```

Check input:

```sh
▶ cat  RUNS/RUNS/arriba.cff | cut -f 8 | sort | uniq -c | wc -l
    3225
▶ cat  RUNS/RUNS/star_fusion.cff | cut -f 8 | sort | uniq -c | wc -l 
    2235
```

> star_fusion 大概 990 个样本检测结果是空，是对上的

> 需要注意下，避免坚果云本地同步导致的结果文件差异（后面直接从云端把数据集全部下载到服务器上去跑）

Test:

```sh
# docker run --rm -it  -v $PWD:/MetaFusion mapostolides/metafusion:latest
# docker run --rm
# docker run -v $PWD:/MetaFusion mapostolides/metafusion:latest sh -c "cd MetaFusion/scripts && bash RUN_MetaFusion.Docker.target.sh"

# Just run the script in docker with one line command
#
# https://stackoverflow.com/questions/68652198/docker-run-entrypoint-with-multiple-commands
docker run --rm --entrypoint /bin/sh -v $PWD:/MetaFusion mapostolides/metafusion:latest -c "cd /MetaFusion/scripts && bash RUN_MetaFusion.Docker.target.sh"
```

## 测试本身的代码和测试数据有没有问题

```sh
docker run --rm --entrypoint /bin/sh -v $PWD:/MetaFusion mapostolides/metafusion:latest -c "cd /MetaFusion/scripts && bash RUN_MetaFusion.Docker.sh"
# --rm 容器停止后自动删除
```

官方bug：·`MetaFusion.sh`  82 行有个变量`outdir`的使用不对

```sh
awk 'BEGIN {OFS="\t"; FS="\t";} ($1 ~ /^(([1-9])|(1[0-9])|(2[0-2])|[XY])$/ && $4 ~ /^(([1-9])|(1[0-9])|(2[0-2])|[XY])$/) { if($3 !~ /^[-+]$/) {$3="NA"}; if($6 !~ /^[-+]$/) {$6="NA"}; print}' $cff > $outdir/$(basename $cff).reformat 
```

## Docker commands

```sh
# https://colobu.com/2018/05/15/Stop-and-remove-all-docker-containers-and-images/
docker ps -aq
docker stop $(docker ps -aq)  # 关闭所有容器
docker rm $(docker ps -aq)  # 删除所有容器
```

## 问题

1. 目前 block 文件的格式还不一致：

```sh
~/NutstoreCloudBridge/MetaFusion                                                                
▶ tail reference_files/blocklist_breakpoints.bedpe 
Y       9160482 9160483 Y       9162338 9162339
Y       9458325 9458326 Y       9460764 9460765
Y       9458325 9458326 Y       9469918 9469919
Y       9633026 9633027 Y       9642383 9642384
Y       9637159 9637160 Y       9649001 9649002
Y       9737733 9737734 Y       9744349 9744350
Y       9745238 9745239 Y       9747321 9747322
Y       9745238 9745239 Y       9747969 9747970
Y       9745238 9745239 Y       9753172 9753173
Y       9834077 9834078 Y       9860804 9860805

~/NutstoreCloudBridge/MetaFusion                                                                
▶ tail reference_files_hg38/blacklist_hg38_GRCh38_v2.4.0.tsv 
Y:9795417       Y:9804774
Y:9799550       Y:9811392
Y:9832210       Y:9835847
Y:9835996       Y:9858069
Y:9839429       Y:9963317
Y:9843675       Y:9944528
Y:9900124       Y:9906740
Y:9907629       Y:9909712
Y:9907629       Y:9910360
Y:9907629       Y:9915563
```

处理：

```sh
$ cd reference_files
$ cat blacklist_hg38_GRCh38_v2.4.0.tsv | grep -v "^#" > blacklist_hg38_clean.tsv

> R
block = data.table::fread("blacklist_hg38_clean.tsv", header = FALSE)
block[, V2 := ifelse(grepl(":", V2), V2, V1)]  
r$> nrow(block)
[1] 16731166
r$> block = unique(block) 

r$> nrow(block)
[1] 16727198

# convert to bedpe format
# r$> block_pe = block |> tidyr::separate("V1", c("chr", "start", "end"), sep = ":|-") |> tidyr::separate("V2", c("chr2", "start2", "end2"), sep = ":|-")

# # transform start pos to 0-based
# block_pe$start = block_pe$start - 1
# block_pe$start2 = block_pe$start2 - 1

# data.table::fwrite(block_pe, file = "blacklist_hg38.bedpe", sep = "\t")


block |> head(10000) |>  tidyr::separate("V1", c("chr", "start", "end"), sep = ":|-") |> tidyr::separate("V2", c("chr2", "start2", "end2"), sep = ":|-")

nrow(block) / 10000

nchunks = ceiling(nrow(block) / 10000)

step = 10000
for (i in 1:nchunks) {
    message("handling chunk", i)
    data = block[((i-1)*step+1):min(i*step, nrow(block)), ]
    chunk_pe = data |> tidyr::separate("V1", c("chr", "start", "end"), sep = ":|-") |> tidyr::separate("V2", c("chr2", "start2", "end2"), sep = ":|-")
    chunk_pe = chunk_pe |>
      dplyr::mutate(end = ifelse(is.na(end), start, end), end2 = ifelse(is.na(end2), start2, end2)) |>
      dplyr::mutate(start = as.integer(start), end = as.integer(end), start2 = as.integer(start2), end2 = as.integer(end2)) |>
      dplyr::mutate(start = start - 1L, start2 = start2 - 1L) |> data.table::as.data.table()

    print(sum(is.na(c(chunk_pe$start, chunk_pe$start2, chunk_pe$end, chunk_pe$end2))))

    data.table::fwrite(chunk_pe, file = "blacklist_hg38.bedpe", sep = "\t", col.names = FALSE, append = TRUE)
}
```


2. 输入 cff 文件有问题：

```r
> table(x$V9)

Normal\r  Tumor\r 
     410    17628 
```

~~用R移除了额外的符号得到`star_fusion2.cff`用于后续测试~~。

彦昆的处理在所有的 cff 文件中都引入了 \r (回车)符号（后续要回溯）

```sh
▶ head -n 2 star_fusion.cff | cat -t 
16^I21352138^I+^I16^I21338921^I+^IRNA^ITARGET-00-RO02176-14A.star_fusion.rna_fusion.tsv^INormal^M^INBL^Istar_fusion^I4^I5^IAC008740.1^INA^ICRYM-AS1^INA
X^I38770575^I+^IX^I38666121^I+^IRNA^ITARGET-00-RO02205-14A.star_fusion.rna_fusion.tsv^INormal^M^INBL^Istar_fusion^I6^I9^IAF241728.1^INA^ITSPAN7^INA
```

> `^M`

处理：

```sh
# linux: sed -i 's/^M//g' arriba.cff
# I use macos
sed -i 's/\r//g' arriba.cff
sed -i '' $'s/\r//g' star_fusion.cff
sed -i '' $'s/\r//g' merged.cff
```

3. 同 symbol 多位点问题

```sh
Warning: Input gene annotations include multiple chr, strand, or regions (5Mb away). Skipping current gene annotation.
set([('7SK', 'chr1', 'r'), ('7SK', 'chr17', 'f'), ('7SK', 'chr3', 'r'), ('7SK', 'chr6', 'f'), ('7SK', 'chr1', 'f'), ('7SK', 'chr11', 'f')])
```

4. cff 转换问题？

输入文件中出现了移位的情况，文件名移动到了 sample_type 列

```sh
Traceback (most recent call last):
  File "rename_cff_file_genes.MetaFusion.py", line 64, in <module>
    fusion = pygeneann.CffFusion(line)
  File "/MetaFusion/scripts/pygeneann_MetaFusion.py", line 526, in __init__
    raise ValueError("sample_type value '" + tmp[8] + "' must be Tumor or Normal\nInvalid entry: " + cff_line) 
ValueError: sample_type value 'TARGET-30-PAPUAR-01A.star_fusion.rna_fusion.tsv' must be Tumor or Normal
Invalid entry: 5        40832545        -       6       70098481        -       RNA   TARGET-30-PANKFE-01A, TARGET-30-PAPUAR-01A.star_fusion.rna_fusion.tsv    Tumor   NBL   star_fusion      126     0       RPL37   NA      RPL37P15        NA
```

```sh
cat star_fusion.cff| grep "TARGET-30-PAPUAR-01A.star_fusion.rna_fusion.tsv"
5       40832545        -       6       70098481        -       RNA     TARGET-30-PANKFE-01A, TARGET-30-PAPUAR-01A.star_fusion.rna_fusion.tsv     Tumor   NBL     star_fusion     126       0       RPL37   NA      RPL37P15        NA
5       40832545        -       6       70098481        -       RNA     TARGET-30-PANKFE-01A, TARGET-30-PAPUAR-01A.star_fusion.rna_fusion.tsv     Tumor   NBL     star_fusion     126       0       RPL37   NA      RPL37P15        NA
```

这个问题完全是 多样本 id 引入的，后面预处理队列生成 cff注意解决即可。

这样测试debug：

```sh
docker run --rm -it  -v $PWD:/MetaFusion mapostolides/metafusion:latest
cd MetaFusion/scripts
# 因为是 python2，测试用下面的方式输出临时信息测试
#        print >> sys.stderr, tmp
# 注意只能临时添加，测试后要去除
```

> cff 格式文件中没有值的列一定要用 NA 填充

5. 其他后续问题

<https://gitea.zhoulab.ac.cn/OncoHarmony-Network/GeneFusionIO/issues/38#issuecomment-1876>
似乎是本地计算资源不足引入的

[随机抽样数据](https://cloud.tencent.com/developer/article/1633612)：

```r
sort -R merged.cff | head -n 100000 > merged_1e5.cff

▶ wc -l merged_1e5.cff 
  100000 merged_1e5.cff
```

```sh
docker run --rm -it  -v $PWD:/MetaFusion mapostolides/metafusion:latest
cd /MetaFusion/scripts

cff=/MetaFusion/NAIVE_MERGE_RUNS/target.Jun-2024-WSX/star_fusion.cff.reformat.renamed
gene_bed=/MetaFusion/reference_files/hg38_gene.bed
genome_fasta=/MetaFusion/reference_files/ref_genome.fa

python reann_cff_fusion.py --cff $cff --gene_bed $gene_bed --ref_fa $genome_fasta
```


6. Run all 的问题

```
Merge cff by genes and breakpoints
create BedTools object with appropriate column names
/opt/conda/lib/python2.7/site-packages/pandas/core/indexing.py:543: SettingWithCopyWarning: 
A value is trying to be set on a copy of a slice from a DataFrame.
Try using .loc[row_indexer,col_indexer] = value instead

See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
  self.obj[item] = s
Intersect fusions: NOTE: rdn=False, keeps self-intersections
/opt/conda/lib/python2.7/site-packages/pybedtools/bedtool.py:3681: UserWarning: Default names for filetype bed are:
['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
but file has 17 fields; you can supply custom names with the `names` kwarg
  "`names` kwarg" % (self.file_type, _names, self.field_count())
RUN_cluster_genes_breakpoints.sh: line 16:   195 Killed                  Rscript $fusiontools_dir/cluster_intersections.R $fid_intersection_file $fid_clusters_file
Traceback (most recent call last):
  File "/MetaFusion/scripts/generate_cluster_file.py", line 98, in <module>
    FID_clusters = [line for line in open(FIDs, "r")]
IOError: [Errno 2] No such file or directory: '/MetaFusion/NAIVE_MERGE_RUNS/target.Jul-2024-WSX-all/FID.clusters.tsv'
output cis-sage.cluster file
ReadThrough, callerfilter 2
python callerfilter_num.py --cluster /MetaFusion/NAIVE_MERGE_RUNS/target.Jul-2024-WSX-all/merged.cff.reformat.renamed.reann.WITH_SEQ.cluster --num_tools 2 | grep -v ReadThrough > /MetaFusion/NAIVE_MERGE_RUNS/target.Jul-2024-WSX-all/merged.cff.reformat.renamed.reann.WITH_SEQ.cluster.RT_filter.callerfilter.2
blocklist filter
ANC adjacent noncoding filter
Rank and generate final.cluster
Traceback (most recent call last):
  File "rank_cluster_file.py", line 27, in <module>
    i=max([len(fusion.tools) for fusion in fusion_list])
ValueError: max() arg is an empty sequence
```

> line 117-134

这样测试debug：

```sh
docker run --rm -it  -v $PWD:/MetaFusion mapostolides/metafusion:latest
# docker run --rm -it  -v $PWD:/MetaFusion metafusion2:latest
cd MetaFusion/scripts
```

```sh
fusiontools=/MetaFusion/scripts
gene_bed=/MetaFusion/reference_files/hg38_gene.bed
genome_fasta=/MetaFusion/reference_files/ref_genome.fa
outdir=/MetaFusion/NAIVE_MERGE_RUNS/target.Jul-2024-WSX-all

cff=$outdir/merged.cff.reformat.renamed.reann.WITH_SEQ

python $fusiontools/extract_closest_exons.py $cff $gene_bed $genome_fasta  > $outdir/$(basename $cff).exons
```


```sh
outdir=/MetaFusion/NAIVE_MERGE_RUNS/target.Jul-2024-WSX-all
fid_intersection_file=$outdir/FID.intersections.tsv
fid_clusters_file=$outdir/FID.clusters.tsv
fusiontools=/MetaFusion/scripts

Rscript $fusiontools/cluster_intersections.R $fid_intersection_file $fid_clusters_file

# R Console
fid_intersection_file = "/MetaFusion/NAIVE_MERGE_RUNS/target.Jul-2024-WSX-all/FID.intersections.tsv"

library(RBGL)
library(data.table)

message("Reading FID.intersections.tsv file")
FIDs=fread(fid_intersection_file, header = TRUE, showProgress = TRUE)
head(FIDs)
FIDs$V1 = NULL

# > nrow(FIDs)
# [1] 849132624
# > nrow(FIDs.nodup)
# [1] 695143960

# For some reason pairToPair output has duplicates, need to remove them
FIDs.nodup <- FIDs[!duplicated(FIDs ), ]
message(nrow(FIDs.nodup), " records to build graph")

# Build graph
message("Building graph")
g <- ftM2graphNEL(as.matrix(FIDs.nodup), W=NULL, V=NULL, edgemode="directed")
edgemode(g) <- "undirected"

# "connections" is the object containing FID clusters
connections <- connectedComp(g)
connections <- lapply(connections, function(x){paste0(x, collapse=",")})
connections <- data.frame(FIDs=matrix(unlist(connections)),stringsAsFactors=FALSE)

#write tsv
message("Writing results")
fwrite(connections, file=fid_clusters_file, quote=FALSE, sep='\t', row.names=FALSE)
```

- 观察到超过 200 G的内存占用

TARGET 3000 多例样本很难得到结果，降到500例左右测试下，测试40万行数据。

```sh
sort -R merged.cff | head -n 400000 > merged_test500.cff
```

可行。

7. 转换为 singularity 的问题。

尝试了 https://www.nas.nasa.gov/hecc/support/kb/converting-docker-images-to-singularity-for-use-on-pleiades_643.html 提到的办法，都会出现下面问题

```
FATAL:   While performing build: packer failed to pack: while unpacking tmpfs: error unpacking rootfs: unpack entry: opt/conda/pkgs/asn1crypto-1.0.1-py27_0/lib/python2.7/site-packages/asn1crypto-1.0.1.dist-info/INSTALLER: link: unpriv.link: unpriv.wrap target: no such file or directory
```

> 这个文件是存在的，但不清楚为何报错

猜测可能是以下原因之一，但目前没有处理手段了：

- docker 镜像的问题：文件权限？用户权限？构建镜像方式？
- singularity 转换的问题：singularity 在哪里出现不兼容问题？

