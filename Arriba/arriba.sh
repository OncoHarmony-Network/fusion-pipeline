# conda create -n arriba
# conda create -n arriba python=3.7
# conda activate arriba
# conda install -c bioconda star=2.7.6a
# conda install -c bioconda arriba=2.0.0
# nohup bash /home/data//EGA/OAK/code/arriba.sh &
# rm -rf ${dir}/arriba

# dir=/home/data #目录
# cd ${dir}
# mkdir -p soft
# cd /home/zhou2/soft
# wget https://github.com/suhrig/arriba/releases/download/v2.0.0/arriba_v2.0.0.tar.gz
# nohup https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz &
# tar -xzf arriba_v2.0.0.tar.gz
# cd arriba_v2.0.0 && make # or use precompiled binaries
##直接用原始fastq
dir=/home/data/EGA/OAK #目录
cd ${dir}
mkdir -p arriba
cd  ${dir}/arriba
# mv /home/data/EGA/OAK/raw/*/*.fastq.gz /home/data/EGA/OAK/raw
# mv /home/data/EGA/OAK/raw/*/*.fastq.gz.md5 /home/data/EGA/OAK/raw
SEND_THREAD_NUM=4 ##并行进程数
tmp_fifofile="/tmp/$$.fifo" # 脚本运行的当前进程ID号作为文件名 
mkfifo "$tmp_fifofile" # 新建一个随机fifo管道文件 
exec 6<>"$tmp_fifofile" # 定义文件描述符6指向这个fifo管道文件 
rm $tmp_fifofile 

for ((i=0;i<$SEND_THREAD_NUM;i++));do 
  echo # for循环 往 fifo管道文件中写入16个空行 
done >&6 

for i in `cat /home/data/EGA/OAK/code/OAK_RNA.txt`;do # 循环 开始
  echo $i # 打印 i 
  read -u6 # 从文件描述符6中读取行（实际指向fifo管道) 
  { 
    if [  ! -f  ${dir}/arriba/${i}.output/fusions.tsv ]; then
    mkdir -p ${dir}/arriba/${i}.output 
    cd ${dir}/arriba/${i}.output
    nohup /home/zhou2/soft/arriba_v2.0.0/run_arriba.sh /home/data/reference/hg38_ek12/STAR_index_arriba \
    /home/data/reference/hg38_ek12/gencode.v34.annotation.gtf \
   /home/data/reference/hg38_ek12/GRCh38.primary_assembly.genome.fa \
    /home/zhou2/soft/arriba_v2.0.0/database/blacklist_hg38_GRCh38_v2.0.0.tsv.gz \
    /home/zhou2/soft/arriba_v2.0.0/database/known_fusions_hg38_GRCh38_v2.0.0.tsv.gz \
    /home/zhou2/soft/arriba_v2.0.0/database/protein_domains_hg38_GRCh38_v2.0.0.gff3 \
    30 \
    ${dir}/raw/${i}_1.fastq.gz ${dir}/raw/${i}_2.fastq.gz
    rm *.bam
    fi
    i=$((i+1))
    sleep 3s #休息三秒
    echo >&6 # 再次往fifo管道文件中写入一个空行。 
  } & 
  
  # {} 这部分语句被放入后台作为一个子进程执行，所以不必每次等待3秒后执行 
  #下一个,这部分的echo $i几乎是同时完成的，当fifo中16个空行读完后 for循环 
  # 继续等待 read 中读取fifo数据，当后台的16个子进程等待3秒后，按次序 
  # 排队往fifo输入空行，这样fifo中又有了数据，for语句继续执行 
  
  pid=$! #打印最后一个进入后台的子进程id 
  echo $pid 

done

wait 
exec 6>&- #删除文件描述符6 
  
exit 0
# nohup STAR --runThreadN 80 \
# --runMode genomeGenerate \
# --genomeDir /home/data/reference/hg38_ek12/STAR_index_arriba \
# --genomeFastaFiles /home/data/reference/hg38_ek12/GRCh38.primary_assembly.genome.fa \
# --sjdbGTFfile /home/data/reference/hg38_ek12/gencode.v34.annotation.gtf \
# --sjdbOverhang 149
###
# STAR \
#    --runThreadN 8 \
#    --genomeDir /path/to/STAR_index --genomeLoad NoSharedMemory \
#    --readFilesIn read1.fastq.gz read2.fastq.gz --readFilesCommand zcat \
#    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
#    --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
#    --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
#    --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
# arriba \
#    -x /dev/stdin \
#    -o fusions.tsv -O fusions.discarded.tsv \
#    -a /path/to/assembly.fa -g /path/to/annotation.gtf \
#    -b /path/to/blacklist.tsv.gz -k /path/to/known_fusions.tsv.gz -t /path/to/known_fusions.tsv.gz -p /path/to/protein_domains.gff3