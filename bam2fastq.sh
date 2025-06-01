## STEP 1: CONVERTING BAMS TO FASTQS WITH BIOBAMBAM - BIOBAMBAM2 2.0.54 
# 首先安装 
# conda create --name bam2fastq
# ## 激活环境
# conda activate bam2fastq
# 需要安装如下：Running autoreconf requires a complete set of tools including autoconf, automake, autoheader, aclocal and libtool.
# conda install -c bioconda libmaus2 #依赖包
# conda install -c dranew biobambam2
# cd ~/soft
# git clone https://gitlab.com/german.tischler/biobambam2.git
# git clone https://gitlab.com/german.tischler/libmaus2.git
# https://gitlab.com/german.tischler/libmaus2

# ~/rclone copy zhou-onedrive-2:/singularity/biobambam2.sif ~/soft/biobambam2
# conda activate trim
# cd /home/data/TCGA/TCGA_LIHC/RNA_seq/raw
# fastqc -t 90 -o . *1.fq.gz
# fastqc -t 90 :140924_UNC15-SN850_0394_AC5DG0ACXX_ACAGTG_L008TCGA_XF_AAMR_01A_1.fq.gz
# fastqc -t 90 :140924_UNC15-SN850_0394_AC5DG0ACXX_ACAGTG_L008TCGA_XF_AAMR_01A_2.fq.gz

# fastqc -t 90 -o . *2.fq.gz
# cd /home/zhou/soft/biobambam2
#!/bin/bash
# cd ${dir}/raw
# nohup bash ${dir}/code/bam2fastq.sh &
########################################################
## the path of each software
# module load singularity
########################################################
# https://manpages.debian.org/unstable/biobambam2/bamcollate2.1.en.html

dir=/home/data/EGA/EGAD00001005500_RNA
cd ${dir}/raw
SEND_THREAD_NUM=10 ##并行进程数
tmp_fifofile="/tmp/$$.fifo" # 脚本运行的当前进程ID号作为文件名 
mkfifo "$tmp_fifofile" # 新建一个随机fifo管道文件 
exec 6<>"$tmp_fifofile" # 定义文件描述符6指向这个fifo管道文件 
rm $tmp_fifofile 

for ((i=0;i<$SEND_THREAD_NUM;i++));do 
  echo # for循环 往 fifo管道文件中写入16个空行 
done >&6 

for i in `cat ${dir}/code/id.txt`;do # 循环 开始
  echo $i # 打印 i 
  read -u6 # 从文件描述符6中读取行（实际指向fifo管道) 
  { 
    if [  ! -f  ${dir}/raw/*${i}_1.fq.gz ]; then
    singularity exec \
    --bind /home/data \
    /home/data/singularity/biobambam2.sif \
    bamtofastq \
    collate=1 \
    exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
    filename=${dir}/bam/${i}.Aligned.sortedByCoord.out.bam \
    gz=1 \
    inputformat=bam \
    level=5 \
    outputdir=${dir}/raw \
    outputperreadgroup=1 \
    outputperreadgroupsuffixF=${i}_1.fq.gz \
    outputperreadgroupsuffixF2=${i}_2.fq.gz \
    outputperreadgroupsuffixO=${i}_o1.fq.gz \
    outputperreadgroupsuffixO2=${i}_o2.fq.gz \
    outputperreadgroupsuffixS=${i}_s.fq.gz \
    tryoq=1
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
