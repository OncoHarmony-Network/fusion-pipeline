#### 运行STAR fusion#####
# conda activate fus
# 链接：https://www.jianshu.com/p/7092a2eb5727
# nohup bash /home/data/EGA/OAK/code/star_fusion.sh &
##直接用原始fastq
# nohup /home/zhou2/rclone copy zhou-onedrive-2:/singularity/star_fusion.v1.9.0.simg /home/data/singularity &
# cd /home/data/reference/hg38.fusion
# tar -xzf GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play.tar.gz
dir=/home/data/EGA/OAK #目录
cd ${dir}
mkdir -p star_fusion_outdir
mkdir -p star_fusion_outdir1
# conda activate fus
cd ${dir}/star_fusion_outdir
SEND_THREAD_NUM=9 ##并行进程数
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
    mkdir -p ${i}
    cd ./${i}
    if [  ! -f   ${dir}/star_fusion_outdir1/${i}/star-fusion.fusion_predictions.tsv ]; then
  # if [  ! -f   ${dir}/star_fusion_outdir/${i}/star-fusion.fusion_predictions.tsv ]; then
  echo STAR-Fusion `date`
  nohup singularity exec \
  --bind /home/data/reference/hg38.fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/:/home/data/reference/hg38.fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ \
  --bind /home/zhou2/:/home/zhou2/ \
  --bind /home/data/:/home/data/ \
  --bind ${dir}/star_fusion_outdir/:${dir}/star_fusion_outdir/ \
  --bind ${dir}/star_fusion_outdir1/:${dir}/star_fusion_outdir1/ \
  /home/data/singularity/star_fusion/star_fusion.v1.9.0.simg \
  /usr/local/src/STAR-Fusion/STAR-Fusion \
  --left_fq ${dir}/raw/${i}_1.fastq.gz \
  --right_fq ${dir}/raw/${i}_2.fastq.gz \
  --genome_lib_dir /home/data/reference/hg38.fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir \
  --FusionInspector validate \
  --examine_coding_effect \
  --denovo_reconstruct \
  --CPU 22 \
  --output_dir ${dir}/star_fusion_outdir/${i} 
  rm ${dir}/star_fusion_outdir/${i}/*.bam
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

# mv ####

# dir=/home/data/EGA/OAK #目录
# cd ${dir}
# mkdir -p star_fusion_outdir1
# cd ${dir}/star_fusion_outdir1
# for i in `cat /home/data/EGA/OAK/code/OAK_RNA.txt`;do # 循环 开始
# mkdir -p ./${i}
# mv ${dir}/star_fusion_outdir/${i}/star-fusion.fusion_predictions.tsv ${dir}/star_fusion_outdir1/${i}
# mv ${dir}/star_fusion_outdir/${i}/star-fusion.fusion_predictions.abridged.tsv ${dir}/star_fusion_outdir1/${i}
# mv ${dir}/star_fusion_outdir/${i}/star-fusion.fusion_predictions.abridged.coding_effect.tsv ${dir}/star_fusion_outdir1/${i}
# done
# 
#  rm -rf ${dir}/star_fusion_outdir

