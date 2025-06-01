#!/usr/bin/env bash
# nohup bash call_metafusion_phs001038.sh >> logs/all.log 2>&1 &

convert_path() {
    local path="$1"
    # 去除文件名和扩展名
    local dirname="${path%/*}"
    local filename="${path##*/}"
    # 去除文件扩展名 .cff
    local basename="${filename%.cff}"

    # 检查目录深度
    IFS='/' read -ra ADDR <<< "$dirname"
    local len=${#ADDR[@]}

    if [[ $len -eq 2 ]]; then
        # 二级目录
        echo "${ADDR[1]}"
    elif [[ $len -eq 3 && "$basename" == "${ADDR[2]}" ]]; then
        # 三级目录，且文件名与第三级目录名相同
        echo "${ADDR[1]}-${ADDR[2]}"
    else
        # 三级目录，文件名与第三级目录名不同，或者文件名包含额外信息
        echo "${ADDR[1]}-${ADDR[2]}_${basename}"
    fi
}

# # 测试路径
# paths=("RUNS/OAK/OAK.cff" "RUNS/TCGA/THCA/THCA.cff" "RUNS/TARGET/AML/Blood_Derived_Normal.cff")

# # 转换路径并打印结果
# for path in "${paths[@]}"; do
#     echo "Original: $path -> Converted: $(convert_path "$path")"
# done
# ----------

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

# cff 的 染色体位置统一要是前缀无chr的
cff_list=`find RUNS -type f -name "*phs001038.cff"`
for cff in ${cff_list}
do
    run=$(convert_path "$cff")
    echo "Running MetaFusion2 container for run ${run}"
    docker run --cpus="20" --memory=200g --rm --entrypoint /bin/sh -v $PWD:/MetaFusion oncoharmony/metafusion2:latest -c "cd /MetaFusion/scripts && bash RUN_MetaFusion.Docker.base.sh ${cff} ${run}" > logs/${run}.log 2>&1
done
