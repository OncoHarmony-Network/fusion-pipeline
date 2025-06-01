#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
# docker run --rm --entrypoint /bin/sh -v $PWD:/MetaFusion mapostolides/metafusion:latest -c "cd /MetaFusion/scripts && bash RUN_MetaFusion.Docker.target.sh"
# MetaFusion/NAIVE_MERGE_RUNS/target.Jul-2024-WSX-t3/final.n2.cluster MetaFusion/NAIVE_MERGE_RUNS/target.Jun-2024-WSX-t2/final.n2.cluster 一致

# Run all
nohup docker run --cpus="10" --memory=50g --rm --entrypoint /bin/sh -v $PWD:/MetaFusion oncoharmony/metafusion2:latest -c "cd /MetaFusion/scripts && bash RUN_MetaFusion.Docker.target.sh" > ../nohup_test500_final.out 2>&1 &
