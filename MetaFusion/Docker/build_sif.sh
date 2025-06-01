# https://www.nas.nasa.gov/hecc/support/kb/converting-docker-images-to-singularity-for-use-on-pleiades_643.html

# 1. Using an Existing Docker Image on Docker Hub
# https://doc.nju.edu.cn/books/e1654/page/ghcr
#singularity build --sandbox metafusion2 docker://ghcr.io/shixiangwang/metafusion2
singularity build --sandbox metafusion2 docker://ghcr.nju.edu.cn/shixiangwang/metafusion2

# 2. Using an Existing Docker Image on Your Local Machine
#docker save bd660d2a5b8f -o metafusion2.tar 
#singularity build --tmpdir=/tmp --sandbox metafusion2 docker-archive://metafusion2.tar
