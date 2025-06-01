#docker tag oncoharmony/metafusion2:latest shixiangwang/metafusion2:latest
docker tag oncoharmony/metafusion2:latest ghcr.io/shixiangwang/metafusion2:latest

export DOCKER_CLIENT_TIMEOUT=900
for i in {1..10}
do
   echo "这是第 $i 次循环"
   #docker push shixiangwang/metafusion2:latest
   docker push ghcr.io/shixiangwang/metafusion2:latest
done

# docker login ghcr.io -u shixiangwang #--password-stdin
