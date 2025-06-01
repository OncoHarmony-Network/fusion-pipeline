docker build -t oncoharmony/metafusion2 -f Dockerfile .
# clean middle layers
# docker rmi $(docker images --filter "dangling=true" -q --no-trunc)
