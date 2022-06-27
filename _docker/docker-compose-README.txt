How to use docker-compose.yml
  
build frontend, backend, and postgres-db first (all in folder container), then use swarm

$ docker swarm init // init manager node
# docker swarm init --advertise-addr <myvm1 ip> //with ip
$ docker stack deploy -c docker-compose.yml modules

# scale replicas and run docker stack deploy again

cheat sheet:
docker stack ls                                            # List stacks or apps
docker stack deploy -c <composefile> <appname>  # Run the specified Compose file
docker service ls                 # List running services associated with an app
docker service ps <service>                  # List tasks associated with an app
docker inspect <task or container>                   # Inspect task or container
docker container ls -q                                      # List container IDs
docker stack rm <appname>                             # Tear down an application
docker node ls # show nodes
docker swarm leave # node leaves swarm
docker stack ps <appname> # show running containers on nodes
docker image prune -a #remove all images which are not used by existing containers
docker container prune #When you stop a container, it is not automatically removed, so remove it
