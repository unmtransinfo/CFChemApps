## Installation with Docker
To run the application on your system, you need docker installed. You can run the app directly using docker by using the following command to do so

    docker-compose up -d --build 
    
You can leave out `-d` flag if you do not wish to run in detached mode

> On successfully running the server, go to http://127.0.0.1:5000/depict/ to check out the portal.

## Depict Portal
In the depict portal, we can enter SMILE representation of the molecules and can get the structure of the molecules as the output. 


## Deployment on Chiltepin
The Chiltepin server is hosted on the `http://206.192.180.166/`. To deploy the application on the Chiltepin, add the following lines in ***/etc/apache2/apache2.conf***

    ProxyPass /depict  [http://localhost:5000](http://localhost:5000/)
    ProxyPassReverse /depict  [http://localhost:5000](http://localhost:5000/)

Finally, run the docker container using the following command in detached mode

    docker-compose -f docker-compose.yml up -d --build
In the above command, we mention the name of the docker compose file that needs to be used while running the container. This way we can separate docker compose files for different environments. 
