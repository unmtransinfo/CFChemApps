## Installation with Docker
To run the application on your system, you need docker installed. You can run the app directly using docker by using the following command to do so

    docker-compose up -d --build 
    
You can leave out `-d` flag if you do not wish to run in detached mode

> On successfully running the server, go to http://127.0.0.1:8000/depict/ to check out the portal.

## Depict Portal
In the depict portal, we can enter SMILE representation of the molecules and can get the structure of the molecules as the output. 


## Deployment on Chiltepin
The current Chiltepin server is hosted on the `http://206.192.180.166/`. To deploy the application on the Chiltepin, add the following lines in ***/etc/apache2/apache2.conf***

    ProxyPass /depict http://localhost:8000/depict
    ProxyPassReverse /depict http://localhost:8000/depict

Finally, run the docker container using the following command in detached mode

    docker-compose -f docker-compose.yml up -d --build
In the above command, we mention the name of the docker compose file that needs to be used while running the container. This way we can separate docker compose files for different environments. 

## Deployment on AWS
With the current AWS structure, we need to create an EC2 instance to run the docker container, a single RDS instance to store the data in Relational-DBMS, and also Route 53 to use DNS service.

 - With EC2 instance, the Github repository can be connected to AWS to
   maintain CI/CD (Amazon Linux Image). 
 - RDS instance needs to be connected to the EC2 instance, which is done
   by adding information about environment variables (will also be
   handled by CI/CD pipelines).
 - The Route 53 needs to be connected to the existing (or new?) domain name, and needs some DNS information as input, which is provided by domain name provider (Varies for providers).

After the EC2 is created, connect using the SSH client option in the AWS console
On successful connection, install docker, docker-compose and git on the instance using the following commands

### Installing Docker

    sudo yum install docker 
    sudo service docker start 
    sudo usermod -a -G docker ec2-user

Use the command below to auto start docker

    sudo chkconfig docker on

### Installing Git

    sudo yum install -y git

### Installing Docker Compose

    sudo curl -L [https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname](https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname) -s)-$(uname -m) -o /usr/local/bin/docker-compose
    sudo chmod +x /usr/local/bin/docker-compose

Finally, after the required installations, reboot the instance by using

    sudo reboot
This will reboot the EC2 instance, thereby closing the SSH connection that was established. 
In order to perform the further steps, connect to the instance again using SSH.

Now that Git and Docker are installed on the EC2 instance, clone the repository using 

    git clone https://github.com/unmtransinfo/CFChemApps.git
Change the current directory to access the files in the repository cloned in the step above. 

To run the server on production, simply run the bash script called ***deploy.sh*** using

    bash deploy.sh
    
The above command will spin up the docker containers and attach them to the ports exposed by the container. 

The server can now be accessed using the public IP of the EC2 instance.

> An environment file is required to be created (.env.prod), and must contain the necessary values required to run the server in production environment. If no values are defined, containers run on default values
> The values that are required to be defined in the environment file are 
		

		DATABASE_ENGINE=<DB_ENGINE>
    		DATABASE_NAME=<DB_NAME>
    		DATABASE_USER=<DB_USERNAME>
    		DATABASE_PASSWORD=<DB_PASSWORD>
    		DATABASE_HOST=<DB_HOST>
    		DATABASE_PORT=<DB_PORT>
