# Depict
The depict portal allows a user to enter SMILES/SDF representations of molecules and get the structure of the molecules as the output. SMARTS-highlighting is also supported. 


The app is publicly available [here](http://3.145.25.193/depict/). 

Depict can also be set up locally by following the instructions below.

## Installation with Docker
To run the application on your system, you need docker installed. You can run the app directly using docker by using the following command to do so
```
docker-compose up -d --build 
```  
You can leave out `-d` flag if you do not wish to run in detached mode

> On successfully running the server, go to http://127.0.0.1:8000/depict/ to check out the portal.

## Deployment on AWS
With the current AWS structure, we need to create an EC2 instance to run the docker container, a single RDS instance to store the data in Relational-DBMS, and also Route 53 to use DNS service.

 - With EC2 instance, the Github repository can be connected to AWS to
   maintain CI/CD (Amazon Linux Image). 
 - RDS instance needs to be connected to the EC2 instance, which is done
   by adding information about environment variables (will also be
   handled by CI/CD pipelines).
 - The Route 53 needs to be connected to the existing (or new?) domain name, and needs some DNS information as input, which is provided by domain name provider (Varies for providers).

> Note: If you want to allow traffic over HTTP and HTTPS, make sure to check the options under **Network Settings** while creating EC2 instance.

After the EC2 is created, connect using the SSH client option in the AWS console.
On successful connection, install docker, docker-compose and git on the instance using the following commands:

```
sudo yum install -y git

sudo curl -L https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname -s)-$(uname -m) -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
```

Finally, after the required installations, reboot the instance by using
```
sudo reboot
```
This will reboot the EC2 instance, thereby closing the SSH connection that was established.

In order to perform the further steps, connect to the instance again using SSH.

Now that Git and Docker are installed on the EC2 instance, clone the repository using 
```
git clone https://github.com/unmtransinfo/CFChemApps.git
```

Change the current directory to access the files in the repository cloned in the step above. 

Create an environment file at **cfchem/.env** with the following parameters:
```
DEBUG=0
SECRET_KEY="<your_secret_key>"
```

To run the server on production, simply run the bash script called ***deploy.sh*** using

```
bash deploy.sh
```
    
The above command will spin up the docker containers and attach them to the ports exposed by the container. 

The server can now be accessed using the public IP of the EC2 instance.

### Updating on AWS
The steps below outline how one can update the public web app after making changes and pushing to master. 
1. Login to the EC2 instance hosting Depict
2. Go to repo directory: `cd CFChemApps/`
3. Pull the changes: `git pull`
4. Compose down: `docker-compose -f docker-compose-prod.yml down`
5. Re-deploy app: `bash deploy.sh`
6. If changes are not observed then try the steps outlined in Production Troubleshooting below.

### Production Troubleshooting
Below are some common issues that can occur launching to production with potential fixes:
* **Server Error: (500)** when initially connecting to site:
    1. Run `docker ps` and find the container with image `cfchemapps-cfchemapps`. Note the container id.
    2. Run `docker exec -it <id> /bin/bash`
    3. Inside the docker container, run `python manage.py migrate`
* CSS changes not showing up because the `static/` directory is not updated inside the docker container.
    1. Stop current docker containers: `docker stop $(docker ps -a -q)`
    2. Remove all docker containers: `docker rm -f $(docker ps -a -q)`
    3. Remove all volumes: `docker volume rm $(docker volume ls -q)`
    4. Re-launch docker containers: `bash deploy.sh`
  * Note that if you've used Depict previously you'll want to clear your browser files cache for CSS-related changes to appear
