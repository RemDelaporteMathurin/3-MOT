# Builds a Docker image with FEniCS stable version built from
# git sources. The image is at:
#
#    https://quay.io/repository/fenicsproject/stable
#
# Authors:
# Jack S. Hale <jack.hale@uni.lu>
#
# docker build -t remi_fenics .
# docker run -it remi_fenics




FROM quay.io/fenicsproject/stable:latest


RUN apt-get --yes install python3-pip


RUN pip3 install tqdm

RUN git clone https://github.com/RemiTheWarrior/3-MOT.git


