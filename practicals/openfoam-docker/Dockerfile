FROM openfoam/openfoam7-paraview56

MAINTAINER Edoardo Cipriano "edoardo.cipriano@polimi.it"

USER root

RUN apt-get -y update && \
    apt-get -y install sudo git

RUN echo "openfoam:openfoam" | chpasswd && adduser openfoam sudo
RUN chown -R openfoam:openfoam /home/openfoam

USER openfoam

RUN cd ~ && git clone https://github.com/acuoci/ATP.git

RUN echo "\n#~~~\n# My Settings\n#~~~\n\n" >> ~/.bashrc && \
    echo "PS1='\[\e[1;32m\]\u@\h:\[\e[1;34m\]\w\[\e[m\] \[\e[1;32m\]\$\[\e[m\] '" >> ~/.bashrc && \
    . ~/.bashrc

USER openfoam

