FROM ubuntu:22.04
COPY HgvsGo /home/
WORKDIR /home/
CMD ["/home/HgvsGo"]
