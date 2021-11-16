FROM gcc:11.2.0

ENV DEBIAN_FRONTEND noninteractive
ENV PATH /usr/bin:$PATH 
ENV CC  /usr/local/bin/gcc
ENV CXX /usr/local/bin/g++

RUN apt-get -y update && \
    apt-get install -y time cmake libgtest-dev libboost-test-dev && \
    rm -rf /var/lib/apt/lists/* 

#RUN mkdir -p /usr/src/tutorial
#COPY . /usr/src/tutorial

# run tests
# WORKDIR /usr/src/tutorial/test-debug
# RUN cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=g++ ../seqan3/test/unit/ && \
#     make -j 4 && ctest .

# run release
#WORKDIR /usr/src/tutorial/build
#RUN cmake -DCMAKE_BUILD_TYPE=Release ../source && make
# CMD ["./hello_world"]
CMD ["bash"]

