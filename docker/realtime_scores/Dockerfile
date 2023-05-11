FROM ubuntu:latest as base

WORKDIR /NeEDL

ENV DEBIAN_FRONTEND=noninteractive

# copy necessary NeEDL repo files
COPY ext ./ext
COPY src ./src
COPY test ./test
COPY CMakeLists.txt ./
COPY install.py ./

RUN mkdir -p ./uploads/tmp

# install dependencies
RUN apt-get update && \
	apt-get install -y build-essential git cmake autoconf libtool pkg-config python3 python3-dev libjsoncpp-dev uuid-dev openssl libssl-dev zlib1g-dev libbrotli-dev gcc g++ sqlite3 libsqlite3-dev


# add symlink to jsoncpp include files
RUN ln -s /usr/include/jsoncpp/json /usr/include/json


# build drogon web framework
RUN git clone https://github.com/drogonframework/drogon drogon
WORKDIR /NeEDL/drogon
RUN git submodule update --init
RUN mkdir build
WORKDIR /NeEDL/drogon/build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_MYSQL=OFF -DBUILD_POSTGRESQL=OFF -DBUILD_REDIS=OFF -DBUILD_SQLITE=OFF -DBUILD_EXAMPLES=OFF -DBUILD_CTL=OFF ..
RUN make
RUN make install
WORKDIR /NeEDL

# build the realtime_score program
RUN ./install.py --target realtime_scores --no-data-unpacking  --clean
