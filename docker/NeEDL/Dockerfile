# ---- stage 1: build R shiny app
FROM rocker/r2u:22.04 as needl_shiny_app

RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c( \
    'shiny', 'shinyjs', 'shinydashboard', 'shinydashboardPlus', \
    'shinyWidgets', 'shinycssloaders', 'shinyFeedback', 'shinyalert', \
    'bslib', 'rclipboard', 'sortable', 'htmltools', 'httr', \
    'data.table', 'DT', 'tidyr', 'dplyr', 'stringr', 'magrittr', \
    'rjson', 'tidyjson', 'jsonlite', 'RSQLite', 'ggplot2', 'gghighlight', \
    'plotly', 'patchwork', 'visNetwork', 'igraph', 'roxygen2', 'argparse'), \
    repos='http://cran.rstudio.com/')"

COPY shinyapp /app

WORKDIR /app
# ---- end of stage 1

# ---- stage 2: use R shiny container and build all NeEDL related binaries
FROM needl_shiny_app AS build
WORKDIR /NeEDL

ENV DEBIAN_FRONTEND=noninteractive

# install build dependencies
RUN apt-get update -y && apt-get install -y build-essential git cmake autoconf libtool pkg-config python3.10 python3.10-dev gcc g++ sqlite3 libsqlite3-dev python3.10-distutils unzip

# copy necessary NeEDL repo files
COPY ext ./ext
COPY src ./src
COPY test ./test
COPY CMakeLists.txt ./
COPY install.py ./
COPY quepistasis ./quepistasis

# copy files in /data
COPY data/dbSNP/inc_pseudogenes /NeEDL/data/dbSNP/inc_pseudogenes
COPY data/BIOGRID /NeEDL/data/BIOGRID

# initialize everything needed for the build and build all targets
RUN python3.10 install.py --clean && \ 
    python3.10 install.py --target NeEDL && \
    python3.10 install.py --target calculate_scores && \
    python3.10 install.py --target convert_to_binary && \
    python3.10 install.py --target epiJSON

# ---- end of stage 2


# ---- stage 3: use shiny app stage and add all previously compiled binaries
FROM needl_shiny_app AS release
WORKDIR /NeEDL

LABEL org.opencontainers.image.source=https://github.com/biomedbigdata/NeEDL

ENV DEBIAN_FRONTEND=noninteractive

# install production dependencies
RUN apt-get update -y && apt-get install -y python3.10 python3.10-dev sqlite3 g++ gcc

# copy everything important over to the new container
COPY --from=build /NeEDL/ext/boost_1_71_0/stage /NeEDL/ext/boost_1_71_0/stage
COPY --from=build /NeEDL/test/model/bin/* /NeEDL/test/model/bin/
COPY --from=build /NeEDL/data/dbSNP/inc_pseudogenes/snps_restruc_full_inc_pseudo.csv /NeEDL/data/dbSNP/inc_pseudogenes/snps_restruc_full_inc_pseudo.csv
COPY --from=build /NeEDL/data/BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.182.tab2.txt /NeEDL/data/BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-current.tab2.txt
COPY --from=build /NeEDL/ext/plink/plink_linux_x86_64_20230116/ /NeEDL/ext/plink/plink_linux_x86_64_20230116/

# create symlinks for all created binaries
RUN ln -s /NeEDL/test/model/bin/NeEDL /usr/local/bin/NeEDL && \
    ln -s /NeEDL/test/model/bin/calculate_scores /usr/local/bin/calculate_scores && \
    ln -s /NeEDL/test/model/bin/convert_to_binary /usr/local/bin/convert_to_binary && \
    ln -s /NeEDL/test/model/bin/convert_to_json /usr/local/bin/epiJSON
# ---- end of stage 3