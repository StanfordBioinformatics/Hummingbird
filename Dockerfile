FROM broadinstitute/gatk3:3.8-0

LABEL gatk.version="3.8-0"

ENV DEBIAN_FRONTEND noninteractive

RUN echo "deb [check-valid-until=no] http://cdn-fastly.deb.debian.org/debian jessie main" > /etc/apt/sources.list.d/jessie.list
RUN echo "deb [check-valid-until=no] http://archive.debian.org/debian jessie-backports main" > /etc/apt/sources.list.d/jessie-backports.list
RUN sed -i '/deb http:\/\/deb.debian.org\/debian jessie-updates main/d' /etc/apt/sources.list
RUN echo "Acquire::Check-Valid-Until \"false\";" > /etc/apt/apt.conf.d/100disablechecks
RUN apt-get update && apt-get install -y \
    git \
    curl \
    time

RUN curl -sL https://aka.ms/InstallAzureCLIDeb | bash

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
  unzip awscliv2.zip && \
  ./aws/install

RUN git clone https://github.com/lh3/seqtk.git && cd seqtk && make install

COPY scripts/ /usr/local/bin/
WORKDIR /tmp

ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]
