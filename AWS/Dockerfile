FROM broadinstitute/gatk:4.1.3.0

LABEL gatk.version="4.1.3.0"

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

RUN apt-get update && apt-get install -y \
    bwa \
    git

RUN git clone https://github.com/lh3/seqtk.git
WORKDIR seqtk
RUN make install

ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh
WORKDIR /tmp

ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]
