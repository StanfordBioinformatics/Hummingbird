FROM broadinstitute/gatk:4.1.3.0

LABEL gatk.version="4.1.3.0"

RUN apt-get update && apt-get install -y \
    bwa \
    git \
    curl

RUN curl -sL https://aka.ms/InstallAzureCLIDeb | bash

RUN git clone https://github.com/lh3/seqtk.git
WORKDIR seqtk
RUN make install

ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh
WORKDIR /tmp

ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]
