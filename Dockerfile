FROM ubuntu:16.04

RUN apt-get update && apt-get install -y git make build-essential zlib1g-dev libncurses5-dev libncursesw5-dev samtools=0.1.19-1ubuntu1 python-pysam libbam-dev

RUN cd /tmp \
  && git clone https://github.com/ding-lab/msisensor.git msisensor_install_dir \
  && export SAMTOOLS_ROOT=/usr/lib/python2.7/dist-packages/pysam/include/samtools/ \
  && cd msisensor_install_dir \
  && git checkout 0.2 \
  && make \
  && mv /tmp/msisensor_install_dir/msisensor /usr/bin/
