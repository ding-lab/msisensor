FROM ubuntu:16.04

RUN apt-get update && apt-get install -y git make build-essential zlib1g-dev libncurses5-dev libncursesw5-dev

RUN cd /tmp \
  && git clone https://github.com/ding-lab/msisensor.git msisensor_install_dir \
  && cd msisensor_install_dir \
  && git checkout 0.6 \
  && make \
  && mv /tmp/msisensor_install_dir/msisensor /usr/bin/
