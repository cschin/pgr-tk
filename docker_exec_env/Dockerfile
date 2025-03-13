FROM ubuntu:24.04
RUN apt-get update
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
RUN apt-get install -y build-essential git ssh curl clang-14 cmake libssl-dev libssl3 pkg-config libzstd-dev zstd
RUN mkdir -p /opt
ENV RUSTUP_HOME=/opt/rustup
ENV CARGO_HOME=/opt/cargo
RUN RUSTUP_HOME=${RUSTUP_HOME} CARGO_HOME=${CARGO_HOME} bash -c "curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y"
RUN . /opt/cargo/env && rustup default stable
RUN . /opt/cargo/env && cargo install --locked maturin
ENV GIT_SSH_COMMAND="ssh -o StrictHostKeyChecking=no"
RUN . /opt/cargo/env && rustup toolchain list
RUN apt-get install -y zlib1g-dev zlib1g libdeflate-dev
RUN cd /opt/ && git clone --recursive https://github.com/cschin/pgr-tk.git
RUN cd /opt/pgr-tk/ && . /opt/cargo/env && cargo install maturin && bash build.sh
RUN apt-get install -y python3-pip python3.12-venv
RUN cd /opt/pgr-tk/ && python3 -m venv /opt/pgr-tk-py/ && /opt/pgr-tk-py/bin/pip install target/wheels/pgrtk-0.6.0-cp312-cp312-linux_x86_64.whl numpy
RUN cd /opt/pgr-tk/ &&  . /opt/cargo/env && cargo clean
ENV PATH="/opt/cargo/bin:$PATH"
