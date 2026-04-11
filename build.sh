# Use the native toolchain for the current architecture
if [[ "$(uname -m)" == "arm64" ]]; then
    rustup default stable-aarch64-apple-darwin
else
    rustup default stable
fi

## if necessary, install maturin with `cargo install --locked maturin`
# cargo install --locked maturin

cargo build -p agc-rs --release
cargo install --path agc-rs

cargo build -p pgr-bin --release
cargo install --path pgr-bin

pushd pgr-tk/
# Set up a uv-managed venv with Python 3.13 if not already present
if [ ! -f .venv/bin/python ]; then
    uv venv --python 3.13
    uv pip install --python .venv/bin/python maturin numpy
fi
# Build with the venv Python, ignoring any active conda environment
env -u CONDA_PREFIX \
    VIRTUAL_ENV="$(pwd)/.venv" \
    PYO3_PYTHON="$(pwd)/.venv/bin/python" \
    PYTHON_SYS_EXECUTABLE="$(pwd)/.venv/bin/python" \
    .venv/bin/maturin build --release
popd
