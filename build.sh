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
maturin build --release
maturin build --release --skip-auditwheel
popd
