rustup default stable

## if necessary, install maturin with `cargo install --locked maturin`
# cargo install --locked maturin

cargo build -p pgr-bin --release
cargo install --path pgr-bin

pushd pgr-tk/
maturin build --release
maturin build --release --skip-auditwheel
popd
