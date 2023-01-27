rustup default stable
cargo build -p pgr-db --release
cargo build -p pgr-bin --release
cargo install --path pgr-bin

