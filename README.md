# Rosalind 🦀🧬
Solutions to the [Rosalind problem set](https://rosalind.info/problems/list-view/) for the rust-users study group

## Getting started

### Install rust

Full instructions here: [rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)

For Linux and macOS:

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

By default this will install into your home directory. You may need to reopen your terminal for `~/.cargo/bin` to appear in your path.

### Fork this repo

Or join the [https://github.com/Rust-Wellcome](Rust@Wellcome) organisation and clone it:

```
git@github.com:Rust-Wellcome/Rosalind.git
```

### Choose a problem and add it

```
cd Rosalind
cargo init --bin PROBLEM
```

The Rosalind problems come with sample input and output. Use these to set up a unittest in `main.rs`:

```rust
#[test]
fn sample_dataset() {
    let sequence = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";

    assert_eq!(count_nucleotides(&sequence), vec![20, 12, 17, 21]);
}
```

### Merge your results to this repo

Don't forget to add it to the table below!

## Problems and Solutions

| Problem | Description |
| ------- | ----------- |
| [KMER](KMER) | [https://rosalind.info/problems/kmer/](https://rosalind.info/problems/kmer/) |
