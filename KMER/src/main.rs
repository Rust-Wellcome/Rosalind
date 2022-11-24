/// Rosalind problem "KMER": https://rosalind.info/problems/kmer/
///
/// Problem
/// For a fixed positive integer k, order all possible k-mers taken from an underlying alphabet lexicographically.
///
/// Then the k-mer composition of a string s can be represented by an array A for which A[m] denotes the number of times that the mth k-mer (with respect to the lexicographic order) appears in s.
///
/// Given: A DNA string s in FASTA format (having length at most 100 kbp).
///
/// Return: The 4-mer composition of s.
extern crate bio_seq;

use std::env;
use std::fs::File;
use std::io::BufReader;

use itertools::Itertools;

use bio_seq::codec::{dna::Dna, Codec};
use bio_seq::{Seq, SeqSlice};

use bio_streams::fasta::Fasta;

/// Take a sequence of DNA and return the tetramer composition
/// as a histogram represented by an array of integers
fn kmer_histogram<C: Codec, const K: usize>(seq: &SeqSlice<C>) -> Vec<usize> {
    // the WIDTH member of a type implementing Codec tells us
    // how many bits encode each character.
    //
    // For dna::Dna this is 2, so our histogram will need (2^2)^4
    // bins to count every possible 4-mer.
    let mut histo = vec![0; 1 << (C::WIDTH * K as u8)];

    for kmer in seq.kmers::<K>() {
        histo[usize::from(kmer)] += 1;
    }

    histo
}

/// Print the vector of integers as a string of numbers with
/// spaces interspersed
#[allow(unstable_name_collisions)]
fn print_histogram(histo: Vec<usize>) -> String {
    histo
        .iter()
        .map(|n| n.to_string())
        .intersperse(" ".to_string())
        .collect()
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let fasta: Fasta<BufReader<File>, Vec<u8>> =
        Fasta::new(BufReader::new(File::open(&args[1]).unwrap()));

    for contig in fasta {
        let seq: Seq<Dna> = Seq::from_vec(
            contig
                .seq
                .iter()
                .rev()
                .map(|c| Dna::from_char(*c as char).unwrap())
                .collect::<Vec<Dna>>(),
        );
        println!("{}", print_histogram(kmer_histogram::<Dna, 4>(&seq)));
    }
}

#[test]
fn sample_dataset() {
    use bio_seq::dna;
    use bio_seq::FromStr;

    let fasta = "CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGGCCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGTTTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCAAATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCGGGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGACTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTACCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG";

    // bio-seq 0.8.3 packs kmers in little-endian, we can hack around this
    // by reversing the sequence
    let rosalind_6431: Seq<Dna> = dna!(&fasta.chars().rev().collect::<String>());
    let bins: Vec<usize> = kmer_histogram::<Dna, 4>(&rosalind_6431);

    assert_eq!(print_histogram(bins), "4 1 4 3 0 1 1 5 1 3 1 2 2 1 2 0 1 1 3 1 2 1 3 1 1 1 1 2 2 5 1 3 0 2 2 1 1 1 1 3 1 0 0 1 5 5 1 5 0 2 0 2 1 2 1 1 1 2 0 1 0 0 1 1 3 2 1 0 3 2 3 0 0 2 0 8 0 0 1 0 2 1 3 0 0 0 1 4 3 2 1 1 3 1 2 1 3 1 2 1 2 1 1 1 2 3 2 1 1 0 1 1 3 2 1 2 6 2 1 1 1 2 3 3 3 2 3 0 3 2 1 1 0 0 1 4 3 0 1 5 0 2 0 1 2 1 3 0 1 2 2 1 1 0 3 0 0 4 5 0 3 0 2 1 1 3 0 3 2 2 1 1 0 2 1 0 2 2 1 2 0 2 2 5 2 2 1 1 2 1 2 2 2 2 1 1 3 4 0 2 1 1 0 1 2 2 1 1 1 5 2 0 3 2 1 1 2 2 3 0 3 0 1 3 1 2 3 0 2 1 2 2 1 2 3 0 1 2 3 1 1 3 1 0 1 1 3 0 2 1 2 2 0 2 1 1");
}
