/// Rosalind problem "DNA": https://rosalind.info/problems/kmer/

fn count_nucleotides(sequence: &str) -> Vec<u32> {
    unimplemented!()
}

fn main() {
//    count_nucleotides();
    unimplemented!()
}

#[test]
fn sample_dataset() {
    let sequence = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";

    assert_eq!(count_nucleotides(&sequence), vec![20, 12, 17, 21]);
}
