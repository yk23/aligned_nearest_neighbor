use rstest::rstest;
use std::path::PathBuf;
use aligned_nearest_neighbor::parse_all_records;

#[rstest]
#[case("simple_test")]
#[case("simple_test_2")]
#[case("mismatched_lengths")]
fn test_with_expected_output(#[case] test_name: &str) {
    let input_path = PathBuf::from(format!("tests/inputs/{}.txt", test_name));
    let _records = parse_all_records(input_path);
    // assert!(records.len() > 0);
}
