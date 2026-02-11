use std::{
    io::{BufRead, BufReader},
    fs::File,
    path::{Path, PathBuf},
};
use bio::io::fasta::{
    Reader as FastaReader,
    Record,
};

pub mod nearest_neighbor;


#[derive(Debug)]
pub enum FastaParseErrorKind {
    IOError,
    EmptyFile,
    LengthMismatch,
}


#[derive(Debug)]
pub struct FastaParseError {
    pub message: String,
    pub kind: FastaParseErrorKind,
}

impl From<std::io::Error> for FastaParseError {
    fn from(err: std::io::Error) -> Self {
        FastaParseError {
            message: format!("IO error: {}", err),
            kind: FastaParseErrorKind::IOError,
        }
    }
}


pub fn parse_record_ids(fpath: &Path) -> Result<Vec<String>, std::io::Error> {
    let file = File::open(fpath)?;

    let reader = BufReader::new(file);
    let mut id_list: Vec<String> = vec![];
    for line in reader.lines() {
        let line = line?.trim().to_owned();
        if line.len() > 0 {
            id_list.push(line);
        }
    }
    Ok(id_list)
}


pub fn parse_all_records(input_fasta: PathBuf) -> Result<Vec<Record>, FastaParseError> {
    let file = File::open(input_fasta)?;
    let reader = BufReader::new(file);

    let fasta_reader =  FastaReader::new(reader);
    let all_fasta_records: Vec<Record> = fasta_reader
        .records()
        .collect::<Result<Vec<Record>, std::io::Error>>()?;

    if all_fasta_records.len() == 0 {
        return Err(FastaParseError {
            message: "No records found.".to_owned(),
            kind: FastaParseErrorKind::EmptyFile,
        })
    }

    let first_len: usize = all_fasta_records.first().unwrap().seq().len();
    for (record_idx, record) in all_fasta_records.iter().enumerate() {
        if record.seq().len() != first_len {
            return Err(FastaParseError {
                message: format!(
                    "Record lengths don't match! FirstLen={}, got Len={} for record index {}",
                    first_len,
                    record.seq().len(),
                    record_idx
                ),
                kind: FastaParseErrorKind::LengthMismatch
            })
        }
    }

    Ok(all_fasta_records)
}


#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use bio::io::fasta::Record;
    use super::{parse_all_records, parse_record_ids};

    #[test]
    fn test_query_db_match() {
        let test_dir = PathBuf::from("tests/inputs/query_db/");
        let fasta_path = test_dir.join("seqs.fasta");
        let db_txt = test_dir.join("db.txt");
        let query_txt = test_dir.join("query.txt");

        let db_ids = parse_record_ids(&db_txt).unwrap();
        let query_ids = parse_record_ids(&query_txt).unwrap();
        let records = parse_all_records(fasta_path).unwrap();

        let query_records: Vec<&Record> = crate::nearest_neighbor::filter_records(&records, Some(query_ids));
        let db_records: Vec<&Record> = crate::nearest_neighbor::filter_records(&records, Some(db_ids));
        let results = crate::nearest_neighbor::compute_nearest_neighbors(&query_records, &db_records).unwrap();

        assert_eq!(results.len(), 2);
        assert_eq!(results.len(), query_records.len());

        let (res, dist) = results[0];
        assert_eq!(res.id(), "db_1");
        assert_eq!(dist, 13);

        let (res, dist) = results[1];
        assert_eq!(res.id(), "db_2");
        assert_eq!(dist, 12);
    }
}
