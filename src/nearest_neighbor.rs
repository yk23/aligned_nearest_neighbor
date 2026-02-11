use std::{
    path::Path,
    fs::File,
    io::{Write, BufWriter},
    sync::Arc,
    collections::HashSet,
    fmt::{Debug, Display, Formatter},
};
use rayon::{
    prelude::*,
};
use indicatif::{ProgressBar, ProgressStyle, ParallelProgressIterator};
use bio::io::fasta::Record;

// ======== boilerplate code START
type NeighborResult<'a> = Vec<(&'a Record, u64)>;


#[derive(Debug, PartialEq)]  // Add PartialEq here
pub enum NearestNeighborError {
    IOError(String),
    HammingDistanceError(String, String),
}


impl Display for NearestNeighborError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            NearestNeighborError::IOError(msg) => { write!(f, "{}", msg) }
            NearestNeighborError::HammingDistanceError(id1, id2) => {
                write!(f, "Hamming distance computation error between: {} and {}", id1, id2)
            }
        }
    }
}

impl From<std::io::Error> for NearestNeighborError {
    fn from(err: std::io::Error) -> NearestNeighborError {
        NearestNeighborError::IOError(format!("{}", err))
    }
}

// ======== boilerplate code END
pub(super) fn filter_records(records: &[Record], id_arr: Option<Vec<String>>) -> Vec<&Record> {
    match id_arr {
        None => records.iter().collect(),
        Some(id_list) => {
            let id_subset: HashSet<String> = HashSet::from_iter(id_list);
            records.iter()
                .filter(|record| id_subset.contains(record.id()))
                .collect()
        }
    }
}


/// Compute all nearest neighbors, and write each result to a TSV file.
pub fn compute_store_nearest_neighbors(
    records: Vec<Record>,
    out_path: &Path,
    query_ids: Option<Vec<String>>,
    db_ids: Option<Vec<String>>,
) -> Result<(), NearestNeighborError> {
    let query_records: Vec<&Record> = filter_records(&records, query_ids);
    let db_records: Vec<&Record> = filter_records(&records, db_ids);

    let results = compute_nearest_neighbors(&query_records, &db_records)?;
    let file = File::create(out_path)?;
    let mut writer = BufWriter::new(file);

    // Pre-computation is done. Now write the results to file.
    assert_eq!(results.len(), query_records.len(), "Results length should always match query length!");
    for (query_record, (neighbor_record, dist)) in query_records.iter().zip(results.iter()) {
        writeln!(writer, "{}\t{}\t{}", query_record.id(), neighbor_record.id(), dist)?;
    }
    Ok(())
}


/// Compute nearest-neighbors using multiple worker threads.
pub(super) fn compute_nearest_neighbors<'a>(
    query_records: &'a Vec<&'a Record>,
    db_records: &'a Vec<&'a Record>,
) -> Result<NeighborResult<'a>, NearestNeighborError> {
    // Setup the loop, including indicatif progress bar styling.
    let db_records = Arc::new(db_records);
    let pbar = ProgressBar::new(query_records.len() as u64);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-")
    );
    // Enable steady tick to prevent multiple threads from causing line breaks
    pbar.enable_steady_tick(std::time::Duration::from_millis(50));

    // Do the calculation, using rayon's par_iter()'s map-reduce pattern.
    let results: NeighborResult<'a> = query_records.par_iter()
        .progress_with(pbar)
        .map(|query_record| {
            let data_ref = Arc::clone(&db_records);
            compute_nearest_neighbors_single(query_record, data_ref)
        })
        .collect();
    Ok(results)
}


/// Compute the nearest neighbor between query and the collection.
/// Single-worker task, meant to be used for the map-reduce in [`compute_nearest_neighbors`].
///
/// # Arguments
///
/// * `query` - The query Fasta record.
/// * `collection` - An Arc-wrapped vector of Fasta Records.
///
/// # Returns
///
/// The nearest-neighbor Fasta record, and the hamming distance between it and the query.
fn compute_nearest_neighbors_single<'a>(query: &'a Record, collection: Arc<&'a Vec<&'a Record>>) -> (&'a Record, u64) {
    let mut best_dist: u64 = u64::MAX;
    let mut best_neighbor: Option<&Record> = None;

    for other in collection.iter().filter(|other| other.id() != query.id()) {
        // Honestly, panicking here is Ok!
        let dist = hamming_distance(query, other).unwrap_or_else(|_| panic!("Hamming distance failed."));
        if dist <= best_dist {
            best_dist = dist;
            best_neighbor = Some(other);
        }
    }

    // honestly, ok to panic here -- the collection ought to be non-empty.
    (best_neighbor.unwrap(), best_dist)
}


fn hamming_distance(x: &Record, y: &Record) -> Result<u64, NearestNeighborError> {
    if x.seq().len() != y.seq().len() {
        return Err(NearestNeighborError::HammingDistanceError(x.id().to_owned(), y.id().to_owned()));
    }

    let dist = x.seq()
        .iter()
        .zip(y.seq().iter())
        .filter(|(xi, yi)| xi != yi)
        .count() as u64;
    Ok(dist)
}


#[cfg(test)]
mod tests {
    use bio::io::fasta::Record;
    use crate::nearest_neighbor::hamming_distance;

    #[test]
    fn test_hamming() {
        let x = Record::with_attrs("input1", None, b"AAAAAAA");
        let y = Record::with_attrs("input2", None, b"AAAACCA");
        assert_eq!(hamming_distance(&x, &y), Ok(2));

        let x = Record::with_attrs("input1", None, b"AAAA");
        let y = Record::with_attrs("input2", None, b"CCCC");
        assert_eq!(hamming_distance(&x, &y), Ok(4));

        let x = Record::with_attrs("input1", None, b"AAAA");
        let y = Record::with_attrs("input2", None, b"AAAA");
        assert_eq!(hamming_distance(&x, &y), Ok(0));

        let x = Record::with_attrs("input1", None, b"AAAA");
        let y = Record::with_attrs("input2", None, b"AAA");
        assert!(hamming_distance(&x, &y).is_err());
    }
}
