use std::{
    process::exit,
    path::{PathBuf},
};
use clap::Parser;

use aligned_nearest_neighbor::{
    parse_all_records, parse_record_ids,
    nearest_neighbor::compute_store_nearest_neighbors,
};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// The path to the aligned multi-FASTA file.
    #[arg(short, long, value_name = "FILE", required = true)]
    input_fasta: PathBuf,

    /// The path to output the result to. The result is a TSV-formatted table.
    #[arg(short, long, value_name = "FILE", required = true)]
    out_path: PathBuf,

    /// The number of worker threads to use.
    #[arg(short, long, value_name = "NUMBER", required = false, default_value_t = 1)]
    num_workers: usize,

    /// An optional text file, listing out fasta record IDs -- one per line.
    /// If provided, restricts the subset of queries to these IDs.
    #[arg(short, long, value_name = "FILE", required = false)]
    query_id_file: Option<PathBuf>,

    /// An optional text file, listing out fasta record IDs -- one per line.
    /// If provided, restricts the subset of database to these IDs.
    #[arg(short, long, value_name = "FILE", required = false)]
    database_id_file: Option<PathBuf>,
}


fn parse_id_file(id_file_path: Option<PathBuf>, arg_name: &str) -> Option<Vec<String>> {
    match id_file_path {
        None => {
            println!("No file specified for {} -- the entire collection will be used.", arg_name);
            None
        },
        Some(fpath) => {
            println!("Parsing {} from file: {}", arg_name, fpath.display());
            Some(parse_record_ids(&fpath).unwrap_or_else(|e| {
                eprintln!("Error reading file {}: {}", fpath.display(), e);
                exit(1);
            }))
        }
    }
}


/// Read a multi-FASTA file, where all sequences have been pre-aligned (possibly with gaps).
/// For each sequence, report the hamming-distance nearest neighbor, as well as statistics for each entry.
fn main() {
    let args = Args::parse();
    let records = parse_all_records(args.input_fasta)
        .unwrap_or_else(|err| {
            eprintln!("Unable to parse FASTA file. Reason: {}", err.message);
            exit(1)
        });
    if records.len() < 2 {
        eprintln!("There must be at least two Fasta records.");
        exit(1);
    }

    let num_workers = args.num_workers;
    println!("Number of workers = {}", num_workers);
    // Set number of threads globally at the start of your program
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_workers)
        .build_global()
        .unwrap_or_else(|err| {
            eprintln!("Failed to build global thread pool. Reason: {}", err);
            exit(1);
        });

    let query_record_ids: Option<Vec<String>> = parse_id_file(args.query_id_file, "query");
    let db_record_ids: Option<Vec<String>> = parse_id_file(args.database_id_file, "database");
    let out_tsv_path = args.out_path;
    if out_tsv_path.exists() {
        println!("The output file {} already exists. It will be overwritten!", out_tsv_path.display());
    }
    let result = compute_store_nearest_neighbors(
        records,
        &out_tsv_path,
        query_record_ids,
        db_record_ids,
    );
    match result {
        Ok(()) => {
            println!("Successfully computed nearest neighbors to: {}", out_tsv_path.display());
        }
        Err(err) => {
            println!("Error while performing nearest neighbors. Reason: {}", err);
            exit(1);
        }
    }
}
