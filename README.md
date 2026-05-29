# aligned_nearest_neighbor

A command-line utility for quickly computing the nearest neighbors of each sequence.
The required input is a FASTA file, containing the multiple alignment of all sequences one would want to query.

Use the `-q` and `-d` option to restrict the subset of queries and databases, respectively.

```
$ ./aligned_nearest_neighbor --help
Usage: aligned_nearest_neighbor [OPTIONS] --input-fasta <FILE> --out-path <FILE>

Options:
  -i, --input-fasta <FILE>       The path to the aligned multi-FASTA file
  -o, --out-path <FILE>          The path to output the result to. The result is a TSV-formatted table
  -n, --num-workers <NUMBER>     The number of worker threads to use [default: 1]
  -q, --query-id-file <FILE>     An optional text file, listing out fasta record IDs -- one per line. If provided, restricts the subset of queries to these IDs
  -d, --database-id-file <FILE>  An optional text file, listing out fasta record IDs -- one per line. If provided, restricts the subset of database to these IDs
  -h, --help                     Print help
  -V, --version                  Print version
```
