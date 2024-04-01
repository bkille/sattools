use std::path::{Path, PathBuf};
use std::fmt::{self, Formatter, Display};
use std::fs::OpenOptions;
use std::iter::Iterator;

use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

use std::collections::HashSet;
use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};

use std::borrow::Borrow;
use bio::utils::TextSlice;
use clap::{arg, command, value_parser, ArgAction, Command};

use itertools::{Chunk, Itertools};

use bio::io::fasta;
use std::io;

use env_logger;
use log::{debug, error, info, trace, warn, LevelFilter};


/// Checks if the path pointed to by v exists.  It can be
/// any valid entity (e.g. disk file, FIFO, directory, etc.).
/// If there is any issue with permissions or failure to properly
/// resolve symlinks, or if the path is wrong, it returns
/// an Err(String), else Ok(PathBuf).
pub fn pathbuf_file_exists_validator(v: &str) -> Result<PathBuf, String> {
    // NOTE: we explicitly *do not* check `is_file()` here
    // since we want to return true even if the path is to
    // a FIFO/named pipe.
    if !Path::new(v).exists() {
        Err(String::from("No valid file was found at this path."))
    } else {
        Ok(PathBuf::from(v))
    }
}


//fn chain_lc_segments(
    //segment_complexities : Vec<f64>,
    //complexity_threshold : f64)
//{
    //j
    //segment_complexities.enumerate().for_each(|(i, complexity)|
        
    //)
//}

fn get_segment_complexities(
    rec : &bio::io::fasta::Record,
    segment_size : usize,
    kmer_size : usize)
{
    let mut write_file = OpenOptions::new()
        .create(true)
        .truncate(true)
        .write(true)
        .open(format!("/tmp/{}.txt", rec.id()))
        .unwrap();

    info!("Starting for {}", rec.id());

    let complexities : Vec<f64> = 
        rec.seq().chunks(segment_size)
        .map(|segment| HashSet::<u32>::from_iter(segment.windows(kmer_size)
            .map(|kmer| fxhash::hash32(&kmer))).len() as f64 / (segment.len() - kmer_size + 1) as f64)
        .collect();

    complexities.iter().enumerate().for_each(
        |(i, c)| writeln!(write_file, "{}\t{}\t{}\t{}", rec.id(), i*segment_size, std::cmp::min(rec.seq().len(), (i+1)*segment_size - 1), c).unwrap()
    );
}

fn main() 
{
    let matches = Command::new("trtools")
        .arg(
            arg!(<INPUT_FASTA> "Input fasta file")
                .required(true)
                .value_parser(pathbuf_file_exists_validator),
        )
        .arg(
            arg!(-s --segment <VALUE> "Segment size")
                .value_parser(value_parser!(usize))
                .default_value("5000"),
        )
        .arg(
            arg!(-k --kmer <VALUE> "Kmer size")
                .value_parser(value_parser!(usize))
                .default_value("11"),
        )
        //.arg(
            //arg!(-c --complexity <VALUE> "Kmer complexity threshold")
                //.value_parser(value_parser!(f64))
                //.default_value("0.8"),
        //)
        .arg(
            arg!(-t --threads <VALUE> "Number of threads to use")
                .value_parser(value_parser!(usize))
                .default_value("1"),
        )
        .arg(
            arg!(-d --debug ... "Turn debugging information on")
        )
        .get_matches();

    env_logger::Builder::new().filter_level(LevelFilter::max()).init();

    let input_fasta = matches.get_one::<PathBuf>("INPUT_FASTA").unwrap();
    let segment_size = *matches.get_one::<usize>("segment").unwrap();
    let kmer_size = *matches.get_one::<usize>("kmer").unwrap();
    //let complexity_threshold = *matches.get_one::<f64>("complexity").unwrap();
    let thread_count = *matches.get_one::<usize>("threads").unwrap();


    rayon::ThreadPoolBuilder::new().num_threads(thread_count).build_global().unwrap();

    let reader = fasta::Reader::new(File::open(input_fasta).unwrap());


    let mut rec_idx_pairs : Vec<(usize, String)> = reader.records().enumerate().par_bridge().map(
        |(i, result)|
        {
            let rec = result.expect("Error during fasta parsing");
            get_segment_complexities(&rec, segment_size, kmer_size);
            (i, rec.id().to_string())
        }
    ).collect();

    rec_idx_pairs.sort();

    let mut write_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("out.txt")
        .unwrap();

    rec_idx_pairs.iter().for_each(|(_i, id)| {
        let mut input = File::open(format!("/tmp/{}.txt", id)).unwrap();
        let _ = io::copy(&mut input, &mut write_file);
    });
}
