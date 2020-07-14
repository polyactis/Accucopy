extern crate bio;
extern crate byteorder;
extern crate clap;
extern crate csv;
extern crate flate2;
extern crate rust_htslib;

use bio::io::fasta;
use byteorder::*;
use std::cmp;
use std::io::prelude::*;
use std::fs;
use std::fs::File;

// The order matters!
// these macros have to be defined before other modules/functions that are will use them.
#[macro_export]
macro_rules! println_stderr(
    ($($arg:tt)*) => { {
        let r = writeln!(&mut ::std::io::stderr(), $($arg)*);
        r.expect("failed printing to stderr");
    } }
);

#[macro_export]
macro_rules! print_stderr(
    ($($arg:tt)*) => { {
        let r = write!(&mut ::std::io::stderr(), $($arg)*);
        r.expect("failed printing to stderr");
    } }
);

pub fn calc_median_usize(numbers: &mut Vec<usize>) -> usize {

    numbers.sort();

    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        (numbers[mid-1] + numbers[mid]) / 2
    } else {
        numbers[mid]
    }

}

pub fn calc_median_i32(numbers: &mut Vec<i32>) -> i32 {

    numbers.sort();

    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        (numbers[mid-1] + numbers[mid]) / 2
    } else {
        numbers[mid]
    }

}

pub mod select_het_snp;

pub mod normalize;

pub mod recall_precision;

pub fn gc_index(input_filename: &str, output_dir: &str) {
    print_stderr!("Opening file {} ...", input_filename);
    let reader = fasta::Reader::from_file(input_filename).unwrap();
    let records = reader.records();
    println_stderr!("Done.");

    if fs::metadata(output_dir).is_err() {
        // .unwrap().is_dir()
        print_stderr!("Creating directory {} ...", output_dir);
        fs::DirBuilder::new()
            .recursive(true)
            .create(output_dir)
            .unwrap();
        println_stderr!("Done.");
    }

    for (_, r) in records.enumerate() {
        let record = r.ok().expect("failed");
        let chromosome_name = record.id().to_string();
        let seq = record.seq();
        let chromosome_length = seq.len();
        println_stderr!("Working on chromosome {}, length={}", chromosome_name, chromosome_length);
        let window_size_list = vec![1, 5, 25, 125];
        for window_size in window_size_list {
            print_stderr!("\tWindow size {} ...", window_size);
            let f_name = output_dir.to_string() + "/" + &chromosome_name + ".gc" +
                &window_size.to_string() + ".bi";
            let mut outfile = File::create(&f_name).ok().expect("failed to open file");
            let mut vec = vec![];
            for start in 0..chromosome_length {
                let end = cmp::min(start + window_size, chromosome_length);
                let mut num: u8 = 0;
                for l in start..end {
                    if seq[l] == 67u8 || seq[l] == 71u8 || seq[l] == 99u8 || seq[l] == 103u8 {
                        num = num + 1;
                    }

                }
                vec.push(num);
            }
            for ii in vec {
                outfile.write_u8(ii).ok().expect("failed");
            }
            println_stderr!("Done.");
        }
    }
}