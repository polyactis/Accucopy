/*
Author:
 Yu S. Huang, polyactis@gmail.com
 Xinping Fan, 897488736@qq.com
 */
use flate2;
use flate2::Compression;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::cmp;
use std::str;
use std::collections::HashMap;
use std::io::prelude::*;
use std::fs::File;
use std::fs;
use std::path::{Path};
extern crate regex;
use self::regex::Regex;

// from lib.rs
use calc_median_usize;
use calc_median_i32;

struct OneChrData{
    chr: String,
    chr_len: usize,
    coverage_per_window: Vec<usize>,
    no_of_fragments: usize,
    coverage_per_base: f32,
    no_of_windows: usize,
}

impl<'a> OneChrData{
    fn new(chr: String,
           chr_len: usize,
           coverage_per_window: Vec<usize>,
           no_of_fragments: usize,
           coverage_per_base: f32,
           no_of_windows: usize,
    ) -> OneChrData {
        OneChrData{
            chr,
            chr_len,
            coverage_per_window,
            no_of_fragments,
            coverage_per_base,
            no_of_windows,
        }
    }
}

pub struct Normalize<'a> {
    tumor_file_path: &'a Path,
    normal_file_path: &'a Path,
    output_folder: &'a Path,
    chromosome_dict: HashMap<String, usize>,
    window_size: usize,
    //TODO expose max_fragment_len as an commandline argument.
    max_fragment_len: usize,
    float_multiplier: usize,
    max_coverage: usize,
    smooth_window_half_size: usize,
    debug: i32,
}



impl<'a> Normalize<'a> {
    pub fn new(tumor_file_path: &'a str,
           normal_file_path: &'a str,
           output_folder: &'a str,
           genome_dict_path: &'a str,
           window_size: usize,
           max_coverage: usize,
           no_of_autosomes: usize,
           smooth_window_half_size: usize,
           debug: i32,
    ) -> Normalize<'a> {
        let selected_chromosome_name: Vec<String> = (1..=no_of_autosomes).map(
            | i | String::from("chr") + &i.to_string()).collect();
        let genome_dict_file: &Path = Path::new(genome_dict_path);
        let contents = fs::read_to_string(genome_dict_file).expect(format!(
            "Something went wrong in reading {}",
            genome_dict_file.to_str().unwrap()).as_str());
        let re = Regex::new(r"SN:(?P<chr_name>chr\d+)\tLN:(?P<chr_len>\d+)").unwrap();
        let mut chromosome_dict: HashMap<String, usize> = HashMap::new();
        for hit in re.captures_iter(contents.as_str()) {
            // only get data of autosome
            let chr_name = hit["chr_name"].to_string();
            if selected_chromosome_name.contains(&chr_name) {
                chromosome_dict.insert(
                    hit["chr_name"].to_string(),
                    hit["chr_len"].parse::<usize>().unwrap());
            }
        }
        println_stderr!("Found {} chromsomes in {:?}.",
            chromosome_dict.len(), genome_dict_path);
        Normalize {
            tumor_file_path: Path::new(tumor_file_path),
            normal_file_path: Path::new(normal_file_path),
            output_folder: Path::new(output_folder),
            chromosome_dict: chromosome_dict,
            window_size,
            max_fragment_len: 1000,
            float_multiplier: 1000,
            max_coverage,
            smooth_window_half_size,
            debug,
        }
    }

    fn smooth_coverage_of_one_chr(&'a self, chr: String, chr_len: usize,
            no_of_fragments: usize, coverage_per_base: f32,
            coverage_per_window: &Vec<usize>) -> OneChrData {
        // calculate gc_ratio_per_base per coverage/fragment in each window
        let mut coverage_per_window_tmp = coverage_per_window.to_vec();
        let no_of_windows: usize = coverage_per_window.len();
        for window_index in 0..no_of_windows {
            //smooth over neighboring windows.
            let smooth_left_index = cmp::max(0 as i32,
                window_index as i32 - self.smooth_window_half_size as i32) as usize;
            let smooth_right_stop = cmp::min(no_of_windows, window_index+self.smooth_window_half_size+1);
            let mut cov_sub_vec = coverage_per_window[smooth_left_index..smooth_right_stop].to_vec();
            let coverage_smoothed = calc_median_usize(&mut cov_sub_vec);
            coverage_per_window_tmp[window_index] = coverage_smoothed;
            // coverage_raw 0 is unknown (and no gc-ratio there), do not collect.
            // use coverage_raw, not coverage_smoothed, because the latter might correspond to windows with no gc.
        }
        return OneChrData::new(chr, chr_len, coverage_per_window_tmp,
            no_of_fragments,
            coverage_per_base, no_of_windows);
    }

    fn read_in_coverage_of_genome(&'a self, input_file_path: &'a Path) -> HashMap<usize, OneChrData> {
        // let if_single_read = self.check_if_single_read(input_file_path);
        println_stderr!("Reading in genome coverage from {:?} ... ", input_file_path);
        let mut chr_idx2one_chr_data: HashMap<usize, OneChrData> = HashMap::new();
        let mut bam_reader = bam::Reader::from_path(&input_file_path).unwrap();
        let header = bam_reader.header().clone();

        let mut no_of_reads: usize = 0;
        let mut no_of_valid_fragments_chr: usize = 0;
        let mut coverage_per_window = vec![0usize];
        let mut total_insert_len_of_chr: usize = 0;
        let mut prev_chr_idx: i32 = -1;
        let mut current_chr_idx: i32;
        let mut chr: String = "chr0".to_string();
        let mut chr_len = 0usize;
        let mut no_of_windows_in_this_chr = 0usize;
        let mut no_of_unique_chrs = 0usize;
        let mut target_name_map: HashMap<i32, (String, i32)> = HashMap::new();

        // generate tid -> (chr_name, chr_idx) to save time
        for name in header.target_names() {
            let tid = header.tid(name).expect("unparsed name") as i32;
            let chr_name = str::from_utf8(name).unwrap().to_string();
            let mut chr_idx = -1i32;
            if self.chromosome_dict.contains_key(&chr_name) {
                chr_idx = chr_name.trim_start_matches("chr").parse::<i32>().unwrap() - 1;
            }
            target_name_map.insert(tid, (chr_name, chr_idx));
        }
        for r in bam_reader.records() {
            let record = r.unwrap();
            no_of_reads += 1;
            let tid = record.tid();
            if tid == -1 {
                // skip unmapped reads
                continue;
            }
            let value = target_name_map[&tid].clone();
            let chr_name = value.0;
            current_chr_idx = value.1;
            if !self.chromosome_dict.contains_key(&chr_name) {
                //skip all remaining irrelevant chromosomes
                continue;
            }

            if current_chr_idx != prev_chr_idx {
                no_of_unique_chrs += 1;
                if prev_chr_idx!=-1 {
                    let chr_idx: usize = prev_chr_idx as usize;
                    let coverage_per_base = total_insert_len_of_chr as f32 / chr_len as f32;
                    println_stderr!("{} reads so far for {:?}. Chromosome {} \
                        contains {} valid fragments.",
                        no_of_reads, &input_file_path, &chr, no_of_valid_fragments_chr);

                    // handle previous chromosome data
                    let one_chr_data = self.smooth_coverage_of_one_chr(
                        chr.clone(), chr_len,
                        no_of_valid_fragments_chr,
                        coverage_per_base,
                        &coverage_per_window);
                    chr_idx2one_chr_data.insert(chr_idx, one_chr_data);
                }
                prev_chr_idx = current_chr_idx;
                no_of_valid_fragments_chr = 0;
                total_insert_len_of_chr = 0;
                chr = String::from("chr") + &(current_chr_idx + 1).to_string();
                chr_len = self.chromosome_dict[&chr];
                no_of_windows_in_this_chr = chr_len / self.window_size;
                if chr_len % self.window_size != 0 {
                    no_of_windows_in_this_chr += 1;
                }
                coverage_per_window.clear();
                coverage_per_window = vec![0usize; no_of_windows_in_this_chr];

                println_stderr!("New chromosome {}, length={}, window size={}, \
                    no_of_windows={}.",
                    chr, chr_len, self.window_size, no_of_windows_in_this_chr);

            }

            if record.mapq()<30 {
                continue
            }
            if record.is_paired() && ( record.insert_size()<0 || 
                record.insert_size()>self.max_fragment_len as i64 ||
                !record.is_proper_pair() || record.is_mate_unmapped() || 
                !record.is_first_in_template() ||
                record.is_secondary() || record.is_duplicate() || 
                record.is_supplementary() ) {
                continue;
            }

            let mut start_pos = record.pos() as usize;
            if !record.is_paired() {
                let mut window_index: usize = (start_pos + self.window_size / 2)/self.window_size;
                if record.is_reverse(){
                    start_pos = start_pos + record.seq().len();
                    window_index = (start_pos - self.window_size / 2)/self.window_size;
                }
                if window_index >= no_of_windows_in_this_chr {
                    window_index = no_of_windows_in_this_chr-1;
                }
                coverage_per_window[window_index] += 1;
                no_of_valid_fragments_chr += 1;
                total_insert_len_of_chr += record.seq().len();
            } else {
                if record.is_reverse(){
                    start_pos = record.mpos() as usize;
                }
                let fragment_len: usize = record.insert_size() as usize;
                let stop_pos: usize = start_pos + fragment_len;
                //window start and stop index is [). The latter is not included.
                let mut window_index_start: usize = start_pos / self.window_size;
                let left_hanger = (window_index_start + 1) * self.window_size - start_pos;
                if left_hanger < self.window_size / 2 && window_index_start < no_of_windows_in_this_chr - 1 {
                    window_index_start += 1;
                }

                let mut window_index_stop: usize = stop_pos / self.window_size;
                let right_hanger = stop_pos - window_index_stop * self.window_size;
                if right_hanger > self.window_size / 2 && window_index_stop < no_of_windows_in_this_chr - 1 {
                    // cover more than half for the last window. add this window.
                    window_index_stop += 1;
                }
                // Make sure the start and stop indices within the bounds
                if window_index_start >= no_of_windows_in_this_chr {
                    window_index_start = no_of_windows_in_this_chr-1;
                }
                if window_index_stop >= no_of_windows_in_this_chr {
                    window_index_stop = no_of_windows_in_this_chr-1;
                }
                for window_index in window_index_start..window_index_stop {
                    coverage_per_window[window_index] += 1;
                    // add all gc_ratio_per_base to this window, and will average later
                }
                no_of_valid_fragments_chr += 1;
                total_insert_len_of_chr += fragment_len as usize;
            }
        }

        // handle the last chromosome
        if prev_chr_idx != -1 && self.chromosome_dict.contains_key(&chr) {
            let chr_idx: usize = prev_chr_idx as usize;
            let coverage_per_base = total_insert_len_of_chr as f32 / chr_len as f32;
            println_stderr!("{} reads so far for {:?}. Chromosome {} contains \
                {} valid fragments.",
                no_of_reads, &input_file_path, &chr, no_of_valid_fragments_chr);

            // handle previous chromosome data
            let one_chr_data = self.smooth_coverage_of_one_chr(chr.clone(), chr_len,
                                    no_of_valid_fragments_chr,
                                    coverage_per_base,
                                    &coverage_per_window);
            chr_idx2one_chr_data.insert(chr_idx, one_chr_data);
        }
        println_stderr!("Reading and smoothing of coverage from {:?} is Done. \
            {} unique chromosomes, {} reads.",
            input_file_path, no_of_unique_chrs, no_of_reads);

        chr_idx2one_chr_data
    }

    fn output_coverage_ratio_of_one_chr(&self, one_chr_data_tumor: &OneChrData, 
            one_chr_data_normal: &OneChrData,
            coverage_mean_tumor: &f32, coverage_mean_normal: &f32){
        print_stderr!("Outputting normalized coverage ratio of {} ... ", one_chr_data_normal.chr);
        let no_of_windows = one_chr_data_tumor.no_of_windows;
        let output_file_path = self.output_folder.join(format!("{}.ratio.w{}.csv.gz", 
            one_chr_data_tumor.chr, self.window_size));
        // let mut writer = csv::Writer::from_path(self.output_file_path).\
        //  expect("Failed to create a writer.");
        let output_f = File::create(&output_file_path)
            .expect(&format!("Error in creating output file {:?}", &output_file_path));
        let mut gz_writer = flate2::GzBuilder::new()
            .filename(output_file_path.file_stem().unwrap().to_str().unwrap())
            .comment("Comment")
            .write(output_f, Compression::default());
        gz_writer.write_fmt(format_args!("#chr: {}\n", 
            one_chr_data_tumor.chr)).unwrap();
        gz_writer.write_fmt(format_args!("#chromosome_len: {}\n", 
            one_chr_data_tumor.chr_len)).unwrap();
        gz_writer.write_fmt(format_args!("#window_size: {}\n", 
            self.window_size)).unwrap();
        gz_writer.write_fmt(format_args!("#no_of_windows: {}\n", 
            one_chr_data_tumor.no_of_windows)).unwrap();
        gz_writer.write_fmt(format_args!("#no_of_fragments_tumor: {}\n", 
            one_chr_data_tumor.no_of_fragments)).unwrap();
        gz_writer.write_fmt(format_args!("#coverage_per_base_tumor: {}\n", 
            one_chr_data_tumor.coverage_per_base)).unwrap();
        gz_writer.write_fmt(format_args!("#no_of_fragments_normal: {}\n", 
            one_chr_data_normal.no_of_fragments)).unwrap();
        gz_writer.write_fmt(format_args!("#coverage_per_base_normal: {}\n", 
            one_chr_data_normal.coverage_per_base)).unwrap();
        gz_writer.write_fmt(format_args!("#genome-wide-coverage-mean-tumor: {}\n",
            coverage_mean_tumor)).unwrap();
        gz_writer.write_fmt(format_args!("#genome-wide-coverage-mean-normal: {}\n", 
            coverage_mean_normal)).unwrap();

        if self.debug>0 {
            gz_writer.write_fmt(format_args!("start,coverage_ratio,coverage_tumor,\
                coverage_tumor_adj,coverage_normal,coverage_normal_adj\n")).unwrap();
        } else {
            gz_writer.write_fmt(format_args!("start,coverage_ratio,\
                coverage_tumor_adj,coverage_normal_adj\n")).unwrap();
        }

        let coverage_per_window_tumor = &one_chr_data_tumor.coverage_per_window;
        let coverage_per_window_normal = &one_chr_data_normal.coverage_per_window;
        //default coverage ratio is -1*self.float_multiplier (unknown), negative will not be outputted.
        let mut cov_ratio_int_smoothed_vec = vec![self.float_multiplier as i32 * -1; no_of_windows];

        for window_index in 0..no_of_windows {
            let coverage_tumor = coverage_per_window_tumor[window_index] as f32;
            let coverage_normal = coverage_per_window_normal[window_index] as f32;
            if coverage_normal > 0.0 && coverage_normal < self.max_coverage as f32
                && coverage_tumor >0.0 && coverage_tumor < self.max_coverage as f32 {
                //coverage_tumor usually won't be 0 because a deletion => zero \
                //  coverage only if it's 100% pure tumor.
                // its gc_ratio_in is -1. cov=0 (unsequenced => unknown , \
                //    not sure if it's deletion or not sequenced).
                // cov=0 data is not fed into GC-regression. it will cause cov_adj_array_tumor[] out of bounds error.
                let coverage_tumor_adj = coverage_tumor / coverage_mean_tumor;
                let coverage_normal_adj = coverage_normal / coverage_mean_normal;
                let coverage_ratio = coverage_tumor_adj / coverage_normal_adj;
                cov_ratio_int_smoothed_vec[window_index] = 
                    (coverage_ratio * self.float_multiplier as f32) as i32;
            }
        }
        for window_index in 0..no_of_windows {
            //smooth over neighboring windows.
            let smooth_left_index = cmp::max(0 as i32,
                window_index as i32 - self.smooth_window_half_size as i32) as usize;
            let smooth_right_stop = cmp::min(no_of_windows, window_index+self.smooth_window_half_size+1);
            let mut ratio_sub_vec = cov_ratio_int_smoothed_vec[smooth_left_index..smooth_right_stop].to_vec();
            let ratio_median_int = calc_median_i32(&mut ratio_sub_vec);
            let coverage_ratio = ratio_median_int as f32/self.float_multiplier as f32;
            if coverage_ratio>0.0 {
                // coverage_ratio=0 is excluded happen because coverage_tumor=0 \
                //  are not included in smooth calculation.
                // if coverage_ratio=0, it means the neighboring -1 (unknown) ratio has been used.
                let coverage_tumor = coverage_per_window_tumor[window_index] as f32;
                let coverage_tumor_adj = coverage_tumor / coverage_mean_tumor;

                let coverage_normal = coverage_per_window_normal[window_index] as f32;
                let coverage_normal_adj = coverage_normal / coverage_mean_normal;
                if self.debug>0 {
                    gz_writer.write_fmt(format_args!("{},{},{},{},{},{}\n", 
                        window_index * self.window_size + 1,
                        coverage_ratio, coverage_tumor, coverage_tumor_adj,
                        coverage_normal, coverage_normal_adj)
                    ).unwrap();
                } else{
                    gz_writer.write_fmt(format_args!("{},{},{},{}\n", 
                        window_index * self.window_size + 1,
                        coverage_ratio, coverage_tumor_adj, coverage_normal_adj)
                    ).unwrap();

                }
            }
            //gz_writer.write(&record.as_byte_record()[1] + "\n").\
            //  expect("Error in writing record into gzipped file.");
        }

        gz_writer.finish()
            .expect(&format!("ERROR finish() failure for gz_writer of {:?}.", &output_file_path));
        println_stderr!("Output done.");
    }

    fn calculate_genome_wide_cov_mean(&self, chr_idx2one_chr_data: &HashMap<usize, OneChrData>) -> f32{
        let mut genome_len = 0usize;
        let mut total_no_of_bases = 0f32;
        for (_chr_idx, ref one_chr_data) in chr_idx2one_chr_data.iter(){
            total_no_of_bases += one_chr_data.coverage_per_base*one_chr_data.chr_len as f32;
            genome_len += one_chr_data.chr_len;
        }
        let coverage_mean = total_no_of_bases as f32/genome_len as f32;
        println_stderr!("Genome wide mean coverage is {}", coverage_mean);
        coverage_mean

    }

    pub fn run(&self) {
        //let chr_idx2gc_map = self.read_gc_indices();
        //TODO parallel tumor and normal. not easy due to shared references (&self) not allowed in threads
        let chr_idx2one_chr_data_tumor = self.read_in_coverage_of_genome(self.tumor_file_path);
        let coverage_mean_tumor = self.calculate_genome_wide_cov_mean(&chr_idx2one_chr_data_tumor);

        let chr_idx2one_chr_data_normal = self.read_in_coverage_of_genome(self.normal_file_path);
        let coverage_mean_normal = self.calculate_genome_wide_cov_mean(&chr_idx2one_chr_data_normal);
        for chr_idx in 0..self.chromosome_dict.len() {
            self.output_coverage_ratio_of_one_chr(&chr_idx2one_chr_data_tumor[&chr_idx],
                &chr_idx2one_chr_data_normal[&chr_idx],
                &coverage_mean_tumor, &coverage_mean_normal);
        }

    }
}
