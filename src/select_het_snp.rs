/*
Author:
 Yu S. Huang, polyactis@gmail.com
 Xinping Fan, 897488736@qq.com
 */
extern crate time;

use flate2;
use flate2::Compression;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use std::fs::File;
use std::io::prelude::*;
use std::path::{Path};
use std::str;


struct Record {
    chr: String,
    pos: u32,
    tumor_depth: i32,
    tumor_ro: i32,
    tumor_ao: i32,
    normal_depth: i32,
    normal_ro: i32,
    normal_ao: i32,
}   

struct SnpSummary {
    no_of_total_records: u32,
    no_of_good_hets_in_normal: u32,
    no_of_good_hets: u32,
}

pub struct SelectHetSNP<'a> {
    snp_file: &'a Path,
    output_file_path: &'a Path,
    min_coverage: usize,
    max_coverage: usize,
}



impl<'a> SelectHetSNP<'a> {
    pub fn new(snp_file: &'a str,
           output_file_path: &'a str,
           min_coverage: usize,
           max_coverage: usize,
    ) -> SelectHetSNP<'a> {
        SelectHetSNP {
            snp_file: Path::new(snp_file),
            output_file_path: Path::new(output_file_path),
            min_coverage,
            max_coverage,
        }
    }

    fn select_het_snp(&self){
        let mut vcf = bcf::Reader::from_path(&self.snp_file).ok().expect(
            "Error opening SNP file.");
        let mut snp_list: Vec<Record> = Vec::new();
        let vcf_header = vcf.header().clone();
        let mut snp_summary = SnpSummary{no_of_total_records:0, 
                              no_of_good_hets_in_normal: 0,no_of_good_hets: 0};
        for rec in vcf.records() {
            snp_summary.no_of_total_records += 1;
            let mut record = rec.ok().expect("Error reading record.");
            let sample_1_genotype: String;
            let _sample_2_genotype: String;
            {
                let genotypes = record.genotypes().expect("Error reading genotypes");
                sample_1_genotype = format!("{}", genotypes.get(0));
                _sample_2_genotype = format!("{}", genotypes.get(1));
            }
            let sample1_allele_depth: Vec<i32>;
            let sample2_allele_depth: Vec<i32>;
            {
                let allele_depth_vec = record.format(b"AD").integer().unwrap();
                sample1_allele_depth = allele_depth_vec[0].iter().map(|x| x.clone()).collect();
                sample2_allele_depth = allele_depth_vec[1].iter().map(|x| x.clone()).collect();
            }
            let chr = String::from_utf8_lossy(vcf_header.rid2name(
                      record.rid().expect("Error read rid.")).unwrap()).to_string();
            let pos = record.pos() + 1; //convert 0 base to 1 base
            let normal_depth = sample1_allele_depth[0] + sample1_allele_depth[1];
            let tumor_depth = sample2_allele_depth[0] + sample2_allele_depth[1];
            // is a good heterogeneous SNP site in normal sample?
            // multi-sample snp calling maybe phase genotype
            if (sample_1_genotype == "0/1" || sample_1_genotype == "0|1" ||
                sample_1_genotype == "1|0") && normal_depth > (self.min_coverage as i32)
                && normal_depth < (self.max_coverage as i32) {
                snp_summary.no_of_good_hets_in_normal += 1;
            } else {
                continue;
            }
            // is a good SNP site in tumor sample?
            // Note: don't filter homologous SNP sites in tumor sample
            if tumor_depth > (self.min_coverage as i32)&& tumor_depth < (self.max_coverage as i32){
                snp_summary.no_of_good_hets += 1;
            } else {
                continue;
            }
            // Have passed all filters, it is a good heterogeneous SNP site
            let index: usize = snp_list.len();
            snp_list.insert(index, Record{chr: chr.clone(), pos: pos as u32,
                normal_ro: sample1_allele_depth[0],
                normal_ao: sample1_allele_depth[1],
                normal_depth: normal_depth, tumor_ro: sample2_allele_depth[0],
                tumor_ao: sample2_allele_depth[1], tumor_depth: tumor_depth});
            // println!("{}:{} {:?}:{:?}, {:?}:{:?}",chr,pos,sample_1_genotype,sample_2_genotype,
            //           sample_1_AD,sample_2_AD);
        }
        println!("total SNP sites: {}\nGood heterogeneous SNP sites in normal: {}\n\
            Good heterogeneous SNP sites keeped: {}",
            snp_summary.no_of_total_records, snp_summary.no_of_good_hets_in_normal,
            snp_summary.no_of_good_hets);
        let output_f = File::create(&self.output_file_path)
            .expect(&format!("Error in creating output file {:?}", &self.output_file_path));
        let mut gz_writer = flate2::GzBuilder::new()
            .filename(self.output_file_path.file_stem().unwrap().to_str().unwrap())
            .comment("Comment")
            .write(output_f, Compression::default());
        gz_writer.write_fmt(format_args!("#min_coverage={}, max_coverage={}\n",
            self.min_coverage, self.max_coverage)).unwrap();
        gz_writer.write_fmt(format_args!("#two sample snp file: {:?}\n", 
            &self.snp_file)).unwrap();
        gz_writer.write_fmt(format_args!("#no_of_total_records: {}\n", 
            snp_summary.no_of_total_records)).unwrap();
        gz_writer.write_fmt(format_args!("#no_of_good_hets in normal: {}\n", 
            snp_summary.no_of_good_hets_in_normal)).unwrap();
        gz_writer.write_fmt(format_args!("#no_of_good hets in two samples: {}\n", 
            snp_summary.no_of_good_hets)).unwrap();
        gz_writer.write_fmt(format_args!("chr\tpos\ttumor_depth\ttumor_ro\t\
                tumor_ao\tnormal_depth\tnormal_ro\tnormal_ao\n")).unwrap();

        for snp_record in snp_list {
            gz_writer.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                snp_record.chr, snp_record.pos, snp_record.tumor_depth,
                snp_record.tumor_ro, snp_record.tumor_ao,
                snp_record.normal_depth, snp_record.normal_ro,
                snp_record.normal_ao)).unwrap();
        }
        gz_writer.finish()
            .expect(&format!("ERROR finish() failure for gz_writer of {:?}.", 
                &self.output_file_path));
        println_stderr!("{} intersect SNPs.", snp_summary.no_of_good_hets);
    }

    pub fn run(&self){
        self.select_het_snp();
    }
}
