use std::fs::File;
use std::io::prelude::*;

const LEN: usize = 22;

static CHR_LIST: [&str;22] = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
    "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
    "chr21","chr22"];

struct Segment {
    chr_id: u32,
    start: u64,
    end: u64,
    copynumber: f64,
    major_cp: f64,
    is_clonal: bool,
}

struct Chromosomes {
    segment_list: Vec<Segment>,
}

struct GenomicSegments {
    chromosomes_list: Vec<Chromosomes>,
}

impl Segment {
    fn new(chr_id: u32, start: u64, end: u64, copynumber: f64, 
           major_cp: f64, is_clonal: bool) -> Segment {
        Segment { chr_id, start, end, copynumber, major_cp, is_clonal}
    }
}

impl Chromosomes {
    fn add(&mut self, seg: Segment) {
        let mut index = 0;
        for each_seg in self.segment_list.iter() {
            if each_seg.start > seg.start && each_seg.end > seg.end {
                break;
            }
            index += 1;
        }
        self.segment_list.insert(index,seg);
    }

    fn len(&self) -> u32 {
        self.segment_list.len() as u32
    }
}

impl GenomicSegments {
    fn new() -> GenomicSegments {
        let chromosomes_list: Vec<Chromosomes> = vec![];
        let mut genoms = GenomicSegments{chromosomes_list};
        for _ in 0..22 {
            let segment_list: Vec<Segment> = vec![];
            genoms.chromosomes_list.push(Chromosomes { segment_list });
        }
        genoms
    }

    fn add(&mut self, seg: Segment) {
        self.chromosomes_list[seg.chr_id as usize - 1].add(seg);
    }
}

pub struct RecallPrecision {
    resultfile: String,
    accurityfile: String,
    outputfile: String,
}

impl RecallPrecision {
    pub fn new(resultfile: &str, accurityfile: &str, outputfile: &str) -> RecallPrecision {
        RecallPrecision {
            resultfile: resultfile.to_string(),
            accurityfile: accurityfile.to_string(),
            outputfile: outputfile.to_string(),
        }
    }

    pub fn run(&self) {
        let mut actual = GenomicSegments::new();
        let mut accurity = GenomicSegments::new();

        read_from_actual_result(&self.resultfile, &mut actual);
        read_from_accurity_result(&self.accurityfile, &mut accurity);
        let (cp_recall, cp_precision, mcp_recall, mcp_precision) = call_recall_and_precision(&actual, &accurity);
        println!("recall of copy number: {:.*}, precision of copy number: {:.*}", 
                  6, cp_recall, 6, cp_precision);
        println!("recall of major allele copy number: {:.*}, precision of major allele copy number: {:.*}", 
                  6, mcp_recall, 6, mcp_precision);

        let mut outputfile = File::create(&self.outputfile).expect("Cannot create file");
        outputfile.write(b"cp_recall\tcp_precision\tmcp_recall\tmcp_precision\n").unwrap();
        outputfile.write_fmt(format_args!("{:.*}\t{:.*}\t{:.*}\t{:.*}\n", 
                             6, cp_recall, 6, cp_precision, 6, mcp_recall, 6, mcp_precision)).unwrap();
    }
}

fn print_genomic_segment(gs: &GenomicSegments) {
    for chromosome in gs.chromosomes_list.iter() {
        for segment in chromosome.segment_list.iter() {
            println!("chr{}: {}--{} {}_{} {}", segment.chr_id, segment.start, 
                     segment.end, segment.copynumber, segment.major_cp, segment.is_clonal);
        }
    }
}

fn read_from_actual_result(resultfile: &String, actual: &mut GenomicSegments)
{
    let mut f = File::open(resultfile).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents).unwrap();

    let mut isheader = true;
    for line in contents.lines() {
        // skip comments
        if line.trim().starts_with("#") {
            continue;
        }
        // skip header
        if isheader {
            isheader = false;
            continue;
        }
        let mut data_list = line.split('\t');
        let chr = data_list.next().unwrap();
        if CHR_LIST.contains(&&chr) == false {
            continue;
        }
        let chr_id: u32 = chr[3..].parse().unwrap();
        let start: u64 = data_list.next().unwrap().parse().unwrap();
        let end: u64 = data_list.next().unwrap().parse().unwrap();
        let tmp1 = data_list.next().unwrap().parse::<f64>();
        let copynumber: f64;
        match tmp1 {
            Ok(data) => {
                copynumber = data;
            },
            Err(_error) => {
                continue;
            }
        }
        let tmp2 = data_list.next().unwrap().parse::<f64>();
        let major_cp: f64;
        let is_clonal: bool;
        match tmp2 {
            Ok(data) => {
                major_cp = data;
                is_clonal = data_list.next().unwrap().parse::<String>().unwrap() == "T";
            },
            Err(_error) => { 
                // set is_clonal to false, then the segment will not be
                // considered when calculate recall and precsion for
                // major allele copy number
                major_cp = -1_f64;
                is_clonal = false;
            }
        }
        let seg = Segment::new(chr_id, start, end, copynumber, 
                            major_cp, is_clonal);
        actual.add(seg);
    }
}

fn read_from_accurity_result(accurityfile: &String, accurity: &mut GenomicSegments)
{
    let mut f = File::open(accurityfile).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents).unwrap();

    let mut isheader = true;
    for line in contents.lines() {
        if isheader {
            isheader = false;
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let mut data_list = line.split('\t');
        let chr_id: u32 = data_list.next().unwrap().parse().unwrap();
        let start: u64 = data_list.next().unwrap().parse().unwrap();
        let end: u64 = data_list.next().unwrap().parse().unwrap();
        let copynumber: f64 = data_list.next().unwrap().parse().unwrap();
        let major_cp: f64;
        let is_clonal: bool;
        let tmp = data_list.next().unwrap().parse::<f64>();
        match tmp {
            Ok(data) => {
                major_cp = data;
                is_clonal = true;
            },
            Err(_error) => {
                println!("chr{}: {}-{} is a subclonal segment", chr_id, start, end);
                major_cp = -1_f64;
                is_clonal = false;
            }, 
        };
        if copynumber == 2.0 && major_cp == 1.0{
            continue;
        }
        let seg = Segment { chr_id, start, end, copynumber, major_cp, is_clonal};
        accurity.add(seg);
    }
}

fn call_recall_and_precision(actual: &GenomicSegments, accurity: &GenomicSegments) -> (f64,f64,f64,f64) {
    let mut xony = GenomicSegments::new();
    y_segments_map_on_x(&actual, &accurity, &mut xony, LEN);
    let (cp_similarity, mcp_similarity) = similarity_of_x_to_y(&xony, LEN);
    // println!("xony");
    // print_genomic_segment(&xony);
    // println!("actual");
    // print_genomic_segment(&actual);
    // println!("estimate");
    // print_genomic_segment(&accurity);

    // println!("similarity: {}",similarity);
    let (actual_cp_sum, actual_mcp_sum) = sum_of_area(&actual, LEN);
    let (accurity_cp_sum, accurity_mcp_sum) = sum_of_area(&accurity, LEN);
    let cp_recall = cp_similarity / actual_cp_sum;
    let cp_precision = cp_similarity / accurity_cp_sum;
    let mcp_recall = mcp_similarity / actual_mcp_sum;
    let mcp_precision = mcp_similarity / accurity_mcp_sum;
    (cp_recall, cp_precision, mcp_recall, mcp_precision)
}

fn max(x: u64, y: u64) -> u64 {
    if x > y {
        x
    } else {
        y
    }
}

fn min(x: u64, y: u64) -> u64 {
    if x > y {
        y
    } else {
        x
    }
}

fn y_segments_map_on_x(x: &GenomicSegments, y: &GenomicSegments, xony: &mut GenomicSegments, length: usize) {
    let mut copynumber_diff: f64;
    for chr_index in 0..length {
        if x.chromosomes_list[chr_index].len() == 0 || y.chromosomes_list[chr_index].len() == 0 {
            continue;
        }
        for y_segment in y.chromosomes_list[chr_index].segment_list.iter() {
            for x_segment in x.chromosomes_list[chr_index].segment_list.iter() {
                let start = max(y_segment.start, x_segment.start);
                let end = min(y_segment.end, x_segment.end);
                let major_cp: f64;
                let is_clonal: bool;
                if y_segment.is_clonal == true && x_segment.is_clonal == true {
                    major_cp = (x_segment.major_cp - y_segment.major_cp).abs();
                    is_clonal = true;
                } else {
                    major_cp = -1_f64;
                    is_clonal = false;
                }
                //Note: for absolute copy number, cp=2, mcp!=1 is negative result
                //becuase we consider cp only. So we should exclude those segments
                //when calculate inner set. But for simplicy in mcp calculating, we
                //set the copynumber_diff to a large number X, so exp(-X) = 0.
                if x_segment.copynumber == 2.0 {
                    copynumber_diff = 10000f64;
                } else {
                    copynumber_diff = (x_segment.copynumber - y_segment.copynumber).abs();
                }
                if start < end + 1 {
                    let seg = Segment {
                        chr_id: chr_index as u32 + 1,
                        start,
                        end,
                        copynumber: copynumber_diff,
                        major_cp,
                        is_clonal
                    };
                    xony.add(seg);
                }
            }
        }
    }
}

fn similarity_of_x_to_y(xony: &GenomicSegments, length: usize) -> (f64, f64) {
    let mut cp_similarity: f64 = 0.0;
    let mut mcp_similarity: f64 = 0.0;
    for chr_idx in 0..length {
        if xony.chromosomes_list[chr_idx].len() == 0 {
            continue;
        }
        for segment in xony.chromosomes_list[chr_idx].segment_list.iter() {
            cp_similarity += (segment.end - segment.start + 1) as f64 * (0.0 - segment.copynumber).exp();
            if segment.is_clonal == true {
                mcp_similarity += (segment.end - segment.start + 1) as f64 * (0.0 - segment.major_cp).exp();
            }
        }
    }
    return (cp_similarity, mcp_similarity);
}

fn sum_of_area(gs: &GenomicSegments, length: usize) -> (f64, f64) {
    let mut cp_sum: u64 = 0;
    let mut mcp_sum: u64 = 0;
    for chr_idx in 0..length {
        if gs.chromosomes_list[chr_idx].len() == 0 {
            continue;
        }
        for segments in gs.chromosomes_list[chr_idx].segment_list.iter() {
            if segments.copynumber != 2_f64 {
                cp_sum += segments.end - segments.start + 1;
            }
            if segments.is_clonal == true {
                mcp_sum += segments.end - segments.start + 1;
            }
        }
    }
    return (cp_sum as f64, mcp_sum as f64);
}
