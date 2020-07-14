extern crate cc;

/*
use std::env;
use std::fs;
use std::path::PathBuf;
*/


fn main() {
    //let out = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    //fs::remove_dir_all(&out).unwrap();
    //fs::create_dir(&out).unwrap();

    cc::Build::new()
        .cpp(true);
        /*
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-lgsl")
        .flag_if_supported("-lgslcblas")
        .file("src_o/infer_class.cpp")
        .file("src_o/BaseGADA.cc")
        .file("src_o/read_para.cpp")
        .file("src_o/prob.cpp")
        .compile("libaccu.a");
        */
}
