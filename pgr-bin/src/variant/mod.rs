pub mod annotate_bed;
pub mod annotate_vcf;
pub mod diploid_vcf;
pub mod merge_sv;
pub mod sv_analysis;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct VariantCmd {
    #[command(subcommand)]
    pub cmd: VariantSubCmd,
}

#[derive(Subcommand)]
pub enum VariantSubCmd {
    /// Generate a diploid VCF from paired haplotype alignment maps
    DiploidVcf(diploid_vcf::Args),
    /// Analyse SV candidates with principal bundle decomposition
    SvAnalysis(sv_analysis::Args),
    /// Merge SV candidate BED records from multiple haplotypes
    MergeSv(merge_sv::Args),
    /// Annotate BED regions with gene annotation features
    AnnotateBed(annotate_bed::Args),
    /// Annotate VCF variants with gene names from a GTF annotation
    AnnotateVcf(annotate_vcf::Args),
}

pub fn run(cmd: VariantCmd) {
    match cmd.cmd {
        VariantSubCmd::DiploidVcf(args) => diploid_vcf::run(args).expect("diploid-vcf failed"),
        VariantSubCmd::SvAnalysis(args) => sv_analysis::run(args).expect("sv-analysis failed"),
        VariantSubCmd::MergeSv(args) => merge_sv::run(args),
        VariantSubCmd::AnnotateBed(args) => annotate_bed::run(args).expect("annotate-bed failed"),
        VariantSubCmd::AnnotateVcf(args) => annotate_vcf::run(args).expect("annotate-vcf failed"),
    }
}
