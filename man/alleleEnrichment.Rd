% File man/alleleEnrichment.Rd
\name{alleleEnrichment}
\alias{alleleEnrichment}
\title{Select a Subsample with Optimized Allele Frequencies from Genotyped Individuals}

\description{
The function selects a user-defined number of samples from a larger cohort so that predefined alleles are guaranteed to be contained at a maximum frequency. It is also possible to require a balanced selection, matching the frequencies in the full cohort. Takes an arbitrary number of independent SNPs for which alleles are to be maximized, and can optionally also enrich phenotypes (have to be linearly scaled integer values). 
}

\usage{
alleleEnrichment(
  subsample.size,
  preferred.alleles,
  map, 
  ped, 
  sampleid.columnindex = 2,
  mode = "max",
  pheno.file = NULL,
  na.samples = na.omit, 
  na.snps = na.pass,
  write.resultfiles = TRUE, 
  ...
)
}

\arguments{
  \item{preferred.alleles}{list. Names have to match the SNPs that should be enriched with identifiers as in the map file. Values have to be allele identifiers from the ped file, specifying which allele to enrich. Make sure that the allele identifier match those in the pedfile (case sensitive). This is not being checked (but can be seen when the number of enriched alleles is 0 in the results). The nucleotides specified should always be uppercase A,C,G or T.}
  \item{mode}{character or list. Can be either one of "max" or "balanced", or a named list of those, where the names correspond to all SNP or phenotype names to enrich. Such a mixing of constraints is still experimental: 'max' constraints are always correctly maximized but balanced constraints might not be optimally balanced.}
  \item{subsample.size}{integer. The target number of samples to select.}
  \item{ped}{character. Path/filename to a pedfile as specified in \url{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped}. Reads various formats (numeric or character nucleotides, delimiters can be either space only or tab with spaces between alleles, ...). May contain genotypes for more SNPs than needed, and arbitrary additional columns ahead of genotype data (i.e. sex, family data, ...).}
  \item{map}{character. Path/filename to a mapfile as specified in \url{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map}. Reads default linkage map format as well as --map3 (CHR, SNP, BP) format.}
  \item{pheno.file}{character. Path/filename to a phenotype file. Format has to be readable by \code{link{read.table}}. Expects a header line with phentoype names. Expects sample IDs in the column at index \code{sampleid.columnindex} that correspond to those in the pedfile. All columns beyond \code{sampleid.columnindex} are read as phenotype variables to be included in the analysis. Make sure that these phenotype variabes are linear and contain only positive integers (noninteger phenotypes are rounded).}
  \item{sampleid.columnindex}{integer. The column index in the ped file where sample IDs are listed. This is also used for the phenotype file if supplied.}
  \item{na.samples}{function. Action to take for NA values in sample columns of source data (e.g. na.omit will remove that sample column if it contains NAs). When neither \code{na.samples} nor \code{na.snps} do remove all NA values, remaining NAs are replaced with a random number between 0 and 1 (not recommended).}
  \item{na.snps}{function. As above, but applied to SNP/phenotype rows (e.g. na.omit will remove that SNP/phenotype column if it contains NAs).}
  \item{write.resultfiles}{boolean, Determines whether to write resulting samples selection, its genotypes and some statistics to outfiles.}
  \item{...}{can be num.perms, search.seconds and weights - see examples.}
}

\details{
Detection of rare causal variants in GWAS-associated loci can be done by resequencing all significantly associated loci at a time for a discovery subsample of the GWAS cohort. Driven by the hypothesis that multiple sought-for causal rare or de novo mutations occur on the same haplotype as the common risk allele, this function provides means to select a subsample for sequencing with maximized occurences of risk alleles for all loci of interest. This way, it can be made sure to have a greater chance of finding rare variation that resides on the same haplotype as the GWAS risk allele. When an assumption about the risk allele being protective or risk-conferrring can not be made (e.g. unknown function of the underlying gene), it is also possible to generate a sample selection with balanced distribution of alleles in relation to the original frequency. This can then reduce the risk of an accidentially misbalanced selection with regard to the original allele frequencies in the entirety of cases especially for small subsamples. 
}

\value{
A list with the elements 
\describe{
 \item{sample.selection}{A vector of sample IDs that have been selected after allele count optimization}
 \item{allele.counts}{A data frame with the first column listing SNP IDs, the second the achieved allele counts in the subsample.}
\item{enriched.percent}{Allele / phenotype frequency in relation to the frequency in the full cohort (i.e. 100 percent means the mean value in the subsample matches the cohort mean.)}
 \item{genotypes}{a ped-like matrix with allele status for the selected samples for all SNPs}
}
In rare cases, less samples than required can be returned, meaning that the frequency constraints have been fulfilled with fewer samples. It is then recommended to just add an additional random sample.
}

\references{
  \url{http://lpsolve.sourceforge.net/}
}

\examples{

\dontshow{
  ## path to example data
  pedmap.files <- system.file("extdata", c("example.map", "example.ped", "example.pheno"), package = "boostSeq")
  gwas.resultfile <- system.file("extdata", "example.gwas", package = "boostSeq")
}

## select 48 samples for 2 SNPs (maximize risk alleles)
alleleEnrichment(
  preferred.alleles = list(rs7960808 = "T", rs8411 = "C"),
  mode = "max",
  subsample.size = 48, 
  map = pedmap.files[1], 
  ped = pedmap.files[2],
  pheno.file = NULL, 
  write.resultfiles = FALSE
)

# now include phentoypes
# also, mix maximized and balanced constraints 
# (sex balanced, exerything else maxed)
# mixture of modes is a bit experimental but usually works 
# (i.e. in the worst case the balanced variable is in fact not 
# balanced in the results, i.e. as if not included)
alleleEnrichment(
  preferred.alleles = list(rs7960808 = "T", rs8411 = "C"),
  mode = list(rs7960808 = "max", rs8411 = "max", SEX = "balanced", FVIII = "max"),
  subsample.size = 48, 
  map = pedmap.files[1], 
  ped = pedmap.files[2],
  pheno.file = pedmap.files[3], 
  write.resultfiles = FALSE
)


# did not enrich much risk alleles: reduce weight of FVIII phenotype
# remark: weights are not considered for balanced variables 
# in the current implementation (but have to be included)
alleleEnrichment(
  preferred.alleles = list(rs7960808 = "T", rs8411 = "C"),
  mode = list(rs7960808 = "max", rs8411 = "max", SEX = "balanced", FVIII = "max"),
  subsample.size = 48, 
  map = pedmap.files[1], 
  ped = pedmap.files[2],
  pheno.file = pedmap.files[3], 
  write.resultfiles = FALSE,
  weights = list(rs8411 = 1, rs7960808 = 1, SEX = 1, FVIII = 0.5)
)

# define weights for SNPs by assoviation log(p) from GWAS
weights.snps <- getSnpWeights(snps = c("rs8411", "rs7960808"), gwas.resultfile)
alleleEnrichment(
  preferred.alleles = list(rs7960808 = "T", rs8411 = "C"),
  mode = list(rs7960808 = "max", rs8411 = "max", SEX = "balanced", FVIII = "max"),
  subsample.size = 48, 
  map = pedmap.files[1], 
  ped = pedmap.files[2],
  pheno.file = pedmap.files[3], 
  write.resultfiles = FALSE,
  weights = c(weights.snps, SEX = 1, FVIII = 0.5)
)

# compare with balanced enrichment, using 50000 random selections:
alleleEnrichment(
  preferred.alleles = list(rs7960808 = "T", rs8411 = "C"),
  mode = "balanced",
  subsample.size = 48, 
  map = pedmap.files[1], 
  ped = pedmap.files[2],
  pheno.file = pedmap.files[3], 
  write.resultfiles = FALSE,
  num.perms = 50000
)

}
