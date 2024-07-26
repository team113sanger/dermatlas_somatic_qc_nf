# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Todo
- Support for tsv -> .xlsx conversion is still required
- Make a variant release

# 0.2.0
### Added
- Calculate and output tumour mutation burden
- Harmonised script name inside and outside of container (MAF) added here

# 0.1.1
### Changed
- Update readme for reliable Farm singularity

# 0.1.0
### Changed
- BCFtools container changed to `bcftools:1.20` for dockerising steps.

### Added 
Initial version of this pipeline. Support for 
- Basic variant filtering and adding annotations for dbSNP155 common variants (Caveman & Pindel)
- Make MAF files and run variant call QC

