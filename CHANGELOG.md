# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## ToDo
- Creation of sample lists from the manifest?
- Make a variant release? 
- BCFtools docker container reversion on secure-lustre to minimise differences (`bcftools:1.20` ->`bcftools:1.9`) Tabix requirement makes life difficult

# 0.5.0 [September 18th 2024]
### Added
- A CI running on secure lustre for nf-test
### Changed
- Tests updated for the new CI.

# 0.4.0 [September 11th 2024]
### Added
- Fixes to the container conflicts caused by QC/MAF cross calling 
- Added support for alternative canonical transcript labels 
- Corrected to use same module call whilst on farm22 as previously

# 0.3.0 [July 2024]
### Added
- Support for tsv -> .xlsx conversion
- Optional exectuion of independent/one tumor per patient/all sample analyses

# 0.2.0 [July 2024]
### Added
- Calculate and output tumour mutation burden
- Harmonised script name inside and outside of container (MAF) added here

# 0.1.1 [July 2024]
### Changed
- Update readme for reliable Farm singularity

# 0.1.0 [June 2024]
### Changed
- BCFtools container changed to `bcftools:1.20` for dockerising steps.

### Added 
Initial version of this pipeline. Support for 
- Basic variant filtering and adding annotations for dbSNP155 common variants (Caveman & Pindel)
- Make MAF files and run variant call QC

