# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3]
### Added
- Support multi individual input VCFs for calling
### Changed
- Output FORMAT tags to allow sample resolution on multi individual VCFs

## [0.2]
### Added
- Build reference, use prebuilt reference
### Fixed
- EDR for 1/1 variants (have ref, but distance should be counted with only alt)
- Set SD to 0 for GT .

## [0.1.0]
### Added
- Changelog, pyproject.toml
- Simple CLI with version argument
- Simple drop calling, filter on edit distance, and overall counts or case average expression
- Cutoff parameter flags, XY flag
- VCF file writing
### Changed
- Moved calculation and annotation to own package
- Updated README


