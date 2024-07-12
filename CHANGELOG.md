# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.2.0] - 2024-06-25
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797)] - add vidarr labels to outputs (changes to medata only)

## [2.1.1] - 2023-08-18
- Fixing typo for mouse assembly ID

## [2.1.0] - 2023-07-04
- Moved assembly-specific code in olive to wdl

## [2.0.2] - 2020-05-31
- Migrate to Vidarr

## [2.0.0] - Unreleased
- Reimplement workflow in WDL

## [1.3.0] - 2019-10-29
- Configurable index_mem and picard_mem_overhead_mb
- Update to picard 2.21.2
- Update to pipedev 2.5

## [1.2.0] - 2017-07-14
### Added
- [GP-1202](https://jira.oicr.on.ca/browse/GP-1202)] - Support RNAseq samples

## [1.1.0] - 2017-04-17
### Added
- [GP-1106](https://jira.oicr.on.ca/browse/GP-1106)
- Support for performing GATK Indel Realignment together on groups of inputs, but provisioning out by the separate groups
- Support grouping of "chr_sizes" intervals

## [1.0.3] - 2017-03-28
### Added
- [GP-1117](https://jira.oicr.on.ca/browse/GP-1117)] - support for "unmapped" input

## [1.0.2] - 2017-03-13
### Changed
- Update to Seqware 1.1.1-gsi

## [1.0.1] - 2016-01-27
### Changed
- fixing metadata attachment, removed do_filter parameter (redundancy to do_sam_filter)

## [1.0.0] - 2015-12-20
### Added
- [GP-604](https://jira.oicr.on.ca/browse/GP-604) Initial Import

