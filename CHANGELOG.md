# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased] - 2025-05-26
- this is a major reorganiztion based off initial work in the GRD-671 branch to remove co-cleaning and indel realignment, both of which are no longer part of gatk best practices
- [GRD-1000](https://jira.oicr.on.ca/browse/GRD-1000)
- implementation of a revised parallelization strategy to address issues identified in duplicate marking. Chromosome intervals subset only with paired reads on the same contig.  Additional intervals are identified by keywords NC : subsets with all non-conanonical chromosomes in the build ( un,random,alt,fix), with paired reads on the same contig.  SPLIT : subsets paired reads on different contigs UNALIGNED : subsets paired reads both unaligned 
- subsetAndFilter task introduced, which divides by interval and filters, replacing perProcessing, using interval strategy described above
- BQSR run on per interval, per input bam subset (previously on the merged within interval)
- withinInterval, crossInput merges either through markduplicates or meregeBams (as before)
- splitNCigarString requires libType="rna" parameter.  This follows duplicate marking
- BQSR Application on perInterval merged bams
- doBamMetrics option will run additional task bamMetrics to generated samtools stats + flagstats, + chrpair counts for ALL input, intermediate, and final bams.  This is useful for troubleshooting, but should not be used in production
- moved memory coefficient preparation to an initial task, prepareIntervals (run once), instead of running within intervals. 
- included "comment" tag, that will insert a comment into the header of the final merged bam, to indicated process for callready bam creation
- general cleanup and reorganization for consistency, readability and comprehenshion
- NOTE this version is not yet tagged.   It needs to incorporate some of the changes that had vollowed GRD-671

## [2.2.1] - 2025-05-26
- Re-deployment to enable labels for optional outputs
- [GRD-948](https://jira.oicr.on.ca/browse/GRD-948)

## [2.2.0] - 2024-06-25
### Added
- Add vidarr labels to outputs (changes to medata only).
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797).

## [2.1.1] - 2023-08-18
### Fixed
- Fixing typo for mouse assembly ID.

## [2.1.0] - 2023-07-04
### Changed
- Moved assembly-specific code in olive to wdl.

## [2.0.5t] - 2023-06-14
### Added
- bamMergePreprocessing_noMarkDup alias
- bamMergePreprocessing_merge_only alias
 
## [2.0.5] - 2022-11-07
### Added
- (this is probably an extra tag, no changes from the previous version)
- created bmpp_by_sample alias
  
## [2.0.4] - 2022-11-07
### Added
- created bmpp_by_sample alias
 
## [2.0.3] - 2022-07-26
### Changed
- [GP-3049](https://jira.oicr.on.ca/browse/GP-3049) removed preprocessingBamRuntimeAttributes

## [2.0.2] - 2020-05-31
### Changed
- Migrate to Vidarr.

## [2.0.0] - Unreleased
### Changed
- Reimplement workflow in WDL.

## [1.3.0] - 2019-10-29
### Changed
- Configurable index_mem and picard_mem_overhead_mb.
- Update to picard 2.21.2
- Update to pipedev 2.5

## [1.2.0] - 2017-07-14
### Added
- [GP-1202](https://jira.oicr.on.ca/browse/GP-1202) - Support RNAseq samples.

## [1.1.0] - 2017-04-17
### Added
- [GP-1106](https://jira.oicr.on.ca/browse/GP-1106)
- Support for performing GATK Indel Realignment together on groups of inputs, but provisioning out by the separate groups.
- Support grouping of "chr_sizes" intervals.

## [1.0.3] - 2017-03-28
### Added
- Support for "unmapped" input.
- [GP-1117](https://jira.oicr.on.ca/browse/GP-1117)] 

## [1.0.2] - 2017-03-13
### Changed
- Update to Seqware 1.1.1-gsi.

## [1.0.1] - 2016-01-27
### Changed
- Fixing metadata attachment, removed do_filter parameter (redundancy to do_sam_filter).

## [1.0.0] - 2015-12-20
### Added
- Initial import. 
- [GP-604](https://jira.oicr.on.ca/browse/GP-604) 

