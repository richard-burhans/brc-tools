# Low-complexity and tandem repeat soft-masking

This workflow soft-masks (lower-cases) low-complexity sequence, tandem repeats and
over-represented sequence in a list of genome assemblies. It runs four detectors on
each assembly and masks the union of their intervals:

- **NCBI dustmasker** (DUST): low-complexity sequence.
- **tantan**: low-complexity sequence and short tandem repeats.
- **FasTAN**: tandem arrays.
- **NCBI WindowMasker**: over-represented sequence, from k-mer counts of the assembly itself.

The workflow first removes any soft-masking that an assembly already has. Thus the
output mask comes only from these four detectors. Each detector uses its tool defaults.
The workflow has no inputs for detector parameters.

This workflow does not use a repeat library to mask interspersed repeats (transposable
elements). For that, use the IWC workflow
[Repeat masking with RepeatModeler and RepeatMasker](https://github.com/galaxyproject/iwc/tree/main/workflows/repeatmasking).
If you want both kinds of mask, run both workflows.

## Inputs

- **Genome assemblies**: a list collection of genome assembly FASTA files, with one
  element for each assembly. The element identifier (for example, a strain name)
  names each per-assembly output and its row in the percent masked table.
  Every step maps over this collection.

## Outputs

All outputs except the table are list collections with one element for each assembly.

- **Soft-masked assemblies**: FASTA files. In each file, the union of the intervals from
  all four detectors is lower-cased (bedtools maskfasta `-soft`).
- **Soft-masked assembly index**: the samtools faidx index (`.fai`) of each soft-masked
  FASTA.
- **dustmasker intervals**, **WindowMasker intervals**, **tantan intervals**,
  **FasTAN intervals**: BED6 interval tracks, one for each detector. You can load
  them as UCSC hub tracks.
- **Percent masked table**: one row for each assembly. It has one column for each
  detector and one for the union. Each value is the percentage of the assembly that
  the detector or the union masks (bedtools genomecov).

## How it works

1. Upper-case each assembly to remove existing soft-masking.
2. Run dustmasker, WindowMasker, tantan and FasTAN on the upper-cased assembly.
3. Concatenate the four BED files. Sort them (bedtools sort) and merge them
   (bedtools merge) into one set of union intervals.
4. Soft-mask the upper-cased assembly at the union intervals (bedtools maskfasta).
   Index the result (samtools faidx).
5. Compute the genome coverage of each detector's intervals and of the union
   (bedtools genomecov). Pivot the results into the percent masked table.

## Related resources

We know of no tutorial for this workflow yet.
