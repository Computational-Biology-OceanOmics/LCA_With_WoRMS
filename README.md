# Calculate LCA using WoRMS

Most LCA pipelines out there rely on NCBI's Taxonomy database. However, that database is large and fish are notoriously prone to change. Other databases like the World Register of Marine Species (WoRMS) are updated more often, so here is a tool that takes a table of BLAST results, takes the hits for each ASV, and queries the WoRMS API to ask for the 'updated'/'current' lineage of the hit, then uses the WoRMS lineages to calculate LCAs. Sometimes several similar species have wrong, outdated families on NCBI Taxonomy that are correct on WoRMS, and when you use NCBI Taxonomy, you will get an order-level LCA, while WoRMS-based LCA will have a family or genus-level LCA.

Be careful, though: WoRMS includes only *marine* species and will return nothing for non-marine species, so if you have a mix of marine and non-marine species in your results better use a different database. We work with marine eDNA so filtering out non-marine hits is very useful for us, lots of similar looking freshwater fish out there!

Here's what LCA levels look like between NCBI Taxonomy ('normal') and WoRMS for one of our internal libraries, thanks to Dr Georgia Nester for the comparison:

![Rplot](https://github.com/user-attachments/assets/781cb094-a06b-42a5-8c61-b36efb150455)

## LCA calculation

The LCA calculation works the same way as [eDNAFlow](https://github.com/mahsa-mousavi/eDNAFlow)'s LCA calculation. Given a group of potential species for an ASV, take the species with the highest identity, subtract 1 from the identity, and then include all species above that cutoff in the LCA. The LCA itself is just a grouping: if there are several species within the cutoff, then the species is set to 'dropped' and we go up one taxonomic level. There's no LCA voting or similar, though that's not hard to add.

## Input

The input is blast-output, tabular, using this output format:

     -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"

## Usage

```
usage: calculateLCAWithWoRMS.py [-h] -f FILE -o OUTPUT [--cutoff CUTOFF] [--pident PIDENT] [--missing_out MISSING_OUT]

Parses a BLAST-tabular output file and produces LCAs by asking the WoRMS API for each hit's lineage. BLAST formatting assumed is: -outfmt "6
qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue
bitscore qcovs qcovhsp"

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  input file of BLAST results
  -o OUTPUT, --output OUTPUT
                        Output file of LCAs. Tab delimited.
  --cutoff CUTOFF       OPTIONAL: Percentage cutoff between best BLAST hit and followup to be considered in LCA. Only species within this
                        percentage identity cutoff will be included in LCA calculation. Default: 1 in line with eDNAFlow's LCA script.
  --pident PIDENT       OPTIONAL: Percentage cutoff for BLAST hits. Hits below this cutoff will be ignored for LCA calculation. Default:
                        Initially consider all BLAST hits, but then filter in line with --cutoff.
  --missing_out MISSING_OUT
                        OPTIONAL: Filename to write missing species (not in WoRMS) to. Default: missing.csv.
```

--cutoff changes how lenient the LCA calculation it is - the larger the cutoff, the more species are included in an LCA.

--pident changes how the BLAST results are parsed, hits below that cutoff will never make it into the LCA.

--missing_out changes the filename of the file missing species are written to, by default 'missing.csv'. Missing species are species present in the BLAST results and on NCBI Taxonomy, but are not in WoRMS. These species are a mix of non-marine species and species with sp./cf. names.

## Output

Looks like this:

```
ASV_name        Class   Order   Family  Genus   Species PercentageID    Species_In_LCA
ASV_17067       Teleostei       Ophidiiformes   Ophidiidae      dropped dropped 89.60   Ventichthys biospeedoi, Bassozetus zenkevitchi
```

A tab-delimited table, one row per unique query in the BLAST results, showing which WoRMS taxonomic levels were included, and which were dropped. It also shows the average BLAST identity of the species-hits included in the LCA, and the species that were included in the LCA.

## Installation

This depends on pytaxonkit (only on bioconda) and pyworms. Follow the [pytaxonkit](https://github.com/bioforensics/pytaxonkit) instructions to download the NCBI taxonomy database.

    conda install -c bioconda pytaxonkit
    pip install pyworms

## FAQ

- Why does this depend on pytaxonkit? I thought it ignores NCBI Taxonomy?

Yes, but I have not found a nice way to pull out the species name from a row of BLAST hits. I have played with [taxonerd](https://github.com/nleguillarme/taxonerd) but it feels like overkill for this problem. Relying on the NCBI Taxonomy ID, then pulling out the WoRMS lineage for that taxonomy ID is the best middle-ground in my opinion. I *could* base this on the sscinames column produced by BLAST, but when BLAST does not have the NCBI taxonomy database present when running that column will be all NA.

- How long does this run for?

In my tests for a whole eDNA dataset with a few thousand ASVs, about 15 minutes. The WoRMS API accepts only 500 species at a time so it takes most of the running time to pull down all lineages.

- You could make this faster using a dump of the WoRMS database!?

Yes but then I have to distribute that database dump with the code, and that dump will be outdated within a week, but others will find the database dump on Github and base their own scripts on the outdated dump and I'd rather not have that happen.

- Can this take other taxonomic databases?

I have thought about adding AFD (Australian Faunal Directory) as another source. Maybe later.

- I have weird taxonomic sub-levels that I am interested in, why are they dropped?

The script only works with class, order, family, genus, species. That way I can make the LCA calculation fairly lazy: instead of having to build a graph of all taxonomic levels and doing weird graph-based magic, I can just calculate the sets of unique species, unique genera, unique families, unique orders, unique classes, and check set size after filtering. If the set size is > 1, set this taxonomic level's taxonomic label to 'dropped'. Easier than breaking my head trying to come up with recursive tree-walking algorithms for the sake of methodological complexity you can publish in a paper, I'd rather have results.

- My BLAST results have NAs for taxonomy IDs so everything breaks.

Your BLAST database was not made using taxonomy IDs. Make a table of sequence IDs and their associated taxonomy IDs, then run this:

    makeblastdb -dbtype nucl -in your_database.fasta -parse_seqids -taxid_map your_database.taxids.txt

- I have more questions!

Please contact me at pbayer AT minderoo.org
