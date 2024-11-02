from collections import OrderedDict, Counter
import pytaxonkit
import pandas as pd
from argparse import ArgumentParser
import pyworms
import statistics

CUTOFF = 1

def get_lca(entries):
    # given a list of entries -could be several species, or several genera,
    # and their percentage identity,
    # pull out the highest percentage and remove all hits with pident - 1 (CUTOFF)
    # an entry can be a species or a genus, it comes with its percent identity
    sorted_spec = sorted(entries)
    top_spec = sorted_spec[-1]
    bottom_spec = sorted_spec[0]
    top_perc = top_spec[0]
    perc_range = top_perc - cutoff
    new_spec = set()
    all_percentages = []
    for s in sorted_spec:
        perc, entry = s
        #print(perc, entry)
        if perc >= perc_range:
            new_spec.add(entry)
            all_percentages.append(perc)

    lca_perc = statistics.mean(all_percentages)
    # do we have only one species/genus/family/order?
    if len(new_spec) == 1:
        # easy!
        lca = list(new_spec)[0]
    else:
        # we have several - go up one level
        lca = 'dropped'
    return lca_perc, lca, new_spec


parser = ArgumentParser(description = 'Parses a BLAST-tabular output file and produces LCAs by asking the WoRMS API for each hit\'s lineage. BLAST formatting assumed is: -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"')
parser.add_argument('-f', '--file', help = 'input file of BLAST results', required = True)
parser.add_argument('-o', '--output', help = 'Output file of LCAs. Tab delimited.', required = True)
parser.add_argument('--cutoff', help = 'OPTIONAL: Percentage cutoff between best BLAST hit and followup to be considered in LCA. Only species within this percentage identity cutoff will be included in LCA calculation.\nDefault: %(default)s in line with eDNAFlow\'s LCA script.', default = CUTOFF, type = float)
parser.add_argument('--pident', help = 'OPTIONAL: Percentage cutoff for BLAST hits. Hits below this cutoff will be ignored for LCA calculation.\nDefault: Initially consider all BLAST hits, but then filter in line with --cutoff.', default = 0, type = float)
parser.add_argument('--missing_out', help = 'OPTIONAL: Filename to write missing species (not in WoRMS) to.\nDefault: %(default)s.', default = 'missing.csv')
args = parser.parse_args()

cutoff = args.cutoff
pident_cutoff = args.pident

assert pident_cutoff >= 0 and pident_cutoff <= 100, 'ERROR: Parameter --pident_cutoff must be between 0 and 100'

asv_hits = OrderedDict()
missing_c = Counter()

taxids = set()
for line in open(args.file):
    ll = line.rstrip().split('\t')

    pident = float(ll[6])
    if pident < pident_cutoff:
        continue

    taxid = ll[2]
    assert taxid != 'NA', f'ERROR: You have NAs for taxids in at least one case. Here is the row: {line.rstrip()}'
    taxids.add(taxid)

taxid_to_name = pytaxonkit.name(list(taxids))
taxid_to_name_dict = dict(zip(taxid_to_name.TaxID, taxid_to_name.Name)) 

all_species = list(set(taxid_to_name_dict.values()))

results = pyworms.aphiaRecordsByMatchNames(all_species)
assert len(results) == len(all_species), 'ERROR: Somehow different number of WoRMS records ({len(results)}) to species ({len(all_species)})'

look_up = {}

for species, b in zip(all_species, results):
     
    if not b:
        # not in worms!
        #print(f'did not find {a}')
        continue

    
    #  [{'AphiaID': 218512, 'url': 'https://www.marinespecies.org/aphia.php?p=taxdetails&id=218512', 'scientificname': 'Pristipomoides auricilla', 'authority': '(Jordan, Evermann & Tanaka, 1927)', 'status': 'accepted', 'unacceptreason': None, 'taxonRankID': 220, 'rank': 'Species', 'valid_AphiaID': 218512, 'valid_name': 'Pristipomoides auricilla', 'valid_authority': '(Jordan, Evermann & Tanaka, 1927)', 'parentNameUsageID': 159804, 'kingdom': 'Animalia', 'phylum': 'Chordata', 'class': 'Teleostei', 'order': 'Eupercaria incertae sedis', 'family': 'Lutjanidae', 'genus': 'Pristipomoides', 'citation': 'Froese, R. and D. Pauly. Editors. (2024). FishBase. Pristipomoides auricilla (Jordan, Evermann & Tanaka, 1927). Accessed through: World Register of Marine Species at: https://www.marinespecies.org/aphia.php?p=taxdetails&id=218512 on 2024-10-30', 'lsid': 'urn:lsid:marinespecies.org:taxname:218512', 'isMarine': 1, 'isBrackish': 0, 'isFreshwater': 0, 'isTerrestrial': 0, 'isExtinct': None, 'match_type': 'exact', 'modified': '2008-01-15T17:27:08.177Z'}]
    this_hit = b[0]
    thisclass = this_hit['class']
    order = this_hit['order']
    family = this_hit['family']
    genus = this_hit['genus']
    realspecies = this_hit['valid_name']

    #print(species, this_hit)
    # some BLAST hits are not on the species level. so what  WoRMS retusn is not on the species-level, either
    if species:
        look_up[species] = [ ("C", thisclass),
                            ("O", order),
                            ("F", family),
                            ("G", genus),
                            ("S", realspecies) ]
    if genus and not species:
        look_up[genus] = [ ("C", thisclass),
                            ("O", order),
                            ("F", family),
                            ("G", genus),
                            ("S", '')
                            ]
    if family and not genus:
        look_up[family] = [ ("C", thisclass),
                            ("O", order),
                            ("F", family),
                            ("G", genus),
                            ("S", '')
                            ]
    if order and not family:
        look_up[family] = [ ("C", thisclass),
                            ("O", order),
                            ("F", family),
                            ("G", genus),
                            ("S", '')
                            ]


for line in open(args.file):
    ll = line.rstrip().split('\t')

    pident = float(ll[6])
    if pident < pident_cutoff:
        continue

    taxid = ll[2]
    species = taxid_to_name_dict[int(taxid)]
    try:
        lineage = look_up[species]

    except KeyError:
        #print(f'NOT FOUND {species}')
        missing_c[ (species, 'species') ] += 1
        try:
            genus = species.split(' ')[0]
            lineage = look_up[genus]
        except:
            #print(f'NOT FOUND {genus}')
            missing_c[ (genus, 'genus') ] += 1
            continue

    if ll[0] not in asv_hits:
        asv_hits[ll[0]] = [ (pident, lineage) ]
    else:
        asv_hits[ll[0]].append( (pident, lineage) )

with open(args.missing_out, 'w') as out:
    for c in missing_c:
        name = '\t'.join(c)
        out.write(f'{name}\t{missing_c[c]}\n')

with open(args.output, 'w') as out:
    out.write('ASV_name\tClass\tOrder\tFamily\tGenus\tSpecies\tPercentageID\tSpecies_In_LCA\n')
    for asv_name in asv_hits:
        orf_hits = asv_hits[asv_name]

        # get all the species
        species = set()
        genera = set()
        families = set()
        orders = set()
        classes = set()
        for a in orf_hits:
            # this is one HIT it has all the levels
            pident, lineage = a
            thisclass, thisorder, thisfamily, thisgenus, thisspecies = [x[1] for x in lineage]
            classes.add( (pident, thisclass) )
            orders.add( (pident, thisorder) )
            families.add( (pident, thisfamily) )
            genera.add( (pident, thisgenus) )
            species.add( (pident, thisspecies) )


        #[('C', 'Actinopterygii'), ('O', 'Ophidiiformes'), ('F', 'Ophidiidae'), ('G', 'Ventichthys'), ('S', 'Ventichthys biospeedoi')]
        # oK now we have all the species and genera

        # we need to find the specs that has >1% difference in pident
      
        lca_spec_perc, lca_spec, included_spec = get_lca(species)
        lca_genus_perc, lca_genus, included_genera = get_lca(genera)
        lca_fam_perc, lca_fam, _ = get_lca(families)
        lca_order_perc, lca_order, _ = get_lca(orders)
        lca_class_perc, lca_class, _ = get_lca(classes)
        out.write(f'{asv_name}\t{lca_class}\t{lca_order}\t{lca_fam}\t{lca_genus}\t{lca_spec}\t{lca_spec_perc:.2f}\t{", ".join(included_spec)}\n')
