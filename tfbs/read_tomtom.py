from xml.etree import ElementTree
import sys


tree = ElementTree.parse(sys.argv[1])
try:
    motifs = {}
    for num, motif_id in enumerate(tree.getiterator('motif')):
        attribs = motif_id.attrib
        motifs[attribs['id']] = {'name': attribs['name'], 'alt': attribs['alt']}
except KeyError:
    pass
finally:
    for num, match_id in enumerate(tree.getiterator('match')):
        matched_motif = motifs[match_id.attrib['target']]
        print(matched_motif['name'] + '\t' + matched_motif['alt'] + '\t' + str(match_id.attrib['pvalue']))
