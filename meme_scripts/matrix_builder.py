import argparse
import random
import pandas
from collections import OrderedDict
from xml.etree import ElementTree

def to_tree(xml):
    ''' 
    Parses an XML file using in-built Python 3.x libraries.
    @param xml: XML file produced from running MAST.
    @return: ElementTree object referencing the parsed XML file.
    '''
    print('Parsing', xml, '... ', end='')
    et_obj = ElementTree.parse(xml)
    tree = et_obj.getroot()
    print('[OK]')
    return tree

def build_mapping(tree):
    ''' 
    MAST references each PWM via an integer-string key.
    By keeping track of this key/value mapping, referencing 
    mapped PWMs enables efficient identification of which PWMs
    mapped and where.
    @param xml: Input XML file.  
    '''
    print('Mapping MAST IDs to PWMs ... ', end='')
    elem_motifs = tree.find("./motifs")
    mapping = {} # key => MAST motif ID, value => PWM name
    for m in elem_motifs:
        if m.tag == 'motif': # process only motifs
            k, v = m.attrib['id'], m.attrib['name']
            mapping[k] = v # map MAST ID-to-PWM data-structure
    print(len(mapping), 'PWMs parsed [OK]')
    return mapping

def to_file(m, f):
    ''' 
    Writes the pre-computed count-matrix to a user-provided
    output file.
    @param m: Count-matrix where headers are PWM names.
    @param f: User-provided output CSV file.
    '''
    df = pandas.DataFrame(m)
    df.to_csv(f)

def build_matrix(tree_cont, tree_query):
    '''
    Builds a count-matrix given PWM mappings in both control
    and query datasets.
    @param tree_cont: Control PWM counts of type ElementTree.
    @param tree_query: Query PWM counts of type ElementTree.
    '''
    # map MAST PWM ID to its PWM; choose control or query 
    dict_pwms = build_mapping(random.choice([tree_cont, tree_query]))
    pwms = list(dict_pwms.values())
    
    # get control and query sequences
    seqs_cont = [x for x in tree_cont.find('./sequences') if x.tag == 'sequence']
    seqs_query = [x for x in tree_query.find('./sequences') if x.tag == 'sequence']
    num_rows = len(seqs_cont) + len(seqs_query) # each sequence is a row
    
    # create an empty skeleton matrix for use in adding PWM counts    
    m = OrderedDict({'Sequence': [''] * num_rows}) # add Sequence vector
    m.update({pwm: [0] * num_rows for pwm in pwms}) # add PWM vectors
    m.update({'Target': [None] * num_rows})
    
    row_num = 0 # row-number counter
    for data_num, dataset in enumerate([seqs_cont, seqs_query]):
        for seq in dataset:
            seq_name = seq.attrib['name']
            print(data_num, row_num, seq_name, seq.attrib)
            hits = seq.findall('.//seg/hit') # get all PWM hits
            for hit in hits:
                pwm_id = hit.attrib['motif'] # get MAST PWM ID
                pwm = dict_pwms[pwm_id] # MAST ID references PWM name
                print('\t', pwm, hit.attrib, pwm in m)
                m[pwm][row_num] += 1 # increment count
            m['Target'][row_num] = data_num # either 0 or 1
            m['Sequence'][row_num] = seq_name # set sequence accession
            row_num += 1 # increment row counter
    return m

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-control', metavar='XML', required=True,
                   help='MAST output XML control file [req]')
    p.add_argument('-query', metavar='XML', required=True,
                   help='MAST output XML query file [req]')
    p.add_argument('-csv', metavar='FILE', required=False, 
                   default='out.csv',
                   help='Output file [out.csv]')
    args = vars(p.parse_args())
    
    # parse both XML files.
    tree_cont = to_tree(xml = args['control'])
    tree_query = to_tree(xml = args['query'])
    
    # begin mapping PWM counts to their respective sequence
    m = build_matrix(tree_cont, tree_query)
    to_file(m=m, f=args['csv']) # lastly, save results
    print('Program successful [OK]')
