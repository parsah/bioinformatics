import argparse

def read_file(fname):
    ''' 
    Reads the user-provided SAX file and saves each entry.
    @param fname: input filename referencing SAX symbols.
    @return: list of parsed symbols.
    '''
    symbols = []
    for line in open(fname):
        line = line.strip()
        symbols.append(line)
    return symbols

def get_chars(symbols):
    ''' 
    Builds a set of all the individual symbols in the SAX collection.
    @param symbols: list of all the parsed symbols.
    @return: dictionary of each symbol and its column-specific abundance.
    '''
    chars = {}
    for symbol in symbols:
        symbol = list(symbol)
        for s in symbol:
            chars[s] = [0]*len(symbol)
    return chars

def compute_abundance(symbols):
    ''' 
    Computes percentages as to abundance per character in each column.
    @param symbols: List of parsed symbols.
    '''
    chars = get_chars(symbols)
    for i, symbol in enumerate(symbols):
        print(symbol)
        for j, char in enumerate(list(symbol)):
            chars[char][j] += 1 # increment column-specific abundance 
    
    for symbol in chars:
        counts = chars[symbol]
        counts = [round(float(i)/len(symbols)*100, 2) for i in chars[symbol]]
        print(symbol, '\t', counts)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', metavar='FILE', required=True, 
                        help='SAX symbols; one per line [na]')
    args = vars(parser.parse_args())
    symbols = read_file(fname=args['in'])
    compute_abundance(symbols) # sends output to standard-output.
