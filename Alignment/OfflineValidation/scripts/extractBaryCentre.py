#!/usr/bin/env python3

from argparse import ArgumentParser
import logging
import json
import sys
import os
import ROOT

class TFileContext(object):
    def __init__(self, *args):
        self.tfile = ROOT.TFile(*args)
    def __enter__(self):
        return self.tfile
    def __exit__(self, type, value, traceback):
        self.tfile.Close()

def display_results(t_data, style='twiki', out=sys.stdout):
    before = after = None
    if  (style == 'twiki'):
        header_fmt = ' | '.join( ['{:^6s}'] + ['{:^9s}' ] * (len(t_data.keys())-1) )
        row_fmt    = ' | '.join( ['{:6d}' ] + ['{:9.6f}'] * (len(t_data.keys())-1) )
        header_fmt = '| '+header_fmt+' |'
        row_fmt    = '| '+row_fmt   +' |'
    elif(style == 'latex'):
        before     = '\\begin{tabular}{' + 'c'*len(t_data.keys()) + '}'
        after      = '\\end{tabular}'
        header_fmt = ' & '.join( ['{:^6s}'] + ['{:^9s}' ] * (len(t_data.keys())-1) )  + ' \\\\\n\\hline'
        row_fmt    = ' & '.join( ['{:6d}' ] + ['{:9.6f}'] * (len(t_data.keys())-1) ) + ' \\\\'
    elif(style == 'csv'):
        header_fmt = ', '.join(  ['{:s}'  ]                *  len(t_data.keys())    )
        row_fmt    = ', '.join(  ['{:d}'  ] + ['{:f}'    ] * (len(t_data.keys())-1) )
    else:
        raise RuntimeError('Unknown style "%s" for table'%(style))

    if(before is not None): out.write(before+'\n')
    out.write( header_fmt.format(*t_data.keys())+'\n' )
    for i, run in enumerate(t_data['run']):
        out.write(row_fmt.format(run, *[t_data[k][i] for k in t_data.keys() if not k == 'run'])+'\n')
    if(after is not None): out.write(after+'\n')

    return out


def main(args):
    logging.debug('args: %s', args)

    folder    = 'PixelBaryCentreAnalyzer' if not args.quality          else 'PixelBaryCentreAnalyzerWithPixelQuality'
    tree_name = 'PixelBarycentre'         if args.type == 'barycentre' else 'BeamSpot'
    if(args.label is not None):
        tree_name += '_'+args.label
    columns = ['run'] + [(args.partition if args.type=='barycentre' else 'BS')+'.'+coord for coord in ('x', 'y', 'z')]

    with TFileContext(args.fname) as tfile:
        if(args.list_content):
            list_content(tfile)
            return 0

        tfolder = tfile.Get(folder)
        if(not tfolder):
            raise KeyError('Folder "%s" not found in "%s"' %(folder, args.fname))

        logging.debug('Opened folder "%s"', folder)
        tree = tfolder.Get(tree_name)
        if(not tree):
            logging.error('Tree "%s" not found; content of file "%s":', tree_name, args.fname)
            list_content(args.fname)
            raise KeyError(tree_name)

        if(args.list_branches):
            list_branches(tree, folder_name=folder)
            return 0

        rdf = ROOT.RDataFrame(tree)
        logging.info('Reading "%s"', tree_name)
        results = rdf.AsNumpy(columns)

    if(args.config is None):
        display_results(results, style=args.style)
    else:
        # When this script is called from the validation framework, the output is
        # written to a list of files, one for each style, with appropriate extensions
        fname_base = '_'.join([args.type, args.partition] + (['quality'] if args.quality else []))
        fname_path = os.path.join(args.config['output'], fname_base) # without the ext
        for style in args.config["styles"]:
            if  (style == 'twiki'): ext = 'txt'
            elif(style == 'latex'): ext = 'tex'
            else: ext = style
            fname = '.'.join([fname_path, ext])
            logging.info('output in "%s"', fname)

            with open(fname, 'w') as fout:
                display_results(results, style=style, out=fout)

    return 0


def list_content(tfile):
    for tfolder in tfile.GetListOfKeys():
        for key in tfolder.ReadObj().GetListOfKeys():
            obj = key.ReadObj()
            print('%s\t%s/%s\t%d branches, %d entries' %(obj.ClassName(), tfolder.GetName(), key.GetName(), len(obj.GetListOfBranches()), obj.GetEntries()))


def list_branches(tree, folder_name=''):
    branches = [b.GetName() for b in tree.GetListOfBranches()]
    logging.info('Branches (%s/%s): %d', (folder_name, tree.GetName(), len(branches)))
    print('\n'.join(branches))


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('fname'            , metavar='FILE')
    parser.add_argument('-p', '--partition', default='BPIX', help='Tracker partition (e.g. BPIX, FPIX, BPIXLYR1, etc.). Default: %(default)s')
    parser.add_argument('-l', '--list-content' , action='store_true', dest='list_content', help='List the contents of file and exit')
    parser.add_argument(      '--list-branches', action='store_true', help='List the branches of the tree and exit')
    parser.add_argument('-t', '--type'     , default='barycentre', choices=('barycentre', 'beamspot'), type=str.lower, help='Default: %(default)s')
    parser.add_argument(      '--label'    , default=None, help='Additional label that is appended to the folder name (i.e. PixelBaryCentreAnalyzer by default)')
    parser.add_argument(      '--quality'  , action='store_true', help='Read results with the WithPixelQuality flag (default: %(default)s)')
    parser.add_argument('-s', '--style'    , default='twiki', choices=('twiki', 'latex', 'csv'), type=str.lower, help='Table style for the output (default: %(default)s)')
    parser.add_argument(      '--loglevel' , metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args = parser.parse_args()

    # If the script is called from the validation framework, the first argument
    # will be a JSON file with the configuration, instead of a ROOT file
    try:
        with open(args.fname) as f:
            config = json.load(f)
        args.fname = config['input']
        args.config = config
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        args.config = None

    return args


if __name__ == '__main__':
    args = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    exit(main(args))
