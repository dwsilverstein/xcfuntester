#! /usr/bin/env python

from __future__ import print_function, division
import sys, os
import subprocess

def main():
    """\
    Program for parsing information to Fortran to test analytical and 
    numerical derivatives of exchange and correlation functionals.  
    Note that the binary xcfuntester can be run alone, this is just
    a wrapper because command line arguments are easier to document
    in Python. 
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-f', '--functional', help='The exchange or correlation'
                        ' functional the derivatives are needed for.', 
                        default='xslater')
    parser.add_argument('-e', '--errortype', help='The type of error analysis'
                        ' to be done.  This can be either residual or ratio',
                        default='residual')
    parser.add_argument('-i', '--ipol', help='Determines whether the calculation'
                        ' is restricted (res) or unrestricted (unres)',
                        default='res')
    args = parser.parse_args()

    # Table of functionals available for testing.
    tfunc = {'xslater'   : 'Dirac Exchange',
             'cvwn1'     : 'VWN 1 Correlation',
             'cvwn2'     : 'VWN 2 Correlation',
             'cvwn3'     : 'VWN 3 Correlation',
             'cvwn4'     : 'VWN 4 Correlation',
             'cvwn5'     : 'VWN 5 Correlation',
             'cvwn1rpa'  : 'VWN 1 RPA Correlation',
             'cpw91lda'  : 'Perdew-Wang 91 LDA Correlation',
             'xbecke88'  : 'Becke 88 Exchange',
             'cperdew86' : 'Perdew 86 Correlation',
             'clyp'      : 'LYP Correlation',
             'xpbe96'    : 'PBE 96 Exchange', 
             'xrevpbe'   : 'revPBE Exchange', 
             'xrpbe'     : 'RPBE Exchange', 
             'cpbe96'    : 'PBE 96 Correlation',
             'xcamb88'   : 'CAM-Becke 88 Exchange',
             'xbnl07'    : 'BNL 07 Exchange',
             'xcamlsd'   : 'CAM-Slater/Dirac Exchange',
             'xcampbe96' : 'CAM-PBE 96 Exchange',
             'xwpbe'     : 'LC-wPBE Exchange'}

    # Program name (should be xcfuntester unless you tampered with the 
    # Makefile)
    program = 'xcfuntester'

    # Current working directory
    cwd = os.getcwd()

    # Get the functional requested by the user
    functional = args.functional

    # Get the type of error analysis
    errortype = args.errortype

    # Get the calculation type (restricted or unrestricted)
    ipol = args.ipol

    # Format for printing
    frmt = ' {0:9s} => {1:21s}'

    # Check that the functional is in the list of available functionals.
    # If not, print an error.
    if functional not in tfunc.keys():
        print(' ' + '/\\'*21)
        print(' WARNING!')
        print(' Functional ' + functional + ' is not recognized.\n')
        print(' Try one of the following (key => meaning):')
        print(frmt.format('xslater', tfunc['xslater']))
        print(frmt.format('cvwn1', tfunc['cvwn1']))
        print(frmt.format('cvwn2', tfunc['cvwn2']))
        print(frmt.format('cvwn3', tfunc['cvwn3']))
        print(frmt.format('cvwn4', tfunc['cvwn4']))
        print(frmt.format('cvwn5', tfunc['cvwn5']))
        print(frmt.format('cvwn1rpa', tfunc['cvwn1rpa']))
        print(frmt.format('cpw91lda', tfunc['cpw91lda']))
        print(frmt.format('xbecke88', tfunc['xbecke88']))
        print(frmt.format('cperdew86', tfunc['cperdew86']))
        print(frmt.format('clyp', tfunc['clyp']))
        print(frmt.format('xpbe96', tfunc['xpbe96']))
        print(frmt.format('xrevpbe96', tfunc['xrevpbe96']))
        print(frmt.format('xrpbe96', tfunc['xrpbe96']))
        print(frmt.format('cpbe96', tfunc['cpbe96']))
        print(frmt.format('xcamlsd', tfunc['xcamlsd']))
        print(frmt.format('xcamb88', tfunc['xcamb88']))
        print(frmt.format('xbnl07', tfunc['xbnl07']))
        print(frmt.format('xcampbe96', tfunc['xcampbe96']))
        print(frmt.format('xwpbe', tfunc['xwpbe']), end='\n\n')
        print(' Please use one of the keys above')
        print(' ' + '/\\'*21)
        sys.exit()

    # Execute the Fortran program.
    # The method to run the program is (if in the current directory):
    # ./xcfuntester functional, where functional is a functional.
    command = cwd + '/' + program
    subprocess.call([command,functional,errortype,ipol])

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
