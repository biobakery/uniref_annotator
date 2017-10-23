#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import argparse
import csv

from utils import say, die, check_path, which, Hit, translate_fasta

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_min_coverage   = 0.80
#c_diamond_filter = "--max-target-seqs 20"
c_diamond_filter = ""
c_output_format  = "6 qseqid sseqid pident qlen qstart qend slen sstart send evalue"
g_force_search   = False

# ---------------------------------------------------------------
# cli
# ---------------------------------------------------------------

description = """
A script to annotate a fasta file of coding sequence against
HUMAnN2-formatted UniRef90/UniRef50 databases. Approximates the UniRef
clustering conventions in that query and target must exceed 80% mutual coverage
with the corresponding percent identity (e.g. 90 for UniRef90).
"""

def get_args( ):
    parser = argparse.ArgumentParser( )
    parser.add_argument( "fasta", 
                         help="Sequences to annotate",
                         )
    parser.add_argument( "--seqtype",
                         choices=["nuc", "cds", "prot"],
                         metavar="<nuc/cds/prot>",
                         default="nuc",
                         help="Sequence type [default: nuc]",
                         )
    parser.add_argument( "--diamond",
                         metavar="<path>",
                         default="diamond",
                         help="Path to diamond binary [default: in PATH]",
                         )
    parser.add_argument( "--uniref90db",
                         metavar="<path>",
                         required=True,
                         help="Path to HUMAnN2-formatted UniRef90 database",
                         )
    parser.add_argument( "--uniref50db",
                         metavar="<path>",
                         required=True,
                         help="Path to HUMAnN2-formatted UniRef50 database",
                         )
    parser.add_argument( "--transitive-map",
                         metavar="<path>",
                         help="Path to UniRef90->UniRef50 idmapping file (optional)",
                         )
    parser.add_argument( "--temp",
                         metavar="<path>",
                         default=".",
                         help="Path for temp files [default: .]",
                         )
    parser.add_argument( "--out",
                         metavar="<path>",
                         default=None,
                         help="Path for output file [default: <fasta>.annotated]",
                         )
    parser.add_argument( "--force-search",
                         action="store_true",
                         help="Rerun searches, even if expected outputs exist",
                         )
    parser.add_argument( "--flags",
                         metavar="<string>",
                         help="Additional flags to pass to diamond, e.g. --threads",
                         )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utils 
# ---------------------------------------------------------------

def say( *args ):
    print( *args, file=sys.stderr )

def die( *args ):
    args = ["FAILED:"] + list( args )
    print( *args, file=sys.stderr )
    sys.exit( )

def check_path( path ):
    if not os.path.exists( path ):
        die( "The specified path is missing: {}".format( path ) )
    return None

def get_mode( path ):
    mode = None
    for test in "90 50".split( ):
        test = "uniref" + test
        if test in path.lower( ):
            mode = test
    if mode is None:
        die( "Could not infer mode from path: {}".format( path ) )
    return mode

def uniref_search( diamond=None, database=None, query=None, seqtype=None, temp=None, flags=None ):
    if which( diamond ) is None:
        die( "<diamond> is not executable as: {}".format( diamond ) )
    for path in [database, query, temp]:
        check_path( path )
    binary = {"nuc":"blastx", "prot":"blastp"}[seqtype]
    mode = get_mode( database )
    results = os.path.split( query )[1]
    results = os.path.join( temp, results )
    results = ".".join( [results, mode, "hits"] )
    command = [
        diamond,
        binary,
        "--db", database,
        "--query", query,
        "--outfmt", c_output_format,
        "--tmpdir", temp,
        "--out", results,
        "--id", get_mode( results ).replace( "uniref", "" ),
        c_diamond_filter,
        ]
    command = " ".join( [str( k ) for k in command] )
    command += (" " + flags) if flags is not None else ""
    if not os.path.exists( results ) or g_force_search:
        say( "Executing:", command )
        os.system( command )
    else:
        say( "Using existing results file:", results )
    return results

def parse_results( results ):
    say( "Parsing results file:", results )
    mapping = {}
    mode = get_mode( results )
    min_pident = float( mode.replace( "uniref", "" ) )
    with open( results ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            h = Hit( row, config=c_output_format )
            if h.qseqid not in mapping:
                if h.pident >= min_pident and h.mcov >= c_min_coverage:
                    uniref = h.sseqid.split( "|" )[0]
                    mapping[h.qseqid] = uniref
    return mapping

def trans_mapping( uniref90map, p_trans_map ):
    say( "Loading transitive mapping file:", p_trans_map )
    check_path( p_trans_map )
    overrides = {}
    uniref90map_r = {}
    for header, uniref90 in uniref90map.items( ):
        uniref90map_r.setdefault( uniref90, set( ) ).add( header )
    with open( p_trans_map ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            uniref90, uniref50 = row
            headers = uniref90map_r.get( uniref90, set( ) )
            for h in headers:
                overrides[h] = uniref50
    return overrides

def reannotate( query=None, out=None, uniref90map=None, uniref50map=None, overrides=None ):
    say( "Writing new output file:", out )
    oh = open( out, "w" )
    with open( query ) as fh:
        for line in fh:
            line = line.strip( )
            if line[0] != ">":
                print( line, file=oh )
            else:
                header = line[1:]
                uniref90code = uniref90map.get( header, "UniRef90_unknown" )
                uniref50code = uniref50map.get( header, "UniRef50_unknown" )
                uniref50code = overrides.get( header, uniref50code )
                print( "|".join( [line, uniref90code, uniref50code] ), file=oh )
    oh.close( )
    return None

# ---------------------------------------------------------------
# main 
# ---------------------------------------------------------------
    
def main( ):
    args = get_args( )
    # global config
    global g_force_search
    if args.force_search:
        g_force_search = True
    # set defaults
    if args.out is None:
        args.out = args.fasta + ".annotated"
    # translate fasta?
    query = args.fasta
    if args.seqtype == "cds":
        query = os.path.split( query )[1]
        query = os.path.join( args.temp, query )
        query = query + ".translated"
        say( "Translating input fasta to:", query )
        translate_fasta( args.fasta, query )
        args.seqtype = "prot"
    # perform uniref90 search
    uniref90hits = uniref_search( 
        diamond=args.diamond, 
        database=args.uniref90db,
        query=query,
        seqtype=args.seqtype,
        temp=args.temp,
        flags=args.flags, )
    uniref90map = parse_results( uniref90hits )
    # perform uniref50 search
    uniref50hits = uniref_search( 
        diamond=args.diamond, 
        database=args.uniref50db,
        query=query,
        seqtype=args.seqtype,
        temp=args.temp,
        flags=args.flags, )
    uniref50map = parse_results( uniref50hits )
    # override mappings?
    overrides = {}
    if args.transitive_map is not None:
        overrides = trans_mapping( uniref90map, args.transitive_map )
    # reannoate the fasta
    reannotate( 
        query=args.fasta, 
        out=args.out, 
        uniref90map=uniref90map, 
        uniref50map=uniref50map, 
        overrides=overrides, )

if __name__ == "__main__":
    main( )
