#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import argparse
import csv

from util import Hit

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_min_coverage = 0.80
c_output_format = "6 qseqid sseqid pident qlen slen qstart qend start send"

# ---------------------------------------------------------------
# cli
# ---------------------------------------------------------------

def get_args( ):
    parser = argparse.ArgumentParser( )
    parser.add_argument( "fasta", 
                         metavar="<path>",
                         help="Sequences to annotate",
                         )
    parser.add_argument( "--seqtype",
                         choices=["nuc", "prot"],
                         metavar="<nuc/prot>",
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
    parser.add_argument( "--tmp",
                         metavar="<path>",
                         default=".",
                         help="Path for temp files (default: here)",
                         )
    parser.add_argument( "--out",
                         metavar="<path>",
                         default="<fasta>.annotated",
                         help="Path for output file (default: <fasta>.annotated)",
                         )
    parser.add_argument( "--diamond-flags",
                         metavar="<string>",
                         help="Additional flags to pass to diamond, e.g. --threads <#>",
                         )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utils 
# ---------------------------------------------------------------

def check_path( path ):
    if not os.path.exists( path ):
        sys.exit( "The specified path is missing: {}".format( path ) )
    return None

def get_mode( path ):
    mode = None
    for test in "90 50".split( ):
        test = "uniref" + test
        if test in database.lower( ):
            mode = test
    if mode is None:
        sys.exit( "Could not infer mode from path: {}".format( database ) )
    return mode

def uniref_search( diamond=None, database=None, query=None, seqtype=None, tmp=None, flags=None ):
    for path in [diamond, database, query, tmp]:
        check_path( path )
    binary = {"nuc":"blastx", "prot":"blastp"}[seqtype]
    mode = get_mode( database )
    results = os.path.split( query )[1]
    results = os.path.join( tmp, results )
    results = ".".join( [results, mode, "hits"] )    
    command = [
        diamond,
        binary,
        "--db", database,
        "--query", query,
        "--outfmt", c_output_format,
        "--max-target-seqs", 1,
        "--id", mode.replace( "uniref", "" ),
        "--tmpdir", tmp,
        "--out", results,
        ]
    command = " ".join( [str( k ) for k in command] )
    command += flags if flags is not None else ""
    print( "Executing:", command, file=sys.stderr )
    os.system( command )
    return results

def parse_results( results ):
    print( "Parsing results file:", results, file=sys.stderr )
    mapping = {}
    mode = get_mode( results )
    min_pident = float( mode.replace( "uniref", "" ) )
    with open( results ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            h = Hit( row, config=c_output_format )
            if h.qseqid not in mapping:
                if h.pident >= min_pident and h.mcov >= c_min_coverage:
                    mapping[h.qseqid] = h.sseqid
    return mapping

def trans_mapping( mapping, p_trans_map, default=None ):
    check_path( p_tran_map )
    trans_map = {}
    with open( p_trans_map ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            trans_map[row[0]] = row[1]
    overrides = {}
    for k, v in mapping.items( ):
        overrides[k] = trans_map.get( v, default )
    return overrides

def reannotate( query=None, out=None, uniref90map=None, uniref50map=None, overrides=None ):
    if out is None:
        out = query + ".annotated"
    print( "Writing new output file:", out, file=sys.stderr )
    oh = open( out, "w" )
    with open( query ) as fh:
        for line in query:
            line = line.strip( )
            if line[0] != ">":
                print( line, file=oh )
            else:
                header = line[1:]
                uniref90code = uniref90map.get( header, "UniRef90_unknown" )
                uniref50code = uniref90map.get( header, "UniRef50_unknown" )
                uniref50code = overrides.get( header, uniref50code )
                print( "|".join( [line, uniref90code, uniref50code] ), file=oh )
    oh.close( )
    return None

# ---------------------------------------------------------------
# main 
# ---------------------------------------------------------------
    
def main( ):
    args = get_args( )
    # perform uniref90 search
    uniref90hits = uniref_search( 
        diamond=args.diamond, 
        database=args.uniref90db,
        query=args.query,
        seqtype=args.seqtype,
        tmp=args.temp,
        flags=args.flags, )
    uniref90map = parse_results( uniref90hits )
    # perform uniref50 search
    uniref50hits = uniref_search( 
        diamond=args.diamond, 
        database=args.uniref50db,
        query=args.query,
        seqtype=args.seqtype,
        tmp=args.temp,
        flags=args.flags, )
    uniref50map = parse_results( uniref50hits )
    # override some mappings?
    overrides = {}
    if args.transitive_map is not None:
        overrides = trans_mapping( uniref90map, args.transitive_map )
    # reannoate the fasta
    reannotate( 
        query=args.query, 
        out=args.out, 
        uniref90map=uniref90map, 
        uniref50map=uniref50map, 
        overrides=overrides, )

if __name__ == "__main__":
    main( )
