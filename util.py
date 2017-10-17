#!/usr/bin/env python

import os
import sys

from zopy.utils import die

blast_fields = [
    ["qseqid",str,"Query Seq-id"],
    ["qgi",str,"Query GI"],
    ["qacc",str,"Query accesion"],
    ["qaccver",str,"Query accesion.version"],
    ["qlen",int,"Query sequence length"],
    ["sseqid",str,"Subject Seq-id"],
    ["sallseqid",str,"All subject Seq-id(s), separated by a ';'"],
    ["sgi",str,"Subject GI"],
    ["sallgi",str,"All subject GIs"],
    ["sacc",str,"Subject accession"],
    ["saccver",str,"Subject accession.version"],
    ["sallacc",str,"All subject accessions"],
    ["slen",int,"Subject sequence length"],
    ["qstart",int,"Start of alignment in query"],
    ["qend",int,"End of alignment in query"],
    ["sstart",int,"Start of alignment in subject"],
    ["send",int,"End of alignment in subject"],
    ["qseq",str,"Aligned part of query sequence"],
    ["sseq",str,"Aligned part of subject sequence"],
    ["evalue",float,"Expect value"],
    ["bitscore",float,"Bit score"],
    ["score",float,"Raw score"],
    ["length",int,"Alignment length"],
    ["pident",float,"Percentage of identical matches"],
    ["nident",int,"Number of identical matches"],
    ["mismatch",int,"Number of mismatches"],
    ["positive",int,"Number of positive-scoring matches"],
    ["gapopen",int,"Number of gap openings"],
    ["gaps",int,"Total number of gaps"],
    ["ppos",float,"Percentage of positive-scoring matches"],
    ["frames",str,"Query and subject frames separated by a '/'"],
    ["qframe",str,"Query frame"],
    ["sframe",str,"Subject frame"],
    ["btop",str,"Blast traceback operations (BTOP)"],
    ["staxid",str,"Subject Taxonomy ID"],
    ["ssciname",str,"Subject Scientific Name"],
    ["scomname",str,"Subject Common Name"],
    ["sblastname",str,"Subject Blast Name"],
    ["sskingdom",str,"Subject Super Kingdom"],
    ["staxids",str,"unique Subject Taxonomy ID(s), separated by a ';'"],
    ["sscinames",str,"unique Subject Scientific Name(s), separated by a ';'"],
    ["scomnames",str,"unique Subject Common Name(s), separated by a ';'"],
    ["sblastnames",str,"unique Subject Blast Name(s), separated by a ';'"],
    ["sskingdoms",str,"unique Subject Super Kingdom(s), separated by a ';'"],
    ["stitle",str,"Subject Title"],
    ["salltitles",str,"All Subject Title(s), separated by a '<>'"],
    ["sstrand",str,"Subject Strand"],
    ["qcovs",float,"Query Coverage Per Subject"],
    ["qcovhsp",float,"Query Coverage Per HSP"],
    ["qcovus",float,"Query Coverage Per Unique Subject (blastn only)"],
]
format = {}
for n, t, d in blast_fields:
    format[n] = t

def contains( items, collection ):
    ret = True
    for i in items:
        if i not in collection:
            ret = False
            break
    return ret

class Hit:
    def __init__( self,
                  row,
                  config="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                  ):
        if config[0:2] == "6 ":
            config = config[2:]
        config = config.split( " " )
        if len( config ) != len( row ):
            die( "config doesn't match row" )
        self.data = {}
        for value, field in zip( row, config ):
            value = format[field]( value )
            self.data[field] = value
        # qcov
        self.data["qcov"] = None
        if contains( "qstart qend qlen".split( ), self.data ):
            self.data["qcov"] = abs( self.data["qend"] - self.data["qstart"] ) + 1
            self.data["qcov"] /= float( self.data["qlen"] )
        # scov
        self.data["scov"] = None
        if contains( "sstart send slen".split( ), self.data ):
            self.data["scov"] = abs( self.data["send"] - self.data["sstart"] ) + 1
            self.data["scov"] /= float( self.data["slen"] )
        # mcov
        self.data["mcov"] = None
        if self.data["qcov"] is not None and self.data["scov"] is not None:
            self.data["mcov"] = min( self.data["qcov"], self.data["scov"] )
        # score
        self.data["strength"] = None
        if self.data["mcov"] is not None and self.data["pident"] is not None:
            self.data["strength"] = self.data["mcov"] * self.data["pident"]
        # set as attr
        for f, v in self.data.items( ):
            setattr( self, f, v )

def which( program ):
    """
    Adapted from:
    https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    ret = None
    def is_exe( fpath ):
        return os.path.isfile( fpath ) and os.access( fpath, os.X_OK )
    fpath, fname = os.path.split( program )
    if fpath and is_exe( program ):
        ret = program
    else:
        for path in os.environ["PATH"].split( os.pathsep ):
            path = path.strip( '"' )
            exe_file = os.path.join( path, program )
            if is_exe( exe_file ):
                ret = exe_file
    return ret
