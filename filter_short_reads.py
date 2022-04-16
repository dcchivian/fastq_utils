#!/usr/bin/python3
'''
Copyright 2022 Dylan Chivian

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
'''

import sys
import argparse
import gzip
import re


# getargs()
#
def getargs():
    default_len = 25
    parser = argparse.ArgumentParser(description="filter out short reads from fastq file(s)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-s", "--singleend", action="store_true", help="single-end library, use only forward reads")
    group.add_argument("-p", "--pairedend", action="store_true", help="paired-end library, use forward and reverse or --interleaved and forward reads")
    parser.add_argument("-i", "--interleaved", action="store_true", help="for paired-end library, use forward reads file only")
    parser.add_argument("-f", "--forwardreads", help="forward reads library")
    parser.add_argument("-r", "--reversereads", help="reverse reads library")
    parser.add_argument("-o", "--outputfile", help="output reads library basename")
    parser.add_argument("-l", "--length", type=int, default=default_len, help="the read length <= to filter (default={})".format(default_len))
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(-1)
    if not args.outputfile:
        print ("must specify --outputfile\n")
        parser.print_help()
        sys.exit (-1)
    if not args.singleend and not args.pairedend:
        print ("must specify either --singleend or --pairedend\n")
        parser.print_help()
        sys.exit (-1)
    if args.singleend and args.reversereads:
        print ("cannot specify --reversereads and --singleend\n")
        parser.print_help()
        sys.exit (-1)
    if args.pairedend and (not args.interleaved and not args.reversereads):
        print ("for pairedend, must also specify either --interleaved or --reversereads\n")
        parser.print_help()
        sys.exit (-1)
    if args.interleaved and args.reversereads:
        print ("for pairedend, must also specify either --interleaved or --reversereads, but not both\n")
        parser.print_help()
        sys.exit (-1)
        
    return args


# parse_read_id()
#
def parse_read_id (header_line, paired_end_flag):
    read_id = re.sub ("[ \t]+.*$", "", header_line)
    # manage read_id edge case: e.g. @SRR5891520.1.1 (forward) & @SRR5891520.1.2 (reverse)                                                               
    if paired_end_flag:
        read_id = ''.join(read_id.rsplit('.',1)) # replace last '.' with ''   
        read_id = re.sub ("[\/\.\_\-\:\;][012lrLRfrFR53]\'*$", "", read_id)    
    return read_id


# get_skip_ids()
#
def get_skip_ids (readsfile, skip_len, paired_end_flag):
    these_skip_ids = dict()
    print ("reading input file {} ...".format(readsfile))
    if readsfile.lower().endswith('.gz'):
        f = gzip.open(readsfile, 'rt')
    else:
        f = open(readsfile, 'r')

    this_id = None
    counter = 0
    read_cnt = 0
    for line in f:
        line = line.rstrip()
        if counter == 4:
            counter = 0
        if line.startswith('@') and counter == 0:
            this_id = parse_read_id(line, paired_end_flag)
            read_cnt += 1
            if read_cnt % 1000000 == 0:
                print ("READs evaluated {}".format(read_cnt))
        elif this_id != None:
            if len(line) <= skip_len:
                these_skip_ids[this_id] = True
            this_id = None
        counter += 1
    print ("READs evaluated {}".format(read_cnt))
    f.close()

    return these_skip_ids


# write_filtered_output()
#
def write_filtered_output (skip_ids=None,
                           outputfile=None,
                           readsfile=None,
                           paired_end_flag=None,
                           read_direction=None):
    
    gzip_output = False
    if outputfile.lower().endswith(".gz"):
        gzip_output = True

    # set this_outputfile name
    this_outputfile = outputfile
    if read_direction != None:
        if gzip_output:
            this_outputfile = re.sub('.gz', '', this_outputfile, flags=re.IGNORECASE)
        if this_outputfile.lower().endswith(".fq"):
            extension = ".fq"
        else:
            extension = ".fastq"
        this_outputfile = re.sub(extension, '', this_outputfile, flags=re.IGNORECASE)
        
        this_outputfile = this_outputfile + '.'+read_direction + extension
        if gzip_output:
            this_outputfile += ".gz"

    # open this_outputfile
    print ("writing output file {} ...".format(this_outputfile))
    if gzip_output:
        out = gzip.open(this_outputfile, 'wt')
    else:
        out = open(this_outputfile, 'w')
    
    # read the lib and skip marked reads
    if readsfile.lower().endswith('.gz'):
        f = gzip.open(readsfile, 'rt')
    else:
        f = open(readsfile, 'r')
        
    this_id = None
    counter = 0
    read_buf = []
    skip_read = False
    read_cnt = 0
    for line in f:
        line = line.rstrip()
        if counter == 4:
            counter = 0
            if not skip_read:
                out.write("\n".join(read_buf)+"\n")
            read_buf = []
            skip_read = False
        if line.startswith('@') and counter == 0:
            this_id = parse_read_id(line, paired_end_flag)
            if this_id in skip_ids:
                skip_read = True
            read_cnt += 1
            if read_cnt % 1000000 == 0:
                print ("READs processed {}".format(read_cnt))
        read_buf += [line]
        counter += 1
    if not skip_read:  # one more
        out.write("\n".join(read_buf)+"\n")
    print ("READs processed {}".format(read_cnt))
    f.close()
    out.close()

    return this_outputfile


# merge_skip_ids()
#
def merge_skip_ids(base_skip_ids, new_skip_ids):
    merged_skip_ids = dict()
    for k in base_skip_ids.keys():
        merged_skip_ids[k] = base_skip_ids[k]
    for k in new_skip_ids.keys():
        if k in base_skip_ids and base_skip_ids[k] != new_skip_ids[k]:
            print ("mismatch in skip_ids for id {}".format(k))
            sys.exit(-2)
        merged_skip_ids[k] = new_skip_ids[k]
    return merged_skip_ids


# main()
#
def main() -> int:
    args = getargs()

    # get skip ids (note: also skips forward read if reverse short, and vice versa)
    skip_ids = get_skip_ids (args.forwardreads, args.length, args.pairedend)
    if args.pairedend and not args.interleaved:
        skip_ids = merge_skip_ids (skip_ids, get_skip_ids (args.reversereads, args.length, args.pairedend))
        
    if not skip_ids:
        print ("no reads less than {} found.  No reads to filter and not creating output".format(args.length))
    else:
        print ("SKIPPING {} reads".format(len(skip_ids.keys())))

        # write filtered forward output
        read_direction = None
        if args.pairedend and not args.interleaved:
            read_direction = 'R1'
        this_outputfile = write_filtered_output (skip_ids=skip_ids,
                                                 outputfile=args.outputfile,
                                                 readsfile=args.forwardreads,
                                                 paired_end_flag=args.pairedend,
                                                 read_direction=read_direction)

        # write filtered reverse output
        read_direction = 'R2'
        if args.pairedend and args.reversereads:
            this_outputfile = write_filtered_output (skip_ids=skip_ids,
                                                     outputfile=args.outputfile,
                                                     readsfile=args.reversereads,
                                                     paired_end_flag=args.pairedend,
                                                     read_direction=read_direction)
                   
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
