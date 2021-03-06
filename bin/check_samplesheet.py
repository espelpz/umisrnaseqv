#!/usr/bin/env python

import os
import sys
import errno
import argparse

def parse_args(args=None):
    Description = 'Reformat samplesheet file and check its contents.'
    Epilog = """Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet file.")
    #parser.add_argument('FILE_OUT', help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def print_error(error,line):
    print("ERROR: Please check samplesheet -> {}\nLine: '{}'".format(error,line.strip()))
    sys.exit(1)


#def check_samplesheet(FileIn,FileOut):
def check_samplesheet(FileIn):
    ## Check header
    #HEADER = ['sample_id', 'fastq_1', 'fastq_2']
    HEADER = ['sample_id', 'fastq_1', 'fastq_2', 'umis']
    fin = open(FileIn,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print("ERROR: Please check samplesheet header -> {} != {}".format(','.join(header),','.join(HEADER)))
        sys.exit(1)

    sampleRunDict = {}
    while True:
        line = fin.readline()
        if line:
            lspl = [x.strip() for x in line.strip().split(',')]

            ## Check valid number of columns per row
            if len(lspl) != len(header):
                print_error("Invalid number of columns (minimum = {})!".format(len(header)),line)

            numCols = len([x for x in lspl if x])
            if numCols < 4: ##
                print_error("Invalid number of populated columns!",line)

            ## Check sample name entries
            sample,fastQFiles = lspl[0],lspl[1:]
            ##if sample:
            if sample.find(' ') != -1:
                print_error("Sample entry contains spaces!",line)
            ##else: #obligatoriamente debe existir sample porque se requiere que numCols sea 4 (que se rellenen todas las columnas)
            ##    print_error("Sample entry has not been specified!",line)

            ## Check FastQ file extension
            for fastq in fastQFiles:
                #if fastq:
                if fastq.find(' ') != -1:
                    print_error("FastQ file contains spaces!",line)
                if fastq[-9:] != '.fastq.gz' and fastq[-6:] != '.fq.gz':
                    print_error("FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",line)
                

            ## Auto-detect paired-end/single-end
            """
            sample_info = []                                                ## [single_end, is_sra, is_ftp, fastq_1, fastq_2, md5_1, md5_2]
            fastq_1,fastq_2 = fastQFiles
            if sample and fastq_1 and fastq_2:                              ## Paired-end short reads
                sample_info = ['0', '0', '0', fastq_1, fastq_2, '', '']
            elif sample and fastq_1 and not fastq_2:                        ## Single-end short reads
                sample_info = ['1', '0', '0', fastq_1, fastq_2, '', '']
            else:
                print_error("Invalid combination of columns provided!",line)
            """
            if sample not in sampleRunDict:
                sampleRunDict[sample] = [fastQFiles]
            else:
                if fastQFiles in sampleRunDict[sample]:
                    print_error("Samplesheet contains duplicate rows!",line)
                else:
                    sampleRunDict[sample].append(fastQFiles)
        else:
            fin.close()
            break

    ## Write validated samplesheet with appropriate columns
#    if len(sampleRunDict) > 0:
#        OutDir = os.path.dirname(FileOut)
#        make_dir(OutDir)
#        fout = open(FileOut,'w')
        #fout.write(','.join(['sample_id', 'single_end', 'is_sra', 'is_ftp', 'fastq_1', 'fastq_2', 'md5_1', 'md5_2']) + '\n')
#        fout.write(','.join(['sample_id', 'fastq_1', 'fastq_2', 'umis']) + '\n')

#        for sample in sorted(sampleRunDict.keys()):

            ## Check that multiple runs of the same sample are of the same datatype
            ##if not all(x[:2] == sampleRunDict[sample][0][:2] for x in sampleRunDict[sample]):
              ##  print_error("Multiple runs of a sample must be of the same datatype","Sample: {}".format(sample))

#            for idx,val in enumerate(sampleRunDict[sample]):
#                fout.write(','.join(["{}_T{}".format(sample,idx+1)] + val) + '\n')
#        fout.close()


def main(args=None):
    args = parse_args(args)
    #check_samplesheet(args.FILE_IN,args.FILE_OUT)
    check_samplesheet(args.FILE_IN)


if __name__ == '__main__':
    sys.exit(main())