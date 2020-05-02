#!/usr/bin/env python
import sys
import os
from os import path
import argparse
import multiprocessing as mp

def parse():
  parser = argparse.ArgumentParser(description="run metaphlan2 or strainphlan")
  parser.add_argument("--run_id", "-i", required=False, default=None, help="samples run id, one line for each sample")
  parser.add_argument("--file", "-f", required=False, default=None, help="the format is sampleName\tread1\tread2")
  parser.add_argument("--outdir", "-o", required=True, help="output direction of aligning result")
  parser.add_argument("--nproc", "-p", required=False, default=10, help="the number of CPU you wanna use")
  parser.add_argument("--threads", "-t", required=False, default=4, help="the number of threads for metaphlan2")
  parser.add_argument("--strainphlan", action='store_true', help="for strainphlan")
  args = parser.parse_args()
  return args

def parseRunid(runid):
  sampledata = "/k11e/pvdisk/bigbase/kbdata/sampledata"
  p = path.join(sampledata, '/'.join(runid[:i+1] for i in range(len(runid))))
  return (runid, path.join(p,runid+".R1.fastq.gz"), path.join(p, runid+".R2.fastq.gz"))

def fun(sample, outdir, args):
  if args.strainphlan:
    metaphlan2 = "metaphlan2.py --nproc {threads} --input_type fastq {r1},{r2} --samout {outdir}/{sample}.sam.bz2 --bowtie2out {outdir}/{sample}.bowtie2out.bz2 -o {outdir}/{sample}.metagenome.txt"
  else:
    metaphlan2 = "metaphlan2.py --nproc {threads} --input_type fastq {r1},{r2} --bowtie2out {outdir}/{sample}.bowtie2out.bz2 -o {outdir}/{sample}.metagenome.txt --biom {outdir}/{sample}.biom"
  os.system(metaphlan2.format(threads=args.threads, sample=sample[0], r1=sample[1], r2=sample[2],outdir=outdir))
  if args.strainphlan:
    marks = "sample2markers.py --ifn_samples {outdir}/{sample}.sam.bz2 --input_type sam --output_dir {outdir}"
    os.system(marks.format(sample=sample[0], outdir=outdir))

def run():
  args = parse()
  if args.file is None:
    with open(args.run_id, 'r') as f:
      samples = [parseRunid(i.strip()) for i in f]
  else:
    with open(args.file, 'r') as f:
      samples = [i.strip().split('\t') for i in f]
  pool = mp.Pool(int(args.nproc))
  for sample in samples:
    pool.apply_async(fun, args=(sample, os.path.abspath(args.outdir), args))
  pool.close()
  pool.join()

def run_one():
  args = parse()
  if args.file is None:
    with open(args.run_id, 'r') as f:
      samples = [parseRunid(i.strip()) for i in f]
  else:
    with open(args.file, 'r') as f:
      samples = [i.strip().split('\t') for i in f]
  for sample in samples:
    fun(sample, os.path.abspath(args.outdir), args)

if __name__=='__main__':  
  run()


