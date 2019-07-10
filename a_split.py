# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, imp
sys.path.append('/cluster/mshen/')
import numpy as np
from collections import defaultdict
from mylib import util
import pandas as pd
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
import seaborn as sns

# Default params
inp_dir = _config.READS_DIR
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

from _config import FILE_FOLDERS_MESC, FILE_FOLDERS_U2OS

##
# Functions
##
def split(inp_fn, out_nm):
  #print inp_fn
  inp_fn_numlines = util.line_count(inp_fn)

  #print out_nm
  
  #print inp_fn
  num_splits = 15
  split_size = int(inp_fn_numlines / num_splits)
  if num_splits * split_size < inp_fn_numlines:
    split_size += 1
  while split_size % 4 != 0:
    split_size += 1
  #print 'Using split size %s' % (split_size)

  split_num = 0
  for idx in range(1, inp_fn_numlines, split_size):
    start = idx
    end = start + split_size  
    out_fn = out_dir + out_nm + '_%s.fastq' % (split_num)
    command = 'tail -n +%s %s | head -n %s > %s' % (start, inp_fn, end - start, out_fn)
    split_num += 1
    print command

  return


##
# Main
##
#@util.time_dec
def main():
  #print NAME  
  
  # Function calls
  for fn in os.listdir(inp_dir):
    if 'fastq' in fn:
        #print fn 
        files = ["190124Gif_D19-{0}_{1}_sequence.fastq".format(pool , read)
                 for pool in range(554,572)
                 for read in [1,2]]

        #print "files are"
        if fn in files:
             split(os.path.join(inp_dir, fn), fn.replace('.fastq', ''))
      
  return


if __name__ == '__main__':
  main()
