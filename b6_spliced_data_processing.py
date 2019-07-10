# CONFIG AND HELPERS
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, copy
import numpy as np
from collections import defaultdict
sys.path.append('/cluster/mshen/')
from mylib import util
from mylib import compbio
import pickle
import pandas as pd
import re

# Default params
inp_dir = _config.OUT_PLACE + 'a_split/'
NAME = util.get_fn(__file__)


out_place = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_place)
exp_design = pd.read_csv(_config.DATA_DIR + '061318_exonskipping_library.csv')


#currently some of the library members will be mapped ambigiously
#this first pass over the data simply ignores higher #ed instances 
#of a given target sequence

names_targets = dict()
target_names = dict()
seen_targets = set()
for idx, row in sorted(exp_design.iterrows(), key = lambda x:int(x[1]['Identifier number']) ):
      #names_targets[row['Identifier number']] = row['Sequence']
      target=  row['Designed 61-bp target site (37i-24e, AG)']
      if target in seen_targets: 
            continue
      
      names_targets[row['Identifier number']] = row['Designed 61-bp target site (37i-24e, AG)']
      target_names[ row['Designed 61-bp target site (37i-24e, AG)']] = row['Identifier number']
      seen_targets.add(target)
        
target_latternames = dict([(k[-24:],v) for k,v in target_names.iteritems()])

        
##CUSTOM CODE FOR DICTIONARY CREATION
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
def rc(inp):
    d = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    return "".join([d[e] for e in inp.strip()[::-1]])


def quality(line):
  q_1 = line.strip()
  qs = [ord(s)-33 for s in q_1]  
  return np.mean(qs)

i = -1

##
# Alignments
##
def alignment(read, cand_idxs):
  seq_align_tool = '/cluster/mshen/tools/seq-align/bin/needleman_wunsch'
  targets = [names_targets[nm] for nm in cand_idxs]
  aligns = []
  for target_seq in targets:
    try:
      targetse =  target_seq
      
      align = subprocess.check_output(seq_align_tool + ' --match 1 --mismatch -1 --gapopen -5 --gapextend -0 --freestartgap ' + read + ' ' + targetse, shell = True)
      aligns.append(align)
    except:
      pass

  if len(aligns) > 1:
    best_align = pick_best_alignment(aligns)
    best_idx = cand_idxs[aligns.index(best_align)]
  else:
    best_align = aligns[0]
    best_idx = cand_idxs[0]
  return best_idx, best_align

def pick_best_alignment(aligns):
  scores = []
  for align in aligns:
    w = align.split()
    s1, s2 = w[0], w[1]
    score = 0
    for i in range(len(s1)):
      if s1[i] == s2[i] and s1[i] != '-':
        score += 1
    scores.append(score)
  best_idx = scores.index(max(scores))
  return aligns[best_idx]

##
# Locality sensitive hashing
##
def build_targets_better_lsh():
  lsh_dict = defaultdict(list)
  for nm in names_targets:
    target = names_targets[nm]
    kmers = get_lsh_kmers(target)
    for kmer in kmers:
      lsh_dict[kmer].append(nm)
  return lsh_dict

def get_lsh_kmers(target):
  kmer_len = 7
  kmers = []
  for idx in range(len(target) - kmer_len):
    kmer = target[idx : idx + kmer_len]
    kmers.append(kmer)
  return kmers

def find_best_designed_target(read, lsh_dict):
  kmers = get_lsh_kmers(read)
  scores = dict()
  for kmer in kmers:
    for exp in lsh_dict[kmer]:
      if exp not in scores:
        scores[exp] = 0
      scores[exp] += 1

  if len(scores) == 0:
    return []

  sorted_scores = sorted(scores, key = scores.get, reverse = True)
  best_score = scores[sorted_scores[0]]
  # cand_idxs = []
  # for exp in sorted_scores:
  #   if scores[exp] + 5 < best_score:
  #     break
  #   cand_idxs.append(exp)
  cand_idxs = [sorted_scores[0]]
  return cand_idxs

##
# IO
##
def store_alignment(alignment_buffer, idx, align_header, align):
  align_string = '%s\n%s' % (align_header, align)
  alignment_buffer[idx].append(align_string)
  return

def init_umis_alignments_buffer():
  umis_alignments_buffer = defaultdict(list)
  return umis_alignments_buffer

def flush_tuples(umis_alignments_buffer, out_dir):
  print 'Flushing... \n%s' % (datetime.datetime.now())
  for umi in umis_alignments_buffer:
    umi_location_code = umi[:6]
    with open(out_dir + '%s.txt' % (umi_location_code), 'a') as f:
      for output in umis_alignments_buffer[umi]:
        f.write(output+"\n")
  umis_alignments_buffer.clear()
  print 'Done flushing.\n%s' % (datetime.datetime.now())
  
  return

def prepare_outfns(out_dir):
  let ="ATGC"
   
  for umi_short in [l1+l2+l3+l4+l5+l6 for l1 in let for l2 in let for l3 in let for l4 in let for l5 in let for l6 in let]:
    out_fn = out_dir + '%s.txt' % (umi_short)
    util.exists_empty_fn(out_fn)
  return




##
# Main
##



##This version of matchmaker looks through the reads, checking for the constant
#sequence at the start of each "read1". If "read1" constant sequence matches,
#computes an alignment between the target library and the read and stores it with the
#umi to record the pre-spliced context

def matchmaker(nm, split):

  ##CUSTOM CODE FOR DICTIONARY CREATION
  from Bio import pairwise2
  from Bio.pairwise2 import format_alignment
  
  from Bio.Seq import Seq
  from Bio.Alphabet import generic_dna
  
  def rc(inp):
      d = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
      return "".join([d[e] for e in inp.strip()[::-1]])
  
  #UNSPLICED DATA PROCESSING
  READ1_TEMPLATE="NNNtaccagctgccctcgTCGaCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNtgattacacatatagacacgcGAGCAGCCATCTTTTATAGAATGGGtagaacccgtcctaaggactcagattgagcatcgtttgcttctcgagtactacctgg"
  READ2_TEMPLATE="NNNaaccgctgtgttctgcACGCGTNNNNNNNNNNNNNNNNNNACCGGTgcaggtaatgggccttactatcagtctcagtccttgtacagctcgtccatgccgagagtgatcccggcggcggtcacgaactccagcaggaccatgtgatcgcgcttctcgttggggtctttgctca"
  
  r1_seq = Seq(READ1_TEMPLATE,"generic_dna".upper())
  r2_seq = Seq(READ2_TEMPLATE,"generic_dna".upper())
  
  
  def quality(line):
    q_1 = line.strip()
    qs = [ord(s)-33 for s in q_1]  
    return np.mean(qs)
  
  i = -1

  

  qc_rejection_count = 0
  read1_rejection_count = 0 
  read_constant_rejection_count = 0 
  accepted_count = 0
  nolib_rejection_count = 0
  num_bad_aligns = 0
  
  print nm, split

  #fq_unspliced_1 = open("/cluster/bh0085/prj/exons/data/{0}_1_sequence.fastq".format(nm))
  #fq_unspliced_2 = open("/cluster/bh0085/prj/exons/data/{0}_2_sequence.fastq".format(nm))

  
  stdout_fn = _config.SRC_DIR + 'b2_status_%s_%s.out' % (nm, split)
  util.exists_empty_fn(stdout_fn)
  
  out_dir = out_place + nm + '/' + split + '/'
  util.ensure_dir_exists(out_dir)

  inp_fn1 = inp_dir + '%s_1_sequence_%s.fastq' % (nm, split)
  inp_fn2 = inp_dir + '%s_2_sequence_%s.fastq' % (nm, split)

  lsh_dict = build_targets_better_lsh()
  umis_alignments_buffer = init_umis_alignments_buffer()
  short_outputs = []
  
  prepare_outfns(out_dir)

  qf = 0

  tot_reads = util.line_count(inp_fn1)
  
  timer = util.Timer(total = tot_reads)

  i = -1

  print "OPENING FILES"
  with open(inp_fn1) as f1:
   with open(inp_fn2) as f2:
     while 1:
       i+=1

       try:
           r2_l = f2.next()
           r1_l = f1.next()
       except StopIteration as e:
           break
   
       if i % 4 == 1 :
           read1 = r1_l
           read2 = r2_l
           
       if i % 4 == 3:
            if quality(r2_l) < 28 or quality(r1_l) < 28: 
               qc_rejection_count+=1
               continue
               
            #PREV SEQUENCE:
            #r1_library_format = "AGNNNNNNNNNNNNNNNNNNNNNNNN"
            #is just not showing up. no AG....

            r1_library_format = "NNNNNNNNNNNNNNNNNNNNNNNN"
            r2_umi_format = "NNNNNNNNNNNNNNN"

            r1_library_constant = "caagatccgccacaacatcgag".upper()
            r1_library_start = len(r1_library_constant) 
            r1_short_skipped_format = "gacggcagcg".upper()
            r1_post_umi_formatb = "tgattacacatat".upper()
            r1_post_umi_formatc = "acggcagcgtgcag".upper()
        
            
            r2_library_constant = "ggggtgttctgctggtagtggtc".upper()
            r2_umi_start = len(r2_library_constant)
            r2_umi_format = "NNNNNNNNNNNNNNNNNN"
            r2_post_umi_format = "ggcgagctgcacgctgc".upper()
        

            try:
                a1_offset = read1.upper().index(r1_library_constant.upper())
            except Exception, e:
                print "C1 REGIONS FAILED TO EXACTLY MATCH"
                read_constant_rejection_count+=1
                a1_offset = None
                continue
            try:
                a2_offset = read2.upper().index(r2_library_constant.upper())
            except Exception, e:
                print "C2 REGIONS FAILED TO EXACTLY MATCH"
                a2_offset = None
                read_constant_rejection_count+=1
                continue
                #a2_umi_end =None

                
    
            read1_const = read1[a1_offset:a1_offset+r1_library_start]
            read1_library_content = read1[a1_offset+r1_library_start:][:len(r1_library_format)]
            read1_short_skipped_content = read1[a1_offset+r1_library_start:][:len(r1_short_skipped_format)]
            read1_extended_content = read1[a1_offset+r1_library_start:]
            
            read2_const = read1[a2_offset:a2_offset+r2_umi_start]
            read2_umi_content = read2[a2_offset+r2_umi_start:][:len(r2_umi_format)]
            read2_extended_content = read2[a2_offset+r2_umi_start:]

            r1_ag = read1_library_content[0:2]

            post_umi_match1b = re.compile(r1_post_umi_formatb).search(read1_extended_content)
            post_umi_match1c = re.compile(r1_post_umi_formatc).search(read1_extended_content)
            post_umi_match2 = re.compile(r2_post_umi_format).search(read2_extended_content)
            
            
            
            #print read2
            print read1_library_content

            
            read1_type = None
            if not ( post_umi_match1b  or post_umi_match1c):
                read_constant_rejection_count+=1
                read1_type= "cryptic"
            elif post_umi_match1b:
                read1_type ="b"
            else:
                read1_type="c"


            #check to see that the splice acceptor is in the right position 
            #and that the read1 constant sequence aligned

            umi = read2_umi_content
            lib = read1_library_content
            lib_extended = read1_extended_content[:len(r1_library_format)+20]

            if umi.count("N") != 0 or lib.count("N") != 0: continue
            #global target_latternames
            #no longer check for perfect matches. Just align
            exp = target_latternames.get(lib,None)

            


            cand_idxs = find_best_designed_target(lib, lsh_dict)
            if len(cand_idxs) == 0:
               print "rejecting for no good match"
               nolib_rejection_count+= 1
               continue


            best_idx = cand_idxs[0]
            #extends a target alignment region to include an extra 20 bases to anchor the alignment for long r1 deletions
            target_alignment_region = names_targets[best_idx]+"tgattacacatatagacacg".upper()
            align = pairwise2.align.localms(target_alignment_region, read1_extended_content, 2, -1, -5, -.1)[0]

            score = align[2]
            if float(score < 65):
                num_bad_aligns +=1
                #continue
            
            output_complete = """>1\n{0}\n{1}\n{2}\n{3}\n{4}""".format(umi,best_idx,align[2],"\n".join(format_alignment(*align).splitlines()[:3]),"read1_type:{0}".format(read1_type))
            output_short = (umi,best_idx,read1_type,score)

            umis_alignments_buffer[umi].append(output_complete)
            short_outputs.append(output_short)
            accepted_count+=1

            if i % int(tot_reads / 1000) < 4 and i > 1:
             print "FLUSHING!"
             # Flush alignment buffer
             flush_tuples(umis_alignments_buffer, out_dir)
             print len(umis_alignments_buffer.keys())

             # Stats for the curious
             with open(stdout_fn, 'a') as outf:
               outf.write('Time: %s\n' % (datetime.datetime.now()))
               outf.write('Progress: %s\n' % (i / int(tot_reads / 100)) )
               outf.write('Quality filtered pct: %s\n' % (qf / (i/4)))
               outf.write("accepted {0}, rejected {1} bad read1, {2} bad lib\n".format(  accepted_count,read1_rejection_count,nolib_rejection_count   ))

            timer.update()
           
  # Final flush
  flush_tuples(umis_alignments_buffer, out_dir)
   
  # Stats for the curious
  with open(stdout_fn, 'a') as outf:
               outf.write('DONE! sorting short outputs')

  
  short_outputs = sorted(short_outputs, key=lambda x:x[0])
  with open(os.path.join(out_dir,"short.csv"),"w" ) as f:
    for so in short_outputs:
          f.write(", ".join([str(e) for e in so]) + "\n")  

  # Stats for the c`urious
  with open(stdout_fn, 'a') as outf:
               outf.write('COMPLETE \n')

  
  return


##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
  print 'Generating qsub scripts...'
  qsubs_dir = _config.QSUBS_DIR + NAME + '/'
  util.ensure_dir_exists(qsubs_dir)
  qsub_commands = []

  num_scripts = 0

  for _nm in ["190124Gif_D19-{0}".format(i) for i in range(557,560)+range(566,569)]:
    for _split in range(15):
      command = 'python %s.py %s %s' % (NAME, _nm, _split)
      script_id = NAME.split('_')[0]

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, _nm, _split)
      with open(sh_fn, 'w') as f:
        f.write('#!/bin/bash\n%s\n' % (command))
      num_scripts += 1

      # Write qsub commands
      qsub_commands.append('qsub -m e -wd %s %s' % (_config.SRC_DIR, sh_fn))

  # Save commands
  with open(qsubs_dir + '_commands.txt', 'w') as f:
    f.write('\n'.join(qsub_commands))

  print 'Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir)
  return

@util.time_dec
def main(nm = '', split = ''):
  print NAME  

  if nm == '' and split == '':
    gen_qsubs()
    return

  # Function calls
  matchmaker(nm, split) 
  return


if __name__ == '__main__':
  if len(sys.argv) > 2:
    main(nm = sys.argv[1], split = sys.argv[2])
  else:
    main()
