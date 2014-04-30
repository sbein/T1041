T1041
=====

Analysis code for T1041 test beam run

Directories:  
- include, src : inputs for Makefile to generate shared libraries
- python       : python scipts to run code
- rootscript   : misc root code/scripts, recomend using for user code
- dqm          : dqm utilities, for running online
- doc          : various documentation


** First Steps **  
source setupTB.sh  
Run make

Generate DB of WC runs 
======================
To allow association of WC tracks, you will need the file: wcdb.txt  
Either copy this from your data area or run wcdbGenerator.py

# generate wcdb, adding only newer files if wcdb exists  
python python/wcdbGenerator.py   
# forces full regeneration of wcdb  
python python/wcdbGenerator.py --force  
# flag to specify directory fo WC files  
python python/wcdbGenerator.py -f WC_file_dir  
# print help message  
python python/wcdbGenerator.py -h  


Make TTrees from TB data
========================
Usage:  
Usage: python TBNtupleMaker [OPTION] [PADE_FILE] [PADE_FILE] ...
      -n max_events  : Maximum (requested) number of events to read
                       Always read at least 1 spill
      -d DIR         : Process all padefiles in DIR
                       Overrides and files given on command line list
      -r DIR         : Process all padefiles in DIR, and all subdirectories
                       Overrides and files given on command line list
      -k             : Keep existing root files, ony process new ones


Produces an output file with the same basename as the PADE_FILE, replacing .txt(.bz2) with .root


Example of reconstruction tools  
===============================  
python python/readerExample.py file.root
 

Display channel mapping  
======================  
root -l rootscript/testChannelMap.C+


Simple event display  
====================  
Make plots averaged over all events -- saves .gif files  
python dqmPlots.py filename.root  

Make plots for a single event -- saves .gif files  
python python/dqmPlots.py filename.root eventNumber  

Simple wave form scanner  
========================  
python python/waveViewer.py rec_capture_20140413_003028.root


Make png plots and html page from a root file
=============================================
python python/root2html.py file.root  
  png files and an index.html file are produced in directory called "file", 
the root file name w/o the .root


Run Jordan's WC display  
=======================  
python dqm/WC_Reader.py file.root




git cheatsheet
==============

clone a working copy

HTTPS clone: git clone https://github.com/ShashlikTB/T1041.git

SSH clone (do this to allow password-less check in w/ ssh key):  
git clone git@github.com:ShashlikTB/T1041.git


update/add:
* git add <file>
* git commit -m "comment"

check it in:
* git push origin master


