#  MIT License
#  
#  Copyright (c) 2016 Matthew Urtnowski
#  
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#  
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#  
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#
#  This script analyzes a 23andMe raw data file that can be downloaded from
#  https://www.23andme.com/you/download/ and compares it to many of the genes
#  SNPedia flags as related to Alzheimer's diseases.  Genes which SNPedia 
#  would indicate an increased risk of Alzheimer's will be highlighted in red.
#
#  This script makes no guarantee of the accuracy of data provided by 23andMe
#  and SNPedia.  This script also most likely contains errors itself as gene
#  analysis can sometime be ambiguous.  23andMe and SNPedia data often to not
#  line up exactly.  See https://www.snpedia.com/index.php/Orientation
#
#  I created this script because 23andMe no longer provides medical analysis
#  of health risks and even when it did, it was very limited in which gene
#  it examined for early onset Alzheimers.
#
#  This author of this script is not a medical expert.
#
#  DO NOT USE THIS SCRIPT TO MAKE MEDICAL AND HEALTH DECISIONS.
#  PLEASE CONSULT A MEDICAL EXPERT IF YOU HAVE CONCERNS ABOUT YOUR HEALTH OR
#  CONCERNS ABOUT ALZHEIMER'S DISEASE.

import csv
import sys
import os.path

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
if len(sys.argv) > 1:
	filename = sys.argv[1]
else:
	print "filename not set"
	exit(1)
	
if not os.path.isfile(filename):
	print "file " + filename + " not found"
	exit(1)

dna = {
	"rs429358": None,
	"rs7412": None,
	"rs199768005": None,
	"rs449647": None,
	"rs405509": None,
	"rs63750847": None,
	"rs145999145": None,
	"rs75932628": None,
	"rs193922916": None,
	"rs63749885": None,
	"rs63749891": None,
	"rs63750577": None,
	"rs63750599": None,
	"rs63750815": None,
	"rs63750886": None,
	"rs63750900": None,
	"rs63751037": None,
	"rs63751141": None,
	"rs63751144": None,
	"rs63751163": None,
	"rs63751223": None,
	"rs63751229": None,
	"rs63751235": None,
	"rs63751320": None,
	"rs661": None,
	"rs63750218": None,
	"rs63750231": None,
	"rs63750526": None,
	"rs63750590": None,
	"rs63749824": None,
	"rs28936379": None,
	"rs28936380": None,
	"rs63749851": None,
	"rs63750048": None,
	"rs63750110": None,
	"rs63750215": None,
	"rs63749884": None,
	"rs63750666": None,
	"rs61757781": None
}

def simpleCheck(snp, geno):
	if dna[snp] == None:
		print bcolors.WARNING + snp + " unknown" + bcolors.ENDC
	elif dna[snp] == geno:
		print bcolors.FAIL + snp + " mutation" + bcolors.ENDC
		print bcolors.OKBLUE + snp + ":   " + dna[snp] + bcolors.ENDC
		print bcolors.OKBLUE + "http://snpedia.com/index.php/" + snp + bcolors.ENDC
	else:
		print bcolors.OKGREEN + snp + " mutation OK" + bcolors.ENDC

with open(filename) as tsv:
	for line in csv.reader(tsv, dialect="excel-tab"): #You can also use delimiter="\t" rather than giving a dialect.
		if line[0] in dna:
			dna[line[0]] = line[3]

## ApoE-4 Check
if dna["rs429358"] == None or dna["rs7412"] == None:
	if dna["rs429358"] == None:
		print bcolors.WARNING + "rs429358 unknown" + bcolors.ENDC
	if dna["rs7412"] == None:
		print bcolors.WARNING + "rs7412 unknown" + bcolors.ENDC
elif dna["rs429358"] == "CT" and dna["rs7412"] == "CC":
	print bcolors.FAIL + "APOE4 Heterogeneous" + bcolors.ENDC
	print bcolors.OKBLUE + "rs429358: " + dna["rs429358"] + bcolors.ENDC
	print bcolors.OKBLUE + "rs7412:   " + dna["rs7412"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/APOE" + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/Rs429358" + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/Rs7412" + bcolors.ENDC
elif dna["rs429358"] == "CT" and dna["rs7412"] == "CC":
	print bcolors.FAIL + "APOE4 Homogeneous" + bcolors.ENDC
	print bcolors.OKBLUE + "rs429358: " + dna["rs429358"] + bcolors.ENDC
	print bcolors.OKBLUE + "rs7412:   " + dna["rs7412"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/APOE" + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/Rs429358" + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/Rs7412" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "APOE4 OK" + bcolors.ENDC
	
## ApoE rs199768005 mutation
if dna["rs199768005"] == None:
	print bcolors.WARNING + "rs199768005 unknown" + bcolors.ENDC
elif dna["rs199768005"] == "AT":
	print bcolors.OKGREEN + "APOE4 rs199768005 mutation - Reduced risk of Alzheimers" + bcolors.ENDC
	print bcolors.OKBLUE + "rs199768005:   " + dna["rs199768005"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/Rs199768005" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "rs199768005 mutation Normal" + bcolors.ENDC

## rs449647/rs405509 mutation
if dna["rs449647"] == None or dna["rs405509"] == None:
	if dna["rs449647"] == None:
		print bcolors.WARNING + "rs449647 unknown" + bcolors.ENDC
	if dna["rs405509"] == None:
		print bcolors.WARNING + "rs405509 unknown" + bcolors.ENDC
elif dna["rs449647"] == "TT" and dna["rs405509"] == "CC":
	print bcolors.FAIL + "rs449647/rs405509 mutation" + bcolors.ENDC
	print bcolors.OKBLUE + "rs449647: " + dna["rs449647"] + bcolors.ENDC
	print bcolors.OKBLUE + "rs405509: " + dna["rs405509"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs449647" + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs405509" + bcolors.ENDC
elif dna["rs449647"] == "TA" and dna["rs405509"] == "AA":
	print bcolors.FAIL + "rs449647/rs405509 mutation" + bcolors.ENDC
	print bcolors.OKBLUE + "rs449647: " + dna["rs449647"] + bcolors.ENDC
	print bcolors.OKBLUE + "rs405509: " + dna["rs405509"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs449647" + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs405509" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "rs449647/rs405509 mutation OK" + bcolors.ENDC
	
## rs63750847 mutation
if dna["rs63750847"] == None:
	print bcolors.WARNING + "rs63750847 unknown" + bcolors.ENDC
elif dna["rs63750847"] == "TT" or dna["rs63750847"] == "TC":
	print bcolors.OKGREEN + "rs63750847 mutation - Reduced risk of Alzheimers" + bcolors.ENDC
	print bcolors.OKBLUE + "rs63750847:   " + dna["rs63750847"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs63750847" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "rs63750847 mutation Normal" + bcolors.ENDC
	
## rs145999145 mutation
if dna["rs145999145"] == None:
	print bcolors.WARNING + "rs145999145 unknown" + bcolors.ENDC
elif dna["rs145999145"] == "AA" or dna["rs145999145"] == "AG":
	print bcolors.FAIL + "rs145999145 mutation" + bcolors.ENDC
	print bcolors.OKBLUE + "rs145999145:   " + dna["rs145999145"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs145999145" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "rs145999145 mutation OK" + bcolors.ENDC
	
## rs75932628 mutation
if dna["rs75932628"] == None:
	print bcolors.WARNING + "rs75932628 unknown" + bcolors.ENDC
elif dna["rs75932628"] == "CT" or dna["rs75932628"] == "TT":
	print bcolors.FAIL + "rs75932628 mutation" + bcolors.ENDC
	print bcolors.OKBLUE + "rs75932628:   " + dna["rs75932628"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs75932628" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "rs75932628 mutation OK" + bcolors.ENDC
	
simpleCheck("rs193922916", "AA")	
	
## rs63750847 mutation
if dna["rs63750847"] == None:
	print bcolors.WARNING + "rs63750847 unknown" + bcolors.ENDC
elif dna["rs63750847"] == "AA" or dna["rs63750847"] == "AG":
	print bcolors.OKGREEN + "rs63750847 mutation - Reduced risk of Alzheimers" + bcolors.ENDC
	print bcolors.OKBLUE + "rs63750847:   " + dna["rs63750847"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs63750847" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "rs63750847 mutation Normal" + bcolors.ENDC

simpleCheck("rs63749885", "CT")	
simpleCheck("rs63749891", "CG")	
simpleCheck("rs63750577", "CT")
simpleCheck("rs63750599", "CT")
simpleCheck("rs63750815", "GT")
simpleCheck("rs63750886", "CG")
simpleCheck("rs63750900", "AG")
simpleCheck("rs63751037", "AG")
simpleCheck("rs63751141", "CG")
simpleCheck("rs63751144", "AC")
simpleCheck("rs63751163", "CT")
simpleCheck("rs63751223", "CG")
simpleCheck("rs63751229", "CT")
simpleCheck("rs63751235", "CG")
simpleCheck("rs63751320", "AC")

## rs661 mutation
if dna["rs661"] == None:
	print bcolors.WARNING + "rs661 unknown" + bcolors.ENDC
elif dna["rs661"] == "AA" or dna["rs661"] == "AG":
	print bcolors.FAIL + "rs661 mutation" + bcolors.ENDC
	print bcolors.OKBLUE + "rs661:   " + dna["rs661"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs661" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "rs661 mutation OK" + bcolors.ENDC

simpleCheck("rs63750218", "CT")
simpleCheck("rs63750231", "AC")
simpleCheck("rs63750526", "AC")
simpleCheck("rs63750590", "AG")
simpleCheck("rs63749824", "CT")

simpleCheck("rs28936379", "AG")
simpleCheck("rs28936380", "CG")
simpleCheck("rs63749851", "AC")
simpleCheck("rs63750048", "CT")
simpleCheck("rs63750110", "AC")
simpleCheck("rs63750215", "AT")
simpleCheck("rs63749884", "AG")
simpleCheck("rs63750666", "CT")

## rs61757781 mutation
if dna["rs61757781"] == None:
	print bcolors.WARNING + "rs61757781 unknown" + bcolors.ENDC
elif dna["rs61757781"] == "AG" or dna["rs61757781"] == "GG":
	print bcolors.FAIL + "rs61757781 mutation" + bcolors.ENDC
	print bcolors.OKBLUE + "rs61757781:   " + dna["rs61757781"] + bcolors.ENDC
	print bcolors.OKBLUE + "http://snpedia.com/index.php/rs61757781" + bcolors.ENDC
else:
	print bcolors.OKGREEN + "rs61757781 mutation OK" + bcolors.ENDC
