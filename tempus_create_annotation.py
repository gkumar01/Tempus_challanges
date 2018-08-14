#!/usr/bin python

import os
import sys
import vcf
import requests
from requests.auth import HTTPDigestAuth
import json
from optparse import OptionParser


def get_allele_freq(json_data):
	"""
	 Extract allele frequency information from json data dictionary
	"""
	
	freq='.'
	try:
		freq = json_data['variant']['allele_freq']
		
	except KeyError:
		pass
		
	return freq
	
def get_snp_id(json_data):
	"""
	Extract snp id from jason data
	"""
	
	rsid='.'
	try:
		rsid = json_data['variant']['rsid']
	except KeyError:
		pass
		
	return rsid

def get_ensemble_geneid(json_data):
	"""
	Extract ensemble gene from json data	
	"""
	ensemble_id = '.'
	try:
		ensemble_id = '/'.join(json_data['variant']['genes'])
	except KeyError:
		pass
		
	return ensemble_id

def request_variant_info(variant_id):
	"""
	Request variant information
	"""
	url_base = "http://exac.hms.harvard.edu"
	resp = requests.get(url_base+'/rest/variant/'+variant_id)
	
	json_data={}
	
	#If call successful, reponse code will be 200(OK)
	if(resp.ok):
		# Loading the response data into a dict variable
		json_data = json.loads(resp.content)
		
	else:
		#failed call result in http error code with description
		resp.raise_for_status()
		pass
		
	return json_data

def get_variant_type(lst=[]):
	
	"""
	For annotation select variant 'del'  from multiple variant sub types 	
	"""
	
	#print(lst)

	new_list = [ ]
	
	if len(lst) > 1 :
		new_list = [x for x in lst if x == 'del']
	else:
		new_list = lst
	
	if len(new_list) == 0:
		new_list=lst[:1]
	
	return new_list[0]

def get_params():
	"""
	set command line arguments
	"""
	parser = OptionParser()
	parser.add_option("-i", "--input_vcf_file",
					dest = 'input_vcf_file',
					default = None,
					help = 'The input vcf file, e.g. test.vcf'
					)
	(options, args) = parser.parse_args()
	return options, args


def run(options,args):
	"""
	 parse vcf file and extract relevant field for annotation
	"""
	#print('%s %s %s' %('hello',options,args))
	
	vcf_file = options.input_vcf_file
	
	#print(os.path.isfile(vcf_file))
	
	#check if vcf file exist
	if os.path.isfile(vcf_file) :
		vcf_reader = vcf.Reader(open(vcf_file, 'r'))
		
		header = ('Variant','Type of Variation',\
				'site_coverage/depth',\
				'read supporting variants',\
				'percentage_variant_read_coverage|'+
				'percentage_reference_read_coverage',\
				 'Allele_Frequency',\
				 'Ensemble_gene_id',\
				 'SNP_id'
				)
		print (",".join(header))
		
		#count=0
		
		#Loop through lines/records for annotation
		for record in vcf_reader:
			
			#create variant-id as per ExAC convention
			variant_str = [record.CHROM,str(record.POS),record.REF]
			#variant_str.extend(record.ALT)
			variant_str.append('/'.join(map(str,record.ALT)))
			variant_name = "-".join(map(str,(variant_str)))
			#print(variant_name)
			
			variant_json_info={}
			
			#request API call for ALT (variant) with single entry
			if (len(record.ALT) == 1):
				variant_json_info = request_variant_info(variant_name)
			
			#parse json data to get allele frequency, gene_id and snp_id
			allele_frequency = get_allele_freq(variant_json_info)
			gene_id = get_ensemble_geneid(variant_json_info)
			snp_id = get_snp_id(variant_json_info)
			#print('%s,%s,%s' %(allele_frequency,gene_id,snp_id))
	
			#count +=1
			#if(count > 1000):
			#	sys.exit('T')
			
			#variant_type = record.INFO['TYPE']
			#variant_type = get_variant_type(record.INFO['TYPE'])
			variant_subtype = record.var_subtype
			
			#If variant subtype is unknow check for multiple annotation
			if variant_subtype == "unknown":
				variant_type = get_variant_type(record.INFO['TYPE'])
				variant_subtype = variant_type
			
			#Read depth/coverage for the site/locus
			site_coverage = float(record.INFO['DP'])
						
			#Variant/Alternative allele observation
			variant_coverage_depth = float (record.INFO['AO'][0])
			
			#Reference allele observation
			ref_coverage_depth =  float(record.INFO['RO'])
						
			#calculate percentage of reads covering variant allele
			pct_variant_read_cov = round(100 *  (variant_coverage_depth/site_coverage),3)
			
			pct_reference_read_cov = round(100 *  (ref_coverage_depth/site_coverage),3)
			
			print ("%s,%s,%s,%s,%s|%s,%s,%s,%s"\
					%(variant_name,\
					variant_subtype,\
					site_coverage,\
					variant_coverage_depth,\
					pct_variant_read_cov,\
					pct_reference_read_cov,\
					allele_frequency,\
					gene_id,\
					snp_id
					)
				)
	else:
		text_message ='No such file exists:' + vcf_file + ' check path!!'
		sys.exit(text_message)
	

if __name__ == '__main__':
		
	"""
	USAGE: python tempus_create_annotation.py -i Challenge_data.vcf > tempus_variant_annotation.csv
	"""
	options,args = get_params()
	run(options,args)
