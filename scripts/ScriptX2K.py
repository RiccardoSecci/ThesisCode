# Import modules
import http.client
import json

##### Function to run X2K, the first transcription factor enrichment tool chosen.
### Input: a Python list of gene symbols
### Output: a dictionary containing the results of X2K, ChEA, G2N, KEA.
def run_X2K(input_genes, options={}):
	# Open HTTP connection
	conn = http.client.HTTPConnection("amp.pharm.mssm.edu")  #

	# Get default options
	default_options = {'text-genes': '\n'.join(input_genes), 'included_organisms': 'human',
		       'TF-target gene background database used for enrichment': 'ChEA 2016', 'path_length': 2, 'min_network_size': 50,
		       'min_number_of_articles_supporting_interaction': 2,
		       'max_number_of_interactions_per_protein': 200, 'max_number_of_interactions_per_article': 100,
		       'enable_Biocarta': True, 'enable_BioGRID': True, 'enable_DIP': True, 'enable_InnateDB': True, 'enable_IntAct': True, 'enable_KEGG': True,
		       'enable_MINT': True, 'enable_ppid': True, 'enable_SNAVI': True, 'number_of_results': 50,
		       'sort transcription factors by': 'combined score', 'sort kinases by': 'combined score', 'kinase interactions to include': 'ARCHS4'}

	 # Update options
	for key, value in options.items():
		if key in default_options.keys() and key != 'text-genes':
		    default_options.update({key: value})

	# Get payload
	boundary = "----WebKitFormBoundary7MA4YWxkTrZu0gW"
	payload = ''.join(
	['--' + boundary + '\r\nContent-Disposition: form-data; name=\"{key}\"\r\n\r\n{value}\r\n'.format(**locals())
	 for key, value in default_options.items()]) + '--' + boundary + '--'



	# Get Headers
	headers = {
	'content-type': "multipart/form-data; boundary=" + boundary,
	'cache-control': "no-cache"
	}

	# Initialize connection
	conn.request("POST", "/X2K/api", payload, headers)

	# Get response
	res = conn.getresponse()

	# Read response
	data = res.read().decode('utf-8')

	# Convert to dictionary
	x2k_results = {key: json.loads(value) if key != 'input' else value for key, value in json.loads(data).items()}

	# Clean results
	x2k_results['ChEA'] = x2k_results['ChEA']['tfs']
	x2k_results['G2N'] = x2k_results['G2N']['network']
	x2k_results['KEA'] = x2k_results['KEA']['kinases']
	x2k_results['X2K'] = x2k_results['X2K']['network']

	# Return results
	return x2k_results

