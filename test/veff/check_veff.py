import sys

'''
	Expected number of global events passing
	And summed weight passing
'''
expected_global_pass = 3
expected_total_weight = 1.99067
expected_total_weight_sigma = 0.0001


'''
	Some utility functions
	The first finds (the first instance of) a string inside a file
	The second finds if a string contains a NAN
'''

def search_file(filename, search_string):
	with open(filename, 'r') as read_obj:
		for line in read_obj:
			if search_string in line:
				return line.rstrip()

def contains_nan(search_string):
	if 'nan' in search_string or \
		'NAN' in search_string or \
		'NaN' in search_string:
		return True


'''
	1: check if the output file exist at all
'''
fin = sys.argv[1]
try:
	with open(fin, 'r') as read_obj:
		read_obj.close()
	print('Output log file exists. File exists test will pass.')
except:
	print('Output log file does not exist. File exists test will fail.')
	sys.exit(-1) # fail out



'''
	2: check if AraSim finished running
	We will check for the 'test is 0' line in the AraSim output
	So, we split on 'test is', then cast the third element to an int
'''
internal_test = search_file(fin, 'test is')
if internal_test is None:
	print('AraSim output file is incomplete. AraSim finished test will fail.')
	sys.exit(-1) # fail out
else:
	internal_test = int(internal_test.split(' ')[2])
	if internal_test == 0:
		print('AraSim reports a successful run. AraSim finished test will pass.')
	else:
		print('AraSim output file is incomplete. AraSim finished test will fail.')
		sys.exit(-1)



'''
	3: check how many events passed globally (integer)
	The AraSim output looks like 
		Total_Global_Pass : #
	So, we first search the file for that string
	Split it on the colon, and cast the object after the colon to an int
'''

global_pass_string = search_file(fin, 'Total_Global_Pass')
if global_pass_string is None:
	print('Total_Global_Pass is missing from AraSim output file. Test will fail.')
	sys.exit(-1)
else:
	has_nan = contains_nan(global_pass_string)
	if has_nan:
		print('Global pass string contains a NAN. Total_Global_Pass test will fail.')
		sys.exit(-1)
	else:
		global_pass = int(global_pass_string.split(":")[1])
		if global_pass == expected_global_pass:
			print('Global pass: {}. Expected {}. Total_Global_Pass test will pass.'.format(
				global_pass, expected_global_pass))
		else:
			print('Global pass: {}. Expected {}. Total_Global_Pass test will fail.'.format(
			global_pass, expected_global_pass))
			sys.exit(-1)


'''
	4: check summed triggering weight (float)
	The AraSim output looks like 
		Total_Weight : #
	So, we first search the file for that string
	Split it on the colon, and cast the object after the colon to a float
'''
total_weight_string = search_file(fin, 'Total_Weight')
if total_weight_string is None:
	print('Total_Weight is missing from AraSim output file. Test will fail.')
	sys.exit(-1)
else:
	has_nan = contains_nan(total_weight_string)
	if has_nan:
		print('Total weight string contains a NAN. Test will fail.')
		sys.exit(-1)
	else:
		total_weight = float(total_weight_string.split(":")[1])
		if abs(total_weight - expected_total_weight) > expected_total_weight_sigma:
			print('Total weight: {:.4f}. Expected {:.4f} +- {:.4f}. Total_Weight test will fail.'.format(
				total_weight, expected_total_weight, expected_total_weight_sigma))
			sys.exit(-1)
		else:
			print('Total weight: {:.4f}. Expected {:.4f} +- {:.4f}. Total_Weight test will pass.'.format(
				total_weight, expected_total_weight, expected_total_weight_sigma))

