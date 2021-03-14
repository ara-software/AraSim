import sys


'''
	Keep track of if tests pass or not
'''
pass_file_exists_test = False
pass_arasim_finished_test = False
pass_global_count_test = False 
pass_total_weight_test = False


'''
	Expected number of global events passing
	And summed weight passing
'''
expected_global_pass = 2
expected_total_weight = 0.998305
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
	1: check if the file exists at all
'''
fin = sys.argv[1]
try:
	with open(fin, 'r') as read_obj:
		read_obj.close()
	print('Output log file exists.')
	pass_file_exists_test = True
except:
	print('Output log file does not exist. Test will fail')



'''
	2: did AraSim finish running
	We will check for the 'test is 0' line in the AraSim output
	So, we split on 'test is', then cast the third element to an int
'''
if pass_file_exists_test:
	internal_test = search_file(fin, 'test is')
	if internal_test is None:
		print('AraSim output file is incomplete. Test will fail.')		
	internal_test = int(internal_test.split(' ')[2])
	if internal_test == 0:
		print('AraSim reports a successful run.')
		pass_arasim_finished_test = True
	else:
		print('AraSim output file is incomplete. Test will fail.')



'''
	3: check how many events passed globally (integer)
	The AraSim output looks like 
		Total_Global_Pass : #
	So, we first search the file for that string
	Split it on the colon, and cast the object after the colon to an int
'''

if pass_arasim_finished_test:
	global_pass_string = search_file(fin, 'Total_Global_Pass')
	if global_pass_string is None:
		print('Total_Global_Pass is missing from AraSim output file. Test will fail.')
	has_nan = contains_nan(global_pass_string)
	if has_nan:
		print('Global pass string contains a NAN. Test will fail.')
	else:
		global_pass = int(global_pass_string.split(":")[1])
		if global_pass == expected_global_pass:
			print("Global pass: {}. Expected {}. Test will pass.".format(
				global_pass, expected_global_pass))
			pass_global_count_test = True
		elif global_pass != expected_global_pass:
			print('Global pass: {}. Expected {}. Test will fail.'.format(
			global_pass, expected_global_pass))


'''
	4: check summed weight triggering (double)
	The AraSim output looks like 
		Total_Weight : #
	So, we first search the file for that string
	Split it on the colon, and cast the object after the colon to a float
'''
if pass_arasim_finished_test:
	total_weight_string = search_file(fin, 'Total_Weight')
	if total_weight_string is None:
		print('Total_Weight is missing from AraSim output file. Test will fail.')
	has_nan = contains_nan(total_weight_string)
	if has_nan:
		print('Total weight string contains a NAN. Test will fail.')
	else:
		total_weight = float(total_weight_string.split(":")[1])
		if abs(total_weight - expected_total_weight) > expected_total_weight_sigma:
			print("Total weight: {:.4f}. Expected {:.4f} +- {:.4f}. Test will fail.".format(
				total_weight, expected_total_weight, expected_total_weight_sigma))
		else:
			print("Total weight: {:.4f}. Expected {:.4f} +- {:.4f}. Test will pass.".format(
				total_weight, expected_total_weight, expected_total_weight_sigma))
			pass_total_weight_test = True

'''
	Check all the tests at the end
	Fail if any of them are passing
'''

if not pass_file_exists_test or not pass_arasim_finished_test or \
	not pass_global_count_test or not pass_total_weight_test:
	print("Some tests failing. Overall test will fail.")
	sys.exit(-1) # fail out

