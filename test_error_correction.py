import error_correction_lib as ec
import numpy as np
from file_utils import codes_from_file

# Choose of the codes pool:
#codes = codes_from_file('codes_4000.txt'); n = 4000
codes = codes_from_file('codes_1944.txt'); n = 1944

# Computing the range of rates for given codes
R_range = []
for code in codes:
	R_range.append(code[0])
print 'R range is', np.sort(R_range)

fname = 'output.txt' # file name for the output
f_start = 1.0 # initial efficiency of decoding
qber_start = 0.03; qber_end = 0.04; qber_step = 0.05 # range of QBERs  
n_tries = 2 # number of keys proccessed for each QBER value

file_output = open(fname,'a')
file_output.write('n = %d, n_tries = %d\n'%(n, n_tries))

for qber in np.arange(qber_start, qber_end, qber_step):
	f_mean, com_iters_mean, R, s_n, p_n, p_n_max, k_n, discl_n, FER = ec.test_ec(qber, R_range, codes, n, n_tries, f_start=f_start, show=2, discl_k=1)
	file_output.write('%8.4f%14.8f%14.8f%14.8f%10d%10d%10d%14d%10d%14.8f\n'%(qber, f_mean, com_iters_mean, R, s_n, p_n, p_n_max, k_n, discl_n, FER))
	file_output.close()