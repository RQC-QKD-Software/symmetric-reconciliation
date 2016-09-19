import numpy as np
from scipy.sparse import dok_matrix as sparse_matrix
from scipy.sparse import find as sparse_find
from time import time
import random
from numpy import zeros, ceil, floor, copy, log2, arange, mean, array, sign


def generate_key(length):
	"""
	Generate random key of length 'length'
	"""
	return np.random.randint(0, 2, (1, length))[0]

def generate_key_zeros(length):
	"""
	Generate key with zeors only of length 'length'
	"""
	return np.zeros(length, dtype=np.int)
    
def add_errors(a, error_prob):
    """
    Flip some values (1->0, 0->1) in 'a' with probability 'error_prob'
    """
    error_mask = np.random.choice(2, size=a.shape, p=[1.0-error_prob, error_prob])
    return np.where(error_mask, ~a+2, a)

def add_errors_prec(a, error_prob):
    """
    Add precisely 'error_prob'*length('a') errors in key 'a' 
    """
    len_a = len(a)
    n_er = int(round(len_a*error_prob))
    list1 = range(0,len_a)
    list2 = random.sample(list1, n_er)
    K_cor = a.copy()
    for i in list2:
        K_cor[i]=1-K_cor[i]
    return K_cor

def h_b(x):
    """
    Binary entropy function of 'x'
    """
    if x > 0:
        return -x*np.log2(x)-(1-x)*np.log2(1-x)
    elif x == 0:
        return 0
    else:
        print("Incorrect argument in binary entropy function")
		
def choose_sp(qber, f, R_range, n):
	'''
	Choose appropriate rate and numbers of shortened and punctured bits
	'''
	def get_sigma_pi(qber, f, R):
		pi = (f*h_b(qber)-1+R)/(f*h_b(qber)-1)
		sigma = 1-(1-R)/f/h_b(qber)
		return max(0,sigma), max(0,pi) 

	delta_min = 1
	R_min = None
	for R in R_range:
		sigma_c, pi_c = get_sigma_pi(qber, f, R)
		delta_c = max(sigma_c,pi_c)
		if delta_c < delta_min:
			delta_min = delta_c; sigma_min = sigma_c; pi_min = pi_c; R_min = R
	if R_min is not None:
		return R_min, int(floor(n*sigma_min)), int(ceil(n*pi_min))

def generate_sp(s_n, p_n, k_n, p_list=None):
	'''
	Generates 's_n', 'p_n' and 'k_n' positions of shortened ('s_pos'), punctured ('p_pos') and key ('k_pos') symbols correspondingly. 
	Punctured symbols are taken from 'p_list' if it is not None. 
	If it is 'p_list' is None of 'p_n' is larger than number of elements in 'p_list', then they are token from the whole key.
	'''
	n = k_n+s_n+p_n  # length of total key
	all_pos = range(int(n)) # array of all indices
	if p_list is None:
		punct_list = all_pos
	elif p_n <= len(p_list):
		punct_list = p_list
	else:
		punct_list = all_pos # taking all postions if it not enough elements in p_list
	if p_n>len(punct_list):
		print 'Error with dimensions p_n=',p_n,'but length of punct_list is',len(punct_list)
	p_pos = np.sort(random.sample(punct_list, p_n)) 
	all_pos1 = np.setdiff1d(all_pos,p_pos)
	s_pos = np.sort(random.sample(all_pos1, s_n)) 
	k_pos = np.setdiff1d(all_pos1, s_pos) 
	return s_pos, p_pos, k_pos

def extend_sp(x, s_pos, p_pos, k_pos):
    '''
    Construct extended key 'x' with shortened/punctured/key bits in positions 's_pos'/'p_pos'/'k_pos'
    '''
    k_n = len(k_pos)
    s_n = len(s_pos)
    p_n = len(p_pos)
    if len(x) != len(k_pos):
        print "Error with dimensions in key and k_pos"
    n = k_n+s_n+p_n # length of extended key
    x_ext = generate_key(n)
    if s_n > 0:
        x_ext[s_pos] = 0
    if p_n > 0:
        x_ext[p_pos] = generate_key(p_n)
    x_ext[k_pos] = x
    return x_ext

def encode_syndrome(x, s_y_joins):
	"""
	Encode vector 'x' with sparse matrix, characterized by 's_y_joins'
	"""
	m = len(s_y_joins)
	s = generate_key_zeros(m)
	for k in range(m):
		s[k] = (sum(x[s_y_joins[k]])%2)
	return np.array(s)

def decode_syndrome_minLLR(y, s, s_y_joins, y_s_joins, qber_est, s_pos, p_pos, k_pos, r_start=None, max_iter=300, x=None, show=1, discl_n=20, n_iter_avg_window=5):
	"""
	INPUT
	'y' is decoding vector. 
	's' is syndrome.
	's_y_joins' and 'y_s_joins' is parity-check matrix information.
	'qber_est' is an estimated level of QBER.
	's_pos'/'p_pos'/'k_pos' stands for positions of shortened/punctured/key bits.
	'r_start' is vector of predefined LLRs.
	'max_iter' is maximal number of iterations in the general cycle.
	'x' is true vecotor decoding vector (for comparison and check of decoding convergence).
	'show' is parameter for additional output: 1 -- output of decoding results, 2 -- output for each iteration.
	'discl_n' is number of bits disclosed in each additional communication round.
	'n_iter_avg_window' is number of iterations for mean LLR averaging for the stop of procedure 
	OUTPUT
	'z' -- decoded vector.
	'minLLR_inds' -- indices of symbols with minimal LLRs.
	"""
	def h_func(x, mode=0):
		"""
		Approximation of log(np.abs(np.exp(x)-1)) for 'mode'=0
		"""
		if mode == 0:
			if x < -3:
				return 0
			elif x < -0.68:
				return -0.25*x-0.75
			elif x < -0.27:
				return -2*x-1.94
			elif x < 0:
				return -8*x-3.56
			elif x < 0.15:
				return 16*x-4
			elif x < 0.4:
				return 4*x-2.2
			elif x < 1.3:
				return 2*x-1.4
			else:
				return x-0.1
		else:
			return np.log(np.abs(np.exp(x)-1))

	def core_func(x, y, mode=1):
		'''
		Core function () for computation of LLRs. 
		'x' and 'y' are arguments. 
		'mode' is approximation method: 0 - piecewise, 1 - table, 2 - exact,  
		'''
		def g_func_piecewise(x):
			"""
			Approximation of log(1+exp(-x)) by linear interpolation between points
			"""
			if x < 0.5:
				return -x/2+0.7
			elif x < 1.6:
				return -x/4+0.575
			elif x < 2.2:
				return -x/8+0.375
			elif x < 3.2:
				return -x/16+0.2375
			elif x < 4.4:
				return -x/32+0.1375
			else:
				return 0

		def g_func_table(x):
			"""
			Aproximation of log(1+exp(-x)) by tabulated values
			"""
			if x < 0.196:
				return 0.65
			elif x < 0.433:
				return 0.55
			elif x < 0.71:
				return 0.45
			elif x < 1.105:
				return 0.35
			elif x < 1.508:
				return 0.25
			elif x < 2.252:
				return 0.15
			elif x < 4.5:
				return 0.05
			else:
				return 0
		if mode == 0:
			return np.sign(x)*np.sign(y)*min(abs(x),abs(y))+g_func_piecewise(abs(x+y))-g_func_piecewise(abs(x-y)) # piecewise
		elif mode == 1:
			return sign(x)*sign(y)*min(abs(x),abs(y))+g_func_table(abs(x+y))-g_func_table(abs(x-y)) # table
		else:
			return sign(x)*sign(y)*min(abs(x),abs(y))+log(1+exp(-abs(x+y)))-log(1+exp(-abs(x-y))) # exact
	
	if not qber_est < 0.5: # Adequate QBER check
		raise ValueError('Aprior error probability must be less than 1/2')
	
	m = len(s_y_joins); n = len(y_s_joins)
	p_n = len(p_pos); s_n = len(s_pos); k_n = len(k_pos)
	v_pos = list(set(p_pos)|set(k_pos))
	
	# Zeroing
	M = np.zeros((m, n)) # Array of messages from symbol nodes to check nodes
	sum_E_abs_mean_hist = [] # Array for mean values of LLRs
	n_iter = 0 # Iteration counter
	
	# Setting initial LLRs:
	if r_start is None:
		r = zeros(n);
		if s_n > 0:
			r[s_pos] = (1-2*y[s_pos])*1000
		if p_n > 0:
			r[p_pos] = 0
		r[k_pos] = (1-2*y[k_pos])*np.log((1-qber_est)/qber_est)
	else:
		r = r_start
		if s_n > 0:
			r[s_pos] = (1-2*y[s_pos])*1000
		
	for j in xrange(m): # Setting initial messages from symbol nodes to check nodes
		M[j, :] = r

	while n_iter < max_iter: # Main cycle
		# Part 1: from check nodes to symbol nodes
		E = np.zeros((m, n)) # Array of messages from check nodes to symbol nodes
		for j in xrange(m): # For all check nodes
			M_cur = M[j][s_y_joins[j]]; M_cur_n = len(M_cur) # All symbol nodes that are connected to current check node and their number
			n_zeros = list(M_cur).count(0.0) # number of zero LLRs
			if n_zeros > 1: # If check node is dead
				E[j,s_y_joins[j]] = np.zeros(M_cur_n) # No messages
			elif n_zeros == 1: # If current check node has one punctured symbol
				E_cur = np.zeros(M_cur_n) # All messages are initializrd with zeros
				M_cur = list(M_cur); zero_ind = M_cur.index(0.0); M_cur.pop(zero_ind) # Excluding zero message
				LS = M_cur[0]
				for k in range(1,M_cur_n-1): # Accumulation of the message
					LS = core_func(LS,M_cur[k])
				E_cur[zero_ind] = LS; E[j,s_y_joins[j]] = E_cur # Filling with nonzero message
			elif n_zeros == 0: # all messages are non zero 
				LS = M_cur[0]
				for k in range(1,M_cur_n):
					LS = core_func(LS,M_cur[k])
				E_cur = zeros(M_cur_n)
				for i1 in range(0,M_cur_n):
					E[j][s_y_joins[j][i1]] = (1-2*s[j])*(h_func(M_cur[i1]+LS)-h_func(M_cur[i1]-LS)-LS) # Computation of messages
			
		# Part 2: from symbol nodes to check nodes    
		sum_E = E.sum(axis=0)+r # Array of sums of messages to symbol nodes (LLRs) 
		z = (1-np.sign(sum_E))/2 # Current decoded message
	
		if (s == encode_syndrome(z, s_y_joins)).all(): # If syndrome is correct
			if np.count_nonzero(z == x) != n:
				print "Convergence error, error positions:"
				print '\n', np.nonzero((z+x)%2)
			if show > 0:
				print 'Done in ', n_iter, 'iters, matched bits:', np.count_nonzero(z == x), '/', n
			return z, None
		if show > 1:
			print 'Matched bits:', np.count_nonzero(z == x), '/', n, 'Mean LLR magnitude:',mean(abs(sum_E[v_pos])), \
			'Averaged mean LLR magnitude:',sum(sum_E_abs_mean_hist[max(0,n_iter-n_iter_avg_window):n_iter])/(min(n_iter,n_iter_avg_window)+10**(-10))

		# Check for procedure stop

		sum_E_abs = list(abs(sum_E))
		sum_E_abs_mean_hist.append(mean(list(abs(sum_E[v_pos]))))

		if n_iter == n_iter_avg_window-1:
			sum_E_mean_avg_old = mean(sum_E_abs_mean_hist)
		if n_iter >= n_iter_avg_window:
			sum_E_mean_avg_cur = sum_E_mean_avg_old+(sum_E_abs_mean_hist[n_iter]-sum_E_abs_mean_hist[n_iter-n_iter_avg_window])/n_iter_avg_window
			if sum_E_mean_avg_cur <= sum_E_mean_avg_old:
				minLLR_inds = []; maxLLR = max(sum_E_abs)
				for cnt in range(discl_n):
					ind = sum_E_abs.index(min(sum_E_abs))
					minLLR_inds.append(ind)
					sum_E_abs[ind] += maxLLR
				return None, minLLR_inds
			else:
				sum_E_mean_avg_old = sum_E_mean_avg_cur

		# Calculating messages from symbol nodes to check nodes  
		M = -E+sum_E
			
		n_iter += 1
		
	minLLR_inds = []; maxLLR = max(sum_E_abs)
	for cnt in range(discl_n):
		ind = sum_E_abs.index(min(sum_E_abs))
		minLLR_inds.append(ind)
		sum_E_abs[ind] += maxLLR
	return None, minLLR_inds

	
def perform_ec(x, y, s_y_joins, y_s_joins, qber_est, s_n, p_n, punct_list=None, discl_n=20, show=0):
	n = len(y_s_joins); m = len(s_y_joins)
	
	s_pos, p_pos, k_pos = generate_sp(s_n, p_n, n-s_n-p_n, p_list=punct_list)

	x_ext = extend_sp(x, s_pos, p_pos, k_pos)
	y_ext = extend_sp(y, s_pos, p_pos, k_pos)

	k_pos_in = copy(k_pos); # For final exclusion
	
	s_x = encode_syndrome(x_ext, s_y_joins)
	s_y = encode_syndrome(y_ext, s_y_joins)

	s_d = (s_x+s_y)%2
	key_sum = (x_ext+y_ext)%2

	e_pat_in = generate_key_zeros(n)

	e_pat, minLLR_inds = decode_syndrome_minLLR(e_pat_in, s_d, s_y_joins, y_s_joins, qber_est, s_pos, p_pos, k_pos, max_iter=100500, x=key_sum, show=show, discl_n=discl_n, n_iter_avg_window=5)

	add_info = 0; com_iters = 0

	while e_pat is None:
		if show>1:
			print 'Additional iteration with p_n=', len(p_pos), 's_n=', len(s_pos), 'k_n=', len(k_pos)
		e_pat_in[minLLR_inds] = (x_ext[minLLR_inds]+y_ext[minLLR_inds])%2
		s_pos = list(set(s_pos) | set(minLLR_inds))
		k_pos = list(set(k_pos) - set(minLLR_inds))
		if p_pos is not None:
			p_pos = list(set(p_pos) - set(minLLR_inds))
		e_pat, minLLR_inds = decode_syndrome_minLLR(e_pat_in, s_d, s_y_joins, y_s_joins, qber_est, s_pos, p_pos, k_pos, r_start=None, max_iter=100500, x=key_sum, show=show, discl_n=discl_n, n_iter_avg_window=5)
		add_info += discl_n 
		com_iters += 1
	
	x_dec = (x_ext[k_pos_in]+e_pat[k_pos_in])%2
	
	ver_check = (x_dec == y).all()
	if not ver_check:
		print "VERIFICATION ERROR"
		#print '\nInitial error pattern:\n', np.nonzero((x_ext+y_ext)%2),'\nFinal error pattern:\n', np.nonzero(e_pat)

	return add_info, com_iters, e_pat[k_pos_in], ver_check


def test_ec(qber, R_range, codes, n, n_tries, f_start=1, show=1, discl_k=1):
	R, s_n, p_n = choose_sp(qber, f_start, R_range, n)
	k_n = n-s_n-p_n
	m = (1-R)*n
	code_params = codes[(R, n)]
	s_y_joins = code_params['s_y_joins']; y_s_joins = code_params['y_s_joins']; punct_list = code_params['punct_list']
	p_n_max = len(punct_list)
	discl_n = int(round(n*(0.0280-0.02*R)*discl_k))
	qber_est = qber
	f_rslt = []
	com_iters_rslt = []
	n_incor = 0
	
	print "QBER = ",qber, "R =", R, "s_n =", s_n, "p_n =", p_n, '(', p_n_max, ')', 'discl_n', discl_n

	for i in range(n_tries):
		print i,
		x = generate_key(n-s_n-p_n)
		y = add_errors(x, qber)
		add_info, com_iters, x_dec, ver_check = perform_ec(x, y, s_y_joins, y_s_joins, qber_est, s_n, p_n, punct_list=punct_list, discl_n=discl_n, show=show)
		f_cur = float(m-p_n+add_info)/(n-p_n-s_n)/h_b(qber)
		f_rslt.append(f_cur)
		com_iters_rslt.append(com_iters)
		if not ver_check:
			n_incor += 1
	print 'Mean efficiency:', np.mean(f_rslt), '\nMean additional communication rounds', np.mean(com_iters_rslt)
	return np.mean(f_rslt), np.mean(com_iters_rslt), R, s_n, p_n, p_n_max, k_n, discl_n, float(n_incor)/n_tries