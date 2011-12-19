#!/usr/bin/python

# source object - initial message
# object is broken down into Z >= 1 source blocks - target for a single raptor code application. identified by SBN
# blocks are broken into K source symbols. K is constant for a given object. each symbol has an ESI (encoding symbol identifier)
# size of source symbols: T. (so K*T = block size). 
# sub-blocks are made from EACH block, such that they can be decoded in working memory. N >= 1 subblocks. 
# sublocks have K sub-symbols, of size T'. 

import numpy
import random
import sys
from bitarray import bitarray

# recommended alg values taken from RFC 5053 (basically represents how XOR
# operations are executed, 4 bytes at a time is standard on 32-bit, probably 8
# bytes on 64-bit? if XOR is done on single memory location)
A1 = 4
# min number of symbols per source block, in bytes
Kmin = 1024 
# target size for sub-block, in bytes. let's say 1kb
W = 1024



# one raptor manager is used per object
class RaptorManager:
	def __init__(self, filename, K=1024, debug=True):
		self.debug = debug
		self.f = open(filename, 'rb')
		# number of symbols in each source block
		self.K = int(K)
		# keep a counter of how many blocks get sent out. 
		self.current_block = 0
		# remember how much padding the last block used 
		self.padding_last = None
		self.last_block = False
		# constraint matrix for codes that use pre-coding.
		self.G = None
		if self.debug:
			print "registered " +self.f.name+ "."
	
	def _encode_binary_block(self):
		# encode a chunk of size k from the input file as binary and return

		# the file read() method reads in 1-byte (8 bit) chunks, so make sure K
		# is a multiple of 8, and then divide K by 8 to determine how many bits
		# to read, to keep num symbols (bits) equal to  K. 
		if self.K % 8 != 0:
			sys.stderr.write("error: K must be a byte multiple (ie multiple of 8).")
			return
		n = self.K/8

		block = bitarray()
		start_pos = self.f.tell()
		try:
			block.fromfile(self.f,n)
		except EOFError:
			# if we reach the end of the file, get just the remaining bytes and
			# pad the block as necessary
			end_pos = self.f.tell()
			print "started at post %d, reached EOF at pos %d" % (start_pos, end_pos)
			remaining_bytes = end_pos - start_pos
			self.f.seek(start_pos)
			block = bitarray()
			block.fromfile(self.f,remaining_bytes)
			padding = n - remaining_bytes
			self.padding_last = padding
			print "remaining bytes: %d, of required %d bits or %d bytes. current length = %d, ==> padding = %d" % (remaining_bytes, self.K, self.K/8, len(block), padding)
			self.last_block = True
			if padding*8 == self.K:
				return None
			else:
				block.extend('0'*padding*8)
		print "length of block is %d" % len(block)
		print "self.K = %d" % self.K
		assert len(block) == self.K
		return block
	
	def generate_constraint_matrix(self, c, d):
		# c is the number of constraint symbols.
		# d is the density (remember LDPC is LOW density)
		# construct a k x c matrix. because c is a constant and not from a
		# distribution, this is considered a "regular" LDPC code. 
		assert (0 < d < 1)
		# given the size k x c, calculate the number of bits to make into 1s
		# note the each of the c rows of the matrix G will form the (n-1)
		# coefficients of a constraint symbol, where n is the number of ones
		# present in the column, calculated as the product of the density term
		# d and k. 
		num_ones = int(c*self.K*d)
		# build G as a long vector at first so it's easier to address specific
		# bits. we'll reshape it below. 
		G = numpy.zeros(self.K*c, int)
		print "G matrix:"
		print G
		one_indices = random.sample(xrange(self.K*c), num_ones)
		G[one_indices] = 1
		G.shape = (self.K,c)
		# both the encoder and deocder need to know G
		self.G = G
		return G

	def get_constraint_matrix(self):
		return self.G

	def num_bits(self, block):
		# block is a bitarray object, which packs bits into longs. so number of
		# bits is the size in bytes * 8. 
		info = block.buffer_info()
		return info[1]*8

	def next_block(self):
		# keep track of where we are and return the next block
		if self.last_block:
			return None
		the_block = self._encode_binary_block()
		self.current_block += 1
		return the_block

class RaptorEncoder:
	def __init__(self, block, G=None, symb_size=1, debug=True):
		# precode and distribution are each function variables
		self.debug = debug
		# symbols is a bitarray()
		self.symbols = block
		# generator matrix for pre-code, if any
		self.G = G
		# precoded is a bitarray()
		self.intermediate = None

	def ldpc_precode(self):
		# NOTE: generate_constraint_matrix() must be called first. this method
		# does not call it so that the same generator matrix can be re-used
		# when desired (eg. to create replicable results)

		# convert symbols to a numpy array for matrix multiplication
		G_rows, G_cols = self.G.shape
		k = numpy.array(self.symbols.tolist())
		print "precoding with generator matrix..."
		print self.G
		print "source symbols"
		print self.symbols
		# now here is the key: we must calculate the c redundant symbols z_i
		# such that z_i xor G[:,i] = 0. 
		z = numpy.array([], int)
		for i in xrange(G_cols):
			# sum([a,b,c...]) % 2 is equivalent to a^b^c^... 
			coefficients = self.G[:,i]
			print "coefficients"
			print coefficients
			other_terms = coefficients*k
			print "other terms"
			print other_terms
			xor_other_terms = sum(other_terms) % 2
			print "xor value of other terms: %d" % xor_other_terms
			# we require xor_other_terms ^ zi = 0. this is true when zi has the
			# same value as xor_other_terms. 
			z = numpy.append(z, xor_other_terms)
		# store as bitarray
		self.z = bitarray(z.tolist())
		self.intermediate = self.symbols + self.z
		print "all intermediate symbols:"
		print self.intermediate
		print len(self.intermediate)
		assert len(self.intermediate) == G_cols+len(k)
		return self.intermediate
	
	def distribution_random_LT(self, num_symbols):
		# return a vector of coefficient indices sampled from the contained
		# distribution. example return value: [17,22,238]

		# sample a weight from uniform distribution
		d = numpy.random.random_integers(1, num_symbols)

		# construct a vector of d coefficient indices from k. (sampling without
		# replacement. 
		v = random.sample(range(num_symbols),d)
		return v

	def generate_encoded(self):
		if self.intermediate:
			symbols = self.intermediate
		else:
			symbols = self.symbols

		# the distribution fn must take as argument the number of symbols it is
		# operating over, and return a vector of coefficients indices, sampled
		# according to the distribution therein.  example output for 10
		# precoded symbols: [2,4,9]
		v = self.distribution_random_LT(len(symbols))
		v.sort()

		# grab the symbols at the index positions that have a 1 in the
		# coefficient vector.
		selected_symbols = [symbols[idx] for idx in v]

		# get the first two values to start the xor process (note this actually
		# removes the values from selected_symbols, so after the for loop below
		# selected_symbols will be in a wierd state and should not be used. 
		xorval = sum(selected_symbols) % 2
			
		# return the xor'ed value and the associated coefficients
		return {'val': xorval, 'coefficients': v}

class RaptorGaussDecoder:

	def __init__(self, K, debug=True):
		self.debug = debug
		self.K = K
		self.A = numpy.array([], dtype=bool)
		self.b = numpy.array([], dtype=bool)
		self.blocks_received = 0
		self.symbols_processed = 0

	def add_block(self, encoded):
		# increment number of blocks received either way
		self.blocks_received += 1

		val = bool(encoded['val'])
		coeff = encoded['coefficients']
		
		# create a new row vector and set it to one as indicated by the
		# coefficient vector
		new_row = bitarray('0'*self.K)
		for i in range(self.K):
			if i in coeff:
				new_row[i] = 1
			if i > max(coeff):
				break
		print "encoded block: " + str(new_row)

		# compare it to the ones we've already received:
		duplicate = False
		for row in self.A:
			if numpy.all(row == new_row):
				duplicate = True
				break

		if not duplicate:
			# add the new row to the bottom
			if not len(self.A):
				self.A = numpy.array(new_row.tolist())
			else:
				self.A = numpy.vstack((self.A, numpy.array(new_row.tolist())))
			self.b = numpy.append(self.b, val)

	def is_full_rank(self):
		if self.A.size == 0:
			return False

		b = numpy.array([self.b])
		both = numpy.hstack((self.A,b.T))
		rank_ab = numpy.linalg.matrix_rank(both)
		rank_a = numpy.linalg.matrix_rank(self.A)

		'''
		# print out list versions to copy over to matlab for sanity checking.
		A = self.A.copy()
		A = A.tolist()
		intA = []
		for row in A:
			intArow = [int(a) for a in row]
			intA.append(intArow)
		print intA
		intb = [int(bb) for bb in self.b.tolist()]
		print intb
		'''
		if rank_ab == rank_a and rank_ab == self.K:
			print "rank_ab = rank_a = %d" % rank_ab
			print "matrix is full rank."
			return True
		print "rank_ab = %d, != rank_a = %d" % (rank_ab, rank_a)
		return False

	def num_blocks(self):
		# how many encoded blocks have we received so far?
		# this is equivalent to the numer of rows in A. shape() returns (rows, cols)
		return self.A.shape[0]

	def remove_null_rows(self, mat):
		# empty rows
		rows, cols = mat.shape
		all_false = numpy.array([], int)
		for r in xrange(rows):
			if not mat[r,:].any():
				all_false = numpy.append(all_false, r)
				
		to_keep = [r for r in xrange(rows) if r not in all_false]
		return mat[to_keep,:]

	def remove_duplicate_rows(self, mat):
		# duplicates
		duplicates = numpy.array([], int)
		rows, cols = mat.shape
		for r in xrange(rows):
			this_row = mat[r,:]
			for rr in xrange(r+1,rows):
				if rr in duplicates:
					continue
				test_row = mat[rr,:]
				diff = test_row ^ this_row
				if not numpy.any(diff):
					duplicates = numpy.append(duplicates, rr)
		to_keep = [r for r in xrange(rows) if r not in duplicates]

		return mat[to_keep,:]


	def decode_gauss_base2(self):
		# use tmp matrices in case our solution fails. 
		b = numpy.array([self.b])
		mat = numpy.hstack((self.A,b.T))
		tri, b = self._triangularize(mat)
		return self._backsub(tri, b)

	def _backsub(self, tri, b):
		mat = numpy.hstack((tri, numpy.array([b]).T))
		rows, cols = tri.shape
		soln = numpy.zeros(cols, int)
		# initialize solution vector with RHS of last row
		for i in (xrange(cols)).__reversed__():
			# the second term is always a 1 since the row in question always
			# has a 1 at the diagonal
			inner_term = (numpy.dot(mat[i,i:cols],soln[i:cols]) + mat[i,-1]) % 2
			soln[i] = numpy.logical_and(inner_term, 1)
			#print "soln for x" + str(i)+" = " + str(soln[i])
		
		return soln
		# verify solution:
		#print mat
		#print "computed solution:"
		#print soln

	def _triangularize(self, mat):

		#mat = self.remove_null_rows(mat)
		#mat = self.remove_duplicate_rows(mat)
		rows, cols = mat.shape
		# we tacked the solution vector onto the end of the matrix, so don't
		# count it in terms of the number of columns to iterate over. 
		cols = cols -1
		
		# first, we want to pivot the rows to put A in upper triangular form
		# (get 0's into all columns positions below the given row)
		for c in xrange(0, cols):
			# examine the row values below the diagonal (note that c starts at
			# 1, not 0. so in the first column, go from 1 to rows, in the
			# second column, go from 2.. rows, etc)
			col_vals = mat[c:rows,c]
			if col_vals.max() == 0:
				print "error: all zeros below row/column (%d, %d). multiple solutions." % (c,c)
				print mat[c:rows, c:rows]
				return None
			# find first row with a 1 in the left-most column (non-zero returns
			# a tuple, and we want the 0'th element of the first dimension of
			# the tuple since this is just a row vector)
			max_i = col_vals.nonzero()[0][0]

			# do the 'partial pivot': swap rows max_r and c in A and b (unless
			# the current row already has a one, then we're not going to get
			# any better). 
			if not (max_i+c) == c:
				upper_row = mat[c,:].copy()
				lower_row = mat[c+max_i,:].copy()
				mat[c,:] = lower_row
				mat[c+max_i,:] = upper_row

			# now zero out the 1's remaining in this column below the diagonal.
			# get the c'th row (yes, c is also the column value - this ensures
			# we start at the row below the diagonal)
			cth_row = mat[c,:]

			# now for all rows below this one, xor the c'th row with those that
			# contain a 1 in this column, in order to make it 0. (make sure to
			# do this with the right hand solution vector, too). 
			for r in xrange(c+1,rows):
				if mat[r,c] == 1:
					mat[r,:] = (cth_row ^ mat[r,:])
		# end column iteration

		# now we can get rid of the dangling rows since our solution is
		# uniquely specified by the top square component. 
		mat = mat[0:cols,:]

		return mat[:, 0:cols], mat[:, -1]


	def decode_gauss_base10(self):
		# attempt decode
		print "attempting solution..."
		if not self.is_full_rank():
			print "A is not full rank, sorry."
			return None
		soln, residues, rank, sing = numpy.linalg.lstsq(self.A, self.b)
		self.decoded_values = soln 
		return soln, residues, rank, sing

	def convert(self):
		# convert the values back to strings
		bits = bitarray(self.decoded_values.tolist())
		return bits.tostring()
		

class RaptorBPDecoder:
	
	def __init__(self, K, G=None):
		# actual data symbols per block
		self.K = K
		# constraint matrix
		self.G = G
		# BP decoder bookkeeping
		self.symbols_processed = 0
		# keep track of how many xor operations have been performed (does not
		# account for looping over lists since these could be optimized out in
		# a more legit implementation)
		self.symbol_operations = 0
		# the number of columns of G is the number of constraint symbols, which
		self.known_symbols = {}
		self.waiting_symbols = []
		if self.G != None:
			self.constraint_symbols = self.G.shape[1] - self.K
			self.prime()
		else:
			self.constraint_symbols = 0

	def prime(self):
		# "prime" the decoding pump by filling in the info we already know. 
		# if pre-coding was used, then the symbols decoded by the BP decoder
		# include redundant 'constraint' symbols that give information. 

		# we know that for the constraint symbols, there are redundant blocks
		# such that xor of coeffs*(x1,x2,...xn) = 0. we have those coeffs
		# already, they are the values of the corresponding columns in G

		constraint_symbols = self.G.shape[1]
		print "there are %d constraint symbols" % constraint_symbols
		for i in xrange(constraint_symbols):
			coeffs = self.G[:,i].nonzero()[0]
			coeffs = numpy.append(coeffs, self.K+i)
			zi = {'coeffs': coeffs.tolist(), 
					'xor_val': 0,
					}
			print zi
			self.waiting_symbols.append(zi)
		print "constraint symbols are"
		print self.waiting_symbols


	def bp_decode(self, block):
		self.symbols_processed += 1
		val = bool(block['val'])
		coeffs = block['coefficients']
		resolved = []
		
		# add the symbol either to the known list if it's of length one, or to
		# the waiting list otherwise. then process the known list against the
		# waiting list. 

		if len(coeffs) == 1:
			print "degree one!"
			print block
			if not coeffs[0] in self.known_symbols.keys():
				self.known_symbols[coeffs[0]] = val
		else:
			print "no such luck. processing..."
			self.waiting_symbols.append({'coeffs': coeffs, 'xor_val': val})
		
		# XXX this is VERY inefficient, but simple. we're not dealing with huge
		# files and not measuring decoding speed, but are measuring symbol
		# overhead, so it's more important to get this right than fast. 
		keep_processing = True
		while keep_processing:
			# initialize to false each time through the loop. set to true only
			# if there's a new symbol released. 
			keep_processing = False
			# always iterate over the known symbols
			new_known = []
			for known_c, known_v in self.known_symbols.iteritems():
				resolved = []
				for item in self.waiting_symbols: 
					if known_c in item['coeffs']:
						item['xor_val'] = item['xor_val'] ^ known_v
						self.symbol_operations += 1
						item['coeffs'].remove(known_c)
					if len(item['coeffs']) == 1:
						keep_processing = True
						resolved.append(item)
				self.waiting_symbols = [item for item in self.waiting_symbols if not item in resolved]
				new_known.extend(resolved)
			for item in new_known:
				if item['coeffs'][0] not in self.known_symbols.keys():
					self.known_symbols[item['coeffs'][0]] = item['xor_val']

		print "state of the BP:"
		print "known symbols"
		print self.known_symbols
		print "still waiting symbols"
		print self.waiting_symbols

		# need the known symbols to be the original k, not (just) the
		# constraint symbols.
		recovered_originals = [k for k in self.known_symbols if k in range(self.K)]
		if len(recovered_originals) == self.K:
			print self.known_symbols
			# return symbols as a bitarray
			return bitarray(self.known_symbols.values()[0:self.K])
		else: return None

def run_gauss(filename):
	DEBUG = True

	# if we want everything to go in one block, then use len(data) as the block
	# length
	K = 8
	epsilon = int(0.5*8)
	manager = RaptorManager(filename, K)
	block = manager.next_block()
	decoded_blocks = []
	processed_blocks = 0
	output_blocks = 0
	while block:
		output_blocks += 1
		print "-----------------------------"
		print "block %d" % output_blocks
		if DEBUG: 
			print "encoding original (source) block: " + str(block)
		# this encoder is non-systematic and uses no pre-code. 
		encoder = RaptorEncoder(block)
		decoder = RaptorGaussDecoder(K)

		# grab new symbols and periodically check to see if we've gathered enough
		# to find a solution. when the decoder matrix is full rank, try and solve
		# for the original symbols. 
		while not decoder.is_full_rank():
			for i in xrange(K+epsilon):
				e = encoder.generate_encoded()
				decoder.add_block(e)

		print "attempting to solve after " + str(decoder.blocks_received) + " blocks."
		original_block = decoder.decode_gauss_base2()
		print "decoded block was..."
		print original_block

		decoded_blocks.append(original_block)
		processed_blocks += decoder.symbols_processed
		block = manager.next_block()

	print "decoder processed %d output blocks after %d received blocks. A total of factor of %d average overhead" % (output_blocks, processed_blocks, processed_blocks/float(output_blocks))
	print "decoded blocks:"
	print decoded_blocks
	print "decoded message"
	for d in decoded_blocks:
		sys.stdout.write(d.tostring())
	print ""

def run_bp(filename, precode, K, c, d):
	DEBUG = True

	manager = RaptorManager(filename, K)
	if precode:
		G = manager.generate_constraint_matrix(c,d)
	else:
		G = None

	decoded_blocks = []
	processed_blocks = 0
	source_blocks = 0
	symops = 0
	block = manager.next_block()
	while block:
		source_blocks += 1
		print "next block... "
		print block
		# this encoder is non-systematic
		encoder = RaptorEncoder(block, G)
		decoder = RaptorBPDecoder(K, G)
		if precode:
			# precoding happens only once per block
			intermediate = encoder.ldpc_precode()

		original_symbols = False
		while not original_symbols:
			e = encoder.generate_encoded()
			original_symbols = decoder.bp_decode(e)
		print block
		print "%d blocks processed for this block of %d source symbols." % (decoder.symbols_processed, K)
		decoded_blocks.append(original_symbols)
		symops += decoder.symbol_operations
		processed_blocks += decoder.symbols_processed
		block = manager.next_block()
		print "for this round, 1 source block was processed after receiving %d encoded blocks." % decoder.symbols_processed
	
	overhead = (processed_blocks/float(source_blocks))
	print "%d output blocks" % source_blocks
	print "%d received blocks." % processed_blocks
	print "overhead: %d" % overhead

	print "decoded blocks:"
	print decoded_blocks
	print "decoded message"
	for d in decoded_blocks:
		sys.stdout.write(d.tostring())
	print ""
	return {'K': K, 'precode':precode, 'c':c, 'd':d, 'source':source_blocks, 'processed': processed_blocks, 'overhead': overhead, 'symops': symops}

	
if __name__ == '__main__':
	if len(sys.argv) != 2:
		sys.stderr.write("Usage: ./raptor filename")
		sys.exit(1)
	filename = sys.argv[1]
	#run_gauss(filename)
	precode_results = []
	noprecode_results = []
	for K in xrange(8,40,8):
		noprecode_results.append(run_bp(filename, precode=False, K=K, c=None, d=None))
		precode_results.append(run_bp(filename, precode=True, K=K, c=3, d=0.2))

	print "Non precoded results"
	for r in noprecode_results:
		print "K: %d. Avg Overhead: %f. Symbol Operations: %d" % (r['K'], r['overhead'], r['symops'])
		print ""
	print "\n"
	
	print "Precoded results"
	for r in precode_results:
		print "K: %d. Avg Overhead: %f. Symbol Operations: %d" % (r['K'], r['overhead'], r['symops'])
		print ""





