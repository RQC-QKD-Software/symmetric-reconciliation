# Simulation of symmetric blind reconciliation protocol for quantum key distribution systems

## 1. General description

Repository contains simulation of symmetric blind reconciliation algorithm realized in Python 2.7.

## 2. Files contents

### 2.1. test_error_correction.py
Launches a simulaiton of the symmetric blind reconciliation for the current number of tries and QBER values. Keeps the result in output.txt.

### 2.2. error_correction_lib.py
Contains basic procedures for performing a test of reconciliation protocol. It includes:
 - generation of random keys (bit strings);
 - adding errors in accordance withc given level of quantum bit error rate (QBER);
 - choosing an appropriate rate of LDPC code among a given range and numbers of shortened and punctured bits;
 - generation of positions for shortened and punctured bit;
 - extending key with shortened and punctured bits;
 - syndrome encoding;
 - syndrome decoding;
 - performing symmetric blind reconciliation for given pair keys;
 - testing of the full procedure of information reconciliation, including generation of keys, addings errors, and collection of statistics;

### 2.3. codes_1944.txt
Pool of four [standard LDPC codes](http://ieeexplore.ieee.org/document/5307322/?arnumber=5307322) of block length 1944 together with positions for [untanited puncturing](http://ieeexplore.ieee.org/document/6290312/?arnumber=6290312). The set of code rates is {5/6, 3/4, 2/3, 1/2}.

### 2.4. codes_4000.txt
Pool of nine LDPC codes with block length 4000, constructed with [improved progressive edge growing algorithm](http://ieeexplore.ieee.org/document/5606185/?arnumber=5606185) with [particular distribuition polynomials](http://ieeexplore.ieee.org/document/5205475/?arnumber=5205475). The set of code rates is {0.9, 0.85, ..., 0.5}.

### 2.5. file_utils.py
Contains some auxiliary procedures for reading files with codes.

## 3. Notes about storage of parity-check matrices
The storage of parity-check matrices is based on two variables: `s_y_joins` and `y_s_joins`. They contains positions of nonzero elements for each row and column correspondigly.
For example, for the matrix
```sh
H = 
1 1 0 1
1 0 1 1
0 1 1 0
```
one has
```sh
s_y_joins = [[0,1,3], [0,2,3],[1,2]]
y_s_joins = [[0,1],[0,2], [1,2],[0,1]]
```
