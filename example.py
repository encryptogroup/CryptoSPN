"""
@author: Amos Treiber
"""
import logging
import sys
import numpy as np
from spn.algorithms.Inference import log_likelihood

from spn.structure.Base import Context
from spn.structure.leaves.parametric.Parametric import Bernoulli
from spn.algorithms.LearningWrappers import learn_parametric
from CryptoSPN import spn_to_aby_exec, spn_to_aby_file, Selection

fmt = logging.Formatter(fmt='{asctime} {levelname:8.8} {process} --- [{threadName:12.12}] {name:32.32}: {message}',
                        style='{')

console = logging.StreamHandler(sys.stdout)
console.setFormatter(fmt)

file = logging.FileHandler('example.log', 'w', 'utf-8')
file.setFormatter(fmt)

logging.basicConfig(handlers=[console, file], level=logging.INFO)

# Create example SPN.
np.random.seed(123)

train_data = np.random.binomial(1, [0.1, 0.8, 0.9, 0.1], size=(100, 4))

ds_context = Context(parametric_types=[Bernoulli, Bernoulli, Bernoulli, Bernoulli]).add_domains(train_data)

spn = learn_parametric(train_data, ds_context, min_instances_slice=20)

# Create file for ABY inputs
np.savetxt("all_data.txt", train_data, delimiter="", fmt='%0.18e;')

# Create ABY .cpp of SPN instance. Selection.OSELECTION specifies an oblivious selection network for RV selection.
# If none is needed, use Selection.NONE.
spn_to_aby_file(spn, bitlen=32, filename="example.cpp", sel=Selection.OSELECTION)

# Automatically create ABY Executable of SPN instance. This does not require a previous call to spn_to_aby_file.
# Since a script for compiling your executable runs in the background, you need to give it executable permission:
# chmod +x CryptoSPN/compiling/compile.sh
# Don't forget to build ABY and to add your CryptoSPN and ABY directories to CryptoSPN/Constants.py. Alternatively, you
# can pass them to spn_to_aby_exec via the aby_path and cryptospn_path parameters of spn_to_aby_exec.
spn_to_aby_exec(spn, bitlen=64, name="example", sel=Selection.OSELECTION, aby_path="/path/to/ABY/",
                cryptospn_path="/path/to/CryptoSPN/CryptoSPN/")

# Now, you can run the executable as client (-r 0) and server (-r 1) on different terminals
# with -b 32 or -b 64 bit precision on the first -i lines of client RVs supplied like "V0;V1;V2;V3;V4" in file -f:
# ./example -r 0 -a 127.0.0.1 -b 64 -i 10 -f "all_data.txt"
# ./example -r 1 -a 127.0.0.1 -b 64 -i 10 -f "all_data.txt"
# The server's IP address must be specified via -a
# When running the exec, the floating point circuit files located in CryptoSPN/aby_files/circ need to be available!

# The CryptoSPN executable will store the output in "{filename}_output_data{role}.txt".
# You can compare that output to the plaintext one:
train_ll = log_likelihood(spn, train_data)
np.savetxt("example_out_data.txt", train_ll, fmt='%0.18e')
