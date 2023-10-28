# The codes help generate ranked res_atom.dat for different interaction types based on the counting of each interaction type
The res_atom.dat can be used in RINRUS following the Example 1 step 9. The github link for RINRUS is attached below:
https://github.com/natedey/RINRUS

# The codes also create a heatmap for different interactions for visualization

# How to execute the codes?
For int_dat.py, do python3 int_dat.py -probe 2cht_h_modify.probe -s A:203 -slv HOH

For heatmap.py, do python3 heatmap.py -hm res_intfreq_dict.dat

# See the output example file for the results of the given probe file
