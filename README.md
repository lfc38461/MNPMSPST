# MNPMSPST
The MNPMSPST instance contains a small-scale set and a large-scale set, whcih are generated based on a production case in automobile manufacturing. The small-scale set includes 9 classes cases with K = {2,3,5},J = {10,20,30}. The large-scale set also includes 9 classes cases with K = {10,15,20},J = {50,80,100}. Every case contains the job information on index, processing time, type, pattern, setup time, initial setup time and quantity.

Folder "code": Contains the source code for all models and methods in the paper, and provides corresponding annotations for necessary code blocks and functions.

Folder "data": It includes all the instances, with "1_K2-J10"-"9_K5-J30" being small-scale examples, "10_K10-J50"-"18_K20-J100" being large-scale examples, and "N2_K10-J50"-"N5_K10-J50" being examples for testing different patterns.

Folder "solution": Including specific solutions for all instances.

"Results.xlsx": Organized detailed experimental data for each example on each method, including running time, upper and lower bounds, gap, etc.