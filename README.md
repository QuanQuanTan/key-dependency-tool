# key-dependency-tool

This tool can help in verifying if a differential characteristic has any key/round dependencies based on the values that allow the differences to pass. Currently, we only support AES, DEOXYS-TBC, SKINNY and GIFT block ciphers.
The tools can be found under their respective folders.
Note that CryptominiSAT has to be installed prior to running the program.
To compile the files, just simply run 
"./run_keyDependent <SIZE>" for SKINNY and GIFT (4 for SKINNY-64, GIFT-64 and 8 for SKINNY-128, GIFT-128) and,
"./run_keyDependent" for AES and DEOXYS-TBC
