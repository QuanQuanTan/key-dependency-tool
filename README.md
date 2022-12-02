# key-dependency-tool

This tool can help in verifying if a differential characteristic has any key/round dependencies based on the values that allow the differences to pass.
Currently, we only support SKINNY and GIFT block ciphers. <br>
The tools can be found under their respective folders.<br>
Note that CryptominiSAT has to be installed prior to running the program (for GIFT).<br>
To compile the files, just go into the respective src folder and enter the command <br>
"make SIZE=4" or "make SIZE=8" for the 64-bit or 12-bit version for the respective ciphers <br>

After compiling it, go to the bin folder where the executables can be found. <br>
To run it, simply give the corresponding trail number as a command line argument. For now, the trails can be found in main.cpp. Most of the trails can be found in eprint where the year and serial number is encoded in the header of the repective trails. <br>

