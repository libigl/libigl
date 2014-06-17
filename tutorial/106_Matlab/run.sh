# Make sure that matlab can be launched from the commandline
# if it is not the case use the following commands:
# cd /usr/local/bin
# sudo ln -s /Applications/MATLAB_R2012b.app/bin/matlab

# Tested on MAC only

DYLD_LIBRARY_PATH=/Applications/MATLAB_R2012b.app/bin/maci64/

./build/106_Matlab
