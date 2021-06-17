
import os

os.system("clear")
os.system("rm run_blow_ups.o")


# Uncomment to run files
# os.system("gcc -o Run3Cycle.o 3Cycle.c /usr/local/lib/libdesignspace.dylib")
# os.system("./Run3Cycle.o")

os.system("gcc -o run_blow_ups.o designspacetest.c /Users/mavalder/Dropbox/Projects\ 2018/UC\ Davis/Research/DST3/design-space-toolbox_3/build/Release/libdesignspace.dylib")
os.system("./run_blow_ups.o")
