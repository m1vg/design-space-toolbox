
import os

os.system("clear")
os.system("rm Run3Cycle.o")


# Uncomment to run files
# os.system("gcc -o Run3Cycle.o 3Cycle.c /usr/local/lib/libdesignspace.dylib")
# os.system("./Run3Cycle.o")

os.system("gcc -o Run3Cycle.o serializationtest.c /Users/mavalder/Dropbox/Projects_2018/UC_Davis/Research/DST3/design-space-toolbox_3/build/Release/libdesignspace.dylib")
os.system("./Run3Cycle.o")
