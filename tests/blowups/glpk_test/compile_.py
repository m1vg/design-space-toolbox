
import os

os.system("clear")
os.system("rm glpk_test.o")


# Uncomment to run files
# os.system("gcc -o Run3Cycle.o 3Cycle.c /usr/local/lib/libdesignspace.dylib")
# os.system("./Run3Cycle.o")

os.system("gcc -o glpk_test.o glpk_test_Case4.c  /usr/local/lib/libglpk.dylib")
os.system("./glpk_test.o")
