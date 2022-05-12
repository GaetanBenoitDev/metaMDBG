




import os, sys

def main(argv):
    for i in range(0, 10000):
        execute_command("/home/gats/workspace/tools/MdbgAssembler/build/bin/mdbgAsmMeta multik -o ../../../run/overlap_test_metaflye/5/ -t 8")


def execute_command(command):

    print(command)
    ret = os.system(command)

    if ret != 0:
        print("Error: ", command)
        exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])  
