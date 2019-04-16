#mass, x, y, z(cartiensian coordinates), u, v, w (velocity vector)
import random
import sys
import getopt

def main():

    file_name = "nbody_input"

    swarm = False
    psm = False
    bhs = False
    fmm = False

    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    try:
        opts, args = getopt.getopt(sys.argv[1:], "abfho:ps", ["all", "help", "output=", "swarm"])
    except getopt.GetoptError as err:
        print str(err)  
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("a", "--all"):
            psm = True
            bhs = True
            fmm = True
        elif o == "b":
            bhs = True
        elif o == "f":
            fmm = True
        elif o in ("o", "--output"):
            file_name = a
        elif o == "p":
            psm = True
        elif o in ("-s", "--swarm"):
            swarm = True
        else:
            assert False, "unhandled option"

    num_bodies = int(sys.argv[len(sys.argv) - 1])
    
    if swarm:
        print("Generating {} bodies as a swarm".format(num_bodies))
        file_name += "_swarm"
    else:
        print("Generating {} bodies".format(num_bodies))

    file_name += "_" + str(num_bodies)

    lim = [-10 * num_bodies**(1/3), 10 * num_bodies**(1/3)]

    bodies = []
    
    if not swarm:
        bodies.append([1, 0, 0, 0, 0, 0, 0]) 
        num_bodies -= 1

    for i in range(num_bodies):
        bodies.append([random.uniform(*lim) for i in range(7)])


    if not swarm:
        num_bodies += 1

    if bhs:
        print("\nWriting BHS Input File")
        output_bhs(bodies, file_name)

    if psm:
        print("\nWriting PSM Input File\n")
        output_psm(num_bodies, bodies, file_name)

def output_bhs(bodies, prefix):
    """
    Writes the specified bodies to an output file that can be parsed
    by the the Barnes-Hut Simulation program
    :param: bodies: The bodies to output
    :param: prefix: The prefix for the output file name
    """

    f_name = prefix + "_bhs.in"

    f = open(f_name, "w")
    
    for b in bodies:
        for c in b[0:3]:
            f.write(str(c) + " ")
        f.write("0 0 0 1")
        f.write("\n")

    f.close()
    

def output_fmm(bodies, prefix):
    """
    Writes the specified bodies to an output file that can be parsed
    by the the Fast-Multipole Methon program
    :param: num_bodies: The number of bodies being written
    :param: bodies: The bodies to output
    :param: prefix: The prefix for the output file name
    """
    pass
    

def output_psm(num_bodies, bodies, prefix):
    """
    Writes the specified bodies to an output file that can be parsed
    by the Parker-Sochaki Method Program
    :param: num_bodies: The number of bodies being written
    :param: bodies: The bodies to output
    :param: prefix: The prefix for the output file name
    """

    f_name = prefix + "_psm.in"
    f = open(f_name, "w")

    f.write(str(num_bodies))
    f.write("\n")

    for b in bodies:
        for e in b:
            f.write(str(e))
            f.write("\n")

    f.close()

def usage():
    print("Usage: {} [options] num_bodies".format(sys.argv[0]))
    print("Options:")
    print("\t-a --all: Generates a file for all supported file types. Equivalent to -bfp")
    print("\t-b: outputs a file for the supported Barnes-Hut Simulation Implementation")
    print("\t-b: outputs a file for the supported Fast-Multipole Method Implementation")
    print("\t-h --help: display the usage")
    print("\t-o --output: [filename]: Determine the file name prefix for the output file")
    print("\t-p: outputs a file for the supported Parker-Sochaki Method Implementation")
    print("\t-s --swarm: generates a swarm based system of bodies")


if __name__ == "__main__":
    main()
