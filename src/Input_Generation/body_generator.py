#mass, x, y, z(cartiensian coordinates), u, v, w (velocity vector)
import random
import sys
import getopt

def main():

    file_name = "nbody_input"
    
    mac_order = 28
    psm_max_t = 1

    swarm = False
    psm = False
    bhs = False
    fmm = False
    verbose = False

    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    # Parse command line arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], "abfhm:o:pst:v", ["all", "help", "output=", "swarm", "time=" "verbose"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)
        
    for o, a in opts:
        if o in ('-h', "--help"):
            usage()
            sys.exit()
        elif o in ('-a', "--all"):
            psm = True
            bhs = True
            fmm = True
        elif o == '-b':
            bhs = True
        elif o == '-f':
            fmm = True
        elif o == '-m':
            mac_order = int(a)
        elif o in ('-o', "--output"):
            file_name = a
        elif o == '-p':
            psm = True
        elif o in ('-s', "--swarm"):
            swarm = True
        elif o in ('-t', '--time'):
            psm_max_t = float(a)
        elif o in ('-v', "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    num_bodies = int(sys.argv[len(sys.argv) - 1])

    if swarm:
        if verbose:
            print("Generating {} bodies as a swarm".format(num_bodies))
        file_name += "_swarm"
    else:
        if verbose:
            print("Generating {} bodies".format(num_bodies))

    file_name += "_" + str(num_bodies)

    # Calculate the maximum and minimum values for body mass/velocity
    lim = [-10 * num_bodies**(1/3), 10 * num_bodies**(1/3)]

    bodies = []
    
    # Add a central body to the system if not in swarm mode
    if not swarm:
        bodies.append([1, 0, 0, 0, 0, 0, 0]) 
        num_bodies -= 1

    # Calculate pseudo-random body values
    for _ in range(num_bodies):
        bodies.append([random.uniform(*lim) for i in range(7)])

    if not swarm:
        num_bodies += 1

    # Output BHS input file
    if bhs:
        if verbose:
            print("\nWriting BHS Input File")
        output_bhs(bodies, file_name)

    # Output PSM input file
    if psm:
        if verbose:
            print("\nWriting PSM Input File (Mac Order: {}, Max Time: {}\n".format(mac_order, psm_max_t))
        output_psm(num_bodies, mac_order, psm_max_t, bodies, file_name)

    # Output FMM input file
    if fmm:
        if verbose:
            print("\nWriting BHS Input File")
        output_fmm(bodies, file_name)


def output_bhs(bodies, prefix):
    """
    Writes the specified bodies to an output file that can be parsed
    by the the Barnes-Hut Simulation program
    :param bodies: The bodies to output
    :param prefix: The prefix for the output file name
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
    by the the Fast-Multipole Method program
    :param num_bodies: The number of bodies being written
    :param bodies: The bodies to output
    :param prefix: The prefix for the output file name
    """

    f_name = prefix + "_fmm.in"

    f = open(f_name, "w")
    
    for b in bodies:
        for c in b[0:3]:
            f.write(str(c) + " ")
        f.write("0 0 0 1")
        f.write("\n")

    f.close()

def output_psm(num_bodies, mac_order, max_t, bodies, prefix):
    """
    Writes the specified bodies to an output file that can be parsed
    by the Parker-Sochaki Method Program
    :param num_bodies: The number of bodies being written
    :param bodies: The bodies to output
    :param prefix: The prefix for the output file name
    """

    f_name = prefix + "_psm.dat"
    f = open(f_name, "w")

    f.write(" {} 0\n".format(num_bodies))
    f.write(" {}\n".format(mac_order))
    f.write(" 0.E+0,  {}D0, -.025\n".format(max_t))
    f.write(" -1.0E-14, .F.\n")

    # Input format:
    # body    mass     x1     x2     x3     v1     v2     v3 
    for b in bodies:
        line = ""
        for i in range (7):
            if b[i] >= 0:
                line += " "
            line += " {0:0.6f}".format(b[i])

        f.write("{}\n".format(line))

    f.close()

def usage():
    print("Usage: {} [options] num_bodies".format(sys.argv[0]))
    print("Options:")
    print("\t-a --all: Generates a file for all supported file types. Equivalent to -bfp")
    print("\t-b: outputs a file for the supported Barnes-Hut Simulation Implementation")
    print("\t-b: outputs a file for the supported Fast-Multipole Method Implementation")
    print("\t-h --help: display the usage")
    print("\t-m: [mac_value]: set the maximum Maclaurin polynomial order (PSM Only)")
    print("\t-o --output: [filename]: Determine the file name prefix for the output file")
    print("\t-p: outputs a file for the supported Parker-Sochaki Method Implementation")
    print("\t-s: --swarm: generates a swarm based system of bodies")
    print("\t-t --time [max_time]: set the maximum time value for calculation (PSM Only)")
    print("\t-v --verbose: run with verbose output")


if __name__ == "__main__":
    main()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
