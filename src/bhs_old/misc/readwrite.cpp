#include "readwrite.h"

#include <string>
#include <fstream>
#include <sstream>

#include <sys/stat.h>

bool file_exists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

std::pair<std::vector<Body>, int> read_bodies(const char * filename, MPI_Comm comm){
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<Body> bodies;
    /* Read bodies from file */
    if(rank != 0){
        MPI_Status status;
        int x;
        MPI_Recv(&x, 1, MPI_INT, rank - 1, 0,
                MPI_COMM_WORLD, &status);
    }

    std::ifstream infile;
    infile.open(filename);

    std::string line;
    
    /*
     * The following code was modified by Richard Bimmer to use the same input
     * format as the rest of our research
     */

    int n = 0;
    int current_body = 0;

    // get total number of n-bodies 
    if (std::getline(infile, line)) {
	 std::string::size_type sz;
	 n = std::stoi(line, &sz);
         int current_line = 0;
         printf("Reading %d n-bodies from file...\n", n);

         double x, y, z, vx, vy, vz, m;
         while (std::getline(infile, line) && current_body < n){
            std::istringstream iss(line);

	    double val = 0;
	    try {
               std::string::size_type sd;
	       val = std::stod(line, &sd);
	       current_line++;
	    } catch (const std::invalid_argument) {
	        continue;
	    }

            if (current_line == 1) m = val;
	    if (current_line == 2) x = val;
	    if (current_line == 3) y = val;
	    if (current_line == 4) z = val;
	    if (current_line == 5) vx = val;
	    if (current_line == 6) vy = val;
            if (current_line == 7) {
                vz = val;

		// add body
		if ((current_body % size) == rank) {
		    bodies.push_back(Body{{x, y, z}, {vx, vy, vz}, m, 1});
                }
 
		// reset variables
		current_line = 0;
                current_body++;
		//printf("Read in n-body %d.\n", current_body++);
            }
         }

         printf("Read %ld nbodies from file.\n", current_body);
    } else {
        printf("There was an error reading from the input file. Please check that it is formatted correctly.\n");
        // TODO: Quit?	
    }	    
    infile.close();

    if(rank != size - 1) {
        int x = 1;
        MPI_Send(&x, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }

    return std::make_pair(bodies, current_body);
}

void write_bodies(const char * filename, const std::vector<Body> & bodies, MPI_Comm comm, bool overwrite){
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    std::ofstream myfile;
    if(rank != 0){
        MPI_Status status;
        int x;
        MPI_Recv(&x, 1, MPI_INT, rank - 1, 0,
                MPI_COMM_WORLD, &status);
        myfile.open(filename, std::ios::app);
    }else{
        if(!file_exists(filename) or overwrite){
            myfile.open(filename);
        }else{
            myfile.open(filename, std::ios::app);
        }
    }

    for(const Body & b : bodies){
        myfile << b.pos[0] << " " << b.pos[1] << " " << b.pos[2] << std::endl;
    }

    if(rank != size - 1){
        myfile.close();
        int x = 1;
        MPI_Send(&x, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }else{
        myfile << std::endl;
        myfile.close();
    }
}

void write_to_file(const char * filename, const double x, bool overwrite){
    std::ofstream myfile;
    if(!file_exists(filename) or overwrite){
        myfile.open(filename);
    }else{
        myfile.open(filename, std::ios::app);
    }
    myfile << x << std::endl;;
    myfile.close();
}

void write_tree(const char * filename, const Tree & tree, bool fulltree, bool overwrite){
    std::ofstream myfile;
    if(!file_exists(filename) or overwrite){
        myfile.open(filename);
    }else{
        myfile.open(filename, std::ios::app);
    }
    myfile << tree.to_string(fulltree) << std::endl;
    myfile.close();
}

void write_summary(InputParser ip, int N, int P){

    std::ofstream file;
    file.open(ip.out_sum_file());

    file << "number of processes: " << P << std::endl;
    file << "number of bodies: " << N << std::endl;
    if(ip.read_bodies()){
        file << "read bodies from file: " << ip.in_file() << std::endl;
    }
    else{
        file << "read bodies from file: none" << std::endl;
    }
    file << "barnes hut approximaiton constant (theta): " << ip.bh_approx_constant() << std::endl;
    file << "gravitational constant (G): " << ip.grav_constant() << std::endl;
    file << "number of time steps: " << ip.n_steps() << std::endl;
    file << "time step: " << ip.time_step() << std::endl;
    file << "sampling interval: " << ip.sampling_interval() << std::endl;
    if(ip.write_positions()){
        file << "output file with positions: " << ip.out_file() << std::endl;
    }
    else{
        file << "output file with positions: none" << std::endl;
    }
    if(ip.write_tree()){
        file << "output file with tree: " << ip.out_tree_file() << std::endl;
    }
    else{
        file << "output file with tree: none" << std::endl;
    }
    if(ip.clock_run()){
        file << "output file with running times: " << ip.out_time_file() << std::endl;
    }
    else{
        file << "output file with running_times: none" << std::endl;
    }

}