#include "input_parser.h"
#include "delaunay_helper.h"
#include "conley_complex.h"

#include <fstream>
#include <limits>
#include <time.h>


int main(int argc, char* argv[])
{
    InputParser input_parser(argc, argv);

    if(input_parser.cmdOptionExists("--help") || input_parser.cmdOptionExists("-h") || argc <= 1){
        std::clog << "Usage: main_conley object.off [-E] [-o output_file] [-h]" << std::endl
        << "Compute flow complex of the 3D object object.off using the delaunay "
        << "triangulation." << std::endl
        << "-o output_file   : write the result in output_file." << std::endl
        << "-h, --help       : display this message." << std::endl;
        exit(EXIT_SUCCESS);
    }

    const char* filename = argv[1];
    // open the mesh file and import into a Polyhedron
    Polyhedron poly;
    std::ifstream input(filename);
    if (!input)
    {
        std::cerr << std::string(filename) + " is not a valid input file." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!(input >> poly)) {
        std::cerr << "Cannot load model" << std::endl;
        exit(EXIT_FAILURE);
    }
    if(poly.empty()) {
        std::cerr << "Empty model" << std::endl;
        exit(EXIT_FAILURE);
    }
    if(!CGAL::is_triangle_mesh(poly)) {
        std::cerr << "Cannot load model" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string output_filename(filename);
    const std::string& output_option = input_parser.getCmdOption("-o");
    if (!output_option.empty()){
        output_filename = output_option;
    }
    else{
        output_filename.erase(output_filename.end()-4, output_filename.end()); // remove ".off" from filename
        output_filename = output_filename + ".flow";
    }
    
    /* Now perform the method */

    // time evaluation variables
    clock_t time_last = clock();
    double time_final = 0.0;


    // Write poly to off file
    std::ofstream out ( "tmp/poly.off", std::ios::out | std::ios::trunc);

    if ( ! out . good () ) {
        std::cerr << "write poly. Fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }
    out << poly ;

    out.close();

    // create flow complex structure
    ConleyComplex conley_cplx(poly, .1);

    // Extract the polyhedron from the Delaunay mesh (and check the inclusion)
//    std::map<Delaunay::Simplex, Polyhedron::Facet> poly_simplices;
    std::set<Delaunay::Simplex> poly_simplices;
    bool included(conley_cplx.extract_shape_from_delaunay(poly, poly_simplices));

    std::cout << "==== INCLUSION" << std::endl;
    std::cout << "poly included in Delaunay: " << included << std::endl;

    std::cout << "===== POSET" << std::endl;
    conley_cplx.print_poset(std::cout);
    std::cout << "===== CONLEY" << std::endl;
    conley_cplx.print_infos();

    conley_cplx.write_vtk("tmp/conley.vtk");

    return 0;
}

// TODO : merge must be bi-directional "edges of poset can be merged in any direction ...
