#include "input_parser.h"
#include "flow_complex.h"

#include <fstream>
#include <limits>
#include <time.h>


int main(int argc, char* argv[])
{
    InputParser input_parser(argc, argv);

    if(input_parser.cmdOptionExists("--help") || input_parser.cmdOptionExists("-h") || argc <= 1){
        std::clog << "Usage: main_flow object.off [-E] [-o output_file] [-h]" << std::endl
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
    if (!input || !(input >> poly) || poly.empty() || !CGAL::is_triangle_mesh(poly))
    {
        std::cerr << std::string(filename) + " is not a valid input file." << std::endl;
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

    // do stuff
    FlowBuilder FB(poly);
    // Delaunay::Finite_vertices_iterator vit = FB.m_dela.finite_vertices_begin();
    // simplices.insert(Simplex(vit));
    std::clog << std::endl << "0D" << std::endl;
    for (Delaunay::Finite_vertices_iterator vit = FB.m_dela.finite_vertices_begin(); 
    vit != FB.m_dela.finite_vertices_end(); vit++) {
        Delaunay::Simplex s = Delaunay::Simplex(vit);
        if (DelaunayHelper::get_critical_info(FB.m_dela, s).c == CriticalType::Critical){
            FlowCell fc = FB.flowcell_from_critical_cell(s);
            std::clog << FB.get_flowcell_id(fc) << "\t" << fc << std::endl;
        }
    }

    // Delaunay::Finite_edges_iterator eit = FB.m_dela.finite_edges_begin();
    // simplices.insert(Simplex(*eit));
    std::clog << std::endl << "1D" << std::endl;
    for (Delaunay::Finite_edges_iterator eit = FB.m_dela.finite_edges_begin(); 
    eit != FB.m_dela.finite_edges_end(); eit++) {
        Delaunay::Simplex s = Delaunay::Simplex(*eit);
        if (DelaunayHelper::get_critical_info(FB.m_dela, s).c == CriticalType::Critical){
            FlowCell fc = FB.flowcell_from_critical_cell(s);
            std::clog << FB.get_flowcell_id(fc) << "\t" << fc << std::endl;
        }
    }

    // Delaunay::Finite_facets_iterator fit = FB.m_dela.finite_facets_begin();
    // simplices.insert(Simplex(*fit));
    std::clog << std::endl << "2D" << std::endl;
    for (Delaunay::Finite_facets_iterator fit = FB.m_dela.finite_facets_begin(); 
    fit != FB.m_dela.finite_facets_end(); fit++) {
        Delaunay::Simplex s = Delaunay::Simplex(*fit);
        if (DelaunayHelper::get_critical_info(FB.m_dela, s).c == CriticalType::Critical){
            FlowCell fc = FB.flowcell_from_critical_cell(s);
            std::clog << FB.get_flowcell_id(fc) << "\t" << fc << std::endl;
        }
    }
    
    // Delaunay::Finite_cells_iterator cit = FB.m_dela.finite_cells_begin();
    // simplices.insert(Simplex(cit));
    std::clog << std::endl << "3D" << std::endl;
    for (Delaunay::Finite_cells_iterator cit = FB.m_dela.finite_cells_begin(); 
    cit != FB.m_dela.finite_cells_end(); cit++) {
        Delaunay::Simplex s = Delaunay::Simplex(cit);
        if (DelaunayHelper::get_critical_info(FB.m_dela, s).c == CriticalType::Critical){
            FlowCell fc = FB.flowcell_from_critical_cell(s);
            std::clog << FB.get_flowcell_id(fc) << "\t" << fc << std::endl;
        }
    }
    std::clog << "\n";
    FB.print_poset(std::clog);
    
    time_final = (double)(clock() - time_last)/CLOCKS_PER_SEC;
    std::clog.precision(3);
    std::clog  << "computed in " << std::fixed << time_final << "s."<< std::endl;

    // eventually saving
    
    return 0;
}
