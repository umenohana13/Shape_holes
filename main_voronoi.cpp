#include "input_parser.h"
#include "filtration_medial.h"
#include "persistence.h"
#include "voronoi_helper.h"

#include <fstream>
#include <limits>
#include <time.h>

// TODO : compute and export (finite) Voronoi


int main(int argc, char* argv[])
{
    const char* filename = argv[1];
    // open the mesh file and import into a Polyhedron
    Polyhedron m_poly;
    std::ifstream input(filename);
    if (!input || !(input >> m_poly) || m_poly.empty() || !CGAL::is_triangle_mesh(m_poly))
    {
        std::cerr << std::string(filename) + " is not a valid input file." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << m_poly << std::endl;

//    LCC lcc_poly;
//    CGAL::polyhedron_3_to_lcc(lcc_poly, m_poly);
//    CGAL::IO::write_VTK("poly.vtk", lcc_poly);

    /* Now perform the method */

    // time evaluation variables
    clock_t time_last = clock();
    double time_final = 0.0;

//    // Create Delaunay triangulation
    Delaunay_index m_dela;
    for (Polyhedron::Vertex_iterator it = m_poly.vertices_begin(); it != m_poly.vertices_end(); ++it)
    {
        m_dela.insert(it->point());
    }

    std::cout << std::endl;

    std::cout << "Delaunay output" << std::endl;
    std::cout << m_dela;

    int cpt(0);
    for(Delaunay_index::Vertex_handle vh : m_dela.finite_vertex_handles()) {
        vh->info().second = cpt++;
    }

    std::cout << std::endl;


    std::cout << "Delaunay to VTK" << std::endl;
    write_VTK(m_dela, "delaunay.vtk", EDGES|CELLS);

    std::cout << "Voronoi" << std::endl;
    Voronoi_medial voronoi_medial(m_dela);

    std::cout << voronoi_medial;

    write_VTK("voronoi.vtk", voronoi_medial);

//    // Output the Delaunay mesh to vtk
//    LCC lcc;
//    CGAL::triangulation_3_to_lcc(lcc, m_dela);
//
////    CGAL::draw(lcc);
//
//    CGAL::IO::write_VTK("delaunay.vtk", lcc);


//    Dart_descriptor dd_copy(dd);
//    display_voronoi(dual_lcc, dd_copy);
//    CGAL::IO::write_VTK("voronoi.vtk", dual_lcc);

    return 0;
}
