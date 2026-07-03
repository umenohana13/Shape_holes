#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3_to_lcc.h>
#include "delaunay_helper.h"

#include <iostream>
#include <fstream>
#include <cassert>

// Code by Yliess Bellargui
// CGAL LCC example

/* // If you want to use exact constructions.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Linear_cell_complex<3,3,
  CGAL::Linear_cell_complex_traits<3, CGAL::Exact_predicates_exact_constructions_kernel> > LCC_3;
*/

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;
typedef LCC_3::Dart_descriptor Dart_descriptor;
typedef LCC_3::Point           Point;

// Function used to display the voronoi diagram.
inline void display_voronoi(LCC_3& alcc, Dart_descriptor adart)
{
  // We remove the infinite volume plus all the volumes adjacent to it.
  // Indeed, we cannot view these volumes since they do not have
  // a "correct geometry".
  std::stack<Dart_descriptor> toremove;
  LCC_3::size_type mark_toremove=alcc.get_new_mark();

  // adart belongs to the infinite volume.
  toremove.push(adart);
  CGAL::mark_cell<LCC_3,3>(alcc, adart, mark_toremove);

  // Now we get all the volumes adjacent to the infinite volume.
  for (LCC_3::Dart_of_cell_range<3>::iterator
         it=alcc.darts_of_cell<3>(adart).begin(),
         itend=alcc.darts_of_cell<3>(adart).end(); it!=itend; ++it)
  {
    if ( !alcc.is_marked(alcc.beta(it,3), mark_toremove) )
    {
      CGAL::mark_cell<LCC_3,3>(alcc, alcc.beta(it,3), mark_toremove);
      toremove.push(alcc.beta(it,3));
    }
  }

  while( !toremove.empty() )
  {
    alcc.remove_cell<3>(toremove.top());
    toremove.pop();
  }

  assert(alcc.is_without_boundary(1) && alcc.is_without_boundary(2));

  std::cout<<"Voronoi subdvision, only finite volumes:"<<std::endl<<"  ";
  alcc.display_characteristics(std::cout) << ", valid="
                                          << alcc.is_valid()
                                          << std::endl;
}

template<typename LCC, typename TR>
void transform_dart_to_their_dual(LCC& alcc, LCC& adual,
                                  std::map<typename TR::Cell_handle,
                                           typename LCC::Dart_descriptor>& assoc)
{
  typename LCC::Dart_range::iterator it1=alcc.darts().begin();
  typename LCC::Dart_range::iterator it2=adual.darts().begin();

  std::map<typename LCC::Dart_descriptor, typename LCC::Dart_descriptor> dual;

  for ( ; it1!=alcc.darts().end(); ++it1, ++it2 )
  {
    dual[it1]=it2;
  }

  for ( typename std::map<typename TR::Cell_handle, typename LCC::Dart_descriptor>
          ::iterator it=assoc.begin(), itend=assoc.end(); it!=itend; ++it)
  {
    assoc[it->first]=dual[it->second];
  }
}

template<typename LCC, typename TR>
void set_geometry_of_dual(LCC& alcc, TR& tr,
                          std::map<typename TR::Cell_handle,
                                   typename LCC::Dart_descriptor>& assoc)
{
  for ( typename std::map<typename TR::Cell_handle, typename LCC::Dart_descriptor>
          ::iterator it=assoc.begin(), itend=assoc.end(); it!=itend; ++it)
  {
    if ( !tr.is_infinite(it->first) )
      alcc.set_vertex_attribute
        (it->second,alcc.create_vertex_attribute(tr.dual(it->first)));
    else
      alcc.set_vertex_attribute(it->second,alcc.create_vertex_attribute());
  }
}

std::ostream& operator<<(std::ostream& out, const std::vector<size_t>& simplex) {
    out << "<";
    for (size_t i : simplex)
        out << i << " ";
    out << "> ";
    return out;
}

template <typename Delaunay>
std::vector<std::list<typename Delaunay::Simplex> > D_sub_faces(const typename Delaunay::Cell_handle ch) {
    std::list<typename Delaunay::Simplex> l;
    l[0].push_back(Delaunay::Simplex(ch->vertex(0)));
    l[0].push_back(Delaunay::Simplex(ch->vertex(1)));
    l[0].push_back(Delaunay::Simplex(ch->vertex(2)));
    l[0].push_back(Delaunay::Simplex(ch->vertex(3)));
    l[1].push_back(Delaunay::Simplex(Delaunay::Edge(ch, 2, 3)));
    l[1].push_back(Delaunay::Simplex(Delaunay::Edge(ch, 1, 3)));
    l[1].push_back(Delaunay::Simplex(Delaunay::Edge(ch, 0, 3)));
    l[1].push_back(Delaunay::Simplex(Delaunay::Edge(ch, 1, 2)));
    l[1].push_back(Delaunay::Simplex(Delaunay::Edge(ch, 0, 2)));
    l[1].push_back(Delaunay::Simplex(Delaunay::Edge(ch, 0, 1)));
    l[2].push_back(Delaunay::Simplex(Delaunay::Facet(ch, 1)));
    l[2].push_back(Delaunay::Simplex(Delaunay::Facet(ch, 0)));
    l[2].push_back(Delaunay::Simplex(Delaunay::Facet(ch, 2)));
    l[2].push_back(Delaunay::Simplex(Delaunay::Facet(ch, 3)));
    l[3].push_back(Delaunay::Simplex(ch));
    return l;
}

// adart belongs to the infinite volume
// Output finite Vornoi cells (any dimension) to VTK
//template<typename LCC, typename TR>

// Assoc maps Delaunay volumes to their dual

template <typename TR>
void finite_voronoi_write_VTK(LCC_3& alcc, Dart_descriptor adart, TR& m_dela, std::map<typename TR::Cell_handle,
                              typename LCC_3::Dart_descriptor>& assoc, std::string filename) {
    std::ofstream out ( filename, std::ios::out | std::ios::trunc);

    if ( ! out . good () ) {
        std::cerr << "finite_voronoi_write_VTK for Delaunay_3. Fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    std::vector<std::vector<std::vector<size_t> > > voronoi_finite_cells;
    voronoi_finite_cells.resize(alcc.ambient_dimension+1);

    // Mark all infinite 0,1,2-cells
    typename LCC_3::size_type mark_infinite=alcc.get_new_mark();
    // Mark all "finite" Voronoi cells
    typename LCC_3::size_type mark_finite=alcc.get_new_mark();

    // Visit Delaunay cells and mark duals of finite non boundary cells
    for ( typename std::map<typename TR::Cell_handle, typename LCC_3::Dart_descriptor>
         ::iterator it=assoc.begin(), itend=assoc.end(); it!=itend; ++it) {
        if( !m_dela.is_infinite(it->first) )
        {
            // Mark the dual Voronoi vertex as finite
            Dart_descriptor dart(it->second);
            alcc.mark(dart, mark_finite);
            // Visit the faces of the Delaunay tetrahedron
            std::vector<std::list<typename Delaunay::Simplex> > faces(TR::D_sub_faces(it));
            // Delaunay triangles
            for (typename TR::Simplex simplex : faces[2]) {
                bool infinite = false;
                // Check if there is an adjacent infinite tetra

                if (!m_dela.is_infinite(simplex))
            }
        }
    }

    // Create points indices
    std::unordered_map<typename LCC_3::Vertex_attribute_const_descriptor, std::size_t>
        index;
    size_t nbpts(0);
    for(auto itv=alcc.vertex_attributes().begin(),
         itvend=alcc.vertex_attributes().end(); itv!=itvend; ++itv)
    {
        index[itv]=nbpts++;
    }

    // Save non marked cells to the vectors of cells voronoi_finite_cells
    typename LCC_3::Dart_const_descriptor sd;
    // Vertices
    for(typename LCC_3::template One_dart_per_cell_range<0>::const_iterator
        itvert=alcc.template one_dart_per_cell<0>().begin(),
         itvertend=alcc.template one_dart_per_cell<0>().end();
        itvert!=itvertend; ++itvert)
    {
        if (!alcc.is_marked(itvert, mark_infinite))
        {
            std::vector<size_t> simplex;
            simplex.push_back(index[alcc.vertex_attribute(itvert)]);
            voronoi_finite_cells[0].push_back(simplex);
        }
    }
    std::cout << "Vertices" << std::endl;
    for (auto simplex : voronoi_finite_cells[0])
        std::cout << simplex << std::endl;

    // Edges
    for(typename LCC_3::template One_dart_per_cell_range<1>::const_iterator
        itvert=alcc.template one_dart_per_cell<1>().begin(),
         itvertend=alcc.template one_dart_per_cell<1>().end();
        itvert!=itvertend; ++itvert)
    {
        if (!alcc.is_marked(itvert, mark_infinite))
        {
            std::vector<size_t> simplex;
            simplex.push_back(index[alcc.vertex_attribute(itvert)]);
            simplex.push_back(index[alcc.vertex_attribute(alcc.beta(itvert,0))]);
            voronoi_finite_cells[1].push_back(simplex);
        }
    }
    std::cout << "Edges" << std::endl;
    for (auto simplex : voronoi_finite_cells[1])
        std::cout << simplex << std::endl;

    // Polygons
    for(typename LCC_3::template One_dart_per_cell_range<2>::const_iterator
        itpoly=alcc.template one_dart_per_cell<2>().begin(),
         itpolyend=alcc.template one_dart_per_cell<2>().end();
        itpoly!=itpolyend; ++itpoly)
    {
        if (!alcc.is_marked(itpoly, mark_infinite))
        {
            std::vector<size_t> simplex;
            for (LCC_3::Dart_of_cell_range<2>::const_iterator
                   it=alcc.darts_of_cell<2>(itpoly).begin(),
                 itend=alcc.darts_of_cell<2>(itpoly).end(); it!=itend; ++it) {
                simplex.push_back(index[alcc.vertex_attribute(it)]);
            }
            voronoi_finite_cells[2].push_back(simplex);
        }
    }
    std::cout << "Polygons" << std::endl;
    for (auto simplex : voronoi_finite_cells[2])
        std::cout << simplex << std::endl;


    out.close();
}
