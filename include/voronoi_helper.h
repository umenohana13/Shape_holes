#include <CGAL/Delaunay_triangulation_3.h>
#include "cgal_typedef.h"
#include "delaunay_helper.h"

#include <iostream>
#include <fstream>
#include <cassert>


//typedef Delaunay_index TR;

struct Voronoi_simplex {
    size_t dim;
    std::vector<size_t> vertices;
};

std::ostream& operator<<(std::ostream& out, const Voronoi_simplex& simplex) {
    out << "<" << simplex.vertices.at(0) ;
    for (size_t i = 1; i < simplex.vertices.size(); ++i)
        out << "," << simplex.vertices.at(i);
    out << " | dim " << simplex.dim << ">" ;
    return out;
}

// Test if a given edge/facet is "Voronoi finite" (ie. not adjacent to any infinite cell) and returns its adjacent cofaces

template <typename TR>
bool voronoi_finite(TR& m_dela, const typename TR::Edge& e, std::vector<typename TR::Cell_handle>& neigh) {
    typename TR::Cell_circulator cell_begin =  m_dela.incident_cells(e), cell_end = cell_begin;
    do {
        if (m_dela.is_infinite(cell_begin)) {
            neigh.clear();
            return false;
        }
        else
            neigh.push_back(cell_begin);
        ++cell_begin;
    } while (cell_begin != cell_end);
    return true;
}

template <typename TR>
bool voronoi_finite(TR& m_dela, const typename TR::Facet& f, std::vector<typename TR::Cell_handle>& neigh) {
    if (m_dela.is_infinite(f.first)) {
        neigh.clear();
        return false;
    }
    else {
        neigh.push_back(f.first);

        typename TR::Cell_handle ch_mirror(m_dela.mirror_facet(f).first);
        if (m_dela.is_infinite(ch_mirror)) {
            neigh.clear();
            return false;
        }
        else
            neigh.push_back(ch_mirror);
    }
    return true;
}

template <typename TR>
class Voronoi_medial {
    //protected:
public:
    std::vector<std::vector<Voronoi_simplex> > medial_voronoi;
    std::unordered_map<typename TR::Cell_handle, size_t> cell_handle_to_index;
    std::vector<std::vector<typename TR::Simplex> > voronoi_to_delaunay;
    std::vector<typename TR::Point> points;

    Voronoi_medial(TR& m_dela) {
        medial_voronoi.resize(4);
        voronoi_to_delaunay.resize(4);
        // Visit Delaunay cells and mark duals of finite non boundary cells

        // First: visit Delaunay cells = Voronoi vertices
        size_t cpt(0);
        for (Delaunay_index::Cell_handle ch : m_dela.finite_cell_handles()) {
            if( !m_dela.is_infinite(ch) )
            {
                // Store the corresponding vertex
                Voronoi_simplex v_simplex;
                v_simplex.dim = 0;
                v_simplex.vertices.push_back(cpt);
                points.push_back(m_dela.dual(ch));
                medial_voronoi[0].push_back(v_simplex);
                voronoi_to_delaunay[0].push_back(typename TR::Simplex(ch));
                cell_handle_to_index[ch] = cpt++;
            }
        }

        // Second: visit Delaunay facets = Voronoi edges
        cpt = 0;
        for (Delaunay_index::Finite_facets_iterator it = m_dela.finite_facets_begin(); it != m_dela.finite_facets_end(); ++it) {
            std::vector<typename TR::Cell_handle> neigh;
            if (voronoi_finite(m_dela, *it, neigh)) {
                // Create dual edge
                CGAL_precondition(neigh.size() == 2);
                Voronoi_simplex v_simplex;
                v_simplex.dim = 1;
                for (typename TR::Cell_handle cell : neigh) {
                    v_simplex.vertices.push_back(cell_handle_to_index[cell]);
                }
                std::sort(v_simplex.vertices.begin(), v_simplex.vertices.end());
                medial_voronoi[1].push_back(v_simplex);
            }
            voronoi_to_delaunay[1].push_back(typename TR::Simplex(*it));
        }

        // Third: visit Delaunay edges = Voronoi faces

        for (Delaunay_index::Finite_edges_iterator it = m_dela.finite_edges_begin(); it != m_dela.finite_edges_end(); ++it) {
            std::vector<typename TR::Cell_handle> neigh;
            if (voronoi_finite(m_dela, *it, neigh)) {
                // Create dual face
                Voronoi_simplex v_simplex;
                v_simplex.dim = 2;
                for (typename TR::Cell_handle cell : neigh) {
                    v_simplex.vertices.push_back(cell_handle_to_index[cell]);
                }
                medial_voronoi[2].push_back(v_simplex);
            }
            voronoi_to_delaunay[2].push_back(typename TR::Simplex(*it));
        }
    }

    template <typename _TR>
    friend std::ostream& operator<<(std::ostream& out, const Voronoi_medial& vm);

    template <typename _TR>
    friend void write_VTK(std::string filename, const Voronoi_medial& voronoi, std::vector<std::vector<double> >* flags);
};

template <typename TR>
std::ostream& operator<<(std::ostream& out, const Voronoi_medial<TR>& vm) {
    for (int i = 0; i<4; ++i) {
        for (Voronoi_simplex vs : vm.medial_voronoi.at(i))
            out << vs ;
        out << std::endl;
    }
    return out;
}

template <typename TR>
void write_VTK(std::string filename, const Voronoi_medial<TR>& voronoi, std::vector<std::vector<double> >* flags = NULL) {
    std::ofstream out ( filename, std::ios::out | std::ios::trunc);
    if ( ! out . good () ) {
        std::cerr << "write_VTK for Voronoi in R^3. Fatal Error:\n  " << filename << " not found.\n";
        throw std::runtime_error("File Parsing Error: File not found");
    }

    assert (voronoi.points.size() == voronoi.medial_voronoi[0].size());

    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "Shape of holes" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl << std::endl;


    // Write vertices
    assert(voronoi.points.size() == voronoi.medial_voronoi[0].size());
    out << "POINTS " <<  voronoi.points.size() << " double" << std::endl;
    for (typename TR::Point p : voronoi.points)
        out << p << std::endl;
    out << std::endl;

    // Cells
    std::vector<size_t> ids, dims;
    size_t cpt(0);
    const size_t n_verts(voronoi.medial_voronoi[0].size()), n_edges(voronoi.medial_voronoi[1].size()), n_faces(voronoi.medial_voronoi[2].size());
    size_t n_simplices(0), size_simplices(0);
    // Vertices
    n_simplices += n_verts;
    size_simplices += n_verts*2;
    // Edges
    n_simplices += n_edges;
    size_simplices += n_edges*3;
    // Faces
    n_simplices += n_faces;
    for (auto cell : voronoi.medial_voronoi[2])
        size_simplices += cell.vertices.size()+1;

    out << "CELLS " << n_simplices << " " << size_simplices << std::endl;
    // Export vertices
    // Vertices
    for (int i=0; i<n_verts; ++i) {
        out << "1 " << i << std::endl;
        dims.push_back(0);
        ids.push_back(i); // index of the dual Delaunay cell
    }

    // Export edges
    cpt = 0 ;
    for (auto edge : voronoi.medial_voronoi[1]) {
        out << 2;
        for (size_t i : edge.vertices)
            out << " " << i;
        out << std::endl;
        dims.push_back(1);
        ids.push_back(cpt++); // index of the dual Delaunay cell
    }

    // Export faces
    cpt = 0 ;
    for (auto facet : voronoi.medial_voronoi[2]) {
        out << facet.vertices.size();
        for (size_t i : facet.vertices)
            out << " " << i;
        out << std::endl;
        dims.push_back(2);
        ids.push_back(cpt++); // index of the dual Delaunay cell
    }

    // Cells types
    out << "CELL_TYPES " << n_simplices << std::endl;
    // Vertices
    for (int i = 0; i < n_verts; ++i)
        out << 1 << " ";
    out << std::endl;
    // Edges
    for (int i = 0; i < n_edges; ++i)
        out << 3 << " ";
    out << std::endl;
    // Facets
    for (int i = 0; i < n_faces; ++i)
        out << 7 << " ";
    out << std::endl;

    out << std::endl;

    // Cells indices
    out << "CELL_DATA " << n_simplices << std::endl;
    assert (ids.size() == n_simplices);
    out << "SCALARS CellId int 1" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (size_t id : ids)
        out << id << " ";
    out << std::endl;

    // Cells dimension
    assert (dims.size() == n_simplices);
    out << "SCALARS Dimension int 1" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (size_t q : dims)
        out << q << " ";
    out << std::endl;

    // Flags (if provided)
    if ((flags != NULL) && (flags->size()>=2)) {
        if (((*flags)[0].size() == n_verts) && ((*flags)[1].size() == n_edges) && ((*flags)[2].size() == n_faces)) {
            out << "SCALARS Flags double 1" << std::endl;
            out << "LOOKUP_TABLE default" << std::endl;
            for (double x : (*flags)[0])
                out << x << " ";
            out << std::endl;
            for (double x : (*flags)[1])
                out << x << " ";
            out << std::endl;
            for (double x : (*flags)[2])
                out << x << " ";
            out << std::endl;
        }
    }
    out.close();
    return out;
}

