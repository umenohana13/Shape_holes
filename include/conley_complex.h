#ifndef CONLEY_COMPLEX_H
#define CONLEY_COMPLEX_H


#include "cgal_typedef.h"
// #include "filtration.h"
#include "delaunay_helper.h"
#include "flow_complex.h"
#include "OSM/OSM.h"

#include <fstream>
#include <limits>

class ConleyComplex : public FlowComplex {
protected:
    double threshold;

public:
    // Value associated to a Poset_edge
    typedef double Poset_edge_value;
//    typedef int Poset_edge_value;
    // For OSM sparse matrices
    typedef Poset_edge_value Coefficient_type;
    typedef CGAL::OSM::Sparse_matrix<CGAL::OSM::Sparse_chain> Sparse_matrix_struct;
    typedef Sparse_matrix_struct::template Sparse_matrix_type<Coefficient_type, CGAL::OSM::COLUMN> Column_matrix;
    typedef Sparse_matrix_struct::template Sparse_matrix_type<Coefficient_type, CGAL::OSM::ROW> Row_matrix;
    typedef Sparse_matrix_struct::template Sparse_chain_type<Coefficient_type, CGAL::OSM::COLUMN> Column_chain;
    typedef Sparse_matrix_struct::template Sparse_chain_type<Coefficient_type, CGAL::OSM::ROW> Row_chain;

    // Edge structure
    struct Poset_edge {
        size_t first, second;
        Poset_edge_value value;
        bool operator< (const Poset_edge& other) const { return value < other.value; };
    };

    // IIS structure
    struct IIS {
        std::set<size_t> cells_ids; // Set of cells merged into the IIS
        size_t designee_id; // Id of the "representative" critical cell (lowest id among cells of the IIS)
    };

protected:
    // Store the poset (weighted by distances between cells) in a sparse OSM matrix
    Row_matrix poset;
    std::map<std::pair<size_t, size_t>,Poset_edge_value> values;
    std::map<size_t, IIS> iiss;

public:
    // Returns the difference between critical points

    Poset_edge_value value_diff_critical_points (size_t cell_id1, size_t cell_id2) {
        // Distance between critical points
        Vector diff(this->flowcell_from_id(cell_id2).get_p() - this->flowcell_from_id(cell_id1).get_p());
        Poset_edge_value res = sqrt(diff.squared_length());
        return res;
    }

    // Returns the difference between distance to border of critical cells

    Poset_edge_value value_diff_df (size_t cell_id1, size_t cell_id2) {
        if (value_diff_critical_points(cell_id1, cell_id2) < threshold) {
            double diff(abs(DelaunayHelper::get_critical_info(m_dela,this->flowcell_from_id(cell_id2).get_critical_simplex()).r) - abs(DelaunayHelper::get_critical_info(m_dela,this->flowcell_from_id(cell_id1).get_critical_simplex()).r));
            Poset_edge_value res = abs(diff);
            return res;
        }
    }

    // Determines if a given edge has a valid value for merge

    bool has_valid_value_for_merge (size_t cell_id1, size_t cell_id2) {
        return  (value_diff_critical_points(cell_id1, cell_id2) < threshold);
    }


    // Sort from largest to lower sizes edges such that the distance between critical points is above threshold
    // -> returns the -(sum of sizes) for edges such that the distance between critical points is lower than threshold

//    Poset_edge_value value (size_t cell_id1, size_t cell_id2) {
//        // Distance between critical points
//        Vector diff(this->flowcell_from_id(cell_id2).get_p() - this->flowcell_from_id(cell_id1).get_p());
//        if (sqrt(diff.squared_length()) < threshold)
//            return -(this->flowcell_from_id(cell_id2).get_simplices().size()+this->flowcell_from_id(cell_id1).get_simplices().size());
//        else
//            return 0;
//    }

//    bool extract_shape_from_delaunay(Polyhedron& poly, std::map<Delaunay::Simplex, Polyhedron::Facet>& simplices) {
    bool extract_shape_from_delaunay(Polyhedron& poly, std::set<Delaunay::Simplex>& simplices) {
        // Visit triangles of the polyhedron and check if they belong to the Delaunay mesh
        // Insert the corresponding simplex

        bool res = true;

        std::map<Polyhedron::Vertex_handle,size_t> poly_vertices_indices;
        std::map<std::vector<size_t>,Delaunay::Simplex> delaunay_triangles_by_indices;

        // Build polyhedron vertices indices
        size_t cpt(0);
        for (Polyhedron::Vertex_iterator it = m_poly.vertices_begin(); it != m_poly.vertices_end(); ++it)
        {
            Polyhedron::Vertex_handle v(it);
            poly_vertices_indices[v] = cpt++;
            //        std::cout << it->point() << std::endl;
        }

        // Build Delaunay triangles indices enumerations (sorted vectors)
        for (Delaunay::Finite_facets_iterator fit = m_dela.finite_facets_begin(); fit != m_dela.finite_facets_end(); ++fit) {
            Delaunay::Facet f(*fit);
            std::vector<size_t> f_indices(DelaunayHelper::simplex_to_indices(m_dela, f));
            delaunay_triangles_by_indices[f_indices] = Delaunay::Simplex(f);
        }

        // Build Polyhedron triangles indices and get their
        for (Polyhedron::Facet_iterator fit = poly.faces_begin(); fit != poly.facets_end(); ++fit) {
            // Visit facet vertices (3)
            if (!fit->is_triangle())
                throw std::runtime_error("Non triangle facet in poly");
            else {
                Polyhedron::Halfedge_around_facet_circulator heit = fit->facet_begin();
                std::vector<size_t> f_indices;

                do {
                    Polyhedron::Vertex_handle vh(heit->vertex());
                    f_indices.push_back(poly_vertices_indices[vh]);
                } while (++heit != fit->facet_begin());
                std::sort(f_indices.begin(), f_indices.end());

                // Search the polyhedron triangle in the Delaunay faces
                auto search(delaunay_triangles_by_indices.find(f_indices));
                if (search != delaunay_triangles_by_indices.end()) {
                    Delaunay::Simplex s(search->second);
//                    simplices[s] = *fit;
                    simplices.insert(s);
                }
                else {
                    std::cout << "triangle " ;
                    for (size_t i : f_indices)
                        std::cout << i << " ";
                    std::cout << " not in Delaunay" << std::endl;
                    res = false;
                }
            }
        }
        return res;
    }

    ConleyComplex(Polyhedron& poly, double thresh) : FlowComplex(poly), threshold(thresh) {
        // Init comparison function

        // Compute flow complex
        this->compute_cells();

        // Init the map of edges values and the sparse matrix encoding the poset
        poset = Column_matrix(this->number_of_flow_cells(), this->number_of_flow_cells());
        for (size_t i=0; i<this->flowcell_faces.size(); ++i) {
            for (size_t j : this->flowcell_faces.at(i)) {
                if (has_valid_value_for_merge(i,j)) {
                    auto e(std::make_pair(i, j));
                    values[e] = value_diff_df(i,j);
                    //                std::cout << "### " << i << ", " << j << ": " << values[e] << std::endl;
                    CGAL::OSM::set_coefficient(poset, i, j, values[e]);
                }
            }
        }
        std::cout << poset ;

        // Get sorted edges
        std::set<Poset_edge> sorted_edges(get_sorted_edges());

        std::cout << "M = [";
        for (Poset_edge edge : sorted_edges)
            std::cout << edge.value << " ";
        std::cout << "]" << std::endl;

        // Init IISs (initially, one IIS by flow_cell/critical point/critical simplex)
        for (size_t i=0; i<this->number_of_flow_cells(); ++i) {
            IIS iis;
            iis.cells_ids.insert(i);
            iis.designee_id = i;
            iiss[i] = iis;
        }

        // Merge edges with value < threshold
        Poset_edge edge(*(sorted_edges.begin()));
        while ((sorted_edges.size() > 0) && (edge.value < threshold)) {
            if (!is_valid_for_merge(edge)) {
                std::cout << "NOT merged " << edge.first << " -> " << edge.second << std::endl;
                sorted_edges.erase(edge);
            }
            else {
                std::cout << "===> merge " << edge.first << " -> " << edge.second << "(value: " << edge.value << ")" << std::endl;
                merge(edge);
//                std::cout << poset;
                sorted_edges = get_sorted_edges();
            }
            edge = *(sorted_edges.begin());
        }
    }

    size_t number_of_conley_iis() { return iiss.size(); }

    void print_infos() {
        std::cout << "number of iis: " << iiss.size() << std::endl;
        for (auto it = iiss.begin(); it != iiss.end(); ++it) {
            std::cout << "iis " << it->first << "(" << it->second.cells_ids.size() << ") -> ";
            size_t cpt(0);
            for (size_t i : it->second.cells_ids) {
                std::cout << i << " ";
                cpt += this->flowcell_from_id(i).get_simplices().size();
            }
            std::cout << " (totale size: " << cpt << ")" << std::endl;
        }
    }

    void write_vtk(std::string filename) {

        size_t id_cell, id_flowcell;
        std::vector<std::vector<int> > conley_ids(4), flow_ids(4); // set to iis index
        const size_t end_cells(number_of_conley_iis()+1), end_flowcells(this->number_of_flow_cells()+1);
        conley_ids.at(0).resize(m_dela.number_of_vertices(),end_cells);
        conley_ids.at(1).resize(m_dela.number_of_finite_edges(),end_cells);
        conley_ids.at(2).resize(m_dela.number_of_finite_facets(),end_cells);
        conley_ids.at(3).resize(m_dela.number_of_finite_cells(),end_cells);
        flow_ids.at(0).resize(m_dela.number_of_vertices(),end_flowcells);
        flow_ids.at(1).resize(m_dela.number_of_finite_edges(),end_flowcells);
        flow_ids.at(2).resize(m_dela.number_of_finite_facets(),end_flowcells);
        flow_ids.at(3).resize(m_dela.number_of_finite_cells(),end_flowcells);

        // Build flow_ids
        // Create a map flowcell_id -> iis id
        std::map<size_t,size_t> flowcell_id_to_iis_id;
        for (auto it = iiss.begin(); it != iiss.end(); ++it) {
            const size_t index(it->first);
            for (size_t i : it->second.cells_ids)
                flowcell_id_to_iis_id[i] = index;
        }

        // Vertices
        size_t cpt(0);
        for (typename Delaunay::Finite_vertices_iterator it = m_dela.finite_vertices_begin(); !(it == m_dela.finite_vertices_end()); ++it) {
            Delaunay::Simplex s(it);
            id_cell = flowcell_id_to_iis_id[simplex_to_flowcell_id[s]];
            id_flowcell = simplex_to_flowcell_id[s];
            flow_ids.at(0).at(cpt) = id_flowcell;
            conley_ids.at(0).at(cpt++) = id_cell;
        }

        // Edges
        cpt = 0;
        for (typename Delaunay::Edge edge : m_dela.finite_edges()) {
            Delaunay::Simplex s(edge);
            id_cell = flowcell_id_to_iis_id[simplex_to_flowcell_id[s]];
            id_flowcell = simplex_to_flowcell_id[s];
            flow_ids.at(1).at(cpt) = id_flowcell;
            conley_ids.at(1).at(cpt++) = id_cell;
        }

        // Triangles
        cpt = 0;
        for (typename Delaunay::Facet facet : m_dela.finite_facets()) {
            Delaunay::Simplex s(facet);
            id_cell = flowcell_id_to_iis_id[simplex_to_flowcell_id[s]];
            id_flowcell = simplex_to_flowcell_id[s];
            flow_ids.at(2).at(cpt) = id_flowcell;
            conley_ids.at(2).at(cpt++) = id_cell;
        }

        // Cells
        cpt = 0;
        for (typename Delaunay::Cell_handle cell : m_dela.finite_cell_handles()) {
            Delaunay::Simplex s(cell);
            id_cell = flowcell_id_to_iis_id[simplex_to_flowcell_id[s]];
            id_flowcell = simplex_to_flowcell_id[s];
            flow_ids.at(3).at(cpt) = id_flowcell;
            conley_ids.at(3).at(cpt++) = id_cell;
        }

        // Export the Delaunay mesh with these flags
        write_VTK(m_dela, filename, CELLS|FACETS|EDGES|VERTICES, &conley_ids, &flow_ids);
    }

protected:


    std::set<Poset_edge> get_sorted_edges() {
        std::set<Poset_edge> sorted_edges;
        for (CGAL::OSM::Bitboard::iterator it = poset.begin(); it != poset.end(); ++it) {
            const size_t i(*it);
            const Row_chain& row(CGAL::OSM::cget_row(poset, i));
            for (Row_chain::const_iterator it2 = row.cbegin(); it2 != row.cend(); ++it2) {
                const size_t j(it2->first);
                const auto e(std::make_pair(i,j));
                Poset_edge pe;
                pe.first = i;
                pe.second = j;
                pe.value = CGAL::OSM::get_coefficient(poset,i,j);
                sorted_edges.insert(pe);
            }
        }
        return sorted_edges;
    }

    bool is_valid_for_merge(const Poset_edge& edge) const {
        // Two IIS can be merged if they are connected by a single path (edge)
        // Flood from the first edge to check if a second path exists
        size_t cell_id(edge.first);
        std::vector<bool> flags(this->number_of_flow_cells(), false); // Visited/non visited
        std::queue<size_t> waiting_list;
        waiting_list.push(edge.first);
        flags.at(cell_id) = true;
        while (!waiting_list.empty()) {
            cell_id = waiting_list.front();
            waiting_list.pop();
            // Visit cells such that cell_id < cell in the poset
            const Row_chain& row(CGAL::OSM::cget_row(poset,cell_id));
            for (Row_chain::const_iterator it = row.cbegin(); it != row.cend(); ++it) {
                const size_t cell_id2(it->first);
                if ( !((cell_id == edge.first) && (cell_id2 == edge.second)) ) { // We do not consider the initial edge
                    if (cell_id2 == edge.second)
                        return false;
                    else {
                        if (!flags.at(cell_id2)) {
                            // First visit
                            waiting_list.push(cell_id2);
                            flags.at(cell_id2) = true;
                        }
                    }
                }
            }
        }
        return true;
    }

    template<typename ChainType>
    static ChainType merge_chains_min (const ChainType& c1, const ChainType& c2) {
        ChainType res(c1);
        for(typename ChainType::const_iterator it = c2.cbegin(); it != c2.cend(); ++it) {
            size_t i(it->first);
            if (c1.get_coefficient(i) != 0) {
                // min of both coefficients (in abs value)
                double coef1(c1.get_coefficient(i)), coef2(c2.get_coefficient(i));
                if (coef1 < coef2)
                    res.set_coefficient(i,coef1);
                else
                    res.set_coefficient(i,coef2);
            }
            else {
                // set coefficient of c2
                res.set_coefficient(i,c2.get_coefficient(i));
            }
        }
        return res;
    }

    template<typename ChainType>
    static ChainType merge_chains_sum (const ChainType& c1, const ChainType& c2) {
        ChainType res(c1);
        res += c2;
        return res;
    }

    // Merge two IIS into the union IIS
    // Pre: is_valid_for_merge(edge)
    void merge(const Poset_edge& edge) {
        const size_t v1(edge.first), v2(edge.second);
        const size_t n1(this->flowcell_from_id(v1).get_simplices().size()), n2(this->flowcell_from_id(v2).get_simplices().size());
        size_t v((n1<n2)?v2:v1), vv((n1<n2)?v1:v2); // Keep the index of the larget initial flow cell
//        std::cout << "v: " << v << " - vv: " << vv << std::endl;

        // Get row/columns v1 and v2
        Row_chain r1(CGAL::OSM::get_row(poset, v1));
        Row_chain r2(CGAL::OSM::get_row(poset, v2));
        Column_chain c1(CGAL::OSM::get_column(poset, v1));
        Column_chain c2(CGAL::OSM::get_column(poset, v2));

        // Empty v1/v2 coefficients
        r1 /= std::vector<size_t>({v1,v2});
        r2 /= std::vector<size_t>({v1,v2});
        c1 /= std::vector<size_t>({v1,v2});
        c2 /= std::vector<size_t>({v1,v2});

        // Compute "min" of chains
        Row_chain r3(merge_chains_min(r1,r2));
        Column_chain c3(merge_chains_min(c1,c2));

        // Compute "sum" of chains
//        Row_chain r3(merge_chains_sum(r1,r2));
//        Column_chain c3(merge_chains_sum(c1,c2));

        // Empty vv row/column
        CGAL::OSM::remove_column(poset, vv);
        CGAL::OSM::remove_row(poset, vv);

        // Set v row/column
        CGAL::OSM::set_column(poset, v, c3);
        CGAL::OSM::set_row(poset, v, r3);

        // Update IIS
        // Union of cells in v
        iiss.at(v).cells_ids.insert(iiss.at(vv).cells_ids.begin(), iiss.at(vv).cells_ids.end());

        // Remove vv
        iiss.erase(vv);

    }

};

#endif // !CONLEY_COMPLEX_H
