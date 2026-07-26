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
public:
    // For OSM sparse matrices
    typedef double Coefficient_type;
    typedef CGAL::OSM::Sparse_matrix<CGAL::OSM::Sparse_chain> Sparse_matrix_struct;
    typedef Sparse_matrix_struct::template Sparse_matrix_type<Coefficient_type, CGAL::OSM::COLUMN> Column_matrix;
    typedef Sparse_matrix_struct::template Sparse_matrix_type<Coefficient_type, CGAL::OSM::ROW> Row_matrix;
    typedef Sparse_matrix_struct::template Sparse_chain_type<Coefficient_type, CGAL::OSM::COLUMN> Column_chain;
    typedef Sparse_matrix_struct::template Sparse_chain_type<Coefficient_type, CGAL::OSM::ROW> Row_chain;

    // Edge structure
    struct Poset_edge {
        size_t first, second;
        double value;
        bool operator< (const Poset_edge& other) const { return abs(value) < abs(other.value); };
    };

    // IIS structure
    struct IIS {
        std::set<size_t> cells_ids; // Set of cells merged into the IIS
        size_t designee_id; // Id of the "representative" critical cell (lowest id among cells of the IIS)
    };

protected:
    std::function<bool(size_t i, size_t j)> compare_cells;

    // Store the poset (weighted by distances between cells) in a sparse OSM matrix
    Row_matrix poset;
    std::vector<double> values;
    std::map<size_t, IIS> iiss;
    double threshold;

public:
    // Returns the difference between distance to border of critical cells
    // TODO: adjust...
    double value (size_t cell_id1) {
        return (DelaunayHelper::get_critical_info(m_dela,this->flowcell_from_id(cell_id1).get_critical_simplex())).r;
    }


    ConleyComplex(Polyhedron& poly, double thresh) : FlowComplex(poly), threshold(thresh) {
        // Init comparison function
        compare_cells = [&](size_t i, size_t j) { return abs(values.at(i)) < abs(values.at(j)); } ;

        // Compute flow complex
        this->compute_cells();
        // Init values
        for (size_t i=0; i<this->number_of_flow_cells(); ++i) {
            values.push_back(value(i));
        }
        // Init the sparse matrix encoding the poset
        poset = Column_matrix(this->number_of_flow_cells(), this->number_of_flow_cells());
        for (size_t i=0; i<this->flowcell_faces.size(); ++i) {
            for (size_t j : this->flowcell_faces.at(i)) {
                CGAL::OSM::set_coefficient(poset, i, j, values.at(j)-values.at(i));
            }
        }
        std::cout << poset ;

        // Get sorted edges
        std::set<Poset_edge> sorted_edges(get_sorted_edges());

        std::cout << "M = [";
        for (Poset_edge edge : sorted_edges)
            std::cout << edge.value << " ";
        std::cout << "]" << std::endl;

        // Init IISs (initially, one critical point by IIS)
        for (size_t i=0; i<this->number_of_flow_cells(); ++i) {
            IIS iis;
            iis.cells_ids.insert(i);
            iis.designee_id = i;
            iiss[i] = iis;
        }

        // Merge edges with value < threshold
        Poset_edge edge(*(sorted_edges.begin()));
        while ((sorted_edges.size() > 0) && (abs(edge.value) < threshold)) {
            if (!is_valid_for_merge(edge)) {
//                std::cout << "NOT merged " << edge.first << " -> " << edge.second << std::endl;
                sorted_edges.erase(edge);
            }
            else {
//                std::cout << "merge " << edge.first << " -> " << edge.second << std::endl;
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
            for (size_t i : it->second.cells_ids)
                std::cout << i << " ";
            std::cout << std::endl;
        }
    }

    void write_vtk(std::string filename) {

        size_t id_cell;
        std::vector<std::vector<int> > conley_ids(4); // set to iis index
        const size_t end_cells(number_of_conley_iis()+1);
        conley_ids.at(0).resize(m_dela.number_of_vertices(),end_cells);
        conley_ids.at(1).resize(m_dela.number_of_finite_edges(),end_cells);
        conley_ids.at(2).resize(m_dela.number_of_finite_facets(),end_cells);
        conley_ids.at(3).resize(m_dela.number_of_finite_cells(),end_cells);

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
            conley_ids.at(0).at(cpt++) = id_cell;
        }

        // Edges
        cpt = 0;
        for (typename Delaunay::Edge edge : m_dela.finite_edges()) {
            Delaunay::Simplex s(edge);
            id_cell = flowcell_id_to_iis_id[simplex_to_flowcell_id[s]];
            conley_ids.at(1).at(cpt++) = id_cell;
        }

        // Triangles
        cpt = 0;
        for (typename Delaunay::Facet facet : m_dela.finite_facets()) {
            Delaunay::Simplex s(facet);
            id_cell = flowcell_id_to_iis_id[simplex_to_flowcell_id[s]];
            conley_ids.at(2).at(cpt++) = id_cell;
        }

        // Cells
        cpt = 0;
        for (typename Delaunay::Cell_handle cell : m_dela.finite_cell_handles()) {
            Delaunay::Simplex s(cell);
            id_cell = flowcell_id_to_iis_id[simplex_to_flowcell_id[s]];
            conley_ids.at(3).at(cpt++) = id_cell;
        }

        // Export the Delaunay mesh with these flags
        write_VTK(m_dela, filename, CELLS|FACETS|EDGES|VERTICES, &conley_ids);
    }

protected:


    std::set<Poset_edge> get_sorted_edges() {
        std::set<Poset_edge> sorted_edges;
        for (CGAL::OSM::Bitboard::iterator it = poset.begin(); it != poset.end(); ++it) {
            const Row_chain& col(CGAL::OSM::cget_row(poset, *it));
            const size_t i(*it);
            for (Column_chain::const_iterator it2 = col.begin(); it2 != col.end(); ++it2) {
                const size_t j(it2->first);
                sorted_edges.insert(Poset_edge({i,j,values.at(j) - values.at(i)}));
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
                if ((cell_id != edge.first) && (cell_id2 != edge.second)) { // We do not consider the initial edge
                    if (cell_id2 == edge.second)
                        return false;
                    else {
                        if (!flags.at(cell_id2)) // First visit
                            waiting_list.push(cell_id2);
                    }
                }
            }
        }
        return true;
    }

    template<typename ChainType>
    static ChainType min_chains (const ChainType& c1, const ChainType& c2) {
        ChainType res(c1);
        for(typename ChainType::const_iterator it = c2.cbegin(); it != c2.end(); ++it) {
            const size_t i(c1.get_coefficient(it->first));
            if (c1.get_coefficient(i) != 0) {
                // min of both coefficients (in abs value)
                double coef1(c1.get_coefficient(i)), coef2(c2.get_coefficient(i));
                if (abs(coef1) < abs(coef2))
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

    // Merge two IIS into the union IIS
    // Pre: is_valid_for_merge(edge)
    void merge(const Poset_edge& edge) {
        const size_t v1(edge.first), v2(edge.second);
        size_t v((v1<v2)?v1:v2), vv((v1<v2)?v2:v1);
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
        Row_chain r3(min_chains(r1,r2));
        Column_chain c3(min_chains(c1,c2));

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
