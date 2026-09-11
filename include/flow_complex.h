#ifndef FLOW_COMPLEX_H
#define FLOW_COMPLEX_H


#include "cgal_typedef.h"
// #include "filtration.h"
#include "delaunay_helper.h"

#include <fstream>
#include <limits>

#define STRICT_SIGN(x) (((x)<0.0) ? -1 : 1)
#define SIGN(x) (((x)<-0.000001) ? -1 : (((x)<0.000001)? 0 : 1 ))
#define DEBUG(i) std::clog << "DEBUG " << i << std::endl
inline int sign(double d){ return (d>=0.0) ? 1 : -1;}

//std::ostream& print_simplex(std::ostream& os,
//                            const Delaunay::Simplex& s)
//{
//    switch (s.dimension()) {
//    case 0: {
//        auto vh = Delaunay::Vertex_handle(s);
//        os << "Vertex(" << vh->point() << ")";
//        break;
//    }
//    case 1: {
//        auto e = Delaunay::Edge(s);
//        os << "Edge("
//           << e.first->vertex(e.second)->point() << ", "
//           << e.first->vertex(e.third)->point() << ")";
//        break;
//    }
//    case 2: {
//        auto f = Delaunay::Facet(s);
//        os << "Facet(";
//        Delaunay::Cell_handle ch = f.first;
//        int i = f.second;
//        os << ch->vertex( (i+1)%4 )->point() << ", ";
//        os << ch->vertex( (i+2)%4 )->point() << ", ";
//        os << ch->vertex( (i+3)%4 )->point() << ")";
//        break;
//    }
//    case 3: {
//        auto c = Delaunay::Cell_handle(s);
//        os << "Cell(";
//        os << c->vertex(0)->point() << ", ";
//        os << c->vertex(1)->point() << ", ";
//        os << c->vertex(2)->point() << ", ";
//        os << c->vertex(3)->point() << ")";
//        break;
//    }
//    default:
//        os << "Invalid simplex";
//    }
//    return os;
//}

std::ostream& print_simplex(std::ostream& os,
                            const Delaunay::Simplex& s)
{
    switch (s.dimension()) {
    case 0: {
        auto vh = Delaunay::Vertex_handle(s);
        os << "Vertex(" << vh->info().second << ")";
        break;
    }
    case 1: {
        auto e = Delaunay::Edge(s);
        os << "Edge("
           << e.first->vertex(e.second)->info().second << ", "
           << e.first->vertex(e.third)->info().second << ")";
        break;
    }
    case 2: {
        auto f = Delaunay::Facet(s);
        os << "Facet(";
        Delaunay::Cell_handle ch = f.first;
        int i = f.second;
        os << ch->vertex( (i+1)%4 )->info().second << ", ";
        os << ch->vertex( (i+2)%4 )->info().second << ", ";
        os << ch->vertex( (i+3)%4 )->info().second << ")";
        break;
    }
    case 3: {
        auto c = Delaunay::Cell_handle(s);
        os << "Cell(";
        os << c->vertex(0)->info().second << ", ";
        os << c->vertex(1)->info().second << ", ";
        os << c->vertex(2)->info().second << ", ";
        os << c->vertex(3)->info().second << ")";
        break;
    }
    default:
        os << "Invalid simplex";
    }
    return os;
}

/**
 * @brief The FlowCell class
 * Contains a reference to the representant critical delaunay simplex
 * and a list of the corresponding delaunay simplices that flow towards it 
 * (flow goes toward the boundary).
 */
class FlowCell
{
private:
    // static int currID;
    // int ID;
    const Delaunay::Simplex crit_simplex;
    std::set<Delaunay::Simplex> del_simplices;
    float df;
    Delaunay::Point p;

public:
    FlowCell(const Delaunay::Simplex _crit_simplex) : //ID(currID++),
    crit_simplex(_crit_simplex) {
        del_simplices.insert(crit_simplex);
        df = 0.0;
    }
    // for the hash?
    bool operator<(const FlowCell& other) const
    {
        return get_critical_simplex() < other.get_critical_simplex();
    }
    bool operator==(const FlowCell& other) const
    {
        return !(crit_simplex < other.crit_simplex) &&
               !(other.crit_simplex < crit_simplex);
    }
    // size_t operator()(const FlowCell& f) const
    // {
    //     return std::hash<Delaunay::Simplex>crit_simplex();
    // }
    
    void add_simplex(const Delaunay::Simplex& s)
    {
        del_simplices.insert(s);
    }

    bool contains(const Delaunay::Simplex& s) const
    {
        return del_simplices.find(s) != del_simplices.end();
    }
    
    std::set<Delaunay::Simplex> get_simplices() const {return del_simplices;}
    Delaunay::Simplex get_critical_simplex() const {return crit_simplex;}
    int dimension() const {return crit_simplex.dimension();}
    void set_df(float val) {df = val;}
    float get_df() const { return df; }
    void set_p(Delaunay::Point pt) { p = pt; }
    Delaunay::Point get_p() { return p; }
    // int id() const { return ID; }
    
    std::ostream& print(std::ostream& os) const
    {
        os << "FlowCell {\n"
           << "  dimension: " << dimension() << '\n'
           << "  df: " << df << '\n'
           << "  critical simplex: ";
        print_simplex(os, crit_simplex) << '\n';

        os << "  simplices (" << del_simplices.size() << "):\n";
        for (const auto& s : del_simplices)
            print_simplex(os, s) << ' ';
        os << "}";
        return os;
    }
    friend std::ostream& operator<<(std::ostream& os, const FlowCell& fc);
};

inline std::ostream& operator<<(std::ostream& os, const FlowCell& fc)
{
    os << "FC{" << fc.dimension() << ", ";
    print_simplex(os, fc.crit_simplex) << ", ";
    os << "[" << fc.del_simplices.size() << "]}";
    return os;
}

/*----------------------------------------------------------------------------*/

/**
* @brief The FlowComplex class
* 
*/
class FlowComplex
{
public:
    typedef Delaunay::Simplex Simplex;
    typedef CGAL::Side_of_triangle_mesh<Polyhedron, Epick>                                  SoT;
    struct Simplex_id {
        size_t dim, i;
    };
    // Member data
    Polyhedron& m_poly;
    Delaunay m_dela;
    std::map<Delaunay::Simplex, Simplex_id> simplex_to_simplex_id; // Map for HDVF indexes
protected:
    std::vector<FlowCell> flowcells;
    std::map<Simplex, size_t> simplex_to_flowcell_id;
    std::vector<std::set<size_t>> flowcell_faces; // encoding the poset
    /* WARNING I think there might be problems with the algorithm, for instance if a flowcell is twice adjacent to another cell...
    * In this case, it will be counted only once in flowcell_faces...
    */

public:
    FlowComplex(Polyhedron& poly) : m_poly(poly)
    {

        //CGAL::draw(m_poly);
        // CGAL::Side_of_triangle_mesh<Polyhedron, Epick> inside(m_poly);
        size_t cpt(0);
        for (Polyhedron::Vertex_iterator it = m_poly.vertices_begin(); it != m_poly.vertices_end(); ++it)
        {
            m_dela.insert(it->point()); // m_dela corresponds to the 3D Delaunay of the surface points
            //        std::cout << it->point() << std::endl;
        }

        // Init Delaunay vertices indices info
        cpt=0;
        for(Delaunay::Vertex_handle vh : m_dela.finite_vertex_handles()) {
            vh->info().second = cpt++;
        }

        std::clog << "Computing flow complex. Delaunay triangulation:" << std::endl;
        std::clog << "nb finite vertices: " << m_dela.number_of_vertices() << std::endl;
        std::clog << "nb edges   : " << m_dela.number_of_edges()  << " | finite: " << m_dela.number_of_finite_edges() << std::endl;
        std::clog << "nb facets  : " << m_dela.number_of_facets() << " | finite: " << m_dela.number_of_finite_facets()<< std::endl;
        std::clog << "nb cells   : " << m_dela.number_of_cells()  << " | finite: " << m_dela.number_of_finite_cells()<< std::endl;
        // CGAL::draw(m_dela);



        // Build indexes of vertices
        cpt=0;
        for (Delaunay::All_vertices_iterator vit = m_dela.all_vertices_begin(); vit != m_dela.all_vertices_end(); vit++) {
            Delaunay::Simplex s = Delaunay::Simplex(vit);
            simplex_to_simplex_id[s] = Simplex_id({0,cpt++});
        }
        // Build indexes of edges
        cpt = 0;
        for (Delaunay::All_edges_iterator eit = m_dela.all_edges_begin(); eit != m_dela.all_edges_end(); eit++) {
            Delaunay::Simplex s = Delaunay::Simplex(*eit);
            simplex_to_simplex_id[s] = Simplex_id({1,cpt++});
        }
        // Build indexes of facets
        cpt = 0;
        for (Delaunay::All_facets_iterator fit = m_dela.all_facets_begin(); fit != m_dela.all_facets_end(); fit++) {
            Delaunay::Simplex s = Delaunay::Simplex(*fit);
            simplex_to_simplex_id[s] = Simplex_id({2,cpt++});
        }
        // Build indexes of cells
        cpt = 0;
        for (Delaunay::All_cells_iterator cit = m_dela.all_cells_begin(); cit != m_dela.all_cells_end(); cit++) {
            Delaunay::Simplex s = Delaunay::Simplex(cit);
            simplex_to_simplex_id[s] = Simplex_id({3,cpt++});
        }
        
        /* Finite version */
        // // Build indexes of vertices
        // cpt=0;
        // for (Delaunay::Finite_vertices_iterator vit = m_dela.finite_vertices_begin(); vit != m_dela.finite_vertices_end(); vit++) {
        //     Delaunay::Simplex s = Delaunay::Simplex(vit);
        //     simplex_to_simplex_id[s] = Simplex_id({0,cpt++});
        // }
        // // Build indexes of edges
        // cpt = 0;
        // for (Delaunay::Finite_edges_iterator eit = m_dela.finite_edges_begin(); eit != m_dela.finite_edges_end(); eit++) {
        //     Delaunay::Simplex s = Delaunay::Simplex(*eit);
        //     simplex_to_simplex_id[s] = Simplex_id({1,cpt++});
        // }
        // // Build indexes of facets
        // cpt = 0;
        // for (Delaunay::Finite_facets_iterator fit = m_dela.finite_facets_begin(); fit != m_dela.finite_facets_end(); fit++) {
        //     Delaunay::Simplex s = Delaunay::Simplex(*fit);
        //     simplex_to_simplex_id[s] = Simplex_id({2,cpt++});
        // }
        // // Build indexes of cells
        // cpt = 0;
        // for (Delaunay::Finite_cells_iterator cit = m_dela.finite_cells_begin(); cit != m_dela.finite_cells_end(); cit++) {
        //     Delaunay::Simplex s = Delaunay::Simplex(cit);
        //     simplex_to_simplex_id[s] = Simplex_id({3,cpt++});
        // }
    }

    
    /////////////////////////////
    //         / \             //
    //        / ! \            //
    //       /_____\ WARNING   //
    /////////////////////////////
    /* /!\ In this code there is two definitions of flowcell_id
     * one is just the id of the representative critical cell (for instance get_flowcell_id())
     * one is an index that works with flowcells (e.g. flowcell_from_simplex, and simplex_to_flowcell_id)
    */
    
    FlowCell flowcell_from_id(size_t i) {
        return flowcells.at(i);
    }
    
    /**
     * \brief Returns the index of the flowcell containing the simplex.
     *
     * \pre Flow cells must have been computed.
    */
    size_t flowcell_from_simplex(const Simplex& simplex) {
        return simplex_to_flowcell_id[simplex];
    }

    /**
    * @brief return the flowcell id of a given flowcell f.
    * In fact, it just check the id associated to the critical cell of f.
    */
    size_t get_flowcell_id(const FlowCell& f) const
    {
      if (auto search = simplex_to_flowcell_id.find(f.get_critical_simplex()); search != simplex_to_flowcell_id.end()) {
          return search->second;
      }
      std::cerr << "Error in get_flowcell_id: crit_simplex not in the map crit_simplex_to_flowcell_id." << std::endl;
      return 0;
    }

    const Delaunay& delaunay_mesh() { return m_dela; }
    
    ////////////////////////// COMPUTATIONS ////////////////////////////////////
    
    /** \brief  Main computation for the flowcomplex
     * First, initialize flowwcells from critical cells
     * For each dimension, up flows to build flowcells
     */
    void compute_cells () {
        init_flowcells_from_critical_cells();
        std::clog << std::endl << "0D" << std::endl;
        for (Delaunay::Finite_vertices_iterator vit = m_dela.finite_vertices_begin();
        vit != m_dela.finite_vertices_end(); vit++) {
            Delaunay::Simplex s = Delaunay::Simplex(vit);
            CriticalInfo crit_info(DelaunayHelper::get_critical_info(m_dela, s));
            if (crit_info.c == CriticalType::Critical){
                FlowCell fc = flowcell_from_critical_cell(s, crit_info);
                std::clog << "\n" << get_flowcell_id(fc) << " '---> " << fc << "\n"; // to display the built flowcells
            }
        }

        std::clog << std::endl << "1D" << std::endl;
        for (Delaunay::Finite_edges_iterator eit = m_dela.finite_edges_begin();
        eit != m_dela.finite_edges_end(); eit++) {
            Delaunay::Simplex s = Delaunay::Simplex(*eit);
            CriticalInfo crit_info(DelaunayHelper::get_critical_info(m_dela, s));
            if (crit_info.c == CriticalType::Critical){
                FlowCell fc = flowcell_from_critical_cell(s, crit_info);
                std::clog << "\n" << get_flowcell_id(fc) << " '---> " << fc << "\n"; // to display the built flowcells
            }
        }

        std::clog << std::endl << "2D" << std::endl;
        for (Delaunay::Finite_facets_iterator fit = m_dela.finite_facets_begin();
        fit != m_dela.finite_facets_end(); fit++) {
            Delaunay::Simplex s = Delaunay::Simplex(*fit);
            CriticalInfo crit_info(DelaunayHelper::get_critical_info(m_dela, s));
            if (crit_info.c == CriticalType::Critical){
                FlowCell fc = flowcell_from_critical_cell(s, crit_info);
                std::clog << "\n" << get_flowcell_id(fc) << " '---> " << fc << "\n"; // to display the built flowcells
            }
        }
        
        std::clog << std::endl << "3D" << std::endl;
        FlowCell fc = flowcell_from_infinite_cells();
        std::clog << "\n" << get_flowcell_id(fc) << " '---> " << fc << "    INF FLOW CELL\n"; // to display the built flowcells
        for (Delaunay::Finite_cells_iterator cit = m_dela.finite_cells_begin();
        cit != m_dela.finite_cells_end(); cit++) {
            Delaunay::Simplex s = Delaunay::Simplex(cit);
            CriticalInfo crit_info(DelaunayHelper::get_critical_info(m_dela, s));
            if (crit_info.c == CriticalType::Critical){
                FlowCell fc = flowcell_from_critical_cell(s, crit_info);
                std::clog << "\n" << get_flowcell_id(fc) << " '---> " << fc << "\n"; // to display the built flowcells
            }
        }
    }

    /**
    * @brief this function builds a flow cell and update the corresponding variables.
    * It should update:
    * std::vector<FlowCell> flowcells;
    * std::map<Simplex, size_t> simplex_to_flowcell_id;
    * std::vector<std::set<size_t>> flowcell_faces;
    *
    * The main idea is that we start from the critical cell and add recursively
    * the simplices given by up_flow_cells.
    */

    void init_flowcell_from_critical_cell(const Simplex crit_simplex, const CriticalInfo& crit_info)
    {
        FlowCell fc(crit_simplex);
        size_t id = flowcells.size();

        if (crit_info.c != CriticalType::Critical)
        {
            std::cerr << "Error in flowcell_from_critical_cell: simplex ";
            print_simplex(std::clog, crit_simplex) << " not critical" << std::endl;
        }
        fc.set_df(crit_info.r);
        fc.set_p(crit_info.p);

        flowcells.push_back(fc);
        simplex_to_flowcell_id[crit_simplex] = id;
        std::clog << "Critical cell: " << simplex_to_simplex_id[crit_simplex].i << " - " << simplex_to_simplex_id[crit_simplex].dim << std::endl;
    }
    
    void init_flowcell_from_infinite_cells()
    {
        Delaunay::Vertex_handle inf_v = m_dela.infinite_vertex();
        FlowCell fc(inf_v);
        size_t id = flowcells.size();
        std::clog << "BEGIN init_flowcell_from_infinite_cells: ";
        print_simplex(std::clog, Simplex(inf_v)) << " infinite vertex" << std::endl;
        std::clog << fc << "\n";

        fc.set_df(0.0);// should be INF
        fc.set_p(Delaunay::Point(0.0,0.0,0.0));// should be INF

        flowcells.push_back(fc);
        simplex_to_flowcell_id[Simplex(inf_v)] = id;
        std::clog << "END init_flowcell_from_infinite_cells" << std::endl;
    }
    
    void init_flowcells_from_critical_cells() {
        std::clog << "BEGIN init_flowcells_from_critical_cells" << std::endl;
        init_flowcell_from_infinite_cells();
        
        for (Delaunay::Finite_vertices_iterator vit = m_dela.finite_vertices_begin();
        vit != m_dela.finite_vertices_end(); vit++) {
            Delaunay::Simplex s = Delaunay::Simplex(vit);
            CriticalInfo crit_info(DelaunayHelper::get_critical_info(m_dela, s));
            if (crit_info.c == CriticalType::Critical){
                init_flowcell_from_critical_cell(s, crit_info);
            }
        }

        for (Delaunay::Finite_edges_iterator eit = m_dela.finite_edges_begin();
        eit != m_dela.finite_edges_end(); eit++) {
            Delaunay::Simplex s = Delaunay::Simplex(*eit);
            CriticalInfo crit_info(DelaunayHelper::get_critical_info(m_dela, s));
            if (crit_info.c == CriticalType::Critical){
                init_flowcell_from_critical_cell(s, crit_info);
            }
        }

        for (Delaunay::Finite_facets_iterator fit = m_dela.finite_facets_begin();
        fit != m_dela.finite_facets_end(); fit++) {
            Delaunay::Simplex s = Delaunay::Simplex(*fit);
            CriticalInfo crit_info(DelaunayHelper::get_critical_info(m_dela, s));
            if (crit_info.c == CriticalType::Critical){
                init_flowcell_from_critical_cell(s, crit_info);
            }
        }

        for (Delaunay::Finite_cells_iterator cit = m_dela.finite_cells_begin();
        cit != m_dela.finite_cells_end(); cit++) {
            Delaunay::Simplex s = Delaunay::Simplex(cit);
            CriticalInfo crit_info(DelaunayHelper::get_critical_info(m_dela, s));
            if (crit_info.c == CriticalType::Critical){
                init_flowcell_from_critical_cell(s, crit_info);
            }
        }
        std::clog << "END init_flowcells_from_critical_cells" << std::endl;
    }

    FlowCell flowcell_from_critical_cell(const Simplex crit_simplex, const CriticalInfo& crit_info)
    {
//        CriticalInfo crit_info = DelaunayHelper::get_critical_info(m_dela, crit_simplex);
        if (crit_info.c != CriticalType::Critical)
        {
            std::cerr << "Error in flowcell_from_critical_cell: simplex ";
            print_simplex(std::clog, crit_simplex) << " not critical" << std::endl;
        }
        // Assert if the flowcell has been initiated with its critical simplex
        auto search(simplex_to_flowcell_id.find(crit_simplex));
        if (search == simplex_to_flowcell_id.end()) {
            std::cerr << "Error in flowcell_from_critical_cell: flowcell non initiated";
            throw std::runtime_error("Error in flowcell_from_critical_cell: flowcell non initiated");
        }
        size_t id = search->second;
        FlowCell fc(flowcells.at(id));
        // print_simplex(std::clog, crit_simplex) << std::endl;

        // Init poset
        std::set<size_t> fc_faces = {};
         
        std::queue<Simplex> to_process;
        for (const Simplex& s : DelaunayHelper::D_faces(crit_simplex)) {
          to_process.push(s);
        }
        
        while (!to_process.empty()) {
            const Simplex s = to_process.front();
            to_process.pop();

            if (auto search = simplex_to_flowcell_id.find(s); search != simplex_to_flowcell_id.end()) {
                // AB: possible error here when search->second is critical but not yet visited
                // AB fix: first visit all critical cells and init their flow_cell with them
                if (search->second != id) {
                    // std::cout << "    - poset: " << id << " -> " << search->second << "("; print_simplex(std::clog, s); std::cout << ")" << std::endl;
                    fc_faces.insert(search->second); // update the poset
                }
            }
            else {
                simplex_to_flowcell_id[s] = id;
                fc.add_simplex(s);
                flowcells.at(id).add_simplex(s);
                std::list<Simplex> ls = up_flow_cells(s);
                
                // std::clog << "    - up_flow : s.dimension = " << s.dimension() << " ";
                // print_simplex(std::clog, s) << "\t[ ";
                
                for (const auto& s_up_flow : ls) {
                    print_simplex(std::clog, s_up_flow);
                    std::clog << " ";
                    // std::clog << s_up_flow.dimension()<< " ";
                    if (s_up_flow.dimension() <= crit_simplex.dimension())
                        to_process.push(s_up_flow);
                    // else
                    //     std::clog << "(not added) "; // It seems that this case can happen in 3D: using our algo, if we start from a critical 2D cell we can arrive on a 3D cell... weird.
                }
                // std::clog << "]\n";
            }
        }
        flowcell_faces.push_back(fc_faces);
        return fc;
    }
    
    FlowCell flowcell_from_infinite_cells()
    {// NOT DEBBUGED YET
    std::clog << "BEGIN flowcell_from_infinite_cells" << std::endl;
        
        size_t id = 0; // WARNING, in the current implementation, the infinite flowcell is the 0 one
        FlowCell fc(flowcells.at(id));
        // print_simplex(std::clog, crit_simplex) << std::endl;
        Delaunay::Vertex_handle inf_v = m_dela.infinite_vertex();
        
        // Init poset
        std::set<size_t> fc_faces = {};
        
        std::queue<Simplex> to_process;
        
        std::vector<Delaunay::Cell_handle> inc_cells;
        m_dela.incident_cells(inf_v, std::back_inserter(inc_cells));
        for (const Delaunay::Cell_handle& cell : inc_cells) {
          to_process.push(Simplex(cell));// add all infinite cells
          std::clog << "incident cell of infinity: ";
          print_simplex(std::clog, Simplex(cell)) << "\n";
        }
        
        while (!to_process.empty()) {
            const Simplex s = to_process.front();
            to_process.pop();

            if (auto search = simplex_to_flowcell_id.find(s); search != simplex_to_flowcell_id.end()) {
                // AB: possible error here when search->second is critical but not yet visited
                // AB fix: first visit all critical cells and init their flow_cell with them
                if (search->second != id) {
                    // std::cout << "    - poset: " << id << " -> " << search->second << "("; print_simplex(std::clog, s); std::cout << ")" << std::endl;
                    fc_faces.insert(search->second); // update the poset
                }
            }
            else {
                simplex_to_flowcell_id[s] = id;
                fc.add_simplex(s);// TODO clean, I don't know why we are using two copies of the same fc
                flowcells.at(id).add_simplex(s);
                std::list<Simplex> ls = up_flow_cells(s);
                
                // std::clog << "    - up_flow : s.dimension = " << s.dimension() << " ";
                // print_simplex(std::clog, s) << "\t[ ";
                
                for (const auto& s_up_flow : ls) {
                    print_simplex(std::clog, s_up_flow);
                    std::clog << " ";
                    // std::clog << s_up_flow.dimension()<< " ";
                    to_process.push(s_up_flow);
                    // else
                    //     std::clog << "(not added) "; // It seems that this case can happen in 3D: using our algo, if we start from a critical 2D cell we can arrive on a 3D cell... weird.
                }
                // std::clog << "]\n";
            }
        }
        flowcell_faces.push_back(fc_faces);
        std::clog << "END flowcell_from_infinite_cells" << std::endl;
        return fc;
    }
    
    /**
    * @brief up_flow_cells(s) return the list of simplices that flows towards s.
    * Precisely, the flow is defined on the dual Voronoi, and is oriented 
    * towards going further away from points.
    * For example, if s is a delaunay critical 3-cell (corresponding to a critical
    * Voronoi vertex), then up_flow_cells(s) is the list of all its 2D faces.
    *
    * Note that simplices in up_flow_cells(s) are either of dimension
    * s.dimension()-1 or s.dimension()+1.
    */
    std::list<Simplex> up_flow_cells(const Simplex s)
    {
        std::list<Simplex> ls = {};
        if (s.dimension() == 0) // Delaunay vertex : Voronoi cell
        {
          // pass
        }
        else if (s.dimension() == 1) // Delaunay edge : Voronoi facet
        {
            // HARD one:
            /* 
            The idea is to think of it in the medial plane in between the two 
            vertices.
            The dual Voronoi is a convex polygon in this plane, and the driver is 
            the mid point in between the two vertices.
            The simplices that flow towards s are the ones whose dual is visible
            from the driver point (I believe).
            */
            const Delaunay::Edge e(s);
            
            if (m_dela.is_infinite(e)) {
                // if it is infinite, go to its finite faces
                std::list<Simplex> e_faces = DelaunayHelper::D_faces(s);
                for (const Simplex face : e_faces){
                    if (!DelaunayHelper::is_infinite(m_dela, face)) 
                        ls.push_back(face);
                }
                return ls;
            }
            
            //std::clog << "2-hole critical" << std::endl;
            const Delaunay::Point p0 = m_dela.point(e.first->vertex(e.second));
            const Delaunay::Point p1 = m_dela.point(e.first->vertex(e.third ));
            ls.push_back(Simplex(e.first->vertex(e.second)));
            ls.push_back(Simplex(e.first->vertex(e.third)));
            const Delaunay::Point driver = p0 + 0.5*(p1-p0);
            
            /* Previously used for test about scalar products: */
            // std::vector<Delaunay::Cell_handle> finite_incident_cells;
            // std::vector<Delaunay::Point> finite_dual_pts;
            // DelaunayHelper::get_finite_incident_cells(m_dela, finite_incident_cells, e);
            // for (Delaunay::Cell_handle cell : finite_incident_cells){
            //     finite_dual_pts.push_back(m_dela.dual(cell));
            // }

            Delaunay::Facet_circulator circ = m_dela.incident_facets(e), past_end(circ);
            do
            {
                Delaunay::Facet f = *circ;
                if (!m_dela.is_infinite(f)){
                  const Delaunay::Cell_handle ch = f.first;
                  const int i = f.second;
                  const Delaunay::Facet f_mirror = m_dela.mirror_facet(f);
                  const Delaunay::Cell_handle ch_mirror = f_mirror.first;
                  const int i_mirror = f_mirror.second; // not used
                  
                  if (m_dela.is_infinite(ch) && m_dela.is_infinite(ch_mirror)){
                    std::cerr << "Error in flowcell_from_critical_cell: in edge case, both incident cell_handle are infinite ";
                  }
                  else if (m_dela.is_infinite(ch) || m_dela.is_infinite(ch_mirror)){
                      // do nothing
                  }
                  else {
                    const Delaunay::Point dual_ch = m_dela.dual(ch);
                    const Delaunay::Point dual_ch_mirror = m_dela.dual(ch_mirror);
                    assert(!CGAL::collinear(dual_ch, dual_ch_mirror, dual_ch+(p0-p1)));// I don't know what to do if they are colinear...
                    Vector normal_edge = CGAL::unit_normal(dual_ch, dual_ch_mirror, dual_ch+(p0-p1));
                    
                    /*
                    /!\ from tests: the scalar product between normal_edge and
                    (pt - dual_ch) (with pt the dual of a finite incident cell)
                    is positive, or eventually zero if pt is dual_ch or dual_ch_mirror.
                    Hence, it points inward of the convex polygon dual of e.
                    The test are commented below:
                    */
                    // std::clog << "FACET ->\t";
                    // int count = 0;
                    // for (Delaunay::Point pt : finite_dual_pts){
                    //     int sign = SIGN(CGAL::scalar_product(normal_edge, pt-dual_ch));
                    //     std::clog << sign << "\t";
                    //     if (sign==0)
                    //         count ++;
                    //     if (sign == -1)
                    //         std::clog << " #####################################################################";
                    // }
                    // std::clog << "  ->  " << count;
                    // if (count != 2){
                    //     std::clog << " WRONG --------------------------------- ";
                    // }
                    // std::clog <<"\n";
                    
                    if (CGAL::scalar_product(normal_edge, dual_ch-driver)>= 0.0) {
                        ls.push_back(Simplex(f));
                    }
                  }
                }
                else {
                    std::cerr << "Warning in up_flow_cells(): infinite FACET\n";
                }      
            }while (++circ != past_end);
                    
        }
        else if (s.dimension() == 2)  // Delaunay facet : Voronoi edge
        {
            const Delaunay::Facet f(s);
            if (m_dela.is_infinite(f)) {
                // if it is infinite, go to its finite faces
                std::list<Simplex> f_faces = DelaunayHelper::D_faces(s);
                for (const Simplex face : f_faces){
                    if (!DelaunayHelper::is_infinite(m_dela, face)) 
                        ls.push_back(face);
                }
                return ls;
            }
            const Delaunay::Cell_handle ch = f.first;
            const int i = f.second;
            const Delaunay::Facet f_mirror = m_dela.mirror_facet(f);
            const Delaunay::Cell_handle ch_mirror = f_mirror.first;
            const int i_mirror = f_mirror.second; // not used

            const Delaunay::Point a = m_dela.point(ch->vertex( (i+1)%4 ));
            const Delaunay::Point b = m_dela.point(ch->vertex( (i+2)%4 ));
            const Delaunay::Point c = m_dela.point(ch->vertex( (i+3)%4 ));
            
            // this check if the triangle is acute (using pythagore) :
            const double acsl = (c-a).squared_length();
            const double absl = (b-a).squared_length();
            const double bcsl = (c-b).squared_length();
            const double maxsl = max3(acsl,absl,bcsl);
            bool is_acute = (2*maxsl <= acsl+absl+bcsl);
            assert(!CGAL::collinear(a,b,c));// I don't know what to do if they are colinear...
            const Vector normal = CGAL::unit_normal(a,b,c);
            const Delaunay::Point p = CGAL::circumcenter(a,b,c);

            if (m_dela.is_infinite(ch) && m_dela.is_infinite(ch_mirror)){
                std::cerr << "Error in flowcell_from_critical_cell: both incident cell_handle are infinite ";
                return ls;
            }
            else if (m_dela.is_infinite(ch)){
                const Delaunay::Point dual_ch_mirror = m_dela.dual(ch_mirror);
                const Delaunay::Point ch_mirror_opposite_pt = m_dela.point(ch_mirror->vertex( (i_mirror)%4 ));
                if(STRICT_SIGN(CGAL::scalar_product(ch_mirror_opposite_pt-p, normal)) != STRICT_SIGN(CGAL::scalar_product(dual_ch_mirror-p, normal)) ){
                  // if opposite pt and dual pt of the finite incident cell are not on the same size,
                  // add the cell
                  ls.push_back(Simplex(ch_mirror));
                }
            }
            else if (m_dela.is_infinite(ch_mirror)){
                const Delaunay::Point dual_ch = m_dela.dual(ch);
                const Delaunay::Point ch_opposite_pt = m_dela.point(ch->vertex( (i)%4 ));
                if(STRICT_SIGN(CGAL::scalar_product(ch_opposite_pt-p, normal)) != STRICT_SIGN(CGAL::scalar_product(dual_ch-p, normal)) ){
                  // if opposite pt and dual pt of the finite incident cell are not on the same size,
                  // add the cell
                  ls.push_back(Simplex(ch));
                }
            }
            else {
                const Delaunay::Point ch_opposite_pt = m_dela.point(ch->vertex( (i)%4 ));
                const int ch_sign = STRICT_SIGN(CGAL::scalar_product(ch_opposite_pt-p, normal));
                const Delaunay::Point dual_ch = m_dela.dual(ch);
                const Delaunay::Point dual_ch_mirror = m_dela.dual(ch_mirror);
                const int sign_dual_ch = STRICT_SIGN(CGAL::scalar_product(dual_ch-p, normal)) ;
                const int sign_dual_ch_mirror = STRICT_SIGN(CGAL::scalar_product(dual_ch_mirror-p, normal)) ;
                if (sign_dual_ch == sign_dual_ch_mirror){
                  // same side 
                  if (sign_dual_ch == ch_sign) {
                    // dual edge -> side of ch
                    ls.push_back(Simplex(ch_mirror));
                  }
                  else { 
                    // dual edge -> side of ch_mirror
                    ls.push_back(Simplex(ch));
                  }
                }
                else {
                  // interesecting the plane!
                  // nothing to add of dimension 3
                }
            }

            // if not acute -> critical -> go to every faces (or should be a WARNING or an error
            // if acute, go to the two smallest sides!
            if (is_acute){
                if (acsl < absl || acsl < bcsl)
                  ls.push_back( Delaunay::Edge(ch, (i+1)%4, (i+3)%4) );
                if (absl < acsl || absl < bcsl)
                  ls.push_back( Delaunay::Edge(ch, (i+1)%4, (i+2)%4) );
                if (bcsl < absl || bcsl < acsl)
                  ls.push_back( Delaunay::Edge(ch, (i+2)%4, (i+3)%4) );
                }
                else {
                std::list<Simplex> f_faces = DelaunayHelper::D_faces(s);
                for (const Simplex face : f_faces){
                   ls.push_back(face);
                }
            }
        

        }
        else if (s.dimension() == 3) // Delaunay cell : Voronoi vertex
        {
            // not so easy:
            const Delaunay::Cell_handle ch(s);
            std::list<Simplex> ch_faces = DelaunayHelper::D_faces(s);
            if (m_dela.is_infinite(ch)) {
                // if it is infinite, go to its finite faces
                for (const Simplex face : ch_faces){
                    if (!DelaunayHelper::is_infinite(m_dela, face)) 
                        ls.push_back(face);
                }
                return ls;
            }
            const Delaunay::Point dual_point = m_dela.dual(ch);
            for (size_t i = 0; i < 4; i++) {
                const Delaunay::Facet f_i = Delaunay::Facet(ch, i);
                const Delaunay::Point pt_i = m_dela.point(ch->vertex(i));
                
                const Delaunay::Point a = m_dela.point(ch->vertex( (i+1)%4 ));
                const Delaunay::Point b = m_dela.point(ch->vertex( (i+2)%4 ));
                const Delaunay::Point c = m_dela.point(ch->vertex( (i+3)%4 ));
                const Delaunay::Point f_circum = CGAL::circumcenter(a,b,c);
                assert(!CGAL::collinear(a,b,c));// I don't know what to do if they are colinear...
                Vector normal = CGAL::unit_normal(a,b,c);
                if (CGAL::scalar_product(f_circum-pt_i, normal)<0.0) {
                    normal = -normal; // we eventually reorient the normal
                }
                // at this point, normal should be the normal of the faces oriented outward from the 3D cell
                
                if (CGAL::scalar_product(f_circum-dual_point, normal)>=0.0) {
                    ls.push_back(Simplex(f_i)); // we add faces that are oriented out the dual point.
                }
                // note that this also work when ch is critical, in this case, all the facets are added
            }
        }
        else
        {
          std::clog << "Error in up_flow_cells : s.dimension = " << s.dimension() << std::endl;
        }
        return ls;
    }

    
    /** \brief check unicity and existence of a flowcell, given a simplex */
    bool check_partition_cells () {
        std::map<Delaunay::Simplex, bool> vertices, edges, faces, cells;

        size_t cpt(0);
        // Vertices
        for (Delaunay::All_vertices_iterator vit = m_dela.all_vertices_begin(); vit != m_dela.all_vertices_end(); vit++) {
            const Delaunay::Simplex s(vit);
            cpt = simplex_to_simplex_id[s].i;
            auto it(simplex_to_flowcell_id.find(s));
            if (it != simplex_to_flowcell_id.end()) { // Simplex found in a cell
                if (vertices.find(s) != vertices.end()) { // Simplex already recorded in another cell
                    std::cout << "Vertex " << cpt << " recorded in two flow cells" << std::endl;
                    return false;
                }
                else
                    vertices[s] = true;
            }
            else {
                std::cout << "Vertex " << cpt << " not in a flow cell" << std::endl;
    //                return false;
            }
        }
        // Edges
        for (Delaunay::All_edges_iterator eit = m_dela.all_edges_begin(); eit != m_dela.all_edges_end(); eit++) {
            const Delaunay::Simplex s(*eit);
            cpt = simplex_to_simplex_id[s].i;
            auto it(simplex_to_flowcell_id.find(s));
            if (it != simplex_to_flowcell_id.end()) { // Simplex found in a cell
                if (edges.find(s) != edges.end()) { // Simplex already recorded in another cell
                    std::cout << "Edge " << cpt << " recorded in two flow cells" << std::endl;
                    return false;
                }
                else
                    edges[s] = true;
            }
            else {
                std::cout << "Edge " << cpt << " not in a flow cell" << std::endl;
    //                return false;
            }
        }
        // Facets
        for (Delaunay::All_facets_iterator fit = m_dela.all_facets_begin(); fit != m_dela.all_facets_end(); fit++) {
            const Delaunay::Simplex s(*fit);
            cpt = simplex_to_simplex_id[s].i;
            auto it(simplex_to_flowcell_id.find(s));
            if (it != simplex_to_flowcell_id.end()) { // Simplex found in a cell
                if (faces.find(s) != faces.end()) { // Simplex already recorded in another cell
                    std::cout << "Facet " << cpt << " recorded in two flow cells" << std::endl;
                    return false;
                }
                else
                    faces[s] = true;
            }
            else {
                std::cout << "Facet " << cpt << " not in a flow cell" << std::endl;
    //                return false;
            }
        }
        // Cells
        for (Delaunay::All_cells_iterator cit = m_dela.all_cells_begin(); cit != m_dela.all_cells_end(); cit++) {
            const Delaunay::Simplex s(cit);
            cpt = simplex_to_simplex_id[s].i;
    //            std::cout << "cell : " << simplex_to_simplex_id[s].i << " / dim " << simplex_to_simplex_id[s].dim << std::endl;
            auto it(simplex_to_flowcell_id.find(s));
            if (it != simplex_to_flowcell_id.end()) { // Simplex found in a cell
                if (cells.find(s) != cells.end()) { // Simplex already recorded in another cell
                    std::cout << "Cell " << cpt << " recorded in two flow cells" << std::endl;
                    return false;
                }
                else
                    cells[s] = true;
            }
            else {
                std::cout << "Cell " << cpt << " not in a flow cell" << std::endl;
    //                return false;
            }
        }
    }

    int sign_inside_out(const Delaunay::Cell_handle& ch, const SoT& inside) {
        // Compute the barycenter of the simplex
        Vector sum;
        // Get points of the cell
        for (int j=0; j<4; ++j)
            sum += m_dela.point(ch->vertex(j)) - CGAL::ORIGIN;
        sum /= 4.;
        Delaunay::Point bary((CGAL::ORIGIN+sum));
        // Check the side of the barycenter
        CGAL::Bounded_side res = inside(bary);
        if (res == CGAL::ON_BOUNDED_SIDE)
            return -1;
        else
            return 1;
    }

    /**
     * @brief FlowComplex::is_on_boundary
     * Check if the facet is on the boundary of the object.
     * @Precondition : m_cell_filtration should have been filled.
     */
    bool is_on_boundary(size_t i, const SoT& inside) {
        Delaunay::Simplex s(flowcells.at(i).get_critical_simplex());
        if (s.dimension() == 2) {
            Delaunay::Facet f(s);

            Delaunay::Cell_handle c1 = f.first;
            Delaunay::Cell_handle c2 = m_dela.mirror_facet(f).first;
            bool inf1=m_dela.is_infinite(c1); bool inf2=m_dela.is_infinite(c2);
            if (inf1 && inf2)
            {
                return false;
            }
            else if (inf1 || inf2) // case where both cells are infinite
            {
                if (inf1)
                    return sign_inside_out(c2,inside) == -1;
                else // inf2
                    return sign_inside_out(c1,inside) == -1;
            }
            return sign_inside_out(c1, inside) != sign_inside_out(c2,inside);
        }
        else
            return false;
    }
    
    /**
    * @brief print the flowcells of the flow complex.
    * It displays all the simplices in each flowcell.
    */
    std::ostream& print_flowcells(std::ostream& os) const
    {
        std::cout << "=== flow cells" << std::endl;
        for (int i=0; i<flowcells.size(); ++i) {
            Delaunay::Simplex cs(flowcells.at(i).get_critical_simplex());
            
            if (auto search = simplex_to_simplex_id.find(cs); search != simplex_to_simplex_id.end()) {
                if (DelaunayHelper::is_infinite(m_dela,cs))
                    std::cout << "- cell " << i << ": critical " << "inf(" << search->second.i << ", " << search->second.dim << ") -   ";
                else 
                    std::cout << "- cell " << i << ": critical " << "(" << search->second.i << ", " << search->second.dim << ") -   ";
            }
            else {
                std::cerr << "Error: no simplex ";
                print_simplex(std::cerr, cs) << "in simplex_to_simplex_id" << "\n";
            }
            
            for (Delaunay::Simplex s : flowcells.at(i).get_simplices()) {
                
                if (auto search = simplex_to_simplex_id.find(s); search != simplex_to_simplex_id.end()) {
                    if (DelaunayHelper::is_infinite(m_dela,s))
                        std::cout << "inf(" << search->second.i << ", " << search->second.dim << ") -   ";
                    else 
                        std::cout << "(" << search->second.i << ", " << search->second.dim << ") -   ";
                }
                else {
                    std::cerr << "Error: no finite simplex ";
                    print_simplex(std::cerr, s) << "in simplex_to_simplex_id" << "\n";
                }
            }
            std::cout << std::endl;
        }
        return os;
    }
    
    /**
    * @brief print the poset related to the flow complex.
    * It displays the faces of each flow cell, with their ids.
    */
    std::ostream& print_poset(std::ostream& os) const
    {
        std::cout << "=== poset" << std::endl;
        size_t i = 0;
        for (const auto& set : flowcell_faces){
            os << i  << " (dim: " << flowcells.at(i).get_critical_simplex().dimension() << " / n: " << flowcells.at(i).get_simplices().size() << ") " << "-> { " << std::endl;
            for (const auto& j : set){
                os << "\t" << j << " (dim: " << flowcells.at(j).get_critical_simplex().dimension() << " / n: "<< flowcells.at(j).get_simplices().size() << ") " << std::endl;
            }
            os << "}\n";
            i++;
        }
        return os;
    }

    ////////////////////////// EXPORT //////////////////////////////////////////

    std::ostream& write_sub(std::ostream& out, const FlowCell& fc) {
        for (Delaunay::Simplex s : fc.get_simplices()) {
            Simplex_id sid(simplex_to_simplex_id[s]);
            out << sid.dim << " " << sid.i << std::endl;
        }
        return out;
    }

    void write_sub(const FlowCell& fc, std::string filename, size_t thres_number = 2) {
        if (fc.get_simplices().size() >= thres_number) {
            std::ofstream out ( filename, std::ios::out | std::ios::trunc);

            if ( ! out . good () ) {
                std::cerr << "write_nodes for Delaunay_3. Fatal Error:\n  " << filename << " not found.\n";
                throw std::runtime_error("File Parsing Error: File not found");
            }

            write_sub(out, fc);
            out.close();
        }
    }

    void write_criticals_sub(std::string filename) {
        std::ofstream out ( filename, std::ios::out | std::ios::trunc);

        if ( ! out . good () ) {
            std::cerr << "write_criticals_sub. Fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }

        for (FlowCell fc : flowcells) {
            Delaunay::Simplex s(fc.get_critical_simplex());
            Simplex_id sid(simplex_to_simplex_id[s]);
            out << sid.dim << " " << sid.i << std::endl;
        }

        out.close();
    }

    // There should be a way of exporting flowcells in VTK format? At least the big ones, to see if it is coherent?
    // Maybe export the delaunay (finite), and add a label corresponding to the id of its flowcell, see simplex_to_flowcell_id

    size_t number_of_flow_cells() const { return flowcells.size(); }

    // write_vtk: write ids for cells and export the Delaunay simplicial complex with these flags (flag: id of the cell, flag2: 0 for the crical cell and 2 for others

    virtual void write_vtk(std::string filename) {
        std::vector<std::vector<int> > simplex_ids(4), flow_ids(4), flow_criticals(4); // set to flow_id or -flow_id (for criticals)
        std::vector<std::vector<double> > flow_df(4);
        const size_t end_cells = number_of_flow_cells()+1;
        simplex_ids.at(0).resize(m_dela.number_of_vertices(),-1);
        simplex_ids.at(1).resize(m_dela.number_of_finite_edges(),-1);
        simplex_ids.at(2).resize(m_dela.number_of_finite_facets(),-1);
        simplex_ids.at(3).resize(m_dela.number_of_finite_cells(),-1);
        flow_ids.at(0).resize(m_dela.number_of_vertices(),end_cells);
        flow_ids.at(1).resize(m_dela.number_of_finite_edges(),end_cells);
        flow_ids.at(2).resize(m_dela.number_of_finite_facets(),end_cells);
        flow_ids.at(3).resize(m_dela.number_of_finite_cells(),end_cells);
        flow_criticals.at(0).resize(m_dela.number_of_vertices(),end_cells);
        flow_criticals.at(1).resize(m_dela.number_of_finite_edges(),end_cells);
        flow_criticals.at(2).resize(m_dela.number_of_finite_facets(),end_cells);
        flow_criticals.at(3).resize(m_dela.number_of_finite_cells(),end_cells);
        flow_df.at(0).resize(m_dela.number_of_vertices(),0.);
        flow_df.at(1).resize(m_dela.number_of_finite_edges(),0.);
        flow_df.at(2).resize(m_dela.number_of_finite_facets(),0.);
        flow_df.at(3).resize(m_dela.number_of_finite_cells(),0.);

        // Build flow_ids

        // Vertices
        size_t cpt(0);
        for (typename Delaunay::Finite_vertices_iterator it = m_dela.finite_vertices_begin(); !(it == m_dela.finite_vertices_end()); ++it) {
            Delaunay::Simplex s(it);
            simplex_ids.at(0).at(cpt) = simplex_to_simplex_id[s].i;
            auto search(simplex_to_flowcell_id.find(s));
            if (search != simplex_to_flowcell_id.end()) {
                size_t id_cell = simplex_to_flowcell_id[s];
                flow_ids.at(0).at(cpt) = id_cell;
                flow_df.at(0).at(cpt) = flowcells.at(id_cell).get_df();
                if (s == flowcell_from_id(id_cell).get_critical_simplex())
                    flow_criticals.at(0).at(cpt++) = 0;
                else
                    flow_criticals.at(0).at(cpt++) = 2;
            }
            else {
                flow_ids.at(0).at(cpt)= -1;//simplex_to_simplex_id[s].i;
                flow_df.at(0).at(cpt) = 10000;
                flow_criticals.at(0).at(cpt++) = -1;
            }
        }

        // Edges
        cpt = 0;
        for (typename Delaunay::Edge edge : m_dela.finite_edges()) {
            Delaunay::Simplex s(edge);
            simplex_ids.at(1).at(cpt) = simplex_to_simplex_id[s].i;
            auto search(simplex_to_flowcell_id.find(s));
            if (search != simplex_to_flowcell_id.end()) {
                size_t id_cell = simplex_to_flowcell_id[s];
                flow_ids.at(1).at(cpt) = id_cell;
                flow_df.at(1).at(cpt) = flowcells.at(id_cell).get_df();
                if (s == flowcell_from_id(id_cell).get_critical_simplex())
                    flow_criticals.at(1).at(cpt++) = 0;
                else
                    flow_criticals.at(1).at(cpt++) = 2;
            }
            else {
                flow_ids.at(1).at(cpt)= -1;//simplex_to_simplex_id[s].i;
                flow_df.at(1).at(cpt) = 10000;
                flow_criticals.at(1).at(cpt++) = -1;
            }
        }

        // Triangles
        cpt = 0;
        for (typename Delaunay::Facet facet : m_dela.finite_facets()) {
            Delaunay::Simplex s(facet);
            simplex_ids.at(2).at(cpt) = simplex_to_simplex_id[s].i;
            auto search(simplex_to_flowcell_id.find(s));
            if (search != simplex_to_flowcell_id.end()) {
                size_t id_cell = simplex_to_flowcell_id[s];
                flow_ids.at(2).at(cpt) = id_cell;
                flow_df.at(2).at(cpt) = flowcells.at(id_cell).get_df();
                if (s == flowcell_from_id(id_cell).get_critical_simplex())
                    flow_criticals.at(2).at(cpt++) = 0;
                else
                    flow_criticals.at(2).at(cpt++) = 2;
            }
            else {
                flow_ids.at(2).at(cpt)= -1;//simplex_to_simplex_id[s].i;
                flow_df.at(2).at(cpt) = 10000;
                flow_criticals.at(2).at(cpt++) = -1;
            }
        }

        // Cells
        cpt = 0;
        for (typename Delaunay::Cell_handle cell : m_dela.finite_cell_handles()) {
            Delaunay::Simplex s(cell);
            simplex_ids.at(3).at(cpt) = simplex_to_simplex_id[s].i;
            auto search(simplex_to_flowcell_id.find(s));
            if (search != simplex_to_flowcell_id.end()) {
                size_t id_cell = simplex_to_flowcell_id[s];
                flow_ids.at(3).at(cpt) = id_cell;
                flow_df.at(3).at(cpt) = flowcells.at(id_cell).get_df();
                if (s == flowcell_from_id(id_cell).get_critical_simplex())
                    flow_criticals.at(3).at(cpt++) = 0;
                else
                    flow_criticals.at(3).at(cpt++) = 2;
            }
            else {
                flow_ids.at(3).at(cpt)= -1;//simplex_to_simplex_id[s].i;
                flow_df.at(3).at(cpt) = 10000;
                flow_criticals.at(3).at(cpt++) = -1;
            }
        }

        // Export the Delaunay mesh with these flags
        write_VTK(m_dela, "tmp/criticals.vtk", simplex_ids, CELLS|FACETS|EDGES|VERTICES, &flow_criticals, &flow_df);
        write_VTK(m_dela, filename, simplex_ids, CELLS|FACETS|EDGES|VERTICES, &flow_ids, &flow_df);
    }

 };

#endif // FLOW_COMPLEX_H
