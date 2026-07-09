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

std::ostream& print_simplex(std::ostream& os,
                            const Delaunay::Simplex& s)
{
    switch (s.dimension()) {
    case 0: {
        auto vh = Delaunay::Vertex_handle(s);
        os << "Vertex(" << vh->point() << ")";
        break;
    }
    case 1: {
        auto e = Delaunay::Edge(s);
        os << "Edge("
           << e.first->vertex(e.second)->point() << ", "
           << e.first->vertex(e.third)->point() << ")";
        break;
    }
    case 2: {
        auto f = Delaunay::Facet(s);
        os << "Facet(";
        Delaunay::Cell_handle ch = f.first;
        int i = f.second;
        os << ch->vertex( (i+1)%4 )->point() << ", ";
        os << ch->vertex( (i+2)%4 )->point() << ", ";
        os << ch->vertex( (i+3)%4 )->point() << ")";
        break;
    }
    case 3: {
        auto c = Delaunay::Cell_handle(s);
        os << "Cell(";
        os << c->vertex(0)->point() << ", ";
        os << c->vertex(1)->point() << ", ";
        os << c->vertex(2)->point() << ", ";
        os << c->vertex(3)->point() << ")";
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
* @brief The FlowBuilder class
* 
*/
class FlowBuilder
{
public:
    typedef Delaunay::Simplex Simplex;
    FlowBuilder(Polyhedron poly) : m_poly(poly)
    {
        //CGAL::draw(m_poly);
        // CGAL::Side_of_triangle_mesh<Polyhedron, Epick> inside(m_poly);
        for (Polyhedron::Vertex_iterator it = m_poly.vertices_begin(); it != m_poly.vertices_end(); ++it)
        {
            m_dela.insert(it->point()); // m_dela corresponds to the 3D Delaunay of the surface points
    //        std::cout << it->point() << std::endl;
        }
        std::clog << "Computing flow complex. Delaunay triangulation:" << std::endl;
        std::clog << "nb finite vertices: " << m_dela.number_of_vertices() << std::endl;
        std::clog << "nb edges   : " << m_dela.number_of_edges()  << " | finite: " << m_dela.number_of_finite_edges() << std::endl;
        std::clog << "nb facets  : " << m_dela.number_of_facets() << " | finite: " << m_dela.number_of_finite_facets()<< std::endl;
        std::clog << "nb cells   : " << m_dela.number_of_cells()  << " | finite: " << m_dela.number_of_finite_cells()<< std::endl;
        // CGAL::draw(m_dela);
    }
    
    // Member data
    Polyhedron m_poly;
    Delaunay m_dela;
    
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
            const Delaunay::Cell_handle ch = f.first;
            const int i = f.second;
            const Delaunay::Facet f_mirror = m_dela.mirror_facet(f);
            const Delaunay::Cell_handle ch_mirror = f_mirror.first;
            const int i_mirror = f_mirror.second; // not used

            const Delaunay::Point a = m_dela.point(ch->vertex( (i+1)%4 ));
            const Delaunay::Point b = m_dela.point(ch->vertex( (i+2)%4 ));
            const Delaunay::Point c = m_dela.point(ch->vertex( (i+3)%4 ));
            const Delaunay::Point ch_opposite_pt = m_dela.point(ch->vertex( (i)%4 ));

            // this check if the triangle is acute (using pythagore) :
            const double acsl = (c-a).squared_length();
            const double absl = (b-a).squared_length();
            const double bcsl = (c-b).squared_length();
            const double maxsl = max3(acsl,absl,bcsl);
            bool is_acute = (2*maxsl <= acsl+absl+bcsl);
            assert(!CGAL::collinear(a,b,c));// I don't know what to do if they are colinear...
            const Vector normal = CGAL::unit_normal(a,b,c);
            const Delaunay::Point p = CGAL::circumcenter(a,b,c);
            const int ch_sign = STRICT_SIGN(CGAL::scalar_product(ch_opposite_pt-p, normal));

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
                if(STRICT_SIGN(CGAL::scalar_product(ch_opposite_pt-p, normal)) != STRICT_SIGN(CGAL::scalar_product(dual_ch-p, normal)) ){
                  // if opposite pt and dual pt of the finite incident cell are not on the same size,
                  // add the cell
                  ls.push_back(Simplex(ch));
                }
            }
            else {
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
                    if (DelaunayHelper::is_infinite(m_dela, face)) 
                        ls.push_back(face);
                }
            }
            else {
                const Delaunay::Point dual_point = m_dela.dual(ch);
                for (size_t i = 0; i < 4; i++) {
                    const Delaunay::Facet f_i = Delaunay::Facet(ch, i);
                    const Delaunay::Point pt_i = m_dela.point(ch->vertex(i));
                    
                    const Delaunay::Point a = m_dela.point(ch->vertex( (i+1)%4 ));
                    const Delaunay::Point b = m_dela.point(ch->vertex( (i+2)%4 ));
                    const Delaunay::Point c = m_dela.point(ch->vertex( (i+3)%4 ));
                    const Delaunay::Point f_circum = CGAL::circumcenter(a,b,c);
                    assert(!CGAL::collinear(a,b,c));// I don't know what to do if they are colinear...
                    Vector normal = CGAL::normal(a,b,c);
                    if (CGAL::scalar_product(f_circum-pt_i, normal)<0.0) {
                        normal = -normal; // we eventually reorient the normal
                    }
                    // at this point, normal should be the normal of the faces oriented outward from the 3D cell
                    
                    if (CGAL::scalar_product(f_circum-dual_point, normal)>0.0) {
                        ls.push_back(Simplex(f_i)); // we add faces that are oriented out the dual point.
                    }
                    // note that this also work when ch is critical, in this case, all the facets are added
                }
            }
        }
        else
        {
          std::clog << "Error in up_flow_cells : s.dimension = " << s.dimension() << std::endl;
        }
        return ls;
    }
    
    FlowCell flowcell_from_critical_cell(const Simplex crit_simplex)
    {
        FlowCell fc(crit_simplex);
        // print_simplex(std::clog, crit_simplex) << std::endl;
        size_t id = flowcells.size();
        flowcells.push_back(fc);
        
        CriticalInfo crit_info = DelaunayHelper::get_critical_info(m_dela, crit_simplex);
        if (crit_info.c != CriticalType::Critical)
        {
            std::cerr << "Error in flowcell_from_critical_cell: simplex ";
            print_simplex(std::clog, crit_simplex) << " not critical" << std::endl;
        }
        fc.set_df(crit_info.r);
        simplex_to_flowcell_id[crit_simplex] = id;
        std::set<size_t> fc_faces = {};
         
        std::queue<Simplex> to_process;
        for (const Simplex& s : DelaunayHelper::D_faces(crit_simplex)) {
          to_process.push(s);
        }
        
        while (!to_process.empty()) {
            const Simplex s = to_process.front();
            to_process.pop();

            if (auto search = simplex_to_flowcell_id.find(s); search != simplex_to_flowcell_id.end()) {
                if (search->second != id)
                    fc_faces.insert(search->second); // update the poset
            }
            else {
                simplex_to_flowcell_id[s] = id;
                fc.add_simplex(s);
                std::list<Simplex> ls = up_flow_cells(s);
                
                std::clog << "    - up_flow : s.dimension = " << s.dimension() << " ";
                print_simplex(std::clog, s) << "\t[ ";
                
                for (const auto& s_up_flow : ls) {
                    //print_simplex(std::clog, s_up_flow) << std::endl;
                    std::clog << s_up_flow.dimension()<< " ";
                    if (s_up_flow.dimension() <= crit_simplex.dimension())
                        to_process.push(s_up_flow);
                    else
                        std::clog << "(not added) "; // It seems that this case can happen in 3D: using our algo, if we start from a critical 2D cell we can arrive on a 3D cell... weird.
                }
                std::clog << "]\n";
            }
        }
        flowcell_faces.push_back(fc_faces);
        return fc;
    }
    
    size_t get_flowcell_id(const FlowCell& f) const
    {
      if (auto search = simplex_to_flowcell_id.find(f.get_critical_simplex()); search != simplex_to_flowcell_id.end()) {
          return search->second;
      }
      std::cerr << "Error in get_flowcell_id: crit_simplex not in the map crit_simplex_to_flowcell_id." << std::endl;
      return 0;
    }
    
    std::ostream& print_poset(std::ostream& os) const
    {
        size_t i = 0;
        for (const auto& set : flowcell_faces){
            os << i <<"->\t{ ";
            for (const auto& j : set){
                os << j << " ";
            }
            os << "}\n";
            i++;
        }
        return os;
    }
    
    // There should be a way of exporting flowcells in VTK format? At least the big ones, to see if it is coherent?
    // Maybe export the delaunay (finite), and add a label corresponding to the id of its flowcell, see simplex_to_flowcell_id
    
private:
    std::vector<FlowCell> flowcells; // probably
    std::map<Simplex, size_t> simplex_to_flowcell_id;
    std::vector<std::set<size_t>> flowcell_faces; // encoding the poset
    /* WARNING I think there might be problems with the algorithm, for instance if a flowcell is twice adjacent to another cell...
    * In this case, it will be counted only once in flowcell_faces...
    */
 };

#endif // FLOW_COMPLEX_H
