#ifndef FLOW_COMPLEX_H
#define FLOW_COMPLEX_H


#include "cgal_typedef.h"
// #include "filtration.h"
#include "delaunay_helper.h"

#include <fstream>
#include <limits>

#define STRICT_SIGN(x) (((x)<0.0) ? -1 : 1)
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
    
    friend std::ostream& operator<<(std::ostream& os, const FlowCell& fc);
};

inline std::ostream& operator<<(std::ostream& os, const FlowCell& fc)
{
    os << "FlowCell {\n"
       << "  dimension: " << fc.dimension() << '\n'
       << "  df: " << fc.df << '\n'
       << "  critical simplex: ";
    print_simplex(os, fc.crit_simplex) << '\n';

    os << "  simplices (" << fc.del_simplices.size() << "):\n";
    for (const auto& s : fc.del_simplices)
        print_simplex(os, s) << ' ';
    os << "}";
    return os;
}
/*----------------------------------------------------------------------------*/

/* geometry */
Vector normal_facet_opposite_to_p(const Delaunay &m_dela, const Delaunay::Facet f, const Delaunay::Point p)
{
    const Delaunay::Cell_handle ch = f.first;
    const int i = f.second;
    
    const Delaunay::Point a = m_dela.point(ch->vertex( (i+1)%4 ));
    const Delaunay::Point b = m_dela.point(ch->vertex( (i+2)%4 ));
    const Delaunay::Point c = m_dela.point(ch->vertex( (i+3)%4 ));
    
    // this check if the triangle is acute (using pythagore) :
    Vector normal = CGAL::normal(a,b,c);
    if (CGAL::scalar_product(a-p, normal)<0.0) {
        return -normal;
    }
    return normal;
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
        //CGAL::Side_of_triangle_mesh<Polyhedron, Epick> inside(m_poly);
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
    
    std::list<Simplex> up_flow_cells(const Simplex s)
    {
      std::list<Simplex> ls = {};
      std::clog << "    - up_flow : s.dimension = " << s.dimension() << "   ";
      print_simplex(std::clog, s) << '\n';
      if (s.dimension() == 0) // Delaunay vertex : Voronoi cell
      {
          // pass
      }
      else if (s.dimension() == 1) // Delaunay edge : Voronoi facet
      {
          // HARD one:
          // TODO
          ls = DelaunayHelper::D_faces(s);
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
          const Vector normal = CGAL::normal(a,b,c);
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
                Vector normal = CGAL::normal(a,b,c);
                if (CGAL::scalar_product(f_circum-pt_i, normal)<0.0) {
                    normal = -normal; // we eventually reorient the normal
                }
                // at this point, normal should be the normal of the faces oriented outward from the 3D cell
                
                if (CGAL::scalar_product(f_circum-dual_point, normal)>0.0) {
                    ls.push_back(Simplex(f_i)); // we add faces that are oriented out the dual point.
                }
                // note that this also work when ch is critical, in this case, we should add all the facets
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
        // TODO
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
                fc_faces.insert(search->second); // update the poset
            }
            else {
                simplex_to_flowcell_id[s] = id;
                fc.add_simplex(s);
                std::list<Simplex> ls = up_flow_cells(s);
                for (const auto& s_up_flow : ls) {
                    //print_simplex(std::clog, s_up_flow) << std::endl;
                    to_process.push(s_up_flow);
                }
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
    
private:
    std::vector<FlowCell> flowcells; // probably
    std::map<Simplex, size_t> simplex_to_flowcell_id;
    std::vector<std::set<size_t>> flowcell_faces; // encoding the poset
    /* WARNING I think there might be problems with the algorithm, for instance if a flowcell is twice adjacent to another cell...
    * In this case, it will be counted only once in flowcell_faces...
    */
 };

#endif // FLOW_COMPLEX_H
