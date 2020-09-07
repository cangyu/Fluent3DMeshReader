#include "xf.h"

/// Convert a boundary condition string literal to unified form within the scope of this code.
/// Outcome will be composed of LOWER case letters and '-' only!
static std::string formalize(const std::string &s)
{
    std::string ret(s);
    std::transform(ret.begin(), ret.end(), ret.begin(), ::tolower);
    for (auto &e : ret)
        if (e == '_')
            e = '-';
    return ret;
}

static void eat(std::istream &in, char c)
{
    char tmp = 0;
    do {
        in >> tmp;
    } while (tmp != c);
}

static void skip_white(std::istream &in)
{
    char tmp = 0;
    do {
        in >> tmp;
    } while (tmp == ' ' || tmp == '\t' || tmp == '\n');

    if (!in.eof())
        in.unget();
}

namespace XF
{
    bool BC::isValidIdx(int x)
    {
        static const std::set<int> candidate_set{
            INTERIOR,
            WALL,
            PRESSURE_INLET,
            PRESSURE_OUTLET,
            SYMMETRY,
            PERIODIC_SHADOW,
            PRESSURE_FAR_FIELD,
            VELOCITY_INLET,
            PERIODIC,
            FAN,
            MASS_FLOW_INLET,
            INTERFACE,
            PARENT,
            OUTFLOW,
            AXIS
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    bool BC::isValidStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "interior",
            "wall",
            "pressure-inlet", "inlet-vent", "intake-fan",
            "pressure-outlet", "exhaust-fan", "outlet-vent",
            "symmetry",
            "periodic-shadow",
            "pressure-far-field",
            "velocity-inlet",
            "periodic",
            "fan", "porous-jump", "radiator",
            "mass-flow-inlet",
            "interface",
            "parent",
            "outflow",
            "axis"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &BC::idx2str(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(INTERIOR, "interior"),
            std::pair<int, std::string>(WALL, "wall"),
            std::pair<int, std::string>(PRESSURE_INLET, "pressure-inlet"),
            std::pair<int, std::string>(PRESSURE_OUTLET, "pressure-outlet"),
            std::pair<int, std::string>(SYMMETRY, "symmetry"),
            std::pair<int, std::string>(PERIODIC_SHADOW, "periodic-shadow"),
            std::pair<int, std::string>(PRESSURE_FAR_FIELD, "pressure-far-field"),
            std::pair<int, std::string>(VELOCITY_INLET, "velocity-inlet"),
            std::pair<int, std::string>(PERIODIC, "periodic"),
            std::pair<int, std::string>(FAN, "fan"),
            std::pair<int, std::string>(MASS_FLOW_INLET, "mass-flow-inlet"),
            std::pair<int, std::string>(INTERFACE, "interface"),
            std::pair<int, std::string>(PARENT, "parent"),
            std::pair<int, std::string>(OUTFLOW, "outflow"),
            std::pair<int, std::string>(AXIS, "axis")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_bc_idx(x);
        else
            return it->second;
    }

    int BC::str2idx(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("interior", INTERIOR),
            std::pair<std::string, int>("wall", WALL),
            std::pair<std::string, int>("pressure-inlet", PRESSURE_INLET),
            std::pair<std::string, int>("inlet-vent", INLET_VENT),
            std::pair<std::string, int>("intake-fan", INTAKE_FAN),
            std::pair<std::string, int>("pressure-outlet", PRESSURE_OUTLET),
            std::pair<std::string, int>("exhaust-fan", EXHAUST_FAN),
            std::pair<std::string, int>("outlet-vent", OUTLET_VENT),
            std::pair<std::string, int>("symmetry", SYMMETRY),
            std::pair<std::string, int>("periodic-shadow", PERIODIC_SHADOW),
            std::pair<std::string, int>("pressure-far-field", PRESSURE_FAR_FIELD),
            std::pair<std::string, int>("velocity-inlet", VELOCITY_INLET),
            std::pair<std::string, int>("periodic", PERIODIC),
            std::pair<std::string, int>("fan", FAN),
            std::pair<std::string, int>("porous-jump", POROUS_JUMP),
            std::pair<std::string, int>("radiator", RADIATOR),
            std::pair<std::string, int>("mass-flow-inlet", MASS_FLOW_INLET),
            std::pair<std::string, int>("interface", INTERFACE),
            std::pair<std::string, int>("parent", PARENT),
            std::pair<std::string, int>("outflow", OUTFLOW),
            std::pair<std::string, int>("axis", AXIS)
        };

        auto it = mapping_set.find(formalize(x));
        if (it == mapping_set.end())
            throw invalid_bc_str(x);
        else
            return it->second;
    }

    RANGE::RANGE(int id, size_t zone, size_t first, size_t last) :
        SECTION(id),
        m_zone(zone),
        m_first(first),
        m_last(last)
    {
        if (first > last)
            throw std::invalid_argument("Invalid range in constructor.");
    }

    RANGE::RANGE(const RANGE &rhs) :
        SECTION(rhs.identity()),
        m_zone(rhs.zone()),
        m_first(rhs.first_index()),
        m_last(rhs.last_index())
    {
        if (first_index() > last_index())
            throw std::invalid_argument("Invalid range in copy-constructor.");
    }

    bool NODE::isValidTypeIdx(int x)
    {
        return x == VIRTUAL || x == ANY || x == BOUNDARY;
    }

    bool NODE::isValidTypeStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "virtual",
            "any",
            "boundary"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &NODE::idx2str(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(VIRTUAL, "virtual"),
            std::pair<int, std::string>(ANY, "any"),
            std::pair<int, std::string>(BOUNDARY, "boundary")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_node_type_idx(x);
        else
            return it->second;
    }

    int NODE::str2idx(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("virtual", VIRTUAL),
            std::pair<std::string, int>("any", ANY),
            std::pair<std::string, int>("boundary", BOUNDARY),
        };

        auto it = mapping_set.find(formalize(x));
        if (it == mapping_set.end())
            throw invalid_node_type_str(x);
        else
            return it->second;
    }

    void NODE::repr(std::ostream &out)
    {
        out << "(" << std::dec << identity();
        out << " (" << std::hex << zone() << " " << first_index() << " " << last_index() << " ";
        out << std::dec << type() << " " << ND() << ")(" << std::endl;

        out.precision(12);
        const size_t N = num();
        for (size_t i = 0; i < N; ++i)
        {
            const auto &node = at(i);
            for (int k = 0; k < m_dim; ++k)
                out << " " << node.at(k);
            out << std::endl;
        }
        out << "))" << std::endl;
    }

    bool CELL::isValidTypeIdx(int x)
    {
        static const std::set<int> candidate_set{
            DEAD,
            FLUID,
            SOLID
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    bool CELL::isValidTypeStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "dead",
            "fluid",
            "solid"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &CELL::idx2str_type(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(FLUID, "fluid"),
            std::pair<int, std::string>(SOLID, "solid"),
            std::pair<int, std::string>(DEAD, "dead")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_cell_type_idx(x);
        else
            return it->second;
    }

    int CELL::str2idx_type(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("fluid", FLUID),
            std::pair<std::string, int>("solid", SOLID),
            std::pair<std::string, int>("dead", DEAD),
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_cell_type_str(x);
        else
            return it->second;
    }

    bool CELL::isValidElemIdx(int x)
    {
        static const std::set<int> candidate_set{
            MIXED,
            TRIANGULAR,
            TETRAHEDRAL,
            QUADRILATERAL,
            HEXAHEDRAL,
            PYRAMID,
            WEDGE,
            POLYHEDRAL
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    bool CELL::isValidElemStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "mixed",
            "triangular",
            "tetrahedral",
            "quadrilateral",
            "hexahedral",
            "pyramid",
            "wedge",
            "prism",
            "polyhedral"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &CELL::idx2str_elem(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(MIXED, "mixed"),
            std::pair<int, std::string>(TRIANGULAR, "triangular"),
            std::pair<int, std::string>(TETRAHEDRAL, "tetrahedral"),
            std::pair<int, std::string>(QUADRILATERAL, "quadrilateral"),
            std::pair<int, std::string>(HEXAHEDRAL, "hexahedral"),
            std::pair<int, std::string>(PYRAMID, "pyramid"),
            std::pair<int, std::string>(WEDGE, "wedge"),
            std::pair<int, std::string>(POLYHEDRAL, "polyhedral")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_elem_type_idx(x);
        else
            return it->second;
    }

    int CELL::str2idx_elem(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("mixed", MIXED),
            std::pair<std::string, int>("triangular", TRIANGULAR),
            std::pair<std::string, int>("tri", TRIANGULAR),
            std::pair<std::string, int>("tetrahedral", TETRAHEDRAL),
            std::pair<std::string, int>("tet", TETRAHEDRAL),
            std::pair<std::string, int>("quadrilateral", QUADRILATERAL),
            std::pair<std::string, int>("quad", QUADRILATERAL),
            std::pair<std::string, int>("hexahedral", HEXAHEDRAL),
            std::pair<std::string, int>("hex", HEXAHEDRAL),
            std::pair<std::string, int>("pyramid", PYRAMID),
            std::pair<std::string, int>("wedge", WEDGE),
            std::pair<std::string, int>("prism", WEDGE),
            std::pair<std::string, int>("polyhedral", POLYHEDRAL)
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_elem_type_str(x);
        else
            return it->second;
    }

    void CELL::repr(std::ostream &out)
    {
        static const size_t NumPerLine = 40;

        out << "(" << std::dec << identity() << " (";
        out << std::hex;
        out << zone() << " " << first_index() << " " << last_index() << " ";
        out << m_type << " " << m_elem << ")";

        if (m_elem != CELL::MIXED)
            out << ")" << std::endl;
        else
        {
            out << "(";
            const size_t N = num();
            for (size_t i = 0; i < N; ++i)
            {
                if (i % NumPerLine == 0)
                    out << std::endl;
                out << " " << at(i);
            }
            out << std::endl << "))" << std::endl;
        }
    }

    bool FACE::isValidIdx(int x)
    {
        static const std::set<int> candidate_set{
            MIXED,
            LINEAR,
            TRIANGULAR,
            QUADRILATERAL,
            POLYGONAL
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    bool FACE::isValidStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "mixed",
            "linear",
            "triangular",
            "quadrilateral",
            "polygonal"
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    const std::string &FACE::idx2str(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(MIXED, "mixed"),
            std::pair<int, std::string>(LINEAR, "linear"),
            std::pair<int, std::string>(TRIANGULAR, "triangular"),
            std::pair<int, std::string>(QUADRILATERAL, "quadrilateral"),
            std::pair<int, std::string>(POLYGONAL, "polygonal")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_face_type_idx(x);
        else
            return it->second;
    }

    int FACE::str2idx(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("mixed", MIXED),
            std::pair<std::string, int>("linear", LINEAR),
            std::pair<std::string, int>("line", LINEAR),
            std::pair<std::string, int>("triangular", TRIANGULAR),
            std::pair<std::string, int>("tri", TRIANGULAR),
            std::pair<std::string, int>("quadrilateral", QUADRILATERAL),
            std::pair<std::string, int>("quad", QUADRILATERAL),
            std::pair<std::string, int>("polygonal", POLYGONAL)
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_face_type_str(x);
        else
            return it->second;
    }

    void FACE::repr(std::ostream &out)
    {
        out << "(" << std::dec << identity() << " (";
        out << std::hex;
        out << zone() << " " << first_index() << " " << last_index() << " ";
        out << bc_type() << " " << face_type() << ")(" << std::endl;

        const size_t N = num();
        if (m_face == MIXED)
        {
            for (size_t i = 0; i < N; ++i)
            {
                const auto &loc_cnect = at(i);
                out << " " << loc_cnect.n.size();
                for (auto e : loc_cnect.n)
                    out << " " << e;
                out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
            }
        }
        else
        {
            for (size_t i = 0; i < N; ++i)
            {
                const auto &loc_cnect = at(i);
                for (auto e : loc_cnect.n)
                    out << " " << e;
                out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
            }
        }

        out << "))" << std::endl;
    }

    bool ZONE::isValidIdx(int x)
    {
        const bool ret = ZONE::DEGASSING <= x && x <= ZONE::WRAPPER;
        return ret;
    }

    bool ZONE::isValidStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "degassing",
            "exhaust-fan",
            "fan",
            "fluid",
            "geometry",
            "inlet-vent",
            "intake-fan",
            "interface",
            "interior",
            "internal",
            "mass-flow-inlet",
            "outflow",
            "outlet-vent",
            "parent-face",
            "porous-jump",
            "pressure-far-field",
            "pressure-inlet",
            "pressure-outlet",
            "radiator",
            "solid",
            "symmetry",
            "velocity-inlet",
            "wall",
            "wrapper"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &ZONE::idx2str(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(DEGASSING, "degassing"),
            std::pair<int, std::string>(EXHAUST_FAN, "exhaust-fan"),
            std::pair<int, std::string>(FAN, "fan"),
            std::pair<int, std::string>(FLUID, "fluid"),
            std::pair<int, std::string>(GEOMETRY, "geometry"),
            std::pair<int, std::string>(INLET_VENT, "inlet_vent"),
            std::pair<int, std::string>(INTAKE_FAN, "intake_fan"),
            std::pair<int, std::string>(INTERFACE, "interface"),
            std::pair<int, std::string>(INTERIOR, "interior"),
            std::pair<int, std::string>(INTERNAL, "internal"),
            std::pair<int, std::string>(MASS_FLOW_INLET, "mass_flow_inlet"),
            std::pair<int, std::string>(OUTFLOW, "outflow"),
            std::pair<int, std::string>(OUTLET_VENT, "outlet_vent"),
            std::pair<int, std::string>(PARENT_FACE, "parent_face"),
            std::pair<int, std::string>(POROUS_JUMP, "porous_jump"),
            std::pair<int, std::string>(PRESSURE_FAR_FIELD, "pressure_far_field"),
            std::pair<int, std::string>(PRESSURE_INLET, "pressure_inlet"),
            std::pair<int, std::string>(PRESSURE_OUTLET, "pressure_outlet"),
            std::pair<int, std::string>(RADIATOR, "radiator"),
            std::pair<int, std::string>(SOLID, "solid"),
            std::pair<int, std::string>(SYMMETRY, "symmetry"),
            std::pair<int, std::string>(VELOCITY_INLET, "velocity_inlet"),
            std::pair<int, std::string>(WALL, "wall"),
            std::pair<int, std::string>(WRAPPER, "wrapper"),
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_zone_type_idx(x);
        else
            return it->second;
    }

    int ZONE::str2idx(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("degassing", DEGASSING),
            std::pair<std::string, int>("exhaust_fan", EXHAUST_FAN),
            std::pair<std::string, int>("fan", FAN),
            std::pair<std::string, int>("fluid", FLUID),
            std::pair<std::string, int>("geometry", GEOMETRY),
            std::pair<std::string, int>("inlet_vent", INLET_VENT),
            std::pair<std::string, int>("intake_fan", INTAKE_FAN),
            std::pair<std::string, int>("interface", INTERFACE),
            std::pair<std::string, int>("interior", INTERIOR),
            std::pair<std::string, int>("internal", INTERNAL),
            std::pair<std::string, int>("mass_flow_inlet", MASS_FLOW_INLET),
            std::pair<std::string, int>("outflow", OUTFLOW),
            std::pair<std::string, int>("outlet_vent", OUTLET_VENT),
            std::pair<std::string, int>("parent_face", PARENT_FACE),
            std::pair<std::string, int>("porous_jump", POROUS_JUMP),
            std::pair<std::string, int>("pressure_far_field", PRESSURE_FAR_FIELD),
            std::pair<std::string, int>("pressure_inlet", PRESSURE_INLET),
            std::pair<std::string, int>("pressure_outlet", PRESSURE_OUTLET),
            std::pair<std::string, int>("radiator", RADIATOR),
            std::pair<std::string, int>("solid", SOLID),
            std::pair<std::string, int>("symmetry", SYMMETRY),
            std::pair<std::string, int>("velocity_inlet", VELOCITY_INLET),
            std::pair<std::string, int>("wall", WALL),
            std::pair<std::string, int>("wrapper", WRAPPER)
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_zone_type_str(x);
        else
            return it->second;
    }

    ZONE::ZONE(int zone, const std::string &zt, const std::string &name, int id) :
        SECTION(SECTION::ZONE),
        m_zoneID(zone),
        m_zoneType(formalize(zt)),
        m_zoneName(name),
        m_domainID(id)
    {
        if (!isValidStr(zt))
            throw invalid_zone_type_str(zt);
    }

    void MESH::readFromFile(const std::string &src, std::ostream &fout)
    {
        // Open grid file
        std::ifstream fin(src);
        if (fin.fail())
            throw std::runtime_error("Failed to open input mesh file: \"" + src + "\".");

        // Clear existing records if any.
        clear_entry();

        // Clear existing size if any.
        reset_counting();

        // Read contents
        while (!fin.eof())
        {
            skip_white(fin);
            eat(fin, '(');
            int ti = -1;
            fin >> std::dec >> ti;
            if (ti == SECTION::COMMENT)
            {
                eat(fin, '\"');
                std::string ts;
                char tc;
                while ((tc = fin.get()) != '\"')
                    ts.push_back(tc);
                eat(fin, ')');
                add_entry(new COMMENT(ts));
                skip_white(fin);
            }
            else if (ti == SECTION::HEADER)
            {
                eat(fin, '\"');
                std::string ts;
                char tc;
                while ((tc = fin.get()) != '\"')
                    ts.push_back(tc);
                eat(fin, ')');
                add_entry(new HEADER(ts));
                skip_white(fin);
            }
            else if (ti == SECTION::DIMENSION)
            {
                int nd = 0;
                fin >> std::dec >> nd;
                eat(fin, ')');
                add_entry(new DIMENSION(nd));
                skip_white(fin);
                m_dim = nd;
                m_is3D = (nd == 3);
            }
            else if (ti == SECTION::NODE)
            {
                eat(fin, '(');
                int zone;
                fin >> std::hex >> zone;
                if (zone == 0)
                {
                    // If zone-id is 0, indicating total number of nodes in the mesh.
                    int tmp;
                    fin >> std::hex;
                    fin >> tmp;
                    if (tmp != 1)
                        throw std::runtime_error("Invalid \"first-index\" in NODE declaration!");
                    fin >> m_totalNodeNum;
                    fout << "Total number of nodes: " << m_totalNodeNum << std::endl;
                    fin >> tmp;
                    if (tmp != 0)
                        throw std::runtime_error("Invalid \"type\" in NODE declaration!");
                    char ndc = fin.get();
                    if (ndc != ')')
                    {
                        fin >> tmp;
                        eat(fin, ')');
                    }
                    eat(fin, ')');
                }
                else
                {
                    // If zone-id is positive, it indicates the zone to which the nodes belong.
                    int first, last, tp, nd;
                    fin >> std::hex;
                    fin >> first >> last;
                    fin >> tp >> nd;
                    auto e = new NODE(zone, first, last, tp, nd);
                    eat(fin, ')');
                    eat(fin, '(');
                    fout << "Reading " << e->num() << " nodes in zone " << zone << " (from " << first << " to " << last << "), whose type is \"" << NODE::idx2str(tp) << "\"  ... ";

                    if (nd != dimension())
                        throw std::runtime_error("Inconsistent with previous DIMENSION declaration!");

                    if (nd == 3)
                    {
                        for (int i = first; i <= last; ++i)
                        {
                            auto &ce = e->at(i - first);
                            fin >> ce.x() >> ce.y() >> ce.z();
                        }
                    }
                    else
                    {
                        for (int i = first; i <= last; ++i)
                        {
                            auto &ce = e->at(i - first);
                            fin >> ce.x() >> ce.y();
                        }
                    }
                    eat(fin, ')');
                    eat(fin, ')');
                    fout << "Done!" << std::endl;
                    add_entry(e);
                }
                skip_white(fin);
            }
            else if (ti == SECTION::CELL)
            {
                eat(fin, '(');
                int zone;
                fin >> std::hex >> zone;
                if (zone == 0)
                {
                    // If zone-id is 0, indicating total number of cells in the mesh.
                    int tmp;
                    fin >> std::hex;
                    fin >> tmp;
                    if (tmp != 1)
                        throw std::runtime_error("Invalid \"first-index\" in CELL declaration!");
                    fin >> m_totalCellNum;
                    fout << "Total number of cells: " << m_totalCellNum << std::endl;
                    fin >> tmp;
                    if (tmp != 0)
                        throw std::runtime_error("Invalid \"type\" in CELL declaration!");
                    char ndc = fin.get();
                    if (ndc != ')')
                    {
                        fin >> tmp;
                        eat(fin, ')');
                    }
                    eat(fin, ')');
                }
                else
                {
                    // If zone-id is positive, it indicates the zone to which the cells belong.
                    int first, last, tp, elem;
                    fin >> std::hex;
                    fin >> first >> last;
                    fin >> tp >> elem;
                    auto e = new CELL(zone, first, last, tp, elem);
                    eat(fin, ')');

                    if (elem == 0)
                    {
                        fout << "Reading " << e->num() << " mixed cells in zone " << zone << " (from " << first << " to " << last << ") ... ";
                        eat(fin, '(');
                        for (int i = first; i <= last; ++i)
                        {
                            fin >> elem;
                            if (CELL::isValidElemIdx(elem))
                                e->at(i - first) = elem;
                            else
                                throw std::runtime_error("Invalid CELL-ELEM-TYPE: \"" + std::to_string(elem) + "\"");
                        }
                        eat(fin, ')');
                        fout << "Done!" << std::endl;
                    }
                    else
                        fout << e->num() << " " << CELL::idx2str_elem(elem) << " in zone " << zone << " (from " << first << " to " << last << ")" << std::endl;

                    eat(fin, ')');
                    add_entry(e);
                }
                skip_white(fin);
            }
            else if (ti == SECTION::FACE)
            {
                eat(fin, '(');
                int zone;
                fin >> std::hex >> zone;
                if (zone == 0)
                {
                    // If zone-id is 0, indicating total number of faces in the mesh.
                    int tmp;
                    fin >> tmp;
                    if (tmp != 1)
                        throw std::runtime_error("Invalid \"first-index\" in FACE declaration!");
                    fin >> m_totalFaceNum;
                    fout << "Total number of faces: " << m_totalFaceNum << std::endl;
                    fin >> tmp;
                    char ndc = fin.get();
                    if (ndc != ')')
                    {
                        fin >> tmp;
                        eat(fin, ')');
                    }
                    eat(fin, ')');
                }
                else
                {
                    // If zone-id is positive, it indicates a regular face section and will be
                    // followed by a body containing information about the grid connectivity.
                    size_t first, last;
                    int bc, face;
                    fin >> first >> last;
                    fin >> bc >> face;
                    auto e = new FACE(zone, first, last, bc, face);
                    eat(fin, ')');
                    eat(fin, '(');
                    fout << "Reading " << e->num() << " " << FACE::idx2str(face) << " faces in zone " << zone << " (from " << first << " to " << last << "), whose B.C. is \"" << BC::idx2str(bc) << "\" ... ";

                    std::vector<size_t> tmp_n;
                    size_t tmp_c[2];
                    if (face == FACE::MIXED)
                    {
                        int x = -1;
                        for (size_t i = first; i <= last; ++i)
                        {
                            // Read connectivity record
                            fin >> x;
                            if (x < 2)
                                throw std::invalid_argument("Invalid node num in the mixed face.");
                            tmp_n.resize(x);
                            for (int j = 0; j < x; ++j)
                                fin >> tmp_n[j];
                            fin >> tmp_c[0] >> tmp_c[1];

                            // Store current connectivity info
                            e->at(i - first).update_included_node(tmp_n);
                            e->at(i - first).update_adjacent_cell(tmp_c[0], tmp_c[1]);
                        }
                    }
                    else
                    {
                        tmp_n.resize(face);
                        for (size_t i = first; i <= last; ++i)
                        {
                            // Read connectivity record
                            for (int j = 0; j < face; ++j)
                                fin >> tmp_n[j];
                            fin >> tmp_c[0] >> tmp_c[1];

                            // Store current connectivity info
                            e->at(i - first).update_included_node(tmp_n);
                            e->at(i - first).update_adjacent_cell(tmp_c[0], tmp_c[1]);
                        }
                    }
                    eat(fin, ')');
                    eat(fin, ')');
                    fout << "Done!" << std::endl;
                    add_entry(e);
                }
                skip_white(fin);
            }
            else if (ti == SECTION::ZONE || ti == SECTION::ZONE_MESHING)
            {
                eat(fin, '(');
                int zone;
                fin >> std::dec >> zone;
                std::string ztp;
                fin >> ztp;
                skip_white(fin);
                std::string zname;
                char t0;
                while ((t0 = fin.get()) != ')')
                    zname.push_back(t0);
                eat(fin, '(');
                eat(fin, ')');
                eat(fin, ')');
                auto e = new ZONE(zone, ztp, zname);
                add_entry(e);
                skip_white(fin);
                fout << "ZONE " << e->zone() << ", named " << R"(")" << e->name() << R"(", )" << "is " << R"(")" << e->type() << R"(")" << std::endl;
                ++m_totalZoneNum;
            }
            else
                throw std::runtime_error("Unsupported section index: " + std::to_string(ti));
        }

        // Close grid file
        fin.close();
        fout << "Done!" << std::endl;
    }

    void MESH::writeToFile(const std::string &dst) const
    {
        if (nCell() == 0)
            throw std::runtime_error("Invalid num of cells.");
        if (nFace() == 0)
            throw std::runtime_error("Invalid num of faces.");
        if (nNode() == 0)
            throw std::runtime_error("Invalid num of nodes.");
        if (m_content.empty())
            throw std::runtime_error("Invalid num of contents.");

        /// Open grid file
        std::ofstream fout(dst);
        if (fout.fail())
            throw std::runtime_error("Failed to open output grid file: " + dst);

        /// Write until dimension declaration
        size_t i = 0;
        while (true)
        {
            m_content[i]->repr(fout);
            bool flag = dynamic_cast<DIMENSION*>(m_content[i]) != nullptr;
            ++i;
            if (flag)
                break;
        }

        /// Declaration of NODE, FACE, CELL
        fout << "(" << std::dec << SECTION::NODE << " (";
        fout << std::hex << 0 << " " << 1 << " " << m_totalNodeNum << " ";
        fout << std::dec << 0 << " " << (m_is3D ? 3 : 2) << "))" << std::endl;
        fout << "(" << std::dec << SECTION::CELL << " (";
        fout << std::hex << 0 << " " << 1 << " " << m_totalCellNum << " ";
        fout << std::dec << 0 << " " << 0 << "))" << std::endl;
        fout << "(" << std::dec << SECTION::FACE << " (";
        fout << std::hex << 0 << " " << 1 << " " << m_totalFaceNum << " ";
        fout << std::dec << 0 << " " << 0 << "))" << std::endl;

        /// Contents
        for (; i < m_content.size(); ++i)
            m_content[i]->repr(fout);

        /// Close grid file
        fout.close();
    }
}
