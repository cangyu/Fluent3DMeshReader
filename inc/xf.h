#ifndef XF_H
#define XF_H

#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <cstddef>
#include <utility>
#include <cmath>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include "vec.h"

namespace XF
{
    struct wrong_index : public std::logic_error
    {
        wrong_index(int idx, const std::string& reason) : std::logic_error("\"" + std::to_string(idx) + "\" " + reason + ".") {}

        wrong_index(size_t idx, const std::string& reason) : std::logic_error("\"" + std::to_string(idx) + "\" " + reason + ".") {}
    };

    struct wrong_string : public std::logic_error
    {
        wrong_string(const std::string& str, const std::string& reason) : std::logic_error("\"" + str + "\" " + reason + ".") {}
    };

    class DIM
    {
    private:
        struct wrong_dimension : public wrong_index
        {
            explicit wrong_dimension(int dim) : wrong_index(dim, "is not a valid dimension") {}
        };

    protected:
        bool m_is3D;
        int m_dim;

    public:
        DIM() = delete;

        explicit DIM(int dim)
        {
            if (dim == 2 || dim == 3)
                m_dim = dim;
            else
                throw wrong_dimension(dim);

            m_is3D = dim == 3;
        }

        DIM(int dim, bool is3d)
        {
            if (dim == 2 || dim == 3)
                m_dim = dim;
            else
                throw wrong_dimension(dim);

            m_is3D = is3d;
        }

        DIM(const DIM& rhs) = default;

        virtual ~DIM() = default;

        [[nodiscard]] bool is3D() const
        {
            return m_is3D;
        }

        [[nodiscard]] int dimension() const
        {
            return m_dim;
        }
    };

    class SECTION
    {
    public:
        enum {
            COMMENT = 0,
            HEADER = 1,
            DIMENSION = 2,
            NODE = 10,
            CELL = 12,
            FACE = 13,
            EDGE = 11,
            ZONE = 39,
            ZONE_MESHING = 45
        };

    private:
        int m_identity;

    public:
        SECTION() = delete;

        explicit SECTION(int id)
        {
            m_identity = id;
        }

        SECTION(const SECTION& rhs) = default;

        virtual ~SECTION() = default;

        virtual void repr(std::ostream& out) = 0;

        [[nodiscard]] int identity() const
        {
            return m_identity;
        }
    };

    class BC
    {
    public:
        struct invalid_bc_idx : public wrong_index
        {
            explicit invalid_bc_idx(int x) : wrong_index(x, "is not a valid B.C. index") {}
        };

        struct invalid_bc_str : public wrong_string
        {
            explicit invalid_bc_str(const std::string& s) : wrong_string(s, "is not a valid B.C. string") {}
        };

    public:
        enum {
            INTERIOR = 2,
            WALL = 3,
            PRESSURE_INLET = 4,
            INLET_VENT = 4,
            INTAKE_FAN = 4,
            PRESSURE_OUTLET = 5,
            EXHAUST_FAN = 5,
            OUTLET_VENT = 5,
            SYMMETRY = 7,
            PERIODIC_SHADOW = 8,
            PRESSURE_FAR_FIELD = 9,
            VELOCITY_INLET = 10,
            PERIODIC = 12,
            FAN = 14,
            POROUS_JUMP = 14,
            RADIATOR = 14,
            MASS_FLOW_INLET = 20,
            INTERFACE = 24,
            PARENT = 31,
            OUTFLOW = 36,
            AXIS = 37
        };

        static bool isValidIdx(int x);

        static bool isValidStr(const std::string& x);

        static const std::string& idx2str(int x);

        static int str2idx(const std::string& x);
    };

    class STR : public SECTION
    {
    private:
        std::string m_msg;

    public:
        STR() = delete;

        STR(int id, const std::string& msg) : SECTION(id), m_msg(msg) {}

        STR(const STR& rhs) : SECTION(rhs.identity()), m_msg(rhs.m_msg) {}

        virtual ~STR() = default;

        [[nodiscard]] const std::string& str() const
        {
            return m_msg;
        }

        std::string& str()
        {
            return m_msg;
        }

        void repr(std::ostream& out) override
        {
            out << "(" << std::dec << identity() << " \"" << str() << "\")" << std::endl;
        }
    };

    class COMMENT : public STR
    {
    public:
        COMMENT() = delete;

        explicit COMMENT(const std::string& info) : STR(SECTION::COMMENT, info) {}

        COMMENT(const COMMENT& rhs) : STR(SECTION::COMMENT, rhs.str()) {}

        ~COMMENT() = default;
    };

    class HEADER : public STR
    {
    public:
        HEADER() = delete;

        explicit HEADER(const std::string& info) : STR(SECTION::HEADER, info) {}

        HEADER(const HEADER& rhs) : STR(SECTION::HEADER, rhs.str()) {}

        ~HEADER() = default;
    };

    class DIMENSION :public SECTION, public DIM
    {
    public:
        DIMENSION() = delete;

        explicit DIMENSION(int dim, bool id3d = true) : SECTION(SECTION::DIMENSION), DIM(dim, id3d) {}

        DIMENSION(const DIMENSION& rhs) : SECTION(SECTION::DIMENSION), DIM(rhs.dimension(), rhs.is3D()) {}

        ~DIMENSION() = default;

        [[nodiscard]] int ND() const
        {
            return dimension();
        }

        void repr(std::ostream& out) override
        {
            out << "(" << std::dec << identity() << " " << ND() << ")" << std::endl;
        }
    };

    class RANGE : public SECTION
    {
    protected:
        size_t m_zone;
        size_t m_first, m_last;

    public:
        RANGE() = delete;

        RANGE(int id, size_t zone, size_t first, size_t last);

        RANGE(const RANGE& rhs);

        virtual ~RANGE() = default;

        [[nodiscard]] size_t zone() const
        {
            return m_zone;
        }

        [[nodiscard]] size_t first_index() const
        {
            return m_first;
        }

        [[nodiscard]] size_t last_index() const
        {
            return m_last;
        }

        [[nodiscard]] size_t num() const
        {
            const size_t ret = last_index() - first_index() + 1;
            return ret;
        }
    };

    class NODE : public RANGE, public DIM
    {
    public:
        struct invalid_node_type_idx : public wrong_index
        {
            explicit invalid_node_type_idx(int x) : wrong_index(x, "is not a valid NODE-TYPE index") {}
        };

        struct invalid_node_type_str : public wrong_string
        {
            explicit invalid_node_type_str(const std::string& s) : wrong_string(s, "is not a valid NODE-TYPE string") {}
        };

        enum {
            VIRTUAL = 0,
            ANY = 1,
            BOUNDARY = 2
        };

        static bool isValidTypeIdx(int x);

        static bool isValidTypeStr(const std::string& x);

        static const std::string& idx2str(int x);

        static int str2idx(const std::string& x);

    private:
        int m_type;
        std::vector<VEC> m_coordinate;

    public:
        NODE() = delete;

        NODE(size_t zone, size_t first, size_t last, int tp, int ND) :
            RANGE(SECTION::NODE, zone, first, last),
            DIM(ND),
            m_coordinate(num()),
            m_type(tp)
        {
            if (!isValidTypeIdx(type()))
                throw std::invalid_argument("Invalid description of node type in constructor.");
        }

        NODE(const NODE& rhs) :
            RANGE(SECTION::NODE, rhs.zone(), rhs.first_index(), rhs.last_index()),
            DIM(rhs.ND(), rhs.is3D()),
            m_coordinate(rhs.m_coordinate.size()),
            m_type(rhs.type())
        {
            if (!isValidTypeIdx(type()))
                throw std::invalid_argument("Invalid description of node type in copy-constructor.");

            for (int i = 0; i < m_coordinate.size(); ++i)
            {
                m_coordinate[i] = rhs.m_coordinate[i];
            }
        }

        ~NODE() = default;

        VEC& at(size_t loc_idx)
        {
            return m_coordinate.at(loc_idx);
        }

        [[nodiscard]] const VEC& at(size_t loc_idx) const
        {
            return m_coordinate.at(loc_idx);
        }

        [[nodiscard]] bool is_virtual_node() const
        {
            return type() == VIRTUAL;
        }

        [[nodiscard]] bool is_boundary_node() const
        {
            return type() == BOUNDARY;
        }

        [[nodiscard]] bool is_internal_node() const
        {
            return type() == ANY;
        }

        int& type()
        {
            return m_type;
        }

        [[nodiscard]] int type() const
        {
            return m_type;
        }

        [[nodiscard]] int ND() const
        {
            return dimension();
        }

        void repr(std::ostream& out) override;
    };

    class CELL : public RANGE
    {
    public:
        struct invalid_cell_type_idx : public wrong_index
        {
            explicit invalid_cell_type_idx(int x) : wrong_index(x, "is not a valid CELL-TYPE index") {}
        };

        struct invalid_cell_type_str : public wrong_string
        {
            explicit invalid_cell_type_str(const std::string& s) : wrong_string(s, "is not a valid CELL-TYPE string") {}
        };

        struct invalid_elem_type_idx : public wrong_index
        {
            explicit invalid_elem_type_idx(int x) : wrong_index(x, "is not a valid CELL-ELEM-TYPE index") {}
        };

        struct invalid_elem_type_str : public wrong_string
        {
            explicit invalid_elem_type_str(const std::string& s) : wrong_string(s, "is not a valid CELL-ELEM-TYPE string") {}
        };

        enum {
            DEAD = 0,
            FLUID = 1,
            SOLID = 17
        };

        static bool isValidTypeIdx(int x);

        static bool isValidTypeStr(const std::string& x);

        static const std::string& idx2str_type(int x);

        static int str2idx_type(const std::string& x);

        enum {
            MIXED = 0,
            TRIANGULAR = 1,
            TETRAHEDRAL = 2,
            QUADRILATERAL = 3,
            HEXAHEDRAL = 4,
            PYRAMID = 5,
            WEDGE = 6,
            POLYHEDRAL = 7
        };

        static bool isValidElemIdx(int x);

        static bool isValidElemStr(const std::string& x);

        static const std::string& idx2str_elem(int x);

        static int str2idx_elem(const std::string& x);

    private:
        int m_type;
        int m_elem;
        std::vector<int> m_desc;

    public:
        CELL() = delete;

        CELL(size_t zone, size_t first, size_t last, int type, int elem_type) :
            RANGE(SECTION::CELL, zone, first, last),
            m_type(type),
            m_elem(elem_type)
        {
            if (!isValidTypeIdx(type))
                throw invalid_cell_type_idx(type);

            if (!isValidElemIdx(elem_type))
                throw invalid_elem_type_idx(elem_type);

            if (elem_type != MIXED)
                m_desc.resize(num());
        }

        CELL(const CELL& rhs) :
            RANGE(SECTION::CELL, rhs.zone(), rhs.first_index(), rhs.last_index()),
            m_type(rhs.type()),
            m_elem(rhs.element_type()),
            m_desc(rhs.m_desc.size())
        {
            if (num() != rhs.num())
                throw std::runtime_error("Default copy operation is inconsistent.");

            for (int i = 0; i < m_desc.size(); ++i)
            {
                m_desc[i] = rhs.m_desc[i];
            }
        }

        ~CELL() = default;

        int& at(size_t loc_idx)
        {
            if (m_elem == MIXED)
                return m_desc.at(loc_idx);
            else
                return m_elem;
        }

        [[nodiscard]] int at(size_t loc_idx) const
        {
            if (m_elem == MIXED)
                return m_desc.at(loc_idx);
            else
                return m_elem;
        }

        /// Type of cells within this section: DEAD cell, FLUID cell or SOLID cell.
        [[nodiscard]] int type() const
        {
            return m_type;
        }

        int& type()
        {
            return m_type;
        }

        /// General description of ALL cell elements within this section.
        [[nodiscard]] int element_type() const
        {
            return  m_elem;
        }

        int& element_type()
        {
            return  m_elem;
        }

        void repr(std::ostream& out);
    };

    struct CONNECTIVITY
    {
        /// Nodes within this face.
        /// Ordered according to right-hand convention.
        std::vector<size_t> n;

        /// Adjacent cells.
        size_t c[2];

        CONNECTIVITY() : c{ 0, 0 } {}

        CONNECTIVITY(const CONNECTIVITY& rhs) : n(rhs.n), c{ rhs.c[0], rhs.c[1] } {}

        ~CONNECTIVITY() = default;

        [[nodiscard]] size_t c0() const
        {
            return c[0];
        }

        [[nodiscard]] size_t c1() const
        {
            return c[1];
        }

        void update_included_node(const std::vector<size_t>& src)
        {
            const auto N = src.size();
            if (N < 2)
                throw std::invalid_argument("Invalid num of nodes within a face.");

            n.resize(N);
            for (size_t i = 0; i < N; ++i)
                n[i] = src[i];
        }

        void update_adjacent_cell(size_t c0_idx, size_t c1_idx)
        {
            c[0] = c0_idx;
            c[1] = c1_idx;
        }

        void nodal_adjacency(size_t loc_idx, size_t& loc_forward, size_t& loc_backward) const
        {
            if (loc_idx == n.size() - 1)
                loc_forward = n[0];
            else
                loc_forward = n[loc_idx + 1];

            if (loc_idx == 0)
                loc_backward = n[n.size() - 1];
            else
                loc_backward = n[loc_idx - 1];
        }
    };

    class FACE : public RANGE
    {
    public:
        struct invalid_face_type_idx : public wrong_index
        {
            explicit invalid_face_type_idx(int x) : wrong_index(x, "is not a valid FACE-TYPE index") {}
        };

        struct invalid_face_type_str : public wrong_string
        {
            explicit invalid_face_type_str(const std::string& s) : wrong_string(s, "is not a valid FACE-TYPE string") {}
        };

        enum {
            MIXED = 0,
            LINEAR = 2,
            TRIANGULAR = 3,
            QUADRILATERAL = 4,
            POLYGONAL = 5
        };

        static bool isValidIdx(int x);

        static bool isValidStr(const std::string& x);

        static const std::string& idx2str(int x);

        static int str2idx(const std::string& x);

    private:
        int m_bc;
        int m_face;
        std::vector<CONNECTIVITY> m_desc;

    public:
        FACE() = delete;

        FACE(size_t zone, size_t first, size_t last, int bc, int face) :
            RANGE(SECTION::FACE, zone, first, last),
            m_desc(num()),
            m_bc(bc),
            m_face(face)
        {
            if (!BC::isValidIdx(bc))
                throw BC::invalid_bc_idx(bc);

            if (!isValidIdx(face))
                throw invalid_face_type_idx(face);
        }

        FACE(const FACE& rhs) :
            RANGE(SECTION::FACE, rhs.zone(), rhs.first_index(), rhs.last_index()),
            m_desc(rhs.m_desc.size()),
            m_bc(rhs.bc_type()),
            m_face(rhs.face_type())
        {
            if (!BC::isValidIdx(bc_type()))
                throw std::runtime_error("Invalid B.C. not detected in previous construction.");

            if (!isValidIdx(face_type()))
                throw std::runtime_error("Invalid FACE-TYPE not detected in previous construction.");

            for (int i = 0; i < m_desc.size(); ++i)
            {
                m_desc[i] = rhs.m_desc[i];
            }
        }

        ~FACE() = default;

        CONNECTIVITY& at(size_t loc_idx)
        {
            return m_desc.at(loc_idx);
        }

        [[nodiscard]] const CONNECTIVITY& at(size_t loc_idx) const
        {
            return m_desc.at(loc_idx);
        }

        /// B.C. of faces within this group.
        [[nodiscard]] int bc_type() const
        {
            return m_bc;
        }

        int& bc_type()
        {
            return m_bc;
        }

        /// Shape of faces within this group.
        [[nodiscard]] int face_type() const
        {
            return m_face;
        }

        int& face_type()
        {
            return m_face;
        }

        void repr(std::ostream& out);
    };

    class ZONE :public SECTION
    {
    public:
        struct invalid_zone_type_idx : public wrong_index
        {
            explicit invalid_zone_type_idx(int x) : wrong_index(x, "is not a valid ZONE-TYPE index") {}
        };

        struct invalid_zone_type_str : public wrong_string
        {
            explicit invalid_zone_type_str(const std::string& s) : wrong_string(s, "is not a valid specification of ZONE-TYPE") {}
        };

    public:
        enum {
            DEGASSING,
            EXHAUST_FAN,
            FAN,
            FLUID,
            GEOMETRY,
            INLET_VENT,
            INTAKE_FAN,
            INTERFACE,
            INTERIOR,
            INTERNAL,
            MASS_FLOW_INLET,
            OUTFLOW,
            OUTLET_VENT,
            PARENT_FACE,
            POROUS_JUMP,
            PRESSURE_FAR_FIELD,
            PRESSURE_INLET,
            PRESSURE_OUTLET,
            RADIATOR,
            SOLID,
            SYMMETRY,
            VELOCITY_INLET,
            WALL,
            WRAPPER
        };

        static bool isValidIdx(int x);

        static bool isValidStr(const std::string& x);

        static const std::string& idx2str(int x);

        static int str2idx(const std::string& x);

    private:
        size_t m_zoneID;
        std::string m_zoneType, m_zoneName;
        int m_domainID;

    public:
        ZONE() = delete;

        ZONE(int zone, const std::string& zt, const std::string& name, int id = 0);

        ZONE(const ZONE& rhs) = default;

        ~ZONE() = default;

        /// Index of this zone, may be any non-consecutive positive integer.
        [[nodiscard]] size_t zone() const
        {
            return m_zoneID;
        }

        size_t& zone()
        {
            return m_zoneID;
        }

        /// B.C. string literal.
        [[nodiscard]] const std::string& type() const
        {
            return m_zoneType;
        }

        std::string& type()
        {
            return m_zoneType;
        }

        /// Name of this zone.
        [[nodiscard]] const std::string& name() const
        {
            return m_zoneName;
        }

        std::string& name()
        {
            return m_zoneName;
        }

        /// Domain ID, NOT used.
        [[nodiscard]] int domain() const
        {
            return m_domainID;
        }

        int& domain()
        {
            return m_domainID;
        }

        void repr(std::ostream& out) override
        {
            out << std::dec << "(" << identity() << " (" << zone() << " " << type() << " " << name() << ")())" << std::endl;
        }
    };

    class MESH : public DIM
    {
    protected:
        std::vector<SECTION*> m_content;
        size_t m_totalNodeNum;
        size_t m_totalCellNum;
        size_t m_totalFaceNum;
        size_t m_totalZoneNum;

    public:
        MESH() = delete;

        MESH(const std::string& inp, std::ostream& log_out) : DIM(3)
        {
            readFromFile(inp, log_out);
        }

        MESH(const MESH& rhs) = delete;

        ~MESH()
        {
            clear_entry();
        }

        /// I/O
        void readFromFile(const std::string& src, std::ostream& f_out);

        void writeToFile(const std::string& dst) const;

        /// Num of elements
        [[nodiscard]] size_t nNode() const
        {
            return m_totalNodeNum;
        }

        [[nodiscard]] size_t nFace() const
        {
            return m_totalFaceNum;
        }

        [[nodiscard]] size_t nCell() const
        {
            return m_totalCellNum;
        }

        [[nodiscard]] size_t nZone() const
        {
            return m_totalZoneNum;
        }

        [[nodiscard]] size_t size() const
        {
            return m_content.size();
        }

        SECTION* at(size_t loc_idx)
        {
            return m_content.at(loc_idx);
        }

    private:
        void add_entry(SECTION* e)
        {
            m_content.push_back(e);
        }

        void clear_entry()
        {
            /// Release previous contents.
            for (auto ptr : m_content)
                delete ptr; /// NO side-effect when deleting a nullptr.

            /// Clear container.
            m_content.clear();
        }

        void reset_counting()
        {
            m_totalNodeNum = 0;
            m_totalCellNum = 0;
            m_totalFaceNum = 0;
            m_totalZoneNum = 0;
        }
    };
}

#endif
