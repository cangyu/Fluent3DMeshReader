#ifndef REP_H
#define REP_H

#include "xf.h"

namespace REP
{
    struct internal_error : public std::runtime_error
    {
        explicit internal_error(int err) : std::runtime_error("Internal error occurred with error code: " + std::to_string(err) + ".") {}

        explicit internal_error(const std::string& msg) : std::runtime_error("Internal error occurred with error message: \"" + msg + "\".") {}

        internal_error(int err, const std::string& msg) : std::runtime_error("Internal error occurred with error code: " + std::to_string(err) + " and error message: \"" + msg + "\".") {}

        internal_error(size_t err, const std::string& msg) : std::runtime_error("Internal error occurred with error code: " + std::to_string(err) + " and error message: \"" + msg + "\".") {}
    };

    struct NODE
    {
        size_t index;

        bool at_boundary;

        VEC coordinate;

        std::vector<size_t> adjacentNode;

        std::vector<size_t> dependentFace;

        std::vector<size_t> dependentCell;

        NODE() = delete;

        NODE(bool flag, double x, double y, double z) : index(0), at_boundary(flag), coordinate(x, y, z) {}
    };

    struct FACE
    {
        size_t index;

        bool at_boundary;

        int shape;

        VEC centroid;

        double area;

        /// Included nodes.
        std::vector<size_t> includedNode;

        /// Adjacent cells.
        /// Set to 0 if not exist.
        size_t c0, c1;

        /// Surface unit normal vector.
        VEC n01, n10;

        FACE() = delete;

        FACE(size_t idx, int sp) :
            index(idx),
            at_boundary(false),
            shape(sp),
            centroid(0.0, 0.0, 0.0),
            area(0.0),
            c0(0),
            c1(0),
            n01(0.0, 0.0, 0.0),
            n10(0.0, 0.0, 0.0)
        {
            if (shape == XF::FACE::TRIANGULAR)
                includedNode.resize(3);
            else if (shape == XF::FACE::QUADRILATERAL)
                includedNode.resize(4);
            else
                throw internal_error("Face shape not recognized");
        }
    };

    struct TRIANGLE : public FACE
    {
        TRIANGLE() = delete;

        explicit TRIANGLE(size_t idx) : FACE(idx, XF::FACE::TRIANGULAR) {}
    };

    struct QUAD : public FACE
    {
        QUAD() = delete;

        explicit QUAD(size_t idx) : FACE(idx, XF::FACE::QUADRILATERAL) {}
    };

    struct CELL
    {
        size_t index;

        int shape;

        VEC centroid;

        double volume;

        std::vector<size_t> includedNode;

        std::vector<size_t> includedFace;

        /// Cell connectivity.
        /// Size is equal to that of "includedFace".
        /// If adjacent cell is boundary, corresponding value will be set to 0.
        std::vector<size_t> adjacentCell;

        /// Surface outward unit normal vector.
        /// Size is equal to that of "includedFace".
        std::vector<VEC> n;

        CELL() = delete;

        CELL(size_t idx, int sp) :
            index(idx),
            shape(sp),
            centroid(0.0, 0.0, 0.0),
            volume(0.0)
        {
            if (shape == XF::CELL::TETRAHEDRAL)
            {
                includedNode.resize(4);
                includedFace.resize(4);
                adjacentCell.resize(4);
                n.resize(4);
            }
            else if (shape == XF::CELL::HEXAHEDRAL)
            {
                includedNode.resize(8);
                includedFace.resize(6);
                adjacentCell.resize(6);
                n.resize(6);
            }
            else if (shape == XF::CELL::PYRAMID)
            {
                includedNode.resize(5);
                includedFace.resize(5);
                adjacentCell.resize(5);
                n.resize(5);
            }
            else if (shape == XF::CELL::WEDGE)
            {
                includedNode.resize(6);
                includedFace.resize(5);
                adjacentCell.resize(5);
                n.resize(5);
            }
            else
                throw internal_error("Cell shape not recognized");
        }

        virtual void add_face(const FACE& f) = 0;
    };

    struct TET : public CELL
    {
    private:
        size_t stage;
        bool face_flag[4];

    public:
        TET() = delete;

        explicit TET(size_t idx) : CELL(idx, XF::CELL::TETRAHEDRAL), stage(0), face_flag{ false } {}

        void add_face(const FACE& f) override;

    private:
        void stage0_handler(const FACE& f);

        void stage1_handler(const FACE& f);

        void stage23_handler(const FACE& f);
    };

    struct HEX : public CELL
    {
    private:
        size_t stage;
        bool face_flag[6];

    public:
        HEX() = delete;

        explicit HEX(size_t idx) : CELL(idx, XF::CELL::HEXAHEDRAL), stage(0), face_flag{ false } {}

        void add_face(const FACE& f) override;

    private:
        void stage0_handler(const FACE& f);

        void stage15_handler(const FACE& f);
    };

    struct PYRAMID : public CELL
    {
        PYRAMID() = delete;

        explicit PYRAMID(size_t idx) : CELL(idx, XF::CELL::PYRAMID) {}

        void add_face(const FACE& f) override;
    };

    struct WEDGE : public CELL
    {
        WEDGE() = delete;

        explicit WEDGE(size_t idx) : CELL(idx, XF::CELL::WEDGE) {}

        void add_face(const FACE& f) override;
    };

    struct ZONE
    {
        /// Identifier of this boundary zone.
        std::string name;

        std::vector<size_t> includedFace;

        std::vector<size_t> includedNode;

        ZONE() = delete;

        explicit ZONE(const std::string& n) : name(n) {}
    };

    class Translator
    {
    protected:
        std::vector<NODE*> m_node;
        std::vector<FACE*> m_face;
        std::vector<CELL*> m_cell;
        std::vector<ZONE> m_zone;

    public:
        Translator() = delete;

        Translator(XF::MESH* mesh, std::ostream& operation_log);

        virtual void write(std::ostream& f_out) = 0;

        void dump_cell_connectivity(std::ostream &f_out);

        void dump_face_connectivity(std::ostream &f_out);

        void dump_node_connectivity(std::ostream &f_out);

    private:
        void extract_node_basic_info(XF::NODE* curObj);

        void extract_cell_basic_info(XF::CELL* curObj);

        void extract_face_basic_info(XF::FACE* curObj);

        void calculate_face_geom_var();

        void count_nodal_connectivity_step1(XF::FACE* curObj);

        void calculate_face_unit_normal();

        void calculate_cell_geom_var();
    };
}

#endif
