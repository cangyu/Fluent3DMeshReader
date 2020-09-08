#include "rep.h"

template<typename T>
static bool contains(const std::vector<T>& src, const T& elem, size_t& idx)
{
    const auto N = src.size();
    for (size_t i = 0; i < N; ++i)
    {
        if (src[i] == elem)
        {
            idx = i;
            return true;
        }
    }
    return false;
}

static inline size_t forward_index(size_t i, size_t mod)
{
    if (i == mod - 1)
        return 0;
    else
        return i + 1;
}

static inline size_t backward_index(size_t i, size_t mod)
{
    if (i == 0)
        return mod - 1;
    else
        return i - 1;
}

void REP::TET::stage0_handler(const FACE& f)
{
    includedFace[0] = f.index;
    if (index == f.c0)
    {
        includedNode[3] = f.includedNode[0];
        includedNode[2] = f.includedNode[1];
        includedNode[1] = f.includedNode[2];
        adjacentCell[0] = f.c1;
    }
    else
    {
        includedNode[1] = f.includedNode[0];
        includedNode[2] = f.includedNode[1];
        includedNode[3] = f.includedNode[2];
        adjacentCell[0] = f.c0;
    }

    face_flag[0] = true;
    face_flag[1] = face_flag[2] = face_flag[3] = false;
}

void REP::TET::stage1_handler(const FACE& f)
{
    /// Find vertex 0
    bool ok = false;
    size_t next, prev;
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t target = f.includedNode[i];
        if (target != includedNode[1] && target != includedNode[2] && target != includedNode[3])
        {
            includedNode[0] = target;
            ok = true;
            next = f.includedNode[forward_index(i, 3)];
            prev = f.includedNode[backward_index(i, 3)];
            break;
        }
    }
    if (!ok)
        throw internal_error("Vertex0 not found in stage 1");

    /// Find corresponding face
    size_t loc;
    if (index == f.c0)
    {
        if (next == includedNode[1] && prev == includedNode[2])
            loc = 3;
        else if (next == includedNode[3] && prev == includedNode[1])
            loc = 2;
        else if (next == includedNode[2] && prev == includedNode[3])
            loc = 1;
        else
            throw internal_error("Face not found in stage 1 part 0");
    }
    else
    {
        if (next == includedNode[2] && prev == includedNode[1])
            loc = 3;
        else if (next == includedNode[1] && prev == includedNode[3])
            loc = 2;
        else if (next == includedNode[3] && prev == includedNode[2])
            loc = 1;
        else
            throw internal_error("Face not found in stage 1 part 1");
    }

    if (face_flag[loc])
        throw internal_error("Duplication detected");
    else
    {
        includedFace[loc] = f.index;
        adjacentCell[loc] = index == f.c0 ? f.c1 : f.c0;
        face_flag[loc] = true;
    }
}

void REP::TET::stage23_handler(const FACE& f)
{
    const std::string STAGE_STR = std::to_string(stage);

    /// Find vertex 0
    size_t i;
    size_t next, prev;
    if (contains(f.includedNode, includedNode[0], i))
    {
        next = f.includedNode[forward_index(i, 3)];
        prev = f.includedNode[backward_index(i, 3)];
    }
    else
        throw internal_error("Vertex0 not found in stage " + STAGE_STR);

    /// Find corresponding face
    size_t loc;
    if (index == f.c0)
    {
        if (next == includedNode[1] && prev == includedNode[2])
            loc = 3;
        else if (next == includedNode[3] && prev == includedNode[1])
            loc = 2;
        else if (next == includedNode[2] && prev == includedNode[3])
            loc = 1;
        else
            throw internal_error("Face not found in stage " + STAGE_STR + " part 0");
    }
    else
    {
        if (next == includedNode[2] && prev == includedNode[1])
            loc = 3;
        else if (next == includedNode[1] && prev == includedNode[3])
            loc = 2;
        else if (next == includedNode[3] && prev == includedNode[2])
            loc = 1;
        else
            throw internal_error("Face not found in stage " + STAGE_STR + " part 1");
    }

    if (face_flag[loc])
        throw internal_error("Duplication detected in stage " + STAGE_STR);
    else
    {
        includedFace[loc] = f.index;
        adjacentCell[loc] = index == f.c0 ? f.c1 : f.c0;
        face_flag[loc] = true;
    }
}

void REP::TET::add_face(const FACE& f)
{
    if (f.shape != XF::FACE::TRIANGULAR)
        throw internal_error(index, "Invalid shape");

    if (index != f.c0 && index != f.c1)
        throw internal_error(index, "Invalid connectivity detected");

    if (stage < 0 || stage >= 4)
        throw internal_error(index, "Unexpected stage");

    if (stage == 0)
        stage0_handler(f);
    else if (stage == 1)
        stage1_handler(f);
    else
        stage23_handler(f);

    ++stage;
}

void REP::HEX::stage0_handler(const FACE& f)
{
    includedFace[0] = f.index;
    if (index == f.c0)
    {
        includedNode[0] = f.includedNode[0];
        includedNode[3] = f.includedNode[1];
        includedNode[7] = f.includedNode[2];
        includedNode[4] = f.includedNode[3];
        adjacentCell[0] = f.c1;
    }
    else
    {
        includedNode[0] = f.includedNode[0];
        includedNode[4] = f.includedNode[1];
        includedNode[7] = f.includedNode[2];
        includedNode[3] = f.includedNode[3];
        adjacentCell[0] = f.c0;
    }

    face_flag[0] = true;
    for (int i = 1; i < 6; ++i)
        face_flag[i] = false;

    includedNode[1] = 0;
    includedNode[2] = 0;
    includedNode[6] = 0;
    includedNode[5] = 0;
}

void REP::HEX::stage15_handler(const FACE& f)
{
    int cnt = 0;
    for (auto e : f.includedNode)
    {
        if (e == includedNode[0] || e == includedNode[3] || e == includedNode[7] || e == includedNode[4])
            ++cnt;
    }

    size_t loc, pivot, pivot_adj;
    if (cnt == 2)
    {
        if (contains(f.includedNode, includedNode[0], pivot))
        {
            if (contains(f.includedNode, includedNode[4], pivot_adj))
            {
                loc = 2;
                size_t v5, v1;
                if (index == f.c0)
                {
                    v5 = f.includedNode[forward_index(pivot_adj, 4)];
                    v1 = f.includedNode[backward_index(pivot, 4)];
                }
                else
                {
                    v5 = f.includedNode[backward_index(pivot_adj, 4)];
                    v1 = f.includedNode[forward_index(pivot, 4)];
                }

                if (includedNode[5] == 0)
                    includedNode[5] = v5;
                else
                {
                    if (includedNode[5] != v5)
                        throw internal_error("Inconsistent at vertex5");
                }

                if (includedNode[1] == 0)
                    includedNode[1] = v1;
                else
                {
                    if (includedNode[1] != v1)
                        throw internal_error("Inconsistent at vertex1");
                }
            }
            else if (contains(f.includedNode, includedNode[3], pivot_adj))
            {
                loc = 4;
                size_t v1, v2;
                if (index == f.c0)
                {
                    v1 = f.includedNode[forward_index(pivot, 4)];
                    v2 = f.includedNode[backward_index(pivot_adj, 4)];
                }
                else
                {
                    v1 = f.includedNode[backward_index(pivot, 4)];
                    v2 = f.includedNode[forward_index(pivot_adj, 4)];
                }

                if (includedNode[1] == 0)
                    includedNode[1] = v1;
                else
                {
                    if (includedNode[1] != v1)
                        throw internal_error("Inconsistent at vertex1");
                }

                if (includedNode[2] == 0)
                    includedNode[2] = v2;
                else
                {
                    if (includedNode[2] != v2)
                        throw internal_error("Inconsistent at vertex2");
                }
            }
            else
                throw internal_error("Vertex0 pair not found");
        }
        else if (contains(f.includedNode, includedNode[7], pivot))
        {
            if (contains(f.includedNode, includedNode[4], pivot_adj))
            {
                loc = 5;
                size_t v6, v5;
                if (index == f.c0)
                {
                    v6 = f.includedNode[forward_index(pivot, 4)];
                    v5 = f.includedNode[backward_index(pivot_adj, 4)];
                }
                else
                {
                    v6 = f.includedNode[backward_index(pivot, 4)];
                    v5 = f.includedNode[forward_index(pivot_adj, 4)];
                }

                if (includedNode[6] == 0)
                    includedNode[6] = v6;
                else
                {
                    if (includedNode[6] != v6)
                        throw internal_error("Inconsistent at vertex6");
                }

                if (includedNode[5] == 0)
                    includedNode[5] = v5;
                else
                {
                    if (includedNode[5] != v5)
                        throw internal_error("Inconsistent at vertex5");
                }
            }
            else if (contains(f.includedNode, includedNode[3], pivot_adj))
            {
                loc = 3;
                size_t v2, v6;
                if (index == f.c0)
                {
                    v2 = f.includedNode[forward_index(pivot_adj, 4)];
                    v6 = f.includedNode[backward_index(pivot, 4)];
                }
                else
                {
                    v2 = f.includedNode[backward_index(pivot_adj, 4)];
                    v6 = f.includedNode[forward_index(pivot, 4)];
                }

                if (includedNode[2] == 0)
                    includedNode[2] = v2;
                else
                {
                    if (includedNode[2] != v2)
                        throw internal_error("Inconsistent at vertex2");
                }

                if (includedNode[6] == 0)
                    includedNode[6] = v6;
                else
                {
                    if (includedNode[6] != v6)
                        throw internal_error("Inconsistent at vertex6");
                }
            }
            else
                throw internal_error("Vertex7 pair not found");
        }
        else
            throw internal_error(index, "Inconsistent connectivity");
    }
    else if (cnt == 0)
        loc = 1;
    else
        throw internal_error(index, "Face does not match");

    if (face_flag[loc])
        throw internal_error("Duplication detected");
    else
    {
        includedFace[loc] = f.index;
        adjacentCell[loc] = index == f.c0 ? f.c1 : f.c0;
        face_flag[loc] = true;
    }
}

void REP::HEX::add_face(const FACE& f)
{
    if (f.shape != XF::FACE::QUADRILATERAL)
        throw internal_error(index, "Invalid shape");

    if (index != f.c0 && index != f.c1)
        throw internal_error(index, "Invalid connectivity detected");

    if (stage < 0 || stage >= 6)
        throw internal_error(index, "Unexpected stage");

    if (stage == 0)
        stage0_handler(f);
    else
        stage15_handler(f);

    ++stage;
}

void REP::PYRAMID::add_face(const FACE& f)
{
    if (index != f.c0 && index != f.c1)
        throw internal_error(index, "Invalid connectivity detected");

    /// TODO
}

void REP::WEDGE::add_face(const FACE& f)
{
    if (index != f.c0 && index != f.c1)
        throw internal_error(index, "Invalid connectivity detected");

    /// TODO
}

void REP::Translator::extract_node_basic_info(XF::NODE* curObj)
{
    /// Node type within this zone
    const bool tp = curObj->is_boundary_node();

    /// 1-based global node index
    const size_t cur_first = curObj->first_index();
    const size_t cur_last = curObj->last_index();
    for (size_t i = cur_first; i <= cur_last; ++i)
    {
        const auto& c = curObj->at(i - cur_first);
        m_node.at(i - 1) = new NODE(tp, c.x(), c.y(), c.z());
        m_node.at(i - 1)->index = i;
    }
}

void REP::Translator::extract_cell_basic_info(XF::CELL* curObj)
{
    /// 1-based global face index
    const size_t cur_first = curObj->first_index();
    const size_t cur_last = curObj->last_index();
    for (size_t i = cur_first; i <= cur_last; ++i)
    {
        /// Element type of current cell in this zone
        const auto elem = curObj->at(i - cur_first);
        switch (elem)
        {
        case XF::CELL::TETRAHEDRAL:
            m_cell.at(i - 1) = new TET(i);
            break;
        case XF::CELL::HEXAHEDRAL:
            m_cell.at(i - 1) = new HEX(i);
            break;
        case XF::CELL::PYRAMID:
            m_cell.at(i - 1) = new PYRAMID(i);
            break;
        case XF::CELL::WEDGE:
            m_cell.at(i - 1) = new WEDGE(i);
            break;
        default:
            throw internal_error("Unexpected cell shape");
        }
    }
}

void REP::Translator::extract_face_basic_info(XF::FACE* curObj)
{
    /// 1-based global face index
    const size_t cur_first = curObj->first_index();
    const size_t cur_last = curObj->last_index();
    for (size_t i = cur_first; i <= cur_last; ++i)
    {
        const auto& c = curObj->at(i - cur_first);

        const auto nn = c.n.size();
        switch (nn)
        {
        case 2:
            throw internal_error("This is a 3D program");
        case 3:
            m_face.at(i - 1) = new TRIANGLE(i);
            break;
        case 4:
            m_face.at(i - 1) = new QUAD(i);
            break;
        default:
            throw internal_error("Unexpected face shape");
        }

        auto dst = m_face.at(i - 1);

        /// Nodes within this face.
        /// 1-based node index are stored.
        /// Right-hand convention is preserved.
        for (int j = 0; j < nn; ++j)
            dst->includedNode[j] = c.n[j];

        /// Adjacent cells.
        /// 1-based cell index are stored, 0 stands for boundary.
        /// Right-hand convention is preserved.
        dst->c0 = c.c0();
        dst->c1 = c.c1();
        if (dst->c0 == 0 && dst->c1 == 0)
            throw internal_error("Surface mesh is not expected currently");

        /// Boundary flag.
        dst->at_boundary = (dst->c0 == 0 || dst->c1 == 0);

        /// Cell construction.
        if (dst->c0 != 0)
            m_cell.at(dst->c0 - 1)->add_face(*dst);
        if (dst->c1 != 0)
            m_cell.at(dst->c1 - 1)->add_face(*dst);
    }
}

void REP::Translator::calculate_face_geom_var()
{
    const size_t N = m_face.size();
    for (size_t i = 0; i < N; ++i)
    {
        auto f = m_face[i];

        if (f->shape == XF::FACE::TRIANGULAR)
        {
            const auto& p1 = m_node.at(f->includedNode[0] - 1)->coordinate;
            const auto& p2 = m_node.at(f->includedNode[1] - 1)->coordinate;
            const auto& p3 = m_node.at(f->includedNode[2] - 1)->coordinate;

            f->area = triangle_area(p1, p2, p3);
            triangle_center(p1, p2, p3, f->centroid);
            triangle_normal(p1, p2, p3, f->n01);
        }
        else if (f->shape == XF::FACE::QUADRILATERAL)
        {
            const auto& p1 = m_node.at(f->includedNode[0] - 1)->coordinate;
            const auto& p2 = m_node.at(f->includedNode[1] - 1)->coordinate;
            const auto& p3 = m_node.at(f->includedNode[2] - 1)->coordinate;
            const auto& p4 = m_node.at(f->includedNode[3] - 1)->coordinate;

            f->area = quadrilateral_area(p1, p2, p3, p4);
            quadrilateral_center(p1, p2, p3, p4, f->centroid);
            quadrilateral_normal(p1, p2, p3, p4, f->n01);
        }
        else if (f->shape == XF::FACE::LINEAR)
            throw internal_error("This is a 3D program");
        else
            throw internal_error("Shape not recognized");
    }
}

void REP::Translator::count_nodal_connectivity_step1(XF::FACE* curObj)
{
    /// 1-based index
    const size_t cur_first = curObj->first_index();
    const size_t cur_last = curObj->last_index();
    for (size_t i = cur_first; i <= cur_last; ++i)
    {
        const auto& c = curObj->at(i - cur_first);

        const size_t loc_c0 = c.c0();
        const size_t loc_c1 = c.c1();

        for (size_t j = 0; j < c.n.size(); ++j)
        {
            auto curNode = m_node.at(c.n[j] - 1);

            size_t loc_backwardNode, loc_forwardNode;
            c.nodal_adjacency(j, loc_forwardNode, loc_backwardNode);

            /// Adjacent nodes
            curNode->adjacentNode.push_back(loc_backwardNode);
            curNode->adjacentNode.push_back(loc_forwardNode);

            /// Dependent faces
            curNode->dependentFace.push_back(i);

            /// Dependent cells
            if (loc_c0 != 0)
                curNode->dependentCell.push_back(loc_c0);
            if (loc_c1 != 0)
                curNode->dependentCell.push_back(loc_c1);
        }
    }
}

void REP::Translator::calculate_face_unit_normal()
{
    const size_t NC = m_cell.size();

    for (size_t i = 0; i < NC; ++i)
    {
        auto curCell = m_cell.at(i);
        VEC approximate_centroid(0.0, 0.0, 0.0);
        for (auto j : curCell->includedNode)
        {
            auto curNode = m_node.at(j - 1);
            approximate_centroid += curNode->coordinate;
        }
        approximate_centroid /= static_cast<double>(curCell->includedNode.size());

        const auto NCF = curCell->includedFace.size();
        for (size_t j = 0; j < NCF; ++j)
        {
            const auto curFaceIdx = curCell->includedFace.at(j);
            auto curFace = m_face.at(curFaceIdx - 1);

            VEC approximate_r = curFace->centroid;
            approximate_r -= approximate_centroid;

            auto& curSN = curCell->n.at(j);
            curSN = curFace->n01;
            if (approximate_r.dot(curSN) < 0.0)
                curSN *= -1.0;
        }
    }

    for (size_t i = 1; i <= NC; ++i)
    {
        auto curCell = m_cell.at(i - 1);

        const auto NCF = curCell->includedFace.size();
        for (size_t j = 0; j < NCF; ++j)
        {
            const auto curFaceIdx = curCell->includedFace.at(j);
            auto curFace = m_face.at(curFaceIdx - 1);

            const auto& curSN = curCell->n.at(j);

            if (curFace->c0 == i)
            {
                curFace->n01 = curSN;
                if (curFace->c1 == 0)
                {
                    curFace->n10 = curSN;
                    curFace->n10 *= -1.0;
                }
            }
            else if (curFace->c1 == i)
            {
                curFace->n10 = curSN;
                if (curFace->c0 == 0)
                {
                    curFace->n01 = curSN;
                    curFace->n01 *= -1.0;
                }
            }
            else
                throw internal_error(-6, "Empty connectivity detected");
        }
    }

    for (auto f : m_face)
    {
        static const auto EPS = std::numeric_limits<double>::epsilon();
        const auto ans = f->n01.dot(f->n10);
        if (std::fabs(ans + 1.0) > EPS)
            throw internal_error(f->index, "Face unit norm vectors not match");
    }
}

void REP::Translator::calculate_cell_geom_var()
{
    const size_t NC = m_cell.size();

    for (size_t i = 0; i < NC; ++i)
    {
        auto curCell = m_cell.at(i);

        curCell->volume = 0.0;
        curCell->centroid.x() = 0.0;
        curCell->centroid.y() = 0.0;
        curCell->centroid.z() = 0.0;

        const auto NCF = curCell->includedFace.size();
        for (size_t j = 0; j < NCF; ++j)
        {
            const auto curFaceIdx = curCell->includedFace.at(j);
            auto curFace = m_face.at(curFaceIdx - 1);

            /// Based on the divergence theorem.
            /// See (5.15) and (5.17) of Jiri Blazek's CFD book.
            const auto& cf_c = curFace->centroid;
            const auto& cf_n = curCell->n.at(j);
            const auto w = cf_c.dot(cf_n) * curFace->area;
            curCell->volume += w;

            curCell->centroid.x() += w * cf_c.x();
            curCell->centroid.y() += w * cf_c.y();
            curCell->centroid.z() += w * cf_c.z();
        }
        curCell->volume /= 3.0;
        curCell->centroid /= (4.0 * curCell->volume);
    }
}

REP::Translator::Translator(XF::MESH* mesh, std::ostream& operation_log)
{
    operation_log << "======================================================================" << std::endl;
    operation_log << "                          Fluent3DMeshReader                          " << std::endl;
    operation_log << "     A package converting 3D FLUENT mesh file into custom format.     " << std::endl;
    operation_log << "======================================================================" << std::endl;

    operation_log << "Counting ..." << std::endl;
    const size_t num_of_node = mesh->nNode();
    const size_t num_of_face = mesh->nFace();
    const size_t num_of_cell = mesh->nCell();
    const size_t num_of_section = mesh->size();
    operation_log << num_of_node << " nodes" << std::endl;
    operation_log << num_of_face << " faces" << std::endl;
    operation_log << num_of_cell << " cells" << std::endl;
    operation_log << num_of_section << " sections" << std::endl;

    operation_log << "Allocating storage ..." << std::endl;
    m_node.resize(num_of_node, nullptr);
    m_face.resize(num_of_face, nullptr);
    m_cell.resize(num_of_cell, nullptr);

    operation_log << "Extracting description from ZONE sections ..." << std::endl;
    std::vector<size_t> zone_idx;
    m_zone.clear();
    for (size_t i0 = 0; i0 < num_of_section; ++i0)
    {
        auto curObj = dynamic_cast<XF::ZONE*>(mesh->at(i0));
        if (curObj == nullptr)
            continue;

        if (XF::ZONE::str2idx(curObj->type()) == XF::ZONE::INTERIOR)
            continue;

        const auto curZoneIdx = curObj->zone();

        for (size_t j0 = 0; j0 < num_of_section; ++j0)
        {
            auto curObj_aux = dynamic_cast<XF::FACE*>(mesh->at(j0));
            if (curObj_aux == nullptr)
                continue;

            if (curObj_aux->zone() == curZoneIdx)
            {
                m_zone.emplace_back(curObj->name());
                zone_idx.push_back(curZoneIdx);
                break;
            }
        }
    }
    for (size_t i = 0; i < zone_idx.size(); ++i)
    {
        const auto curZoneIdx = zone_idx[i];
        auto& cz = m_zone.at(i);

        for (size_t j = 0; j < num_of_section; ++j)
        {
            auto curObj = dynamic_cast<XF::FACE*>(mesh->at(j));
            if (curObj == nullptr)
                continue;

            if (curObj->zone() == curZoneIdx)
            {
                cz.includedFace.resize(curObj->num());
                for (size_t k = 0; k < cz.includedFace.size(); ++k)
                    cz.includedFace[k] = curObj->first_index() + k;

                std::set<size_t> nl;
                for (size_t k = 0; k < cz.includedFace.size(); ++k)
                {
                    const auto& connection = curObj->at(k);
                    for (const auto& e : connection.n)
                        nl.insert(e);
                }
                cz.includedNode.assign(nl.begin(), nl.end());
            }
        }
    }

    operation_log << "Extracting basic description from NODE sections ..." << std::endl;
    for (size_t i0 = 0; i0 < num_of_section; ++i0)
    {
        auto curPtr = mesh->at(i0);
        if (curPtr->identity() == XF::SECTION::NODE)
        {
            auto curObj = dynamic_cast<XF::NODE*>(curPtr);
            if (curObj == nullptr)
                throw internal_error(-1);
            else
                extract_node_basic_info(curObj);
        }
    }

    operation_log << "Extracting basic description from CELL sections ..." << std::endl;
    for (size_t i0 = 0; i0 < num_of_section; ++i0)
    {
        auto curPtr = mesh->at(i0);
        if (curPtr->identity() == XF::SECTION::CELL)
        {
            auto curObj = dynamic_cast<XF::CELL*>(curPtr);
            if (curObj == nullptr)
                throw internal_error(-2);
            else
                extract_cell_basic_info(curObj);
        }
    }

    operation_log << "Extracting basic description from FACE sections ..." << std::endl;
    for (size_t i0 = 0; i0 < num_of_section; ++i0)
    {
        auto curPtr = mesh->at(i0);
        if (curPtr->identity() == XF::SECTION::FACE)
        {
            auto curObj = dynamic_cast<XF::FACE*>(curPtr);
            if (curObj == nullptr)
                throw internal_error(-3);
            else
                extract_face_basic_info(curObj);
        }
    }

    operation_log << "Calculating geometric attributes of faces ..." << std::endl;
    calculate_face_geom_var();

    operation_log << "Counting adjacent nodes, dependent faces, and dependent cells of each node ..." << std::endl;
    /// Step1: Count all occurrence
    for (size_t i0 = 0; i0 < num_of_section; ++i0)
    {
        auto curPtr = mesh->at(i0);
        if (curPtr->identity() == XF::SECTION::FACE)
        {
            auto curObj = dynamic_cast<XF::FACE*>(curPtr);
            if (curObj == nullptr)
                throw internal_error(-4);
            else
                count_nodal_connectivity_step1(curObj);
        }
    }
    /// Step2: Remove duplication
    for (size_t i = 0; i < num_of_node; ++i)
    {
        auto curNode = m_node.at(i);

        const std::set<size_t> st1(curNode->adjacentNode.begin(), curNode->adjacentNode.end());
        curNode->adjacentNode.assign(st1.begin(), st1.end());

        const std::set<size_t> st2(curNode->dependentCell.begin(), curNode->dependentCell.end());
        curNode->dependentCell.assign(st2.begin(), st2.end());
    }

    operation_log << "Calculating surface unit normal ..." << std::endl;
    calculate_face_unit_normal();

    operation_log << "Calculating CELL volume and centroid ..." << std::endl;
    calculate_cell_geom_var();

    operation_log << "Finished!" << std::endl;
}

void REP::Translator::dump_cell_connectivity(std::ostream& f_out)
{
    static const char SEP = ' ';

    f_out << m_cell.size() << std::endl;

    for (auto c : m_cell)
    {
        for (const auto& e : c->includedNode)
        {
            f_out << e << SEP;
        }
        f_out << std::endl;
    }
}

void REP::Translator::dump_face_connectivity(std::ostream& f_out)
{
    static const char SEP = ' ';

    f_out << m_face.size() << std::endl;

    for (auto f : m_face)
    {
        f_out << f->includedNode.size() << SEP;
        for (const auto& e : f->includedNode)
        {
            f_out << e << SEP;
        }
        f_out << f->c0 << SEP << f->c1 << std::endl;
    }
}

void REP::Translator::dump_node_connectivity(std::ostream& f_out)
{
    static const char SEP = ' ';

    f_out << m_node.size() << std::endl;

    for (auto n : m_node)
    {
        f_out << n->adjacentNode.size() << SEP;
        for (const auto& e : n->adjacentNode)
            f_out << e << SEP;

        f_out << n->dependentFace.size() << SEP;
        for (const auto& e : n->dependentFace)
            f_out << e << SEP;

        f_out << n->dependentCell.size() << SEP;
        for (const auto& e : n->dependentCell)
            f_out << e << SEP;

        f_out << std::endl;
    }
}
