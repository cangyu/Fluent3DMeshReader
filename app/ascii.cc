#include <cstring>
#include "rep.h"

class AsciiWriter : public REP::Translator
{
public:
    AsciiWriter(XF::MESH* mesh, std::ostream& f_out) : Translator(mesh, f_out) {}

    void write(std::ostream& f_out) override
    {
        static const char SEP = ' ';

        const size_t N_NODE = m_node.size();
        const size_t N_FACE = m_face.size();
        const size_t N_CELL = m_cell.size();
        const size_t N_ZONE = m_zone.size();

        f_out << N_NODE << SEP << N_FACE << SEP << N_CELL << SEP << N_ZONE << std::endl;

        for (size_t i = 0; i < N_NODE; ++i)
        {
            auto cn = m_node.at(i);

            f_out << (cn->at_boundary ? 1 : 0) << SEP;

            f_out << cn->coordinate.x() << SEP;
            f_out << cn->coordinate.y() << SEP;
            f_out << cn->coordinate.z() << SEP;

            f_out << cn->adjacentNode.size() << SEP;
            for (const auto& e : cn->adjacentNode)
                f_out << e << SEP;

            f_out << cn->dependentFace.size() << SEP;
            for (const auto& e : cn->dependentFace)
                f_out << e << SEP;

            f_out << cn->dependentCell.size() << SEP;
            for (const auto& e : cn->dependentCell)
                f_out << e << SEP;

            f_out << std::endl;
        }

        for (size_t i = 0; i < N_FACE; ++i)
        {
            auto cf = m_face.at(i);

            f_out << (cf->at_boundary ? 1 : 0) << SEP;

            f_out << cf->shape << SEP;

            f_out << cf->centroid.x() << SEP;
            f_out << cf->centroid.y() << SEP;
            f_out << cf->centroid.z() << SEP;

            f_out << cf->area << SEP;

            for (const auto& e : cf->includedNode)
                f_out << e << SEP;

            f_out << cf->c0 << SEP;
            f_out << cf->c1 << SEP;

            f_out << cf->n01.x() << SEP;
            f_out << cf->n01.y() << SEP;
            f_out << cf->n01.z() << SEP;

            f_out << cf->n10.x() << SEP;
            f_out << cf->n10.y() << SEP;
            f_out << cf->n10.z() << SEP;

            f_out << std::endl;
        }

        for (size_t i = 0; i < N_CELL; ++i)
        {
            auto cc = m_cell.at(i);

            f_out << cc->shape << SEP;

            f_out << cc->centroid.x() << SEP;
            f_out << cc->centroid.y() << SEP;
            f_out << cc->centroid.z() << SEP;

            f_out << cc->volume << SEP;

            for (const auto& e : cc->includedNode)
                f_out << e << SEP;

            for (const auto& e : cc->includedFace)
                f_out << e << SEP;

            for (const auto& e : cc->adjacentCell)
                f_out << e << SEP;

            for (const auto& e : cc->n)
            {
                f_out << e.x() << SEP;
                f_out << e.y() << SEP;
                f_out << e.z() << SEP;
            }

            f_out << std::endl;
        }

        for (size_t i = 0; i < N_ZONE; ++i)
        {
            const auto& cz = m_zone.at(i);

            f_out << cz.name << SEP << cz.includedFace.size() << SEP << cz.includedNode.size() << std::endl;

            for (const auto& e : cz.includedFace)
                f_out << e << SEP;
            f_out << std::endl;

            for (const auto& e : cz.includedNode)
                f_out << e << SEP;
            f_out << std::endl;
        }
    }
};

int main(int argc, char* argv[])
{
    std::string input_mesh_path;
    std::string translated_mesh_path;
    std::string cell_connectivity_path;
    std::string face_connectivity_path;
    std::string node_connectivity_path;

    int cnt = 1;
    while (cnt < argc)
    {
        if (!std::strcmp(argv[cnt], "--mesh"))
            input_mesh_path = argv[cnt + 1];
        else if (!std::strcmp(argv[cnt], "--translate-to"))
            translated_mesh_path = argv[cnt + 1];
        else if (!std::strcmp(argv[cnt], "--cell-connectivity"))
            cell_connectivity_path = argv[cnt + 1];
        else if (!std::strcmp(argv[cnt], "--face-connectivity"))
            face_connectivity_path = argv[cnt + 1];
        else if (!std::strcmp(argv[cnt], "--node-connectivity"))
            node_connectivity_path = argv[cnt + 1];
        else
            throw std::invalid_argument("Unrecognized option: \"" + std::string(argv[cnt]) + "\".");

        cnt += 2;
    }

    auto* mesh = new XF::MESH(input_mesh_path, std::cout);

    auto* trans = new AsciiWriter(mesh, std::cout);

    std::ofstream f_out(translated_mesh_path);
    if (f_out.fail())
        throw std::runtime_error("Failed to open target file.");
    f_out.setf(std::ios::scientific, std::ios::floatfield);
    f_out.precision(15);
    trans->write(f_out);
    f_out.close();

    if (!cell_connectivity_path.empty())
    {
        std::ofstream f_out0(cell_connectivity_path);
        if (f_out0.fail())
            throw std::runtime_error("Failed to open target file.");
        trans->dump_cell_connectivity(f_out0);
        f_out0.close();
    }

    if (!face_connectivity_path.empty())
    {
        std::ofstream f_out0(face_connectivity_path);
        if (f_out0.fail())
            throw std::runtime_error("Failed to open target file.");
        trans->dump_face_connectivity(f_out0);
        f_out0.close();
    }

    if (!node_connectivity_path.empty())
    {
        std::ofstream f_out0(node_connectivity_path);
        if (f_out0.fail())
            throw std::runtime_error("Failed to open target file.");
        trans->dump_node_connectivity(f_out0);
        f_out0.close();
    }

    delete mesh;
    delete trans;

    return 0;
}
