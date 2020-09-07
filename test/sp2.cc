#include "rep.h"

int main(int argc, char *argv[])
{
    const std::string input_mesh = "../mesh/hex32.msh";
    const std::string output_mesh = "../mesh/hex32_blessed.msh";

    auto *mesh = new XF::MESH(input_mesh, std::cout);
    mesh->writeToFile(output_mesh);

    auto *trans = new REP::Translator(mesh, std::cout);

    const std::string translated_mesh = "../mesh/hex32_translated.txt";
    std::ofstream f_out(translated_mesh);
    if(f_out.fail())
        throw std::runtime_error("Failed to open target binary mesh file.");

    f_out.setf(std::ios::scientific, std::ios::floatfield);
    f_out.precision(15);
    trans->write(f_out);

    f_out.close();

    delete mesh;
    delete trans;

    return 0;
}
