#include "xf.h"

int main(int argc, char* argv[])
{
    const std::string input_mesh = "../mesh/hex32.msh";
    const std::string output_mesh = "../mesh/hex32_blessed.msh";

    auto* obj = new XF::MESH(input_mesh, std::cout);
    obj->writeToFile(output_mesh);
    delete obj;

    return 0;
}
