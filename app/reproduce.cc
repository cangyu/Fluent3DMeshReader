#include "rep.h"

int main(int argc, char* argv[])
{
    std::string input_mesh;
    std::string output_mesh;

    if(argc == 3)
    {
        input_mesh = argv[1];
        output_mesh = argv[2];
    }
    else
        return -1;

    auto* mesh = new XF::MESH(input_mesh, std::cout);
    mesh->writeToFile(output_mesh);

    delete mesh;

    return 0;
}
