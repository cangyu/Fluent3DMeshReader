#include "rep.h"

class BinaryWriter : public REP::Translator
{
public:
    BinaryWriter(XF::MESH* mesh, std::ostream& f_out) : Translator(mesh, f_out) {}

    void write(std::ostream& f_out) override;
};

void BinaryWriter::write(std::ostream& f_out)
{
    /// TODO
}

int main(int argc, char* argv[])
{
    /// TODO

    return 0;
}
