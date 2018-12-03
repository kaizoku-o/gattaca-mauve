#include "iostream"
#include "fstream"
#include "string"

using namespace std;

int main(int argc, char* argv[] ) {
    string in_file;
    in_file = argv[1];

    cout << "in_file: " << in_file << endl;
    ifstream inp_file(in_file);

    std::string str;
    string elapsed_time;
    string memAlloc;

    while (inp_file >> str) {
	if (str == "ElapsedTime:") {
	    inp_file >> str;
	    elapsed_time = str;
	}
	if (str == "MemAlloc:") {
	    inp_file >> str;
	    memAlloc = str;
	}
    }
    string fin_out = in_file + "," + elapsed_time + "," + memAlloc + "\n";
    cout << "Got: " << fin_out << endl;
    std::ofstream outfile;
    outfile.open("l_out.csv", std::ios_base::app);
    cout << "Wrote: " << fin_out << endl;
    outfile << fin_out;
    return 0;
}
