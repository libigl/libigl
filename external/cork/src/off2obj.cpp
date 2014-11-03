
#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;

int main(int argc, const char* argv[])
{
    string input;
    if(!(cin >> input) || input != "OFF")   return 1;
    
    int numvertices;
    int numfaces;
    int numedges;
    
    if(!(cin >> numvertices >> numfaces >> numedges) || numedges > 0)
        return 1;
    
    vector<double> vertices(numvertices*3);
    vector< vector<int> > faces(numfaces);
    
    for(int i=0; i<numvertices; i++)
    {
        cout << "v";
        for(int v=0; v<3; v++)
        {
            double coord;
            cin >> coord;
            cout << ' ' << coord;
        }
        cout << endl;
    }
    
    for(int i=0; i<numfaces; i++)
    {
        int vcount;
        if(!(cin >> vcount) || vcount < 3) return 1;
        
        cout << "f ";
        for(int v=0; v<vcount; v++)
        {
            int index;
            cin >> index;
            cout << ' ' << index+1;
        }
        cout << endl;
    }
    
    return 0;
}


