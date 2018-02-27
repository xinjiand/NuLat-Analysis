#include <iostream>
#include <fstream>
using namespace std;
int main()
{
        int minfile, maxfile;
        cout << "input the nubmer of the first file needs to be analysised" << endl;
        cin >> minfile;
        cout << "input the nubmer of the last file needs to be analysised" << endl;
        cin >> maxfile;
        ofstream fout;
        fout.open("event-list");
        for (int i=0; i<(maxfile-minfile+1); i++)
        {     
                           
                fout << minfile+i<< endl;
        }
        fout.clear();
        fout.close();
        return 0;
}

