#include <iostream>
#include <fstream>
using namespace std;
int main()
{
        /*NuLatAB-PSD exp_0043_run_0703_subtracted.glenn 2-2-1-Event-B-2-2-1.txt test*/

        int minfile, maxfile, channel;
        cout << "input the nubmer of the first file needs to be analysised" << endl;
        cin >> minfile;
        cout << "input the nubmer of the last file needs to be analysised" << endl;
        cin >> maxfile;
        cout << "input the channel of the last file needs to be analysised" << endl;
        cin >> channel;
        string eventfile;
        cout << "input event file needs to be used" << endl;
        cin >> eventfile;
        string anafile;
        cout << "input the name of analysised file" << endl;
        cin >> anafile;
       
        ofstream fout;
        fout.open("PSD-shell");
        fout << "NuLatAB-PSD ";
        for (int i=0; i<(maxfile-minfile+1); i++)
        {     
                if ((minfile+i)<1000)          
                fout << "exp_0043_run_0"<<minfile+i<<"_subtracted.glenn ";
                else
                fout << "exp_0043_run_"<<minfile+i<<"_subtracted.glenn ";                
        }
        fout << channel << " " << eventfile << " " << anafile << endl;
        fout.clear();
        fout.close();
        return 0;
}

