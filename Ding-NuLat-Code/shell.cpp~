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
        string pedfile;
        cout << "input pedistal file needs to be used" << endl;
        cin >> pedfile;
        string anafile;
        cout << "input the name of analysised file" << endl;
        cin >> anafile;
       
        ofstream fout;
        fout.open("analysis-shell");
        for (int i=0; i<(maxfile-minfile+1); i++)
        {                
                if ((minfile+i)<1000)                
                {                
                        fout << "scrod2subtract exp_0043_run_0" <<minfile+i<<".dat "<<pedfile << endl;
                        fout << "scrod2glenn exp_0043_run_0"<<minfile+i<<"_subtracted.dat"<<endl;
                       // fout << "rm exp_0043_run_0"<<minfile+i<<".dat"<<endl;
                }
                
                else
                {
                        fout << "scrod2subtract exp_0043_run_" <<minfile+i<<".dat "<<pedfile << endl;
                        fout << "scrod2glenn exp_0043_run_"<<minfile+i<<"_subtracted.dat"<<endl;
                        //fout << "rm exp_0043_run_"<<minfile+i<<".dat"<<endl;
             
                }
                
        }
        fout << "NuLatMap" << endl;
        for (int i=0; i<(maxfile-minfile+1); i++)
        {     
                
                if ((minfile+i)<1000)                
                {                
                     fout << "NuLatAB exp_0043_run_0"<<minfile+i<<"_subtracted.glenn "<<anafile << endl;
                }
                //fout << "rm exp_0043_run_0"<<minfile+i<<".dat"<<endl;
                else
                {
                     fout << "NuLatAB exp_0043_run_"<<minfile+i<<"_subtracted.glenn "<<anafile << endl;
                }
           
           
        }
        fout.clear();
        fout.close();
        return 0;
}

