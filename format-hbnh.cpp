#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>

// This code corresponds to depositing HBNH (4-atom molecule) into the system
using namespace std;

const int N = 1000;	// Maximum number of possible atoms

int main(int argc, const char* argv[])
{
    string line;
    ifstream infile;
    ofstream ofile1,ofile2;
    
    char element[N];
    int i,num=0,num_final;
    double x[N],y[N],z[N],vx[N],vy[N],vz[N], cmx=0, cmy=0, cmz=0;
    double vel_HBNH,theta,phi,theta2,phi2;
    
  // Usage hint
    if (argc < 2) {
        cerr<<"Usage:"<<argv[0]<<" [# of hit] "<<endl;
        cerr<<"It will generate input file named bnnt.gen and veloc.dat"<<endl;
        return 1;
    }

    infile.open("./geo_stp.xyz");
    ofile1.open("./bnnt.gen");
    ofile2.open("./veloc.dat");
    
  // Check the last frame of last run
    do {
        getline(infile,line);
    }while(line!="MD iter: 4000");
  // Read in the coordinates and velocities
    while(infile){
        getline(infile, line);
        stringstream ss(line);
        ss >> element[num] >> x[num] >> y[num] >> z[num] >> vx[num] >> vy[num] >> vz[num];
        ++num;
    };
    --num; // This is a correction of the total number of atoms. (some side effect of while loop)
  // Calculate the center of mass
    int num_B=0,num_N=0,num_H=0;
    num_final=num;
    for (i=0;i<num;i++){
      // Check the isolated atoms; make a mark on the data and correct the actual atom number
        if ((x[i]*x[i]+y[i]*y[i]+z[i]*z[i])>1300) {
            element[i]='X';
            --num_final;
        }
        switch (element[i]){
            case 'B':
                cmx+=10.811*x[i];
                cmy+=10.811*y[i];
                cmz+=10.811*z[i];
                num_B++;
                break;
            case 'N':
                cmx+=14.007*x[i];
                cmy+=14.007*y[i];
                cmz+=14.007*z[i];
                num_N++;
                break;
            case 'H':
                cmx+=1.008*x[i];
                cmy+=1.008*y[i];
                cmz+=1.008*z[i];
                num_H++;
                break;
            default:
                break;
        }
    }
    cmx=cmx/(num_B*10.811+num_N*14.007+num_H*1.008);
    cmy=cmy/(num_B*10.811+num_N*14.007+num_H*1.008);
    cmz=cmz/(num_B*10.811+num_N*14.007+num_H*1.008);
  // Make a shift of the stucture to make sure cm is at (0,0,0)
    for (i=0;i<num;i++){
        x[i]-=cmx;
        y[i]-=cmy;
        z[i]-=cmz;
    }
   
  // Generate an N atom with equal distribution 
    theta=(atoi(argv[1])%11+1)*M_PI/12;
    phi=(atoi(argv[1])%5)*2*M_PI/5;
    theta2=(atoi(argv[1])%7+1)*M_PI/8;
    phi2=(atoi(argv[1])%13)*2*M_PI/13;
    vel_HBNH=3.0488e-1*sqrt(2000);  // The velocity corresponds to 2000K
    
  // format the output file
    ofile1 << num_final+4 << "  " << "C" << endl; // Add a 2-atom molecule
    ofile1 << "B  N  H" << endl;
    for (i=0; i<num; i++) {
        switch (element[i]) {
            case 'B':
                ofile1 << setw(4) << i+1 << setw(4)<< "1" << "\t" << setw(12) << x[i] << "\t" << setw(12) << y[i] << "\t" << setw(12) << z[i] << endl;
                ofile2 << setw(12) << vx[i] << "\t" << setw(12) << vy[i] << "\t" << setw(12) << vz[i] <<endl;
                break;
            case 'N':
                ofile1 << setw(4) << i+1 << setw(4)<< "2" << "\t" << setw(12) << x[i] << "\t" << setw(12) << y[i] << "\t" << setw(12) << z[i] << endl;
                ofile2 << setw(12) << vx[i] << "\t" << setw(12) << vy[i] << "\t" << setw(12) << vz[i] <<endl;
                break;
            case 'H':
                ofile1 << setw(4) << i+1 << setw(4)<< "3" << "\t" << setw(12) << x[i] << "\t" << setw(12) << y[i] << "\t" << setw(12) << z[i] << endl;
                ofile2 << setw(12) << vx[i] << "\t" << setw(12) << vy[i] << "\t" << setw(12) << vz[i] <<endl;
                break;
            default:
                break;
        }
    }
    ofile1 << setw(4) << num_final+1 << setw(4)<< "1" << "\t" << setw(12) << 35*sin(theta)*cos(phi) << "\t" << setw(12) << 35*sin(theta)*sin(phi) << "\t" << setw(12) << 35*cos(theta) << endl;
    ofile1 << setw(4) << num_final+2 << setw(4)<< "2" << "\t" << setw(12) << 35*sin(theta)*cos(phi)+1.2509*sin(theta2)*cos(phi2) << "\t" << setw(12) << 35*sin(theta)*sin(phi)+1.2509*sin(theta2)*sin(phi2) << "\t" << setw(12) << 35*cos(theta)+1.2509*cos(theta2) << endl;
    ofile1 << setw(4) << num_final+3 << setw(4)<< "3" << "\t" << setw(12) << 35*sin(theta)*cos(phi)-1.1982*sin(theta2)*cos(phi2) << "\t" << setw(12) << 35*sin(theta)*sin(phi)-1.1982*sin(theta2)*sin(phi2) << "\t" << setw(12) << 35*cos(theta)-1.1982*cos(theta2) << endl;
    ofile1 << setw(4) << num_final+4 << setw(4)<< "3" << "\t" << setw(12) << 35*sin(theta)*cos(phi)+2.2905*sin(theta2)*cos(phi2) << "\t" << setw(12) << 35*sin(theta)*sin(phi)+2.2905*sin(theta2)*sin(phi2) << "\t" << setw(12) << 35*cos(theta)+2.2905*cos(theta2) << endl;
    ofile2 << setw(12) << -vel_HBNH*sin(theta)*cos(phi) << "\t" << setw(12) << -vel_HBNH*sin(theta)*sin(phi) << "\t" << setw(12) << -vel_HBNH*cos(theta) <<endl;
    ofile2 << setw(12) << -vel_HBNH*sin(theta)*cos(phi) << "\t" << setw(12) << -vel_HBNH*sin(theta)*sin(phi) << "\t" << setw(12) << -vel_HBNH*cos(theta) <<endl;
    ofile2 << setw(12) << -vel_HBNH*sin(theta)*cos(phi) << "\t" << setw(12) << -vel_HBNH*sin(theta)*sin(phi) << "\t" << setw(12) << -vel_HBNH*cos(theta) <<endl;
    ofile2 << setw(12) << -vel_HBNH*sin(theta)*cos(phi) << "\t" << setw(12) << -vel_HBNH*sin(theta)*sin(phi) << "\t" << setw(12) << -vel_HBNH*cos(theta) <<endl;
    
    infile.close();
    ofile1.close();
    ofile2.close();
    return 0;
}
