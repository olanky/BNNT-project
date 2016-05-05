#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

const int N = 500;	// Maximum number of possible atoms
const int HITS=187;	// # of hits
const int FRAME=20; // # of frames for each hit (Don't include the frame at T=0)
const int SKIP_FRAME=10;
const int CLUSTER_CUTOFF=10;	// Cutoff radius of cluster [A]

int main(int argc, const char* argv[])
{
    string line;
    ifstream infile;
    ofstream ofile;
    
    char element[N];
    int i,num=0,num_final,num_r,num_atom[HITS+1];
    double x[N],y[N],z[N],lindemann[HITS+1];
    
    infile.open("./traj.xyz");		// Set the name of output (xyz) file
    ofile.open("./lindemann.txt");
    
    for(int hit=1;hit<=HITS;hit++){
		// Skip the first few frame where atoms are impacting
		for(i=0;i<=SKIP_FRAME;i++){
			getline(infile,line);
			stringstream ss(line);
			ss >> num;
			getline(infile,line);
			for (int i=0;i<num;i++){
				getline(infile, line);
			}
		}
		double rij,ave_rij[N*(N-1)/2]={0},ave_rij_sq[N*(N-1)/2]={0};
		for(i=SKIP_FRAME+1;i<=FRAME;i++){
			// Read in data
			getline(infile,line);
			stringstream ss(line);
			ss >> num;				// read in number of atoms in this frame
			num_final=num;
			getline(infile,line);
			// Within each frame
            int num_B=0,num_N=0,num_H=0;
            double cmx=0,cmy=0,cmz=0;
			for (int i=0;i<num;i++){
				getline(infile, line);
				stringstream ss(line);
				ss >> element[i] >> x[i] >> y[i] >> z[i];
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
            for (int i=0;i<num;i++){
                x[i]-=cmx;
                y[i]-=cmy;
                z[i]-=cmz;
				if ((x[i]*x[i]+y[i]*y[i]+z[i]*z[i])>CLUSTER_CUTOFF*CLUSTER_CUTOFF) {
				element[i]='X'; 
				--num_final;
				}
			}
			// Average rij calculation
			num_r=0;
			for (int i=0;i<num-1;i++){
				for (int j=i+1,ind;j<num;j++,num_r++){
					ind=(element[i]!='X'&&element[j]!='X');
					rij=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
					ave_rij[num_r-1]=ind*(ave_rij[num_r-1]+rij);
					ave_rij_sq[num_r-1]=ind*(ave_rij_sq[num_r-1]+rij*rij);
				}
			}
		}
		for(i=0;i<num_r;i++){
			if(ave_rij[i]==0||ave_rij_sq[i]==0) continue;
			ave_rij[i]/=(FRAME-SKIP_FRAME);
			ave_rij_sq[i]/=(FRAME-SKIP_FRAME);
			lindemann[hit]+=sqrt(fabs(ave_rij_sq[i]-ave_rij[i]*ave_rij[i]))/ave_rij[i];
		}
		lindemann[hit]*=2.0/(num_final*(num_final-1));
		num_atom[hit]=num_final;
	}
	
	// Output
	ofile << "Lindemann index:" << endl;
	ofile << "# of_hits" << " # of_atom" << "   Index" << endl;
	for (i=1;i<=HITS;i++){
		ofile << setw(9) << i << setw(10) << num_atom[i] << setw(8) << setprecision(4) << lindemann[i] << endl;
	}
	ofile << endl;
	
	infile.close();
    ofile.close();
    return 0;
}
	
					
					
			
			
