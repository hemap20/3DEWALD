#include "libinclude.h" 
#include "potdec.h"
#include "mathdec.h"
#include "fundec.h"
#include "const.h"
#include "fftw3.h"
#include "complex" 
//modif_1
int main(int argc, char **argv){

	string filename;
	string posfile;
	string garbage, garbage1;
	char garbage2[200];
	int tmp;
	int i, j, MC, MD;
	int  nMD, md_print;
	float Temp;
	int  MDtime, MDeq;
	float dt;
	int genpairlist;
	int randnum;
	int hyb, NMCMD, rdvel;
	float tmp2;
	double unitzer=(AVG*COL*CHARELEC*CHARELEC/rtoa)/1000/4.184;

	//get the name of the file for input variables, default is input.in 
	if(argc <2)
		filename = "input.in";
	else
		filename = argv[1];

	///////////////////////////////////////////////
	//////initialize variables////////////////////

	ifstream InputIn(filename.c_str(),ios::in);
	if(!InputIn){
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	InputIn>>garbage>>garbage;
	getline(InputIn, garbage);

	cout<<"Project Name: "<<garbage<<endl;

	InputIn>>garbage>>garbage;
	InputIn>>posfile;

	InputIn>>garbage>>garbage;
	InputIn>>genpairlist;

	InputIn>>garbage>>garbage;
	InputIn>>Temp;

	// cout<<"Temperature of the system:\t\t"<<Temp<<" K"<<endl;

	InputIn>>garbage>>garbage;
	InputIn>>hyb;

	InputIn>>garbage>>garbage;
	InputIn>>NMCMD;

	InputIn>>garbage>>garbage;
	InputIn>>rdvel;

	bool comdens;

	InputIn>>garbage>>garbage;
	InputIn>>comdens;

	// cout<<"Compute Density: "<<comdens<<"\n";

	double lim1, lim2, dx;

	InputIn>>garbage>>garbage;
	InputIn>>lim1>>lim2>>dx;

	// cout<<"Range of x-values considered for computing densities: "<<lim1<<"\t"<<lim2<<endl;

	//Parameters for simulated annealing		
	bool simann;

	InputIn>>garbage>>garbage;
	InputIn>>simann;

	float TempB, TempE;

	InputIn>>garbage>>garbage;
	InputIn>>TempB>>TempE;

	int rampsize, rampstep;

	InputIn>>garbage>>garbage;
	InputIn>>rampsize>>rampstep;

	int ncycle;

	InputIn>>garbage>>garbage;
	InputIn>>ncycle;

	if(simann==1){	
		cout<<"Perform Simulated Annealing"<<endl;
		cout<<"Start Temperature: "<<TempB<<" K"<<endl;
		cout<<"End Temperature: "<<TempE<<" K"<<endl;
		cout<<"Delta T:"<<rampsize<<" K"<<endl;
		cout<<"Increase temperature after: "<<rampstep<<" MD steps"<<endl;
		cout<<"ncycles: "<<ncycle<<endl;
	}

	bool minimize;

	InputIn>>garbage>>garbage;
	InputIn>>minimize;                      //whether to minimize or not
	if(minimize)cout<<"****Perform minimization****"<<endl;

	InputIn.close();

	/////////////////////////////// ewald input file /////////

	// ifstream EWALDIn("ewald.in",ios::in);
	// if(!EWALDIn){
	// 	cerr << "File EWALDIn could not be opened" <<endl;
	// 	exit(1);
	// }

	// int gx, gy, gz, nd;

	// EWALDIn>>garbage>>garbage;
	// EWALDIn>>gx;

	// EWALDIn>>garbage>>garbage;
	// EWALDIn>>gy;                                    

	// EWALDIn>>garbage>>garbage;
	// EWALDIn>>gz;                             

	// EWALDIn>>garbage>>garbage;
	// EWALDIn>>nd;                             

	// EWALDIn.close();

	// cout<<"Initial positions read from the file: "<<posfile<<endl;

	//////////////////////get positions from the POSCAR file//////////////////////////////////////////////////////////

	ifstream PosIn(posfile.c_str(),ios::in);
	if(!PosIn){
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	getline(PosIn, garbage);
	istringstream StrStream(garbage);

	int n_atomtype=0;
	// Counting through the types of the atoms present in the cell and storing in the n_atomtype variable
	while(StrStream){
		getline(StrStream, garbage1, ' ');
		if(garbage1.compare("") != 0)
			n_atomtype = n_atomtype + 1;
	}
	
	string *atomtype;
	atomtype=new string [n_atomtype];

	istringstream strstream(garbage);

	tmp = 0;

	while(strstream){
		getline(strstream, garbage1, ' ');
		if(garbage1.compare("") != 0){
			atomtype[tmp]=garbage1;
			tmp = tmp + 1;
		}
	}
	
	int *natoms_type, natoms;  //number of atoms for each type and total atoms in the unit cell

	natoms_type=new int [n_atomtype];
	natoms=0;

	//Creating 3x3 boxcell array
	float **boxcell;
	boxcell=new float * [3];
	for(i = 0; i<3;i++){
		boxcell[i]=new float [3];
	}

	getline(PosIn, garbage); //read overall scaling factor

	//get box cell vectors
	for(i = 0; i<3;i++){
		for(j=0;j<3;j++){
			PosIn>>boxcell[i][j];
		}
	}

	//get number of atoms for each type
	for(i=0;i<n_atomtype;i++){
		PosIn>>natoms_type[i];
		natoms=natoms+natoms_type[i];
	}

	getline(PosIn, garbage);
	getline(PosIn, garbage);

	double **PosIons, **ForceIons, **vel;
	int **fixatoms;
	float *mass;

	mass=new float [natoms];
	PosIons=new double * [natoms];
	ForceIons=new double *[natoms];
	vel=new double *[natoms];
	fixatoms=new int *[natoms];

	for(i=0;i<natoms;i++){
		PosIons[i]=new double [3];
		ForceIons[i]=new double [3];
		vel[i]=new double [3];		
		fixatoms[i]=new int [3];
	}

	//getting positions of each atom
	for(i=0;i<natoms;i++){
		PosIn>>PosIons[i][0]>>PosIons[i][1]>>PosIons[i][2];
	}

	PosIn.close();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// charge input file ////////////////////////////////////////////////////////////

	ifstream CHARGEIn("charge.in",ios::in);
	if(!CHARGEIn){
		cerr << "File CHARGEIn could not be opened" <<endl;
		exit(1);
	}

	float *chg;
	chg=new float [n_atomtype];

	for(int i=0; i<n_atomtype; i++){
		CHARGEIn>>chg[i];
	}

	CHARGEIn.close();

	//creating the charge array for each atom present in the unit cell ////////////////////////////////////
	float *ion_charges;
	ion_charges=new float [natoms];
	int c=0;
	for (int i = 0; i < n_atomtype; i++){
		for (int j = 0; j < natoms_type[i]; j++){
			ion_charges[c]=chg[i];
			c++;
		}		
	}

	// print_carcoor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR");
	// print_coor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "COOR");
	// print_lammps_input_file(PosIons, chg, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "poscar7.data");
	float a=5.42/boxcell[0][0];

	// chrono::time_point<std::chrono::system_clock> start, end;
	// start = chrono::system_clock::now();
	double selfenergy=selfe(n_atomtype, natoms_type, chg, a)*unitzer;
	cout<<fixed<<setprecision(5)<<"Self Energy: "<<selfenergy<<" Kcal/mol"<<"\n\n";
	// end = chrono::system_clock::now();
	// chrono::duration<double> elapsed_seconds = end - start;
    // time_t end_time = std::chrono::system_clock::to_time_t(end);
	// cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds.count() << " sec\n\n";
	
	chrono::time_point<std::chrono::system_clock> start1, end1;
	start1 = chrono::system_clock::now();
	double recienergy=reci_energy(PosIons, ion_charges, natoms, a, boxcell,6)*unitzer;
	cout<<fixed<<setprecision(5)<<"Reciprocal Energy: "<<recienergy<<" Kcal/mol"<<"\n";
	end1 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds1 = end1- start1;
    time_t end_time1 = std::chrono::system_clock::to_time_t(end1);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds1.count() << " sec\n\n";

	chrono::time_point<std::chrono::system_clock> start2, end2;
	start2 = chrono::system_clock::now();
	double realenergy=real_energy(PosIons, ion_charges, natoms, a, boxcell)*unitzer;
	cout<<fixed<<setprecision(5)<<"Real Energy: "<<realenergy<<" Kcal/mol"<<"\n";
	end2 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds2 = end2 - start2;
    time_t end_time2 = std::chrono::system_clock::to_time_t(end2);
	cout<< "Elapsed time: " << elapsed_seconds2.count() << " sec\n\n";

	chrono::time_point<std::chrono::system_clock> start3, end3;
	start3 = chrono::system_clock::now();
	double recienergy_bs=bspline(PosIons, ion_charges, natoms, a, boxcell,GRID_SIZE,6,BSPLINE_ORDER)*unitzer;
	cout<<fixed<<setprecision(5)<<"Reciprocal Energy FFTW: "<<recienergy_bs<<" Kcal/mol"<<"\n";
	end3 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds3 = end3 - start3;
    time_t end_time3 = std::chrono::system_clock::to_time_t(end3);
	cout<< "Elapsed time: " << elapsed_seconds3.count() << " sec\n\n";

	// cout<<fixed<<setprecision(10)<<","<< error(recienergy,recienergy_bs);
	cout<<fixed<<setprecision(10)<<"Relative Error with FFTW: "<< error(recienergy,recienergy_bs)<<"\n";
/* using std::chrono::duration_cast; */
/* using HR = std::chrono::high_resolution_clock; */
/* using HRTimer = HR::time_point; */
/* using std::chrono::microseconds; */
/* using std::chrono::seconds; */

/*   HRTimer start = HR::now(); */
/*  recienergy_bs=bspline(PosIons, ion_charges, natoms, a, boxcell,60,6,5)*unitzer; */
/* 	cout<<fixed<<setprecision(5)<<"Reciprocal Energy FFTW: "<<recienergy_bs<<" Kcal/mol"<<"\n"; */
/*   HRTimer end = HR::now(); */
/*   auto duration = duration_cast<microseconds>(end - start).count(); */
/* 	cout<< "Elapsed time: " << duration << " usec"; */

	// for benchmarking using the lammps data
	// cout<<fixed<<setprecision(5)<<" "<<realenergy<<" ";
	// cout<<fixed<<setprecision(5)<<selfenergy+recienergy<<"\n";
	// cout<<fixed<<setprecision(5)<<selfenergy+recienergy+realenergy<<" ";

	// delete dynamic variables 
	for(i=0;i<3;i++){
		delete [] boxcell[i]; 
	}

	for(i=0;i<natoms;i++){
		delete [] PosIons[i];
		delete [] ForceIons[i];
		delete [] vel[i];
		delete [] fixatoms[i];
	}


	delete [] PosIons;
	delete [] boxcell;
	delete [] atomtype;
	delete [] natoms_type;
	delete [] ForceIons;
	delete [] vel;
	delete [] fixatoms;
    delete [] chg;

	return 0;
}
//end of main loop
