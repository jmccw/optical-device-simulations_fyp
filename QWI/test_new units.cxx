/*
 * FYP_CXX.png
 * PY4115 MAJOR RESEARCH PROJECT
 * Copyright 2024 jmccw <jmcc0@DESKTOP-65IEBIH>
 */

#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <cmath>     // for std::abs
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <complex.h>
#include <vector>
#include <Eigen/Core>
//#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
//#include <unsupported/Eigen/FFT>
using namespace Eigen;
#include <algorithm>
#include <fstream>
#include "gnuplot.cxx"
//#include <fftw3.h>
#define complex _Complex
const double pi = M_PI;

//Global Variables
bool over_ride_offsets = false;
vector<double> CBO_override, VBO_override;
double electron_charge = 1.0; //we use eV, so just let this = 1.
double electric_field = 0.0;
double mass_electron = 1.0; //9.109534e-31; // 1.0 # [kg]
//double hbar = 6.582119569e-16; // [eV s] source: wikipedia
double h_c = 1240; // eV nm
double hbar_sqaured_2m = 3.81; // eV Amstrong squared

const int number_steps = 500;
const int max_matrix = 500;
const int ND = 2048;

//Global Vectors - fast solving solution

// Input (destroyed during calculation of eigenvalues and eigenvectors)
double trimatrix_diag [max_matrix];
double trimatrix_subdiag [max_matrix];
// Output
double trimatrix_eigenvalue [max_matrix];
double trimatrix_eigenvector [max_matrix][max_matrix];
// only used during calculation of eigenvalues and eigenvectors
double trimatrix_result[max_matrix][max_matrix];
double* pointer_matrix [max_matrix];

class Material {
	private:
		string name;
		double electron_affinity; // eV
		double band_gap; // eV
		double e_eff_mass; // kg
		double lh_eff_mass; // kg
		double hh_eff_mass; // kg

	public:
		// Constructor
		Material(string name_, double electron_affinity_, double band_gap_, double e_eff_mass_, double lh_eff_mass_, double hh_eff_mass_) : name(name_), electron_affinity(electron_affinity_), band_gap(band_gap_), e_eff_mass(e_eff_mass_), lh_eff_mass(lh_eff_mass_), hh_eff_mass(hh_eff_mass_) {}

		// Getter methods
		string getName() const { return name; }
		double getBG() const { return band_gap; } //get band gap // eV
		double getEF() const { return electron_affinity; } //get electron affinity // eV
		double getEffectiveMass(int p) const {
			if (p == 0){ // electron
				return e_eff_mass;
			} 
			if (p == 1){ // light hole
				return lh_eff_mass;
			} 
			if (p == 2){ // heavy hole
				return hh_eff_mass;
			} 
			else return 0.0;
		}
		
		void display() const {
			cout << "Material: " << name << std::endl;
			cout << "Electron Affinity: " << electron_affinity << " eV" << endl;
			cout << "Band Gap: " << band_gap << " eV" << endl;
			cout << "Effective Masses" << endl;
			for(int i=0; i<=2; i++){
				cout << i << " : " << this->getEffectiveMass(i) << " kg" << endl;
			}
		}
};

struct Layer {
	private:
		Material material;
		double thickness; // in Amstrong
		
	public:
		// Constructor 
		Layer(Material material_, double thickness_) : material(material_), thickness(thickness_) {}
		
		double getThickness() const { return thickness; }
		Material getMaterial() const { return material; }
		
		void display() {
			cout << "Thickness: " << this->getThickness() << " A" << endl;
			cout << "Material: " << this->getMaterial().getName() << endl;
		}
};

//Default band offset calculations | These functions utilise ^ Materials ^ and are utilised by Heterostructure (next class)
double CBO(Material material1, Material material2){ //conduction band offset default
	return material1.getEF()-material2.getEF();
}
double VBO(Material material1, Material material2){ //valence band offset default
	double calc1 = material1.getEF() + material1.getBG();
    double calc2 = material2.getEF() + material2.getBG();
    return calc2 - calc1;
}

class Heterostructure {
	private:
		vector<Layer> layers;
		double heterostructure_thickness;
	public:
		// Constructor
		Heterostructure(vector<Layer> layers_) : layers(layers_) {
			double total_thickness = 0;
			for(Layer layer : layers) {
				total_thickness += layer.getThickness();
			}
			this->heterostructure_thickness = total_thickness;
		}
		
		vector<Layer> getLayers() const {
			return this->layers;
		}
		
		// Getters / Setters
		double getThickness() const {
			return this->heterostructure_thickness;
		}
		
		double eff_mass(int particle, double x) const { //
			double material_threshold = 0.0;
			for(Layer layer : this->layers) {
				material_threshold += layer.getThickness();		//Material threshold determines whether or not we need to add a band offset to the potential at x
				if(x <= material_threshold) {			
					return layer.getMaterial().getEffectiveMass(particle);
				}
			}
			throw std::runtime_error("Error: Unable to find effective mass for the given position."+to_string(x));
			return 0.0; //this will never happen under normal circumstances
		}
		
		double potential(int particle, double x) { // particle = 0,1,2 -> elecron, lh, hh
			double material_threshold = 0.0;
			double U = 0.0, V = 0.0;
			int i = 0;
			for(Layer layer : this->layers) {
				material_threshold += layer.getThickness();		//Material threshold determines whether or not we need to add a band offset to the potential at x
				if(x >= material_threshold) {
					if(over_ride_offsets == false) {				// so if x is not within material of layer[i] add/subtract relevant band offset
						U += CBO(this->layers[i].getMaterial(), this->layers[i+1].getMaterial());
						V += VBO(this->layers[i].getMaterial(), this->layers[i+1].getMaterial());						
					} else if (over_ride_offsets == true) {
						U += CBO_override[i]; //TEMPORARY SOLUTION ??
						V += VBO_override[i]; //TEMPORARY SOLUTION ??
					}
					i++;
				}
				else { // x within new material threshold, return relevant potential - if used correctly, there should never be an error here (one might occur when an x input is out of analysis bounds)
					if (particle == 1 || particle == 2){ // if particle is hole
						return V - electron_charge*electric_field*(x-this->getThickness()); // not as efficient as defining a strict function globally but gives ALOT of freedom (only adds +(number of layers)*3 computations at most  
					}
					else if (particle == 0){             // else it is electron
						return U + electron_charge*electric_field*x;
					}
				}
			}
			throw std::runtime_error("Error: Unable to find valid potential for the given position."+to_string(x));
			return 0.0; //this will never happen under normal circumstances
		}
		
		void display() { // print heterostructure details
			cout << "Layers: " << endl;
			int i = 1;
			for(Layer layer : layers){
				cout << i << " : " << layer.getMaterial().getName() << " : " << layer.getThickness() << " A"<< endl;
				i++;
			}	
			cout << "Total Thickness : " << this->getThickness() << " A" << endl;		
		}
};

double relative_energy(double energy, int p, Heterostructure& QW) { // accept a 2 paramters [energy relative to well of..][particle type]
    // corrects an (inital) calculated energy from QW solution for respective band of type particle type p
    double EF_offset = electron_charge*electric_field*QW.getThickness();
    double E_REL;
    
    // this should be okay with and without override on
    double BG = abs(QW.getLayers()[0].getMaterial().getBG());
	for (int i = 1; i < (int)QW.getLayers().size() - 1; ++i) {
		double band_gap = abs(QW.getLayers()[i].getMaterial().getBG());
		BG = min(BG, band_gap);
	}
	
    if(p==0){
        if (over_ride_offsets == true) {
			double max_abs_CBO_override = std::abs(CBO_override[0]);
			for (int i = 1; i < (int)QW.getLayers().size() - 1; ++i) {
				double abs_CBO_override = std::abs(CBO_override[i]);
				max_abs_CBO_override = std::max(max_abs_CBO_override, abs_CBO_override);
			}
            E_REL = energy + BG + max_abs_CBO_override;
        } else {
			double max_abs_CBO = std::abs(CBO(QW.getLayers()[0].getMaterial(), QW.getLayers()[1].getMaterial()));
			for (int i = 1; i < (int)QW.getLayers().size() - 1; ++i) {
				double abs_CBO = std::abs(CBO(QW.getLayers()[i].getMaterial(), QW.getLayers()[i+1].getMaterial()));
				max_abs_CBO = std::max(max_abs_CBO, abs_CBO);
			}
            E_REL = energy + BG + max_abs_CBO;
		}
    } else { // i.e. p==1 or p==2
        if (over_ride_offsets == true) {
			double max_abs_neg_VBO_override = std::abs(-VBO_override[0]);
			for (int i = 1; i < (int)QW.getLayers().size() - 1; ++i) {
				double abs_neg_VBO_override = std::abs(-VBO_override[i]);
				max_abs_neg_VBO_override = std::max(max_abs_neg_VBO_override, abs_neg_VBO_override);
			}
            E_REL = -max_abs_neg_VBO_override - energy + EF_offset;
        } else {
			double max_abs_VBO = std::abs(VBO(QW.getLayers()[0].getMaterial(), QW.getLayers()[1].getMaterial()));
			for (int i = 1; i < (int)QW.getLayers().size() - 1; ++i) {
				double abs_VBO = std::abs(VBO(QW.getLayers()[i].getMaterial(), QW.getLayers()[i+1].getMaterial()));
				max_abs_VBO = std::max(max_abs_VBO, abs_VBO);
			}
            E_REL = -max_abs_VBO - energy + EF_offset;
		}
	}
    return E_REL;
}

// Define a function to solve the heterostructure
void solve(Heterostructure& heterostructure, std::vector<double>& x_out,
           std::vector<std::vector<double>>& energies,
           std::vector<std::vector<std::vector<double>>>& eigenVectors) {
			    
	double delta_x = heterostructure.getThickness() / (number_steps + 1.0); 
	double x[number_steps];
	
	x[0] = delta_x; 
	x_out[0] = x[0]; 
	
	for(int i = 1; i<number_steps; i++){ //initialise x (once)
		x[i] = x[0] + i*delta_x;
		x_out[i] = x[i];
	}

	//implementation of Kowano & Kito : 'Optical Waveguide Analysis' : solution of SE with effective mass approximation for all bands/particles
    for(int p = 0; p<=2; p++){ //for all particle types
		//initialise solution matrix M
		double M[number_steps][number_steps];
		for(int i = 0; i<number_steps; i++){
			for(int j = 0; j<number_steps; j++){
				M[i][j] = 0;
			}
		}
		
		double alpha_w[number_steps], alpha_e[number_steps], alpha_x[number_steps];
		alpha_w[0] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[0])+heterostructure.eff_mass(p,x[0]));
        alpha_e[0] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[0])+heterostructure.eff_mass(p,x[1]));
        alpha_x[0] = -alpha_w[0]-alpha_e[0];
        
        M[0][0] = alpha_x[0] + heterostructure.potential(p,x[0]);
        M[0][1] = alpha_e[0];
        
        for(int nr = 1; nr < number_steps-1; nr++){
            alpha_w[nr] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[nr])+heterostructure.eff_mass(p,x[nr-1]));
            alpha_e[nr] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[nr])+heterostructure.eff_mass(p,x[nr+1]));
            alpha_x[nr] = -alpha_w[nr]-alpha_e[nr];

            M[nr][nr-1] = alpha_w[nr];    //sub-diagonal
            M[nr][nr] = alpha_x[nr] + heterostructure.potential(p,x[nr]); //diagonal
            M[nr][nr+1] = alpha_e[nr];   //upper diagonal   
		}
		cout << "test0 ";
        alpha_w[number_steps-1] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[number_steps-1])+heterostructure.eff_mass(p,x[number_steps-1-1]));
        alpha_e[number_steps-1] = -hbar_sqaured_2m * 1.0/(delta_x*delta_x) * 2.0/(heterostructure.eff_mass(p,x[number_steps-1])+heterostructure.eff_mass(p,x[number_steps-1])); // assuming m(x_edge-dx) = m(x_edge) as boundary condition
        alpha_x[number_steps-1] = -alpha_w[number_steps-1]-alpha_e[number_steps-1];
        M[number_steps-1][number_steps-2] = alpha_w[number_steps-1];
        M[number_steps-1][number_steps-1] = alpha_x[number_steps-1] + heterostructure.potential(p,x[number_steps-1]);
		
		cout << "test0 ";
		//solve Matrix (using Eigen)
		Eigen::MatrixXd M_eigen(number_steps, number_steps);
		cout << "test1 ";
		for (int i = 0; i < number_steps; ++i) { // really not a speed issue here
			for (int j = 0; j < number_steps; ++j) {
				M_eigen(i, j) = M[i][j];
			}
		}
		cout << "test2 ";
		
		// Solve the matrix using Eigen's EigenSolver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M_eigen);
		// Retrieve the eigenvalues and eigenvectors
		//~ cout << "TEST1";
		//~ Eigen::VectorXcd eigenvalues_com = solver.eigenvalues();
		//~ cout << "TEST2";
		//~ Eigen::MatrixXcd eigenvectors_com = solver.eigenvectors();
		//~ cout << "TEST3";
		
		// Convert eigenvalues and eigenvectors to real-valued types
		Eigen::VectorXd eigenvalues = solver.eigenvalues();
		cout << "TEST4";
		Eigen::MatrixXd eigenvectors = solver.eigenvectors();
		cout << "TEST5";

		// Print the eigenvalues
		cout << "test3 ";
		//double test[number_steps];
		for (int i = 0; i < number_steps; ++i) {
			// Assign eigenvalue to energies
			energies[p][i] = eigenvalues(i);
			// Assign eigenvector to eigenVectors
			for (int j = 0; j < number_steps; ++j) {
				eigenVectors[p][i][j] = eigenvectors(j, i);
			}
		}
		
		//~ double x_array_2[number_steps], y_array_2[number_steps];
		//~ for (int i = 0; i < number_steps; i++) {
			//~ y_array_2[i] = eigenVectors[p][0][i];
			//~ x_array_2[i] = x[i];
		//~ }
		//~ cout << "VOIGT RESULT:\n";
		//~ gnuplot_two_functions ("Numerical Eigenfunctions", "linespoints", "x", "Eigenfunction |E_1> (x)",
			//~ x_array_2, y_array_2, number_steps, "v_0 = 10.00", x_array_2, y_array_2, number_steps, "v_0 = 400.00");

		//for (int i = 0; i < (int)energies[0].size(); i++) cout << eigenvectors(i, 0) << endl;
		cout << "test4 " << endl;
	}
	//~ double x_array_2[number_steps], y_array_1[number_steps], y_array_2[number_steps];
	//~ for (int i = 0; i < number_steps; i++) {
		//~ y_array_2[i] = eigenVectors[0][5][i];
		//~ y_array_1[i] = eigenVectors[0][0][i];
		//~ x_array_2[i] = x[i];
	//~ }
	//~ cout << "VOIGT RESULT:\n";
	//~ gnuplot_two_functions ("Numerical Eigenfunctions"+to_string(electric_field), "linespoints", "x", "Eigenfunction |E_1> (x)",
		//~ x_array_2, y_array_1, number_steps, "Ground state", x_array_2, y_array_2, number_steps, "5th excited");

}

vector<vector<vector<double>>> findTransitions(vector<vector<double>> energies){
	vector<vector<vector<double>>> E_gap(number_steps, vector<vector<double>>(number_steps, vector<double>(2))); //Gap between [electron state] and [hole state] of [hole type]
	for (int k=0; k<2; k++){ // for 2 hole types
        for (int i=0; i<(int)energies[0].size(); i++){ // for all electrons
            for (int j=0; j<(int)energies[0].size(); j++){ // and all hole states
                E_gap[i][j][k] = abs(energies[0][i]-energies[k+1][j]);
			}
		}
	}
    return E_gap; // a matrix/array of energy transtions indexed as [electron state][hole state][hole type] 
}

vector<vector<double>> findEnergiesRelative(vector<vector<double>> energies, Heterostructure& heterostructure) {
	vector<vector<double>> energies_relative(3, vector<double>(number_steps));
	for (int k=0; k<3; k++){ // for 3 particle types
        for (int i=0; i<(int)energies[0].size(); i++){ //for all states
            energies_relative[k][i] = relative_energy(energies[k][i], k, heterostructure);
		}
	}
	return energies_relative; //returns [3][number_solutions] matrix of energies relative to their ectual value in well structure
}

vector<double> overlapIntegral(vector<double> vector1, vector<double> vector2){
	if((int)vector1.size() != (int)vector2.size()) throw std::runtime_error("vector sizes not equal, "+to_string((int)vector1.size())+" and "+to_string((int)vector2.size()));
	vector<double> overlap((int)vector1.size()), vector1_dummy((int)vector1.size()), vector2_dummy((int)vector2.size());
	double N1 = 0.0, N2 = 0.0;
	// possibly might need to declare some dummy vectors ?? - this was a python issue
	for(int i = 0; i < (int)overlap.size(); i++){
		N1 += fabs(vector1[i])*fabs(vector1[i]);
		N2 += fabs(vector2[i])*fabs(vector2[i]);
	}
	for(int i = 0; i < (int)overlap.size(); i++){
		vector1_dummy[i] = vector1[i]/N1;
		vector2_dummy[i] = vector2[i]/N2;
		overlap[i] = vector1_dummy[i]*vector2_dummy[i];
	}
	return overlap;
}

double I_squared(vector<double> vector1, vector<double> vector2){
	vector<double> overlap = overlapIntegral(vector1, vector2);
	double I_squared = 0;
	for(int i = 0; i < (int)overlap.size(); i++) I_squared += fabs(overlap[i]);
	I_squared *= I_squared; // square result
	return I_squared;
}

vector<vector<vector<double>>> findOverlapsAll(vector<vector<vector<double>>> wavefunctions) {
	vector<vector<vector<double>>> I_squared_matrix(number_steps, vector<vector<double>>(number_steps, vector<double>(2)));
	// [electron state][hole state][hole type]
	for (int k=0; k<2; k++) { // for 2 hole types
        for (int i=0; i<(int)wavefunctions[0].size(); i++) { // for all electrons
			vector<double> state1 = wavefunctions[0][i];
            for (int j=0; j<(int)wavefunctions[0].size(); j++) { // and all hole states
				vector<double> state2 = wavefunctions[1+k][j];
                I_squared_matrix[i][j][k] = I_squared(state1, state2);
			}
		}
	}
	return I_squared_matrix;
}

// - //

vector<double> pad_func_zeros(vector<double> func){
	vector<double> func_padded(2*func.size());
	for(int i = 0; i<(int)(0.25*func_padded.size()); i++) func_padded[i] = 0.0;
	for(int i = (int)(0.75*func_padded.size()); i<(int)(func_padded.size()); i++) func_padded[i] = 0.0;
	int j = 0;
	for(int i = (int)(0.25*func_padded.size()); i<(int)(0.75*func_padded.size()); i++){
		func_padded[i] = func[j];
		j++;
	}
	return func_padded;
}

std::vector<double> shiftEdgeToCenter(const std::vector<double>& input) {
    std::vector<double> shiftedVector = input;
    int halfSize = input.size() / 2;
    
    // Calculate the number of rotations needed to center the edge values
    int rotations = halfSize % input.size();

    // Perform cyclic rotations
    std::rotate(shiftedVector.begin(), shiftedVector.begin() + rotations, shiftedVector.end());

    return shiftedVector;
}

// - //

vector<vector<vector<double>>> wavelengthTransformation(vector<vector<vector<double>>> data_in) {
	//The intention here is that the user passes in E_GAP to find an adjacent matrix in wavelength terms.
	vector<vector<vector<double>>> data_out(number_steps, vector<vector<double>>(number_steps, vector<double>(2)));
	for (std::size_t i = 0; i < data_out.size(); ++i) {
		for (std::size_t j = 0; j < data_out[i].size(); ++j) {
			for (std::size_t k = 0; k < data_out[i][j].size(); ++k) {
				data_out[i][j][k] = h_c / data_in[i][j][k]; //h_c = 1240 eV nm :: so [eV nm / eV] = [nm]
			}
		}
	}
	return data_out;
}

// Ioffe Materials
double x_ratio = 0.3;
double BG_GaAs = 1.424; // [eV]
double EF_GaAs = 4.07; // [eV] @ 300K !
double EF_AlGaAs = 4.07 - 1.1*x_ratio; // [eV] (x<0.45) @ 300K !
double BG_AlGaAs = 1.424 + 1.247*x_ratio; // [eV] (x<0.45)

// Decleration: Material(EF, BG, e_eff_mass, lh_eff_mass, hh_eff_mass) @ 300K
Material InAs("InAs", 4.9, 0.354, 0.023*mass_electron, 0.026*mass_electron, 0.41*mass_electron);
Material GaInAs("GaInAs", 4.5, (0.36+0.63*x_ratio+0.43*x_ratio*x_ratio), 0.041*mass_electron, 0.052*mass_electron, 0.45*mass_electron);
Material GaAs("GaAs", EF_GaAs, BG_GaAs, 0.063*mass_electron, 0.082*mass_electron, 0.51*mass_electron);
Material AlGaAs("AlGaAs", EF_AlGaAs, BG_AlGaAs, (0.063+x_ratio*0.083)*mass_electron, (0.082+x_ratio*0.068)*mass_electron, (0.51+x_ratio*0.25)*mass_electron);
Material InP("InP", 4.38, 1.344, 0.077*mass_electron, 0.14*mass_electron, 0.6*mass_electron);
Material AlAs("AlAs", 0.0, 2.12, 0.15*mass_electron, 0.16*mass_electron, 0.79*mass_electron);


// Alloy Functions
double InGaAlAs_BG(double x, double y){
	double result = 0.36 + 2.093*y + 0.629*x + 0.577*y*y + 0.436*x*x + 1.013*x*y - 2.0*x*y*(1-x-y);
	return result;
}
	// P(In{1-x-y}Ga{x}Al{y}As)=P(InAs{1-x-y})+P(GaAs{x})+P(AlAs{y}) where P is any material parameter
		// using this for effective mass
double eff_mass_InGaAlAs(int particle, double x, double y){
	double result = InAs.getEffectiveMass(particle)*(1-x-y) + GaAs.getEffectiveMass(particle)*(x) + AlAs.getEffectiveMass(particle)*(y);
	return result;
}

// AlGaInAs Alloys
Material AlGaInAs("Al_y Ga_x In_(1-x-y) (Ref: Mondry Babic) BG@1.4K", 0.0, 1.140, 0.24*mass_electron, 0.265*mass_electron, 1.69*mass_electron);
double x = 0.124, y = 0.245;
Material AlGaInAs_two("Al_y Ga_x In_(1-x-y) (Ref: Mondry Babic) BG@1.4K", 0.0, InGaAlAs_BG(x,y), eff_mass_InGaAlAs(0,x,y), eff_mass_InGaAlAs(1,x,y), eff_mass_InGaAlAs(2,x,y));

void plot_potential(Heterostructure& QW){
	double potential1[number_steps], potential2[number_steps];
	double x_test[number_steps];
	for(int i = 0; i<number_steps; i++){
		x_test[i] = (QW.getThickness()/(number_steps+1)) + i*(QW.getThickness()/(number_steps+1));
		potential1[i] = relative_energy(QW.potential(2, x_test[i]), 2, QW);
		potential2[i] = relative_energy(QW.potential(0, x_test[i]), 0, QW);
	}
    gnuplot_two_functions ("Potentials", "linespoints", "x [nm]", "Effective Potentials",
			x_test, potential1, number_steps, "Conduction", x_test, potential2, number_steps, "Valence");
}

void script(){
	//double length = 30.0;
	cout << InGaAlAs_BG(0.24,0.2) << endl;
	cout << InGaAlAs_BG(0.24,0.42) << endl;
	cout << InGaAlAs_BG(0.2412,0.01) << endl;
	AlGaInAs_two.display();
	cout << "Layer Thicknesses: " << endl;
	cout << "1 (well) [A] > ";
	double t_1, t_2;
	cin >> t_1;
	cout << "2 (barrier) [A] > ";
	cin >> t_2;
	//Layer layer2(GaAs, t_1); // Material, latyer thickness [A]
	//Layer layer1(AlGaAs, t_2); // Material, latyer thickness [A]
	
	//Layer layer2(GaAs, 100); // Material, latyer thickness [A]
	//Layer layer1(AlGaAs, 100); // Material, latyer thickness [A]
	
	Layer layer2(AlGaInAs_two, t_1); // Material, latyer thickness [A]
	Layer layer1(InP, t_2); // Material, latyer thickness [A]
	
	layer1.display();
	layer2.display();
	
	vector<Layer> layers;
	layers = {layer1, layer2, layer1};
	electric_field = 0;
	
	cout << "Set a linear value to iterate the electric field.\ndelta_field [V/um] = ";
	double delta_field;
	cin >> delta_field;
	delta_field *= 0.0001; //conversion to V/Amstrong
	
	over_ride_offsets = true;
	double E_gap_diff = fabs(-layer2.getMaterial().getBG() + layer1.getMaterial().getBG());
	double conduction_band_offset_ratio = 0.7;
	double valence_band_offset_ratio = 0.3;
	VBO_override = {-E_gap_diff*valence_band_offset_ratio, E_gap_diff*valence_band_offset_ratio};
	CBO_override = {-E_gap_diff*conduction_band_offset_ratio, E_gap_diff*conduction_band_offset_ratio};
	
	Heterostructure QW(layers);
	QW.display();
	plot_potential(QW);
		
	for(int i = 0; i < 11; i++) {
		electric_field = i*delta_field;
		cout << "E-Field = " << electric_field << endl;
		cout << "Run " << i << endl;
		
		vector<vector<vector<double>>> RESULT; //ISSUE HERE
		
		//cout << "test ";
		
		// solve QW with Electric field
		vector<double> x(number_steps); //length element
		vector<vector<double>> energies(3, vector<double>(number_steps)); //particle, energy_level
		vector<vector<vector<double>>> eigenVectors(3, vector<vector<double>>(number_steps, vector<double>(number_steps)));
		cout << "test1 ";
		solve(QW, x, energies, eigenVectors);
		if (i==0){
			//PLOTTING FINAL
			//plot_potential(QW);
			double x_array_2[number_steps], y_array_2[number_steps], y_array_1[number_steps];
			for (int i = 0; i < number_steps; i++) {
				y_array_2[i] = eigenVectors[0][0][i];
				y_array_1[i] = eigenVectors[0][5][i];
				x_array_2[i] = x[i];
			}
			gnuplot_two_functions ("Numerical Eigenfunctions", "linespoints", "x", "Eigenfunction |E_1> (x)",
				x_array_2, y_array_2, number_steps, "Ground", x_array_2, y_array_1, number_steps, "5th excited state");

		}

		//plot_potential(QW);
		// sort results
		vector<vector<double>> energies_relative(3, vector<double>(number_steps));
		energies_relative = findEnergiesRelative(energies, QW);
		// Absorption Routine
		cout << "Calculating overlaps.\n";
		vector<vector<vector<double>>> I_squared_matrix = findOverlapsAll(eigenVectors);
		cout << "Sorting transitions.\n";
		vector<vector<vector<double>>> E_GAP = findTransitions(energies_relative); //Gap between [electron state] and [hole state] of [hole type]
		cout << "Converting to wavelength.\n";
		vector<vector<vector<double>>> E_GAP_WL = wavelengthTransformation(E_GAP);
		
		// Write I_squared_matrix to a file
		ofstream I_squared_file("I_squared_matrix_"+to_string(i)+".txt");
		if(I_squared_file.is_open()){
			for(const auto& row : I_squared_matrix){
				for(const auto& col : row){
					for(double val : col){
						I_squared_file << val << " ";
					}
					I_squared_file << "\n";
				}
				I_squared_file << "\n";
			}
			I_squared_file.close();
		}
		else{
			cerr << "Unable to open I_squared_matrix file for writing.\n";
		}

		// Write E_GAP_WL to a file
		ofstream E_GAP_WL_file("E_GAP_WL_"+to_string(i)+".txt");
		if(E_GAP_WL_file.is_open()){
			for(const auto& row : E_GAP_WL){
				for(const auto& col : row){
					for(double val : col){
						E_GAP_WL_file << val << " ";
					}
					E_GAP_WL_file << "\n";
				}
				E_GAP_WL_file << "\n";
			}
			E_GAP_WL_file.close();
		}
		else{
			cerr << "Unable to open E_GAP_WL file for writing.\n";
		}
		if (i==10){
			//PLOTTING FINAL
			plot_potential(QW);
			double x_array_2[number_steps], y_array_2[number_steps], y_array_1[number_steps];
			for (int i = 0; i < number_steps; i++) {
				y_array_2[i] = eigenVectors[0][0][i];
				y_array_1[i] = eigenVectors[0][5][i];
				x_array_2[i] = x[i];
			}
			gnuplot_two_functions ("Numerical Eigenfunctions", "linespoints", "x", "Eigenfunction |E_1> (x)",
				x_array_2, y_array_2, number_steps, "Ground", x_array_2, y_array_1, number_steps, "5th excited state");

		}
	}
}

int main(int argc, char **argv)
{
	script();
	return 0;
}
