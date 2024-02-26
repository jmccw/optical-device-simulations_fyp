/*
 * InGaAlAs_BG.cxx
 * 
 * Copyright 2024 jmcc0 <jmcc0@DESKTOP-65IEBIH>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
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

double InGaAlAs_BG(double x, double y){
	double result = 0.36 + 2.093*y + 0.629*x + 0.577*y*y + 0.436*x*x + 1.013*x*y - 2.0*x*y*(1-x-y);
	return result;
}

int main(int argc, char **argv)
{
	while(1){
		double x,y;
		cout << "x > ";
		cin >> x;
		cout << "y > ";
		cin >> y;
		cout << "Bandgap: " << InGaAlAs_BG(x,y);
	}
	
	return 0;
}

