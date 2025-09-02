#ifndef SHOWER_MAX_RESPONSE_LOOKUP_H
#define SHOWER_MAX_RESPONSE_LOOKUP_H

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include "RtypesCore.h"
#include "TF2.h"


using std::cout, std::cerr, std::endl;
using std::vector, std::string;
using std::ifstream, std::ofstream;
using std::stringstream, std::getline;
using std::pair;

const int nParticles = 5; //(e-, gamma, mu-, pi-, neutron) 
inline string particleList[] = {"e-", "gamma", "mu-", "pi-", "neutron"};
inline const Double_t energyList[] = {5, 10, 50, 100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000};

inline std::pair<Double_t, Double_t> get_energy_bounds(Double_t energy) {
    Double_t lower = 0;
    Double_t upper = 0;
    int energyListSize = sizeof(energyList) / sizeof(energyList[0]);

    for (int i = 0; i < energyListSize; ++i) {
        if (energyList[i] <= energy) {
            lower = energyList[i];
        }
        if (energyList[i] > energy) {
            upper = energyList[i];
            break;
        }
    }

    return {lower, upper};
}

inline vector<vector<Double_t>> retrieve_fit_data(const string& particleName) {
    string csv_file_name = Form("../data/fit_param_xy_%s_ifarm.csv", particleName.c_str());
    ifstream file(csv_file_name);
    vector<vector<Double_t>> fit_params;

    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << csv_file_name << endl;
        exit(1);
    }

    string line;
    getline(file, line);
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') { continue; }

        stringstream ss(line);
        string token;
        vector<Double_t> row;

        while (getline(ss, token, ',')) {
            row.push_back(std::stod(token));
        }
        fit_params.push_back(row);
    }

    file.close();
    return fit_params;
}

inline TF2* create_fit_function(Double_t energy, const vector<vector<Double_t>>& fit_data) {
    Double_t xWeight = 0.50;
    Double_t yWeight = 0.50;
    string fit_function = "[0]*([1] + x*[2]) + [3]*([4] + y*[5] + y*y*[6])";

    TF2* fitFormula = new TF2("fitFormula", fit_function.c_str());

    for (const auto& data : fit_data) {
        if (data[0] == energy) {
            Double_t xp0 = data[1];
            Double_t xp1 = data[2];
            Double_t yp0 = data[3];
            Double_t yp1 = data[4];
            Double_t yp2 = data[5];
            Double_t tp0 = data[6];

            fitFormula->SetParameter(0, xWeight);
            fitFormula->SetParameter(1, xp0);
            fitFormula->SetParameter(2, xp1);
            fitFormula->SetParameter(3, yWeight);
            fitFormula->SetParameter(4, yp0);
            fitFormula->SetParameter(5, yp1);
            fitFormula->SetParameter(6, yp2);
        }
    }
    return fitFormula;
}

inline Double_t interpolate_fit_values(Double_t energy, Double_t x, Double_t y) {
    vector<vector<Double_t>> fit_data = retrieve_fit_data("e-");

    pair<Double_t, Double_t> energy_bounds = get_energy_bounds(energy);
    Double_t lower = energy_bounds.first;
    Double_t upper = energy_bounds.second;

    TF2* fitLower = create_fit_function(lower, fit_data);
    TF2* fitUpper = create_fit_function(upper, fit_data);

    Double_t valueLower = fitLower->Eval(x, y);
    Double_t valueUpper = fitUpper->Eval(x, y);

    Double_t value = valueLower + (valueUpper - valueLower) * (energy - lower) / (upper - lower);
    return value;
}

inline Double_t get_PE_response(vector<vector<Double_t>> fit_data, Double_t energy, Double_t x, Double_t y) {
    // Get the energy bounds
    pair<Double_t, Double_t> energy_bounds = get_energy_bounds(energy);
    Double_t lowerBoundE = energy_bounds.first;
    Double_t upperBoundE = energy_bounds.second;
    // cout << "Energy bounds: " << lowerBoundE << " " << upperBoundE << endl;

    // Fit function
    // cout << "Creating fit function for energy: " << lowerBoundE << " MeV" << endl;
    TF2* fitFnLower_e = create_fit_function(lowerBoundE, fit_data);

    // cout << "Ceating fit function for energy: " << upperBoundE << " MeV" << endl;
    TF2* fitFnUpper = create_fit_function(upperBoundE, fit_data);

    // From fit function, get the value of the function at given x and y
    Double_t valueLower = fitFnLower_e->Eval(x, y);
    Double_t valueUpper = fitFnUpper->Eval(x, y);

    // cout << "Evaluated values (upper and lower): " << valueLower << ", " << valueUpper << endl;

    // Interpolate the values
    Double_t value = (valueUpper - valueLower) / (upperBoundE - lowerBoundE) * (energy - lowerBoundE) + valueLower;

    // cout << "Estimated PEs: " << value << endl;
    return value;
}

#endif  // SHOWER_MAX_RESPONSE_LOOKUP_H

