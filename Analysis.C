
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "TH1F.h"
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include <math.h> 
#include <stdio.h>
#include <string.h>
#include <TSpectrum.h>
#include <TSpectrumTransform.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <regex>

void FancyPlot();
vector<string> GetFileNames();
void SanitizeFile(string file);
std::vector<double> ReadFile_wavelength(string fileName);
std::vector<double> ReadFile_absorption(string fileName);
void FitRayleigh(TGraph* graph, double fitMin, double fitMax);
void findDuplicates(vector<std::string> &filenames, vector<vector <double>> &WaveData, vector<vector <double>> &AbsData, vector<std::string> &filenamesAvg, vector<vector <double>> &WaveDataAvg, vector<vector <double>> &AbsDataAvg);

void Analysis()
{

	FancyPlot();
	system("rm *.png");
 	std::vector<string> fileNames = GetFileNames();
	vector<string> AvgNames;
	vector<string> correctedNames;

	std::vector<std::vector<double>>     Wavelength,     Absorption;
	std::vector<std::vector<double>>  WavelengthAvg,  AbsorptionAvg;
	std::vector<std::vector<double>> WavelengthBcor, AbsorptionBcor;
	vector<double> auxWcorrected, auxAcorrected;

	std::vector<TGraph*> BareDataPlots;
	std::vector<TGraph*> AvgDataPlots;
	std::vector<TGraph*> CorrectedDataPlots;



	cout << endl << "-------------------------------------------------------------------------" << endl;
	cout << "Using the data at filenames.dat to read each data file" << endl << endl;

 	for (int i = 0; i < fileNames.size(); ++i)
 	{
 		cout << endl;
 		cout << "----- " << fileNames[i] << " -----" << endl;
 		SanitizeFile(fileNames[i]);
		Wavelength.push_back(ReadFile_wavelength(fileNames[i]));
		Absorption.push_back(ReadFile_absorption(fileNames[i]));
 	}


// first plot ----------------------------------------------------------------

	TCanvas *Simple = new TCanvas("SimpleDataPlot", "SimpleDataPlot", 1920, 1080);

	TMultiGraph *mg = new TMultiGraph();
 	auto legend = new TLegend(0.65,0.45,0.85,0.85);

	for (int i = 0; i < fileNames.size(); ++i)
	{
		TGraph* auxGraph =	new TGraph(Wavelength[i].size(),&Wavelength[i][0],&Absorption[i][0]);
		auxGraph->SetTitle(fileNames[i].c_str());
		mg->Add(auxGraph);
		legend->AddEntry(auxGraph, fileNames[i].c_str() , "l");
		
		BareDataPlots.push_back(auxGraph); //also stores the TGraph in a vector for possible further use
	}

	mg->SetTitle(";Wavelength (nm);Absorbance");
	mg->Draw("A plc");
 	legend->Draw();
 	Simple->SaveAs("RawSpectrum.png");




// // second plot - averages ----------------------------------------------------------------

	TCanvas *averages = new TCanvas("AverageOverDuplicates", "AverageOverDuplicates", 1920, 1080);

	cout << endl << "-------------------------------------------------------------------------" << endl;
	cout << "looking for duplicates in the files" << endl << endl;
	findDuplicates(fileNames, Wavelength, Absorption, AvgNames, WavelengthAvg, AbsorptionAvg);
	cout << endl << "-------------------------------------------------------------------------" << endl;
	cout << "The now concatenated (averaged) measurements are:" << endl <<  endl;

	for (int i = 0; i < AvgNames.size(); ++i)
	{
		cout << AvgNames[i] << endl;
	}

	TMultiGraph *mgAvg = new TMultiGraph();
 	auto legendAvg = new TLegend(0.65,0.45,0.85,0.85);

	for (int i = 0; i < AvgNames.size(); ++i)
	{
		AvgDataPlots.push_back(new TGraph(WavelengthAvg[i].size(),&WavelengthAvg[i][0],&AbsorptionAvg[i][0]));
		AvgDataPlots[i]->SetTitle(AvgNames[i].c_str());
		mgAvg->Add(AvgDataPlots[i]);
		legendAvg->AddEntry(AvgDataPlots[i], AvgNames[i].c_str() , "l");

	}
	
	mgAvg->SetTitle(";Wavelength (nm);Absorbance");


// next find automatically the baseline index from after and before
// and perform the splines automatically.
	int newBaselineIndex = -1;
	int oldBaselineIndex = -1;

	for (int i = 0; i < AvgNames.size(); ++i)
	{
		if (AvgNames[i].find("baseline_after") != std::string::npos) {
            newBaselineIndex = i; // Found a match
        }
		if (AvgNames[i].find("baseline_before") != std::string::npos) {
            oldBaselineIndex = i; // Found a match
        }
	}

	if (newBaselineIndex == -1 || oldBaselineIndex == -1)
	{
		cout << "baseline not found... " << endl;
		return;
	}

 	// // preparing some info for next section
	TGraph* newBaselineGraph = new TGraph(WavelengthAvg[newBaselineIndex].size(),&WavelengthAvg[newBaselineIndex][0],&AbsorptionAvg[newBaselineIndex][0]);
	TSpline3* newBaselineSpline = new TSpline3("New Baseline Spline", newBaselineGraph);
	
	TGraph* oldBaselineGraph = new TGraph(WavelengthAvg[oldBaselineIndex].size(),&WavelengthAvg[oldBaselineIndex][0],&AbsorptionAvg[oldBaselineIndex][0]);
	TSpline3* oldBaselineSpline = new TSpline3("Old Baseline Spline", oldBaselineGraph);
	
	legendAvg->AddEntry(newBaselineSpline, "newBaseline" , "l");
	legendAvg->AddEntry(oldBaselineSpline, "oldBaseline" , "l");

	newBaselineSpline->SetLineWidth(2);
	oldBaselineSpline->SetLineWidth(1);
	newBaselineSpline->SetLineColor(2);
	oldBaselineSpline->SetLineColor(2);
	newBaselineSpline->SetLineStyle(7);
	oldBaselineSpline->SetLineStyle(2);

	mgAvg->Draw("A plc");
	newBaselineSpline->Draw("same");
	oldBaselineSpline->Draw("same");
 	legendAvg->Draw();
	averages->Update();
 	averages->SaveAs("Averages.png");

// 	// Third plot - New Baseline ----------------------------------------------------------------


 	TCanvas *newBaseline = new TCanvas("New Baseline", "New Baseline plots", 1920, 1080);

 	int MinWave = 280;
 	int MaxWave = 700;
 	cout << " ----->>> Using data " << AvgNames[newBaselineIndex] << " as the new baseline." << endl;
 	cout << " ----->>> Using data " << AvgNames[oldBaselineIndex] << " as the old baseline." << endl;

 	cout << "We will also restric the limits between: ";
	cout << MinWave << " nm and " << MaxWave << " nm" << endl; 

	for (int i = 0; i < AvgNames.size(); ++i)
	{
		correctedNames.push_back(AvgNames[i]);
		for (int j = 0; j < WavelengthAvg[i].size(); ++j)
		{
			if (WavelengthAvg[i][j] >= MinWave && WavelengthAvg[i][j] <= MaxWave)
				{
					auxWcorrected.push_back(WavelengthAvg[i][j]);
					double newbaselineCorrection = newBaselineSpline->Eval(WavelengthAvg[i][j]);
					double oldbaselineCorrection = oldBaselineSpline->Eval(WavelengthAvg[i][j]);
					auxAcorrected.push_back(AbsorptionAvg[i][j] + oldbaselineCorrection - newbaselineCorrection);
				}	
		}

        WavelengthBcor.push_back(auxWcorrected);
        AbsorptionBcor.push_back(auxAcorrected);
        auxWcorrected.clear();
        auxAcorrected.clear();
	}

	TMultiGraph *mgCorrected = new TMultiGraph();
 	auto legendCorrected = new TLegend(0.65,0.45,0.85,0.85);

	cout << "----------------------------" << endl;
	for (int i = 0; i < correctedNames.size(); ++i)
	{
		CorrectedDataPlots.push_back(new TGraph(WavelengthBcor[i].size(),&WavelengthBcor[i][0],&AbsorptionBcor[i][0]));
		CorrectedDataPlots[i]->SetTitle(correctedNames[i].c_str());
		mgCorrected->Add(CorrectedDataPlots[i]);
			FitRayleigh(CorrectedDataPlots[i], 320, 600);
		legendCorrected->AddEntry(CorrectedDataPlots[i], correctedNames[i].c_str() , "l");
		cout << correctedNames[i] << endl;


	}
	
	mgCorrected->SetTitle("Baseline Corrected Spectrum;Wavelength (nm);Absorbance");
	mgCorrected->Draw("A plc");
 	legendCorrected->Draw();

 	newBaseline->SaveAs("BaselineCorrected.png");

	return;
} 


void FancyPlot()
{
	gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.95);
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.2);
	gStyle->SetStatBorderSize(3);
	gStyle->SetLineWidth(3);
}



std::vector<string> GetFileNames()
{
	cout << "Reading the folder for .txt files. " << endl;

	string fileRm = "rm filenames.dat";
	string fileCt = "touch filenames.dat";
	string fileLs = "ls | grep .txt >> filenames.dat";

	cout << "Removing old filenames.dat" << endl;
	int resultRm = system(fileRm.c_str());
	if (resultRm == 0) 
	{
	   cout << "  -> Previous filename.dat removal was successful." << endl;
	} 
	else {
	    std::cerr << "  -> Previous filename.dat removal was NOT successful." << std::endl;
	}


	cout << "Creating a new filenames.dat" << endl;
	int resultCt = system(fileCt.c_str());
	if (resultRm == 0) 
	{
	   cout << "  -> filename.dat creation was successful." << endl;
	} 
	else {
	    std::cerr << "  -> filename.dat creation was NOT successful." << std::endl;
	}

	cout << "Filling filenames.dat with the ls contents" << endl;
	int resultLs = system(fileLs.c_str());
	if (resultLs == 0) {
	    std::cout << "  -> Filling filenames.dat  was successful!" << std::endl;
	} else {
	    std::cerr << "  -> Filling filenames.dat was NOT successful." << std::endl;
	}


	cout << "Now reading this new filenames.dat and storing at the fileNames vector" << endl;
	string filename = "filenames.dat";

	std::vector<std::string> fileNames;
    FILE* file = fopen(filename.c_str(), "r");

    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return fileNames;
    }

    char buffer[256];
    while (fscanf(file, "%255s", buffer) == 1) {
        fileNames.push_back(buffer);
    }

    fclose(file);
    return fileNames;
}



void SanitizeFile(string fileName)
{
	cout << "Starting sanitization for file " << fileName << endl;
	string sedLine = "sed -i '/^[0-9]*\\.[0-9]*,\\s*$/d' " + fileName;

	int result = system(sedLine.c_str());

	if (result == 0) {
	    std::cout << "  -> Sanitization was successful!" << std::endl;
	} else {
	    std::cerr << "  -> Error while Sanitizing" << std::endl;
	}
}




std::vector<double> ReadFile_wavelength(string fileName) {


		double auxW, auxA; //auxiliary variables to read the file
		std::vector<double> WavelengthData; //vector to store the wavelength data

		FILE* ptr = fopen(fileName.c_str(),"r");

		cout << "Trying to read file " << fileName.c_str() << " - wavelength data" << endl;
		
		if (ptr==NULL)
		{
			cout << "------> Problems reading file " << fileName.c_str() << endl;
			return WavelengthData;
		}

		fscanf(ptr, "%*[^\n]\n"); // Skip first  line
		fscanf(ptr, "%*[^\n]\n"); // Skip second line


		while (fscanf(ptr,"%lf,%lf", &auxW, &auxA)!= EOF)
		{
			WavelengthData.push_back(auxW);
		}


		int size   = WavelengthData.size();
		auto PointtoMin = std::min_element(WavelengthData.begin(), WavelengthData.end());
		auto PointtoMax = std::max_element(WavelengthData.begin(), WavelengthData.end());

		double min = *PointtoMin;
		double max = *PointtoMax;

		cout << "  -> There was " << size << " data points.";
		cout << " From " << min << " nm to " << max << " nm." << endl;

		fclose(ptr);

		return WavelengthData;
} 

std::vector<double> ReadFile_absorption(string fileName) {


		double auxW, auxA; //auxiliary variables to read the file
		std::vector<double> AbsorptionData; //vector to store the wavelength data
		
		FILE* ptr = fopen(fileName.c_str(),"r");

		cout << "Trying to read file " << fileName.c_str() << " - absorption data" << endl;
		
		if (ptr==NULL)
		{
			cout << "------> Problems reading file " << fileName.c_str() << endl;
			return AbsorptionData;
		}

		fscanf(ptr, "%*[^\n]\n"); // Skip first  line
		fscanf(ptr, "%*[^\n]\n"); // Skip second line


		while (fscanf(ptr,"%lf,%lf", &auxW, &auxA)!= EOF)
		{
			AbsorptionData.push_back(auxA);
		}


		int size   = AbsorptionData.size();
		auto PointtoMin = std::min_element(AbsorptionData.begin(), AbsorptionData.end());
		auto PointtoMax = std::max_element(AbsorptionData.begin(), AbsorptionData.end());

		double min = *PointtoMin;
		double max = *PointtoMax;

		cout << "  -> There was " << size << " data points.";
		cout << " From " << min << " to " << max << "." << endl;

		fclose(ptr);

		return AbsorptionData;
} 


void findDuplicates(vector<std::string> &filenames, vector<vector <double>> &WaveData, vector<vector <double>> &AbsData, vector<std::string> &filenamesAvg, vector<vector <double>> &WaveDataAvg, vector<vector <double>> &AbsDataAvg) 
{
	// honestly this came from ChatGPT, so I'm not 100% sure how it works
	std::vector<double> auxW, auxA;

    std::unordered_map<std::string, std::vector<int>> fileGroups;
    std::regex pattern(R"((.*)_\d+\.txt$)");  // Regex for "XXX_N.txt"

    for (size_t i = 0; i < filenames.size(); ++i) 
    {
        std::smatch match;
        if (std::regex_match(filenames[i], match, pattern)) 
        {
            fileGroups[match[1]].push_back(i);  // Store base name + index
        }
    }

        // Output duplicates
    for (const auto& entry : fileGroups) 
    {
        if (entry.second.size() > 0) 
        {  

        	// Multiple entries for the same base name
 			std::cout << "Duplicate found for: " << entry.first << "\nIndices: ";
            for (int idx : entry.second)
            {
                std::cout << idx << " ";
            }
            std::cout << "\n";
        
            //now to average all this values

        	filenamesAvg.push_back(entry.first + "_Avg");

            // auto* fist_element_pointer = entry.second.begin();
            int first_element = entry.second.front();

            for (int i = 0; i < WaveData[first_element].size(); ++i)
            { 
            	double auxValueW = 0, auxValueA = 0;
            	int indexer = 0;
	            for (int idx : entry.second) 
	            {
	                auxValueW = auxValueW + WaveData[idx][i];
	                auxValueA = auxValueA +  AbsData[idx][i];
	                indexer++;
	            }
	            auxW.push_back(auxValueW/indexer);
	            auxA.push_back(auxValueA/indexer);
            }

            WaveDataAvg.push_back(auxW);
            AbsDataAvg.push_back(auxA);
            auxW.clear();
            auxA.clear();
        }
    }
}

void FitRayleigh(TGraph* graph, double fitMin, double fitMax) {
    if (!graph) {
        std::cerr << "Error: Null graph pointer." << std::endl;
        return;
    }

    TF1* rayleighFit = new TF1("rayleighFit", "[0]/x^4 + [1]", fitMin, fitMax);
    rayleighFit->SetParameter(0, 1.0);  // Initial guess for amplitude

    graph->Fit(rayleighFit, "R");  // 'R' option restricts fit to range
}
