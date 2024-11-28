#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TPaveText.h"
#include "TDirectory.h"
#include "TClass.h"
#include "TRegexp.h"
#include "TMath.h"
#include <regex>
#include "TH1D.h"
#include "THStack.h"
#include <TH1F.h>
#include <TLegend.h>
#include "TPaveStats.h"
#include <algorithm> // For std::max
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include <set>
#include <iomanip>
#include <vector>
#include <utility>
#include <map>

#include <sstream>


// Data structure to store the interaction compositions
std::vector<std::string> interaction_compositions;


void SaveHistogramAsPNG(TObject* histogram, const std::string& suffix);
void CustomizeStatBox(TH1D* hist, const std::string& suffix, int total_tracks, double weighted_mean, double weighted_stddev);
std::tuple<double, double, double, double> ProcessHistogram(TH1D* hist, const std::string& suffix);
void SaveHistogramToCSV(TH1* histogram, const std::string& suffix);
void SaveStatisticsToCSV(const std::string& suffix, const std::vector<std::pair<TH1*, std::string>>& histograms);



int caf_plotter(std::string input_file_list, std::string output_rootfile, bool dataOnly) {

	int true_interactions = 0;
	int true_interactions_cc = 0;
	int true_interactions_ccLArFV = 0;

	int reco_interactions = 0;
	int reco_interactions_cc = 0;
	int reco_interactions_ccLArFV = 0;
	int reco_interactions_ccLArFV_minerva = 0;

	
	
		

    double mnvOffsetX = 0; 
	double mnvOffsetY = 0;
    if (dataOnly){ mnvOffsetX=-10; mnvOffsetY=5;}
	
	int mult_bins_genie = 20;
    int mult_edge_L_genie = 0;
    int mult_edge_U_genie = 20;
		
    TH1D *reco_mult_prim_total = new TH1D("reco_mult_prim_total", "Reconstructed Track Multiplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);
	TH1D *reco_mult_prim_shower = new TH1D("reco_mult_prim_shower", "Reconstructed Shower Multiplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);

    TH1D *true_mult_genie_total = new TH1D("true_mult_genie_total", "True GENIE Multplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);
    TH1D *true_mult_trackOnly = new TH1D("true_mult_trackOnly", "True Track Multplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);

	TH1D *reco_mult_muon = new TH1D("reco_mult_muon", "Muon Multiplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);
	TH1D *reco_mult_pion = new TH1D("reco_mult_pion", "Charged Pion Multiplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);
	TH1D *reco_mult_proton = new TH1D("reco_mult_proton", "Proton Multiplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);
	TH1D *reco_mult_kaon = new TH1D("reco_mult_kaon", "Charged Kaon Multiplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);


	int length_bins = 50;
	int length_bins_L = 0;
	int length_bins_U = 100;

	TH1F *reco_length_prim_muon = new TH1F("reco_length_prim_muon", "Reconstructed Muons;Track Length (cm);Counts", length_bins, length_bins_L, length_bins_U);	
	TH1F *reco_length_prim_pion = new TH1F("reco_length_prim_pion", "Reconstructed Charged Pions;Track Length (cm);Counts", length_bins, length_bins_L, length_bins_U);
	TH1F *reco_length_prim_proton = new TH1F("reco_length_prim_proton", "Reconstructed Protons;Track Length (cm);Counts", length_bins, length_bins_L, length_bins_U);
	TH1F *reco_length_prim_kaon = new TH1F("reco_length_prim_kaon", "Reconstructed Charged Kaons;Track Length (cm);Counts", length_bins, length_bins_L, length_bins_U);


	int energy_bins = 50;
	int energy_bins_L = 0;
	int energy_bins_U = 1;

	TH1F *reco_energy_prim_muon = new TH1F("reco_energy_prim_muon", "Reconstructed Muons;Energy (GeV);Counts", energy_bins, energy_bins_L, energy_bins_U);
	TH1F *reco_energy_prim_pion = new TH1F("reco_energy_prim_pion", "Reconstructed Charged Pions;Energy (GeV);Counts", energy_bins, energy_bins_L, energy_bins_U);
	TH1F *reco_energy_prim_proton = new TH1F("reco_energy_prim_proton", "Reconstructed Protons;Energy (GeV);Counts", energy_bins, energy_bins_L, energy_bins_U);
	TH1F *reco_energy_prim_kaon = new TH1F("reco_energy_prim_kaon", "Reconstructed Charged Kaons;Energy (GeV);Counts", energy_bins, energy_bins_L, energy_bins_U);

	int cosTheta_bins = 100;
	int cosTheta_bins_L = -1;
	int cosTheta_bins_U = 1;

	TH1F *reco_cosTheta_prim_muon = new TH1F("reco_cosTheta_prim_muon", "Reconstructed Muons;Cosine Theta;Counts", cosTheta_bins, cosTheta_bins_L, cosTheta_bins_U);
	TH1F *reco_cosTheta_prim_pion = new TH1F("reco_cosTheta_prim_pion", "Reconstructed Charged Pions;Cosine Theta;Counts", cosTheta_bins, cosTheta_bins_L, cosTheta_bins_U);
	TH1F *reco_cosTheta_prim_proton = new TH1F("reco_cosTheta_prim_proton", "Reconstructed Protons;Cosine Theta;Counts", cosTheta_bins, cosTheta_bins_L, cosTheta_bins_U);
	TH1F *reco_cosTheta_prim_kaon = new TH1F("reco_cosTheta_prim_kaon", "Reconstructed Charged Kaons;Cosine Theta;Counts", cosTheta_bins, cosTheta_bins_L, cosTheta_bins_U);

    int start_bins = 50;
    int start_bins_L = -70;
    int start_bins_U = 70;
	
	TH1F *reco_start_x_muon = new TH1F("reco_start_x_muon", "Start X of Muons;Start X (cm);Counts", start_bins, start_bins_L, start_bins_U);
	TH1F *reco_start_y_muon = new TH1F("reco_start_y_muon", "Start Y of Muons;Start Y (cm);Counts", start_bins, start_bins_L, start_bins_U);
	TH1F *reco_start_z_muon = new TH1F("reco_start_z_muon", "Start Z of Muons;Start Z (cm);Counts", start_bins, start_bins_L, start_bins_U);

	TH1F *reco_start_x_pion = new TH1F("reco_start_x_pion", "Start X of Charged Pions;Start X (cm);Counts", start_bins, start_bins_L, start_bins_U);
	TH1F *reco_start_y_pion = new TH1F("reco_start_y_pion", "Start Y of Charged Pions;Start Y (cm);Counts", start_bins, start_bins_L, start_bins_U);
	TH1F *reco_start_z_pion = new TH1F("reco_start_z_pion", "Start Z of Charged Pions;Start Z (cm);Counts", start_bins, start_bins_L, start_bins_U);
	
	TH1F *reco_start_x_proton = new TH1F("reco_start_x_proton", "Start X of Protons;Start X (cm);Counts", start_bins, start_bins_L, start_bins_U);
	TH1F *reco_start_y_proton = new TH1F("reco_start_y_proton", "Start Y of Protons;Start Y (cm);Counts", start_bins, start_bins_L, start_bins_U);
	TH1F *reco_start_z_proton = new TH1F("reco_start_z_proton", "Start Z of Protons;Start Z (cm);Counts", start_bins, start_bins_L, start_bins_U);

	TH1F *reco_start_x_kaon = new TH1F("reco_start_x_kaon", "Start X of Charged Kaons;Start X (cm);Counts", start_bins, start_bins_L, start_bins_U);
	TH1F *reco_start_y_kaon = new TH1F("reco_start_y_kaon", "Start Y of Charged Kaons;Start Y (cm);Counts", start_bins, start_bins_L, start_bins_U);
	TH1F *reco_start_z_kaon = new TH1F("reco_start_z_kaon", "Start Z of Charged Kaons;Start Z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1D *reco_nu_vtx_x_CC_LArFV = new TH1D("reco_nu_vtx_x_CC_LArFV", "Neutrino Vertex-x (within LArFV);Vertex-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_y_CC_LArFV = new TH1D("reco_nu_vtx_y_CC_LArFV", "Neutrino Vertex-y (within LArFV);Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_z_CC_LArFV = new TH1D("reco_nu_vtx_z_CC_LArFV", "Neutrino Vertex-z (within LArFV);Vertex-z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1D *true_nu_vtx_x_CC_LArFV = new TH1D("true_nu_vtx_x_CC_LArFV", "Neutrino Vertex-x (within LArFV);Vertex-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *true_nu_vtx_y_CC_LArFV = new TH1D("true_nu_vtx_y_CC_LArFV", "Neutrino Vertex-y (within LArFV);Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *true_nu_vtx_z_CC_LArFV = new TH1D("true_nu_vtx_z_CC_LArFV", "Neutrino Vertex-z (within LArFV);Vertex-z (cm);Counts", start_bins, start_bins_L, start_bins_U);


    std::ifstream caf_list(input_file_list.c_str());
    if (!caf_list.is_open()) {
        std::cerr << "File not found: " << input_file_list << std::endl;
        return 1;
    }

    std::string tmp;
    TChain *caf_chain = new TChain("cafTree");
    while (caf_list >> tmp) {
        caf_chain->Add(tmp.c_str());
        std::cout << "Adding File: " << tmp << std::endl;
    }

    if (!caf_chain) {
        std::cerr << "There is no tree in the input files." << std::endl;
        return 1;
    }

    long Nentries = caf_chain->GetEntries();
    std::cout << "Total number of spills = " << Nentries << std::endl;

    auto sr = new caf::StandardRecord;
    caf_chain->SetBranchAddress("rec", &sr);
	
    const double distance_from_wall = 5.0;
    const double minX = -63.931 + distance_from_wall;
    const double maxX = +63.931 - distance_from_wall;
    const double minY = -62.076 + distance_from_wall;
    const double maxY = +62.076 - distance_from_wall;
    const double minZ = -64.538 + distance_from_wall;
    const double maxZ = +64.538 - distance_from_wall;

    const double X_average = (minX + maxX) / 2;
    const double Z_average = (minZ + maxZ) / 2;
    const double inner_epsilon = 5.0;
    const double Xinner_neg = X_average - inner_epsilon;
    const double Xinner_pos = X_average + inner_epsilon;
    const double Zinner_neg = Z_average - inner_epsilon;
    const double Zinner_pos = Z_average + inner_epsilon;

	// Module 0 bounds
	const double module_0_minX = Xinner_pos;
	const double module_0_maxX = maxX;
	const double module_0_minZ = Zinner_pos;
	const double module_0_maxZ = maxZ;
	const double module_0_minY = -57.076;
    const double module_0_maxY = +57.076;
	
	// Module 1 bounds
	const double module_1_minX = Xinner_pos;
	const double module_1_maxX = maxX;
	const double module_1_minZ = minZ;
	const double module_1_maxZ = Zinner_neg;
	const double module_1_minY = -57.076;
    const double module_1_maxY = +57.076;

	// Module 2 bounds
	const double module_2_minX = minX;
	const double module_2_maxX = Xinner_neg;
	const double module_2_minZ = Zinner_pos;
	const double module_2_maxZ = maxZ;
	const double module_2_minY = -57.076;
    const double module_2_maxY = +57.076;

	// Module 3 bounds
	const double module_3_minX = minX;
	const double module_3_maxX = Xinner_neg;
	const double module_3_minZ = minZ;
	const double module_3_maxZ = Zinner_neg;
	const double module_3_minY = -57.076;
    const double module_3_maxY = +57.076;


	// Global counters for individual particle types
	int total_muons = 0;
	int total_pions = 0;
	int total_protons = 0;
	int total_kaons = 0;
	int total_tracks = 0;


    for (long n = 0; n < Nentries; n++) {
        if (n % 10000 == 0) std::cout << "Processing trigger " << n << " of " << Nentries << std::endl;
        caf_chain->GetEntry(n);


		// Loop over truth interactions
        for (long unsigned ntrue = 0; ntrue < sr->mc.nu.size(); ntrue++) {
            auto vertex = sr->mc.nu[ntrue].vtx;
            double true_nu_vtxX = vertex.x;
            double true_nu_vtxY = vertex.y;
            double true_nu_vtxZ = vertex.z;

            auto truePrimary = sr->mc.nu[ntrue].prim;

            int trueGenieMultiplicity = 0;         // for true_mult_genie_total
            int trueTrackMultiplicity = 0;     // true_mult_trackOnly

			// Count total interactions
			true_interactions++;
			
			// Check for CC inetraction and count
            if (sr->mc.nu[ntrue].iscc == false) continue;
			true_interactions_cc++; 	
			
			// Check for CC interaction in LAr
            if (sr->mc.nu[ntrue].targetPDG != 1000180400) continue;
			
			// Conditions for All Modules			
			if (abs(true_nu_vtxX) > maxX || abs(true_nu_vtxX) < minX) continue;
			if (abs(true_nu_vtxY) > maxY || abs(true_nu_vtxY) < minY) continue;
			if (abs(true_nu_vtxZ) > maxZ || abs(true_nu_vtxZ) < minZ) continue;
			if (true_nu_vtxX > Xinner_neg && true_nu_vtxX < Xinner_pos) continue; 
			if (true_nu_vtxZ > Zinner_neg && true_nu_vtxZ < Zinner_pos) continue; 

			
			// int muon_count = 0;
			// int other_particle_count = 0;
			
			// for (long unsigned primaries = 0; primaries < truePrimary.size(); primaries++) {
				// int pdg = truePrimary[primaries].pdg;

				// Count particles
				// if (pdg == 13) muon_count++; // Muon
				// else other_particle_count++; // Any other particle
			// }

			// Keep interactions with exactly one muon and at least one other particle
			// if (muon_count != 1 || other_particle_count < 1) continue;
			
			
			// Count selected interactions within LArFV
			true_interactions_ccLArFV++;
						
			// Plot these interactions 
			true_nu_vtx_x_CC_LArFV->Fill(true_nu_vtxX);
			true_nu_vtx_y_CC_LArFV->Fill(true_nu_vtxY);
			true_nu_vtx_z_CC_LArFV->Fill(true_nu_vtxZ);	
			

            for (long unsigned primaries = 0; primaries < truePrimary.size(); primaries++) {
                int pdg = truePrimary[primaries].pdg;
                auto start_pos_true = sr->mc.nu[ntrue].prim[primaries].start_pos;
                auto end_pos_true = sr->mc.nu[ntrue].prim[primaries].end_pos;

                double dX_true = (end_pos_true.x - start_pos_true.x);
                double dY_true = (end_pos_true.y - start_pos_true.y);
                double dZ_true = (end_pos_true.z - start_pos_true.z);
                double length_true = TMath::Sqrt(dX_true * dX_true + dY_true * dY_true + dZ_true * dZ_true);
                double cosTheta_true = dZ_true / length_true;    

                double currentX_prim_true = start_pos_true.x;
                double currentY_prim_true = start_pos_true.y;
                double currentZ_prim_true = start_pos_true.z;


                if ((abs(pdg) == 13 || abs(pdg) == 2212 || abs(pdg) == 211 || abs(pdg) == 13 ) && length_true > 0) {
                // if ((abs(pdg) == 2212 || abs(pdg) == 211 || abs(pdg) == 13 ) && length_true > 0) {
                    trueTrackMultiplicity++;
                }

                trueGenieMultiplicity++;
            }

            true_mult_genie_total->Fill(trueGenieMultiplicity);
            true_mult_trackOnly->Fill(trueTrackMultiplicity);
        }


        // for (long unsigned nixn = 1; nixn < 2; nixn++) {
        for (long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++) {
			
			// Count all interactions
			reco_interactions++;
			
			// Declare MINERvA variables (starts)
			bool 	oneContained 	= false; 
			bool 	oneNotContained = false; 
			bool 	goodInteraction = false;
			double	maxDotProductDS = -999; 
			double	maxDotProductUS = -999;
			int 	maxEventPar 	= -999; 
			int		maxEventTyp 	= -999; 
			int 	maxEventIxn 	= -999;
			double 	dirZExiting 	= -999;   
			double 	startZMuonCand 	= -999;
			double	longestTrk 		= -999;
			// Declare MINERvA variables (ends)

            double reco_nu_vtxX = sr->common.ixn.dlp[nixn].vtx.x;
            double reco_nu_vtxY = sr->common.ixn.dlp[nixn].vtx.y;
            double reco_nu_vtxZ = sr->common.ixn.dlp[nixn].vtx.z;


			int muon_count = 0;
			int other_particle_count = 0;

			for (size_t k = 0; k < sr->common.ixn.dlp[nixn].part.dlp.size(); k++) {
				int pdg = std::abs(sr->common.ixn.dlp[nixn].part.dlp[k].pdg);

				// Only process primary particles
				if (!sr->common.ixn.dlp[nixn].part.dlp[k].primary) continue;

				// Count particles
				if (pdg == 13) muon_count++; // Muon
				else other_particle_count++; // Any other particle
			}

			// Keep interactions with exactly one muon and at least one other particle
			// if (muon_count != 1 || other_particle_count < 1) continue;

			// Count selected interactions within LArFV
			reco_interactions_cc++;


			// Conditions for All Modules			
			// if (abs(reco_nu_vtxX) > maxX || abs(reco_nu_vtxX) < minX) continue;
			// if (abs(reco_nu_vtxY) > maxY || abs(reco_nu_vtxY) < minY) continue;
			// if (abs(reco_nu_vtxZ) > maxZ || abs(reco_nu_vtxZ) < minZ) continue;
			// if (reco_nu_vtxX > Xinner_neg && reco_nu_vtxX < Xinner_pos) continue; 
			// if (reco_nu_vtxZ > Zinner_neg && reco_nu_vtxZ < Zinner_pos) continue; 
			
			if (/*abs(abs(sr->common.ixn.dlp[nixn].vtx.x)-33)<1 || */ abs(sr->common.ixn.dlp[nixn].vtx.x)>59 || abs(sr->common.ixn.dlp[nixn].vtx.x)<5 || abs(sr->common.ixn.dlp[nixn].vtx.y)>57 || abs(sr->common.ixn.dlp[nixn].vtx.z)<5 || abs(sr->common.ixn.dlp[nixn].vtx.z)>59.5)  continue;

			
			// Count reco interactions in FV
			// reco_interactions_ccLArFV++;
				

			// Conditions for Module 0 (top-right in the array):
			// if (reco_nu_vtxX < module_0_minX || reco_nu_vtxX > module_0_maxX) continue;  // Module 0 X bounds
			// if (reco_nu_vtxY < module_0_minY || reco_nu_vtxY > module_0_maxY) continue;  // Module 0 Y bounds
			// if (reco_nu_vtxZ < module_0_minZ || reco_nu_vtxZ > module_0_maxZ) continue;  // Module 0 Z bounds

			// Conditions for Module 1 (top-left in the array):
			// if (reco_nu_vtxX < module_1_minX || reco_nu_vtxX > module_1_maxX) continue;  // Module 1 X bounds
			// if (reco_nu_vtxY < module_1_minY || reco_nu_vtxY > module_1_maxY) continue;  // Module 1 Y bounds
			// if (reco_nu_vtxZ < module_1_minZ || reco_nu_vtxZ > module_1_maxZ) continue;  // Module 1 Z bounds

			// Conditions for Module 2 (bottom-right in the array):
			// if (reco_nu_vtxX < module_2_minX || reco_nu_vtxX > module_2_maxX) continue;  // Module 2 X bounds
			// if (reco_nu_vtxY < module_2_minY || reco_nu_vtxY > module_2_maxY) continue;  // Module 2 Y bounds
			// if (reco_nu_vtxZ < module_2_minZ || reco_nu_vtxZ > module_2_maxZ) continue;  // Module 2 Z bounds

			// Conditions for Module 3 (bottom-left in the array):
			// if (reco_nu_vtxX < module_3_minX || reco_nu_vtxX > module_3_maxX) continue;  // Module 3 X bounds
			// if (reco_nu_vtxY < module_3_minY || reco_nu_vtxY > module_3_maxY) continue;  // Module 3 Y bounds
			// if (reco_nu_vtxZ < module_3_minZ || reco_nu_vtxZ > module_3_maxZ) continue;  // Module 3 Z bounds



			int reco_multiplicity = 0;
			int reco_muon_count = 0;
			int reco_pion_count = 0;	
			int reco_proton_count = 0;
			int reco_kaon_count = 0;

			
			int reco_muon_count_interaction = 0;
			int reco_pion_count_interaction = 0;
			int reco_proton_count_interaction = 0;
			int reco_kaon_count_interaction = 0;

			std::ostringstream interaction_details;
			

			for (size_t k = 0; k < sr->common.ixn.dlp[nixn].part.dlp.size(); k++) {

				int pdg = std::abs(sr->common.ixn.dlp[nixn].part.dlp[k].pdg);
				if (abs(pdg)!=321 && abs(pdg)!=211 && abs(pdg)!=2212 && abs(pdg)!=13) continue;
				if (sr->common.ixn.dlp[nixn].part.dlp[k].primary!=true) continue;

				
				double start_x = sr->common.ixn.dlp[nixn].part.dlp[k].start.x;
				double start_y = sr->common.ixn.dlp[nixn].part.dlp[k].start.y;
				double start_z = sr->common.ixn.dlp[nixn].part.dlp[k].start.z;
				double end_x = sr->common.ixn.dlp[nixn].part.dlp[k].end.x;
				double end_y = sr->common.ixn.dlp[nixn].part.dlp[k].end.y;
				double end_z = sr->common.ixn.dlp[nixn].part.dlp[k].end.z;
				double length = TMath::Sqrt(TMath::Power(end_x - start_x, 2) + TMath::Power(end_y - start_y, 2) + TMath::Power(end_z - start_z, 2));				
				double cosTheta = (end_z - start_z, 2) / length;
				double energy = sr->common.ixn.dlp[nixn].part.dlp[k].E;
				
				
				// MINERvA-pairing starts
				auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[k].start;
				auto end_pos = sr->common.ixn.dlp[nixn].part.dlp[k].end;
				
				double diffVertexdZ = abs(start_pos.z - sr->common.ixn.dlp[nixn].vtx.z);
				double diffVertexdX = abs(start_pos.x - sr->common.ixn.dlp[nixn].vtx.x);
				double diffVertexdY = abs(start_pos.y - sr->common.ixn.dlp[nixn].vtx.y);
				double diffVertex = TMath::Sqrt(diffVertexdZ * diffVertexdZ + diffVertexdY * diffVertexdY + diffVertexdX * diffVertexdX);
	
				if (diffVertex > 5) continue;
    	
				double dX = (end_pos.x - start_pos.x);
				double dY = (end_pos.y - start_pos.y);
				double dZ = (end_pos.z - start_pos.z);
				double length_again = TMath::Sqrt(dX*dX + dY*dY + dZ*dZ);
				double dirX = dX/length_again; 
				double dirY = dY/length_again; 
				double dirZ = dZ/length_again;

				if (dirZ < 0) { 
					dirZ = -dirZ; 
					dirX = -dirX; 
					dirY = -dirY; 
					auto temp = start_pos; 
					end_pos = start_pos; 
					end_pos = temp;
				}
	
				if (std::isnan(start_pos.z)) length = -999;
	
				if (length_again > longestTrk) longestTrk = length_again;

				int maxPartMinerva	 = -999; 
				int maxTypeMinerva	 = -999;  
				int maxIxnMinerva 	 = -999;
				int maxPartMinervaUS = -999; 
				int maxTypeMinervaUS = -999;  
    
				if ((abs(start_pos.z) > maxZ || abs(end_pos.z) > maxZ) ) { 

					int minervaPass 		= 0;
					double dotProductDS 	= -999; 
					double deltaExtrapYUS 	= -999; 
					double deltaExtrapY 	= -999; 
					double dotProductUS 	= -999; 
					double deltaExtrapX 	= -999; 
					double deltaExtrapXUS 	= -999;	
			
					for(int i = 0; i < sr->nd.minerva.ixn.size(); i++) {
						for (int j = 0; j < sr->nd.minerva.ixn[i].ntracks; j++) {
							
							double dir_z	= sr->nd.minerva.ixn[i].tracks[j].dir.z;
							double end_z	= sr->nd.minerva.ixn[i].tracks[j].end.z;
							double start_z	= sr->nd.minerva.ixn[i].tracks[j].start.z;
							double end_x	= sr->nd.minerva.ixn[i].tracks[j].end.x;
							double start_x	= sr->nd.minerva.ixn[i].tracks[j].start.x;
							double end_y	= sr->nd.minerva.ixn[i].tracks[j].end.y;
							double start_y	= sr->nd.minerva.ixn[i].tracks[j].start.y;
					

							if (start_z > 0 && ((start_pos.z) > maxZ || (end_pos.z) > maxZ) ) {

								int truthPart			= sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
								double dXMnv			= (sr->nd.minerva.ixn[i].tracks[j].end.x - sr->nd.minerva.ixn[i].tracks[j].start.x);
								double dYMnv			= (sr->nd.minerva.ixn[i].tracks[j].end.y - sr->nd.minerva.ixn[i].tracks[j].start.y);
								double dZMnv			= (sr->nd.minerva.ixn[i].tracks[j].end.z - sr->nd.minerva.ixn[i].tracks[j].start.z);
								double lengthMinerva	= TMath::Sqrt(dXMnv*dXMnv + dYMnv*dYMnv + dZMnv*dZMnv);

								if (lengthMinerva < 10) continue;

								double dirXMinerva	= dXMnv/lengthMinerva;
								double dirYMinerva	= dYMnv/lengthMinerva;
								double dirZMinerva	= dZMnv/lengthMinerva;
								double dotProduct	= dirXMinerva*dirX + dirYMinerva*dirY + dirZ*dirZMinerva;
								double extrapdZ		= start_z - end_pos.z;
								double extrapY		= dirY/dirZ * (extrapdZ) + end_pos.y - start_y;
								double extrapX		= dirX/dirZ * (extrapdZ) + end_pos.x - start_x;
								double diffExtrap	= TMath::Sqrt(TMath::Power(extrapY - start_y, 2));


								if (dotProductDS < dotProduct && abs(extrapY-mnvOffsetY) < 15 && abs(TMath::ATan(dirXMinerva/dirZMinerva) - TMath::ATan(dirX/dirZ)) < 0.06 && abs(TMath::ATan(dirYMinerva/dirZMinerva) - TMath::ATan(dirY/dirZ)) < 0.06 && abs(extrapX-mnvOffsetX) < 15) {
									dotProductDS	= dotProduct;
									deltaExtrapY	= extrapY;
									deltaExtrapX	= extrapX;
									maxPartMinerva	= sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
									maxTypeMinerva	= sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
									maxIxnMinerva	= sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;	
									
									if (end_z > 300) {
										minervaPass = 1; 
										if(dirZExiting < dirZ) dirZExiting = dirZ;
									}	// if (end_z > 300) {	
								}	// if (dotProductDS < dotProduct && abs(extrapY-mnvOffsetY) < 15 && ...	
							}	// if (start_z > 0 && ((start_pos.z) > maxZ || (end_pos.z) > maxZ) ) {			
						}	// for (int j = 0; j < sr->nd.minerva.ixn[i].ntracks; j++) {	
					}	// for(int i = 0; i < sr->nd.minerva.ixn.size(); i++) {		
								
					if (dotProductDS > maxDotProductDS){ 
						maxDotProductDS=dotProductDS;
						maxEventPar=maxPartMinerva;
						maxEventTyp=maxTypeMinerva;
						maxEventIxn=maxIxnMinerva;
					}	// if (dotProductDS > maxDotProductDS){				

					// Count reco interactions in FV
					reco_interactions_ccLArFV_minerva++;

				}	// if ((abs(start_pos.z) > maxZ || abs(end_pos.z) > maxZ) ) {



				// Count tracks contributing to multiplicity (no maxDotProductDS yet)
				if (length > 2) {
					if (pdg == 13) reco_muon_count++;
					if (pdg == 211) reco_pion_count++;
					if (pdg == 2212) reco_proton_count++;
					if (pdg == 321) reco_kaon_count++;
					if (pdg == 13 || pdg == 211 || pdg == 2212 || pdg == 321) reco_multiplicity++;
					
				}
				


				// if (length > 2) {
					
					// if (abs(pdg) == 13 || abs(pdg) == 211 || abs(pdg) == 2212 || abs(pdg) == 321) {
						// reco_multiplicity++;
					// }
					
					// if (length > 2 && pdg == 13) total_muons++;
					// if (length > 2 && pdg == 211) total_pions++; 
					// if (length > 2 && pdg == 2212) total_protons++;
					// if (length > 2 && pdg == 321) total_kaons++;

					// if (abs(pdg) == 13 || abs(pdg) == 211 || abs(pdg) == 2212 || abs(pdg) == 321) total_tracks++;


					 // if (pdg == 13) {
						// reco_muon_count_interaction++;
					// } else if (pdg == 211) {
						// reco_pion_count_interaction++;
					// } else if (pdg == 2212) {
						// reco_proton_count_interaction++;
					// } else if (pdg == 321) {
						// reco_kaon_count_interaction++;
					// }


					// if (pdg == 13) {
						// reco_muon_count++;
					// }
					
					// if (pdg == 211) {
						// reco_pion_count++;
					// }
					
					// if (pdg == 2212) { 
						// reco_proton_count++;
					// }
					
					// if (pdg == 321) {
						// reco_kaon_count++;
					// }

				// }	// if (length > 0) {
					
					

					
			}	// for (size_t k = 0; k < sr->common.ixn.dlp[nixn].part.dlp.size(); k++) {


			// After the loop over particles, fill the proton multiplicity histogram if there was at least one proton in the interaction
			// if (reco_muon_count > 0 && maxDotProductDS>0.99) {				
				// reco_mult_muon->Fill(reco_muon_count);
			// }
			
			// if (reco_proton_count > 0 && maxDotProductDS > 0.99) {
				// reco_mult_proton->Fill(reco_proton_count);
			// } 
			
			// if (reco_pion_count > 0 && maxDotProductDS>0.99) {
				// reco_mult_pion->Fill(reco_pion_count);
			// }

			// if (reco_kaon_count > 0 && maxDotProductDS>0.99) {
				// reco_mult_kaon->Fill(reco_kaon_count);
			// }
			


			// Fill multiplicity histograms if maxDotProductDS condition is satisfied
			if (maxDotProductDS > 0.99) {
				if (reco_muon_count > 0) reco_mult_muon->Fill(reco_muon_count);
				if (reco_proton_count > 0) reco_mult_proton->Fill(reco_proton_count);
				if (reco_pion_count > 0) reco_mult_pion->Fill(reco_pion_count);
				if (reco_kaon_count > 0) reco_mult_kaon->Fill(reco_kaon_count);
				
				if (reco_multiplicity > 0) {
					reco_interactions_ccLArFV++;				
					reco_mult_prim_total->Fill(reco_multiplicity);
					
					reco_nu_vtx_x_CC_LArFV->Fill(reco_nu_vtxX);
					reco_nu_vtx_y_CC_LArFV->Fill(reco_nu_vtxY);
					reco_nu_vtx_z_CC_LArFV->Fill(reco_nu_vtxZ);	
										
					interaction_details << reco_muon_count << " ";
					interaction_details << reco_proton_count << " ";
					interaction_details << reco_pion_count << " ";
					interaction_details << reco_kaon_count << " ";
					
					interaction_compositions.push_back(interaction_details.str());	
				}
				
			}
			
			
			// Second pass: Fill property histograms
			if (maxDotProductDS > 0.99) {
				for (size_t k = 0; k < sr->common.ixn.dlp[nixn].part.dlp.size(); k++) {
					int pdg = std::abs(sr->common.ixn.dlp[nixn].part.dlp[k].pdg);
					if (!sr->common.ixn.dlp[nixn].part.dlp[k].primary) continue;

					double start_x = sr->common.ixn.dlp[nixn].part.dlp[k].start.x;
					double start_y = sr->common.ixn.dlp[nixn].part.dlp[k].start.y;
					double start_z = sr->common.ixn.dlp[nixn].part.dlp[k].start.z;
					double end_x = sr->common.ixn.dlp[nixn].part.dlp[k].end.x;
					double end_y = sr->common.ixn.dlp[nixn].part.dlp[k].end.y;
					double end_z = sr->common.ixn.dlp[nixn].part.dlp[k].end.z;
					double length = TMath::Sqrt(TMath::Power(end_x - start_x, 2) + TMath::Power(end_y - start_y, 2) + TMath::Power(end_z - start_z, 2));
					double cosTheta = (end_z - start_z) / length;
					double energy = sr->common.ixn.dlp[nixn].part.dlp[k].E;

					// Fill property histograms only for particles contributing to counts
					if (length > 2) {
						if (reco_muon_count > 0 && pdg == 13) {
							reco_length_prim_muon->Fill(length);
							reco_energy_prim_muon->Fill(energy);
							reco_cosTheta_prim_muon->Fill(cosTheta);
							reco_start_x_muon->Fill(start_x);
							reco_start_y_muon->Fill(start_y);
							reco_start_z_muon->Fill(start_z);
						} else if (reco_pion_count > 0 && pdg == 211) {
							reco_length_prim_pion->Fill(length);
							reco_energy_prim_pion->Fill(energy);
							reco_cosTheta_prim_pion->Fill(cosTheta);
							reco_start_x_pion->Fill(start_x);
							reco_start_y_pion->Fill(start_y);
							reco_start_z_pion->Fill(start_z);
						} else if (reco_proton_count > 0 && pdg == 2212) {
							reco_length_prim_proton->Fill(length);
							reco_energy_prim_proton->Fill(energy);
							reco_cosTheta_prim_proton->Fill(cosTheta);
							reco_start_x_proton->Fill(start_x);
							reco_start_y_proton->Fill(start_y);
							reco_start_z_proton->Fill(start_z);
						} else if (reco_kaon_count > 0 && pdg == 321) {
							reco_length_prim_kaon->Fill(length);
							reco_energy_prim_kaon->Fill(energy);
							reco_cosTheta_prim_kaon->Fill(cosTheta);
							reco_start_x_kaon->Fill(start_x);
							reco_start_y_kaon->Fill(start_y);
							reco_start_z_kaon->Fill(start_z);
						}
					}
				}
			}




			// if (reco_multiplicity > 0 && maxDotProductDS>0.99) {
				// reco_interactions_ccLArFV++;				
				// reco_mult_prim_total->Fill(reco_multiplicity);
				
				// reco_nu_vtx_x_CC_LArFV->Fill(reco_nu_vtxX);
				// reco_nu_vtx_y_CC_LArFV->Fill(reco_nu_vtxY);
				// reco_nu_vtx_z_CC_LArFV->Fill(reco_nu_vtxZ);	
				
				// Record the interaction details
				// interaction_details << "Interaction " << (nixn + 1) << ": ";
				// interaction_details << reco_muon_count << " muon(s), ";
				// interaction_details << reco_proton_count << " proton(s), ";
				// interaction_details << reco_pion_count << " pion(s), ";
				// interaction_details << reco_kaon_count << " kaon(s)";
				
				// interaction_details << reco_muon_count << " ";
				// interaction_details << reco_proton_count << " ";
				// interaction_details << reco_pion_count << " ";
				// interaction_details << reco_kaon_count << " ";
				
				// interaction_compositions.push_back(interaction_details.str());	
			// }
			


				
			std::cout << "Total tracks contributing to reco_mult_prim_total:" << std::endl;
			std::cout << "  Muons: " << reco_muon_count << std::endl;
			std::cout << "  Pions: " << reco_pion_count << std::endl;
			std::cout << "  Protons: " << reco_proton_count << std::endl;
			std::cout << "  Kaons: " << reco_kaon_count << std::endl;
			std::cout << "  total_tracks: " << total_tracks << std::endl;
				
			// Initialize electron and photon multiplicity
			int reco_multiplicity_shower = 0;

			// New loop for counting electron and photon showers
			for (size_t k = 0; k < sr->common.ixn.dlp[nixn].part.dlp.size(); k++) {
				int pdg = std::abs(sr->common.ixn.dlp[nixn].part.dlp[k].pdg);
				
				if (sr->common.ixn.dlp[nixn].part.dlp[k].primary!=true) continue;
				double start_x = sr->common.ixn.dlp[nixn].part.dlp[k].start.x;
				double start_y = sr->common.ixn.dlp[nixn].part.dlp[k].start.y;
				double start_z = sr->common.ixn.dlp[nixn].part.dlp[k].start.z;
				double end_x = sr->common.ixn.dlp[nixn].part.dlp[k].end.x;
				double end_y = sr->common.ixn.dlp[nixn].part.dlp[k].end.y;
				double end_z = sr->common.ixn.dlp[nixn].part.dlp[k].end.z;
				double length = TMath::Sqrt(TMath::Power(end_x - start_x, 2) + TMath::Power(end_y - start_y, 2) + TMath::Power(end_z - start_z, 2));
				

				if (pdg == 11 || pdg == 22) {
					if (length > 1) {
						reco_multiplicity_shower++;
					}	// if (length > 1) {
				}	// if (pdg == 11 || pdg == 22) {
			}	// for (size_t k = 0; k < sr->common.ixn.dlp[nixn].part.dlp.size(); k++) {

			if (reco_multiplicity_shower > 0 && maxDotProductDS > 0.99) {
				reco_mult_prim_shower->Fill(reco_multiplicity_shower);
			}	// if (reco_multiplicity_shower > 0 && maxDotProductDS > 0.99) {
	
				
				
        }	// for (long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++) {

    }	// for (long n = 0; n < Nentries; n++) {

std::ofstream interaction_log("interaction_compositions.csv");
if (interaction_log.is_open()) {
    // Write header to the file
    interaction_log << "Interaction,Muons,Protons,Pions,Kaons\n";

    // Print header to console
    std::cout << std::left << std::setw(20) << "Interaction"
              << std::setw(10) << "Muons"
              << std::setw(10) << "Protons"
              << std::setw(10) << "Pions"
              << std::setw(10) << "Kaons" << std::endl;

    std::cout << std::string(62, '-') << std::endl;

    int total_muons = 0;
    int total_protons = 0;
    int total_pions = 0;
    int total_kaons = 0;

    // Iterate through interactions
    for (size_t i = 0; i < interaction_compositions.size(); ++i) {
        // Parse composition into separate particle counts
        std::stringstream ss(interaction_compositions[i]);
        int muons = 0, protons = 0, pions = 0, kaons = 0;
        char delim; // For skipping commas
        // ss >> muons >> delim >> protons >> delim >> pions >> delim >> kaons;
		ss >> muons >> protons >> pions >> kaons; // Space-separated parsing


        // Write to the file
        interaction_log << i + 1 << "," << muons << "," << protons << "," << pions << "," << kaons << "\n";

        // Print to the console
        std::cout << std::left << std::setw(12) << (i + 1)
                  << std::setw(10) << muons
                  << std::setw(10) << protons
                  << std::setw(10) << pions
                  << std::setw(10) << kaons << std::endl;

        // Accumulate totals
        total_muons += muons;
        total_protons += protons;
        total_pions += pions;
        total_kaons += kaons;
    }

    // Print totals to the console
    std::cout << std::string(42, '-') << std::endl;
    std::cout << std::left << std::setw(12) << "Total"
              << std::setw(10) << total_muons
              << std::setw(10) << total_protons
              << std::setw(10) << total_pions
              << std::setw(10) << total_kaons << std::endl;

    // Close the file
    interaction_log.close();
    std::cout << "Interaction compositions saved to interaction_compositions.csv" << std::endl;
} else {
    std::cerr << "Failed to open file interaction_compositions.csv for writing." << std::endl;
}





TFile *caf_out_file = new TFile(output_rootfile.c_str(), "recreate");

// std::string fv_modules = "sandboxV4_modules_0123";
// std::string suffix = "sandboxV4_modules_0123";

std::string fv_modules = "MR6p1_modules_0123";
std::string suffix = "MR6p1_modules_0123";


true_mult_genie_total->Write();
ProcessHistogram(true_mult_genie_total, fv_modules);
SaveHistogramToCSV(true_mult_genie_total, suffix);

true_mult_trackOnly->Write();
ProcessHistogram(true_mult_trackOnly, fv_modules);
SaveHistogramToCSV(true_mult_trackOnly, suffix);

reco_mult_prim_total->Write();
// ProcessHistogram(reco_mult_prim_total, fv_modules);
SaveHistogramAsPNG(reco_mult_prim_total, suffix);
SaveHistogramToCSV(reco_mult_prim_total, suffix);

reco_mult_prim_shower->Write();
// ProcessHistogram(reco_mult_prim_shower, fv_modules);
SaveHistogramAsPNG(reco_mult_prim_shower, suffix);
SaveHistogramToCSV(reco_mult_prim_shower, suffix);

reco_mult_muon->Write();
// ProcessHistogram(reco_mult_muon, fv_modules);
SaveHistogramAsPNG(reco_mult_muon, suffix);
SaveHistogramToCSV(reco_mult_muon, suffix);

reco_mult_pion->Write();
// ProcessHistogram(reco_mult_pion, fv_modules);
SaveHistogramAsPNG(reco_mult_pion, suffix);
SaveHistogramToCSV(reco_mult_pion, suffix);

reco_mult_proton->Write();
// ProcessHistogram(reco_mult_proton, fv_modules);
SaveHistogramAsPNG(reco_mult_proton, suffix);
SaveHistogramToCSV(reco_mult_proton, suffix);

reco_mult_kaon->Write();
// ProcessHistogram(reco_mult_kaon, fv_modules);
SaveHistogramAsPNG(reco_mult_kaon, suffix);
SaveHistogramToCSV(reco_mult_kaon, suffix);


reco_length_prim_muon->Write();
SaveHistogramAsPNG(reco_length_prim_muon, fv_modules);
SaveHistogramToCSV(reco_length_prim_muon, fv_modules);

reco_length_prim_pion->Write();
SaveHistogramAsPNG(reco_length_prim_pion, fv_modules);
SaveHistogramToCSV(reco_length_prim_pion, fv_modules);

reco_length_prim_proton->Write();
SaveHistogramAsPNG(reco_length_prim_proton, fv_modules);
SaveHistogramToCSV(reco_length_prim_proton, fv_modules);

reco_length_prim_kaon->Write();
SaveHistogramAsPNG(reco_length_prim_kaon, fv_modules);
SaveHistogramToCSV(reco_length_prim_kaon, fv_modules);


reco_energy_prim_muon->Write();
SaveHistogramAsPNG(reco_energy_prim_muon, fv_modules);
SaveHistogramToCSV(reco_energy_prim_muon, fv_modules);

reco_energy_prim_pion->Write();
SaveHistogramAsPNG(reco_energy_prim_pion, fv_modules);
SaveHistogramToCSV(reco_energy_prim_pion, fv_modules);

reco_energy_prim_proton->Write();
SaveHistogramAsPNG(reco_energy_prim_proton, fv_modules);
SaveHistogramToCSV(reco_energy_prim_proton, fv_modules);

reco_energy_prim_kaon->Write();
SaveHistogramAsPNG(reco_energy_prim_kaon, fv_modules);
SaveHistogramToCSV(reco_energy_prim_kaon, fv_modules);


reco_cosTheta_prim_muon->Write();
SaveHistogramAsPNG(reco_cosTheta_prim_muon, fv_modules);
SaveHistogramToCSV(reco_cosTheta_prim_muon, fv_modules);

reco_cosTheta_prim_pion->Write();
SaveHistogramAsPNG(reco_cosTheta_prim_pion, fv_modules);
SaveHistogramToCSV(reco_cosTheta_prim_pion, fv_modules);

reco_cosTheta_prim_proton->Write();
SaveHistogramAsPNG(reco_cosTheta_prim_proton, fv_modules);
SaveHistogramToCSV(reco_cosTheta_prim_proton, fv_modules);

reco_cosTheta_prim_kaon->Write();
SaveHistogramAsPNG(reco_cosTheta_prim_kaon, fv_modules);
SaveHistogramToCSV(reco_cosTheta_prim_kaon, fv_modules);

reco_start_x_muon->Write();
reco_start_y_muon->Write();
reco_start_z_muon->Write();
SaveHistogramAsPNG(reco_start_x_muon, fv_modules);
SaveHistogramAsPNG(reco_start_y_muon, fv_modules);
SaveHistogramAsPNG(reco_start_z_muon, fv_modules);
SaveHistogramToCSV(reco_start_x_muon, fv_modules);
SaveHistogramToCSV(reco_start_y_muon, fv_modules);
SaveHistogramToCSV(reco_start_z_muon, fv_modules);

reco_start_x_pion->Write();
reco_start_y_pion->Write();
reco_start_z_pion->Write();
SaveHistogramAsPNG(reco_start_x_pion, fv_modules);
SaveHistogramAsPNG(reco_start_y_pion, fv_modules);
SaveHistogramAsPNG(reco_start_z_pion, fv_modules);
SaveHistogramToCSV(reco_start_x_pion, fv_modules);
SaveHistogramToCSV(reco_start_y_pion, fv_modules);
SaveHistogramToCSV(reco_start_z_pion, fv_modules);

reco_start_x_proton->Write();
reco_start_y_proton->Write();
reco_start_z_proton->Write();
SaveHistogramAsPNG(reco_start_x_proton, fv_modules);
SaveHistogramAsPNG(reco_start_y_proton, fv_modules);
SaveHistogramAsPNG(reco_start_z_proton, fv_modules);
SaveHistogramToCSV(reco_start_x_proton, fv_modules);
SaveHistogramToCSV(reco_start_y_proton, fv_modules);
SaveHistogramToCSV(reco_start_z_proton, fv_modules);

reco_start_x_kaon->Write();
reco_start_y_kaon->Write();
reco_start_z_kaon->Write();
SaveHistogramAsPNG(reco_start_x_kaon, fv_modules);
SaveHistogramAsPNG(reco_start_y_kaon, fv_modules);
SaveHistogramAsPNG(reco_start_z_kaon, fv_modules);
SaveHistogramToCSV(reco_start_x_kaon, fv_modules);
SaveHistogramToCSV(reco_start_y_kaon, fv_modules);
SaveHistogramToCSV(reco_start_z_kaon, fv_modules);



reco_nu_vtx_x_CC_LArFV->Write();
reco_nu_vtx_y_CC_LArFV->Write();
reco_nu_vtx_z_CC_LArFV->Write();
SaveHistogramAsPNG(reco_nu_vtx_x_CC_LArFV, fv_modules);
SaveHistogramAsPNG(reco_nu_vtx_y_CC_LArFV, fv_modules);
SaveHistogramAsPNG(reco_nu_vtx_z_CC_LArFV, fv_modules);
SaveHistogramToCSV(reco_nu_vtx_x_CC_LArFV, fv_modules);
SaveHistogramToCSV(reco_nu_vtx_y_CC_LArFV, fv_modules);
SaveHistogramToCSV(reco_nu_vtx_z_CC_LArFV, fv_modules);

true_nu_vtx_x_CC_LArFV->Write();
true_nu_vtx_y_CC_LArFV->Write();
true_nu_vtx_z_CC_LArFV->Write();
SaveHistogramAsPNG(true_nu_vtx_x_CC_LArFV, fv_modules);
SaveHistogramAsPNG(true_nu_vtx_y_CC_LArFV, fv_modules);
SaveHistogramAsPNG(true_nu_vtx_z_CC_LArFV, fv_modules);
SaveHistogramToCSV(true_nu_vtx_x_CC_LArFV, fv_modules);
SaveHistogramToCSV(true_nu_vtx_y_CC_LArFV, fv_modules);
SaveHistogramToCSV(true_nu_vtx_z_CC_LArFV, fv_modules);


std::vector<std::pair<TH1*, std::string>> histograms = {
    {reco_mult_prim_total, "reco_mult_prim_total"},
    {reco_mult_prim_shower, "reco_mult_prim_shower"},
    {reco_mult_muon, "reco_mult_muon"},
    {reco_mult_pion, "reco_mult_pion"},
    {reco_mult_proton, "reco_mult_proton"},
    {reco_mult_kaon, "reco_mult_kaon"},	
    {reco_length_prim_muon, "reco_length_prim_muon"},
    {reco_length_prim_pion, "reco_length_prim_pion"},
    {reco_length_prim_proton, "reco_length_prim_proton"},
    {reco_length_prim_kaon, "reco_length_prim_kaon"},	
    {reco_energy_prim_muon, "reco_energy_prim_muon"},
    {reco_energy_prim_pion, "reco_energy_prim_pion"},
    {reco_energy_prim_proton, "reco_energy_prim_proton"},
    {reco_energy_prim_kaon, "reco_energy_prim_kaon"},	
    {reco_cosTheta_prim_muon, "reco_cosTheta_prim_muon"},
    {reco_cosTheta_prim_pion, "reco_cosTheta_prim_pion"},
    {reco_cosTheta_prim_proton, "reco_cosTheta_prim_proton"},
    {reco_cosTheta_prim_kaon, "reco_cosTheta_prim_kaon"},
	{reco_start_x_muon, "reco_start_x_muon"},
	{reco_start_y_muon, "reco_start_y_muon"},
	{reco_start_z_muon, "reco_start_z_muon"},
	{reco_start_x_pion, "reco_start_x_pion"},
	{reco_start_y_pion, "reco_start_y_pion"},
	{reco_start_z_pion, "reco_start_z_pion"},	
	{reco_start_x_proton, "reco_start_x_proton"},
	{reco_start_y_proton, "reco_start_y_proton"},
	{reco_start_z_proton, "reco_start_z_proton"},
	{reco_start_x_kaon, "reco_start_x_kaon"},
	{reco_start_y_kaon, "reco_start_y_kaon"},
	{reco_start_z_kaon, "reco_start_z_kaon"},
	{reco_nu_vtx_x_CC_LArFV, "reco_nu_vtx_x_CC_LArFV"},
	{reco_nu_vtx_y_CC_LArFV, "reco_nu_vtx_y_CC_LArFV"},
	{reco_nu_vtx_z_CC_LArFV, "reco_nu_vtx_z_CC_LArFV"},	
	{true_nu_vtx_x_CC_LArFV, "true_nu_vtx_x_CC_LArFV"},
	{true_nu_vtx_y_CC_LArFV, "true_nu_vtx_y_CC_LArFV"},
	{true_nu_vtx_z_CC_LArFV, "true_nu_vtx_z_CC_LArFV"}
};

// Call the function with the histogram data
SaveStatisticsToCSV(fv_modules, histograms);




caf_out_file->Close();

std::cout << "Total number of True interactions (All): " << true_interactions << std::endl;
std::cout << "Total number of True interactions (CC): " << true_interactions_cc << std::endl;
std::cout << "Total number of True interactions (CC LArFV): " << true_interactions_ccLArFV << std::endl;

std::cout << "Total number of Reco interactions (All): " << reco_interactions << std::endl;
std::cout << "Total number of Reco interactions (CC): " << reco_interactions_cc << std::endl;
std::cout << "Total number of Reco interactions (CC LArFV): " << reco_interactions_ccLArFV << std::endl;
std::cout << "Total number of Reco interactions (CC LArFV + MINERVA): " << reco_interactions_ccLArFV_minerva << std::endl;



std::cout << "Efficiency of All Interactions (ccLArFV reco to all truth): " 
          << reco_interactions_ccLArFV << " / " << true_interactions 
          << " = " << static_cast<double>(reco_interactions_ccLArFV) / static_cast<double>(true_interactions) << std::endl;

std::cout << "Efficiency of CC Interactions (ccLArFV reco to cc truth): " 
          << reco_interactions_ccLArFV << " / " << true_interactions_cc 
          << " = " << static_cast<double>(reco_interactions_ccLArFV) / static_cast<double>(true_interactions_cc) << std::endl;

std::cout << "Efficiency of CC LAr Interactions (ccLArFV reco to cc truth): " 
          << reco_interactions_ccLArFV << " / " << true_interactions_ccLArFV 
          << " = " << static_cast<double>(reco_interactions_ccLArFV) / static_cast<double>(true_interactions_ccLArFV) << std::endl;

std::cout << "Efficiency of All Interactions (all reco to all truth): " 
          << reco_interactions << " / " << true_interactions 
          << " = " << static_cast<double>(reco_interactions) / static_cast<double>(true_interactions) << std::endl;

std::cout << "Efficiency of CC Interactions (cc reco to cc truth): " 
          << reco_interactions_cc << " / " << true_interactions_cc 
          << " = " << static_cast<double>(reco_interactions_cc) / static_cast<double>(true_interactions_cc) << std::endl;

std::cout << "Efficiency of CC LAr Interactions (ccLArFV reco to ccLArFV truth): " 
          << reco_interactions_ccLArFV << " / " << true_interactions_ccLArFV 
          << " = " << static_cast<double>(reco_interactions_ccLArFV) / static_cast<double>(true_interactions_ccLArFV) << std::endl;

std::cout << "Efficiency of CC LAr MINERvA Interactions (ccLArFV MINERvA reco to all truth): " 
          << reco_interactions_ccLArFV_minerva << " / " << true_interactions 
          << " = " << static_cast<double>(reco_interactions_ccLArFV_minerva) / static_cast<double>(true_interactions) << std::endl;

std::cout << "Efficiency of CC LAr MINERvA Interactions (ccLArFV MINERvA reco to ccLArFV truth): " 
          << reco_interactions_ccLArFV_minerva << " / " << true_interactions_ccLArFV 
          << " = " << static_cast<double>(reco_interactions_ccLArFV_minerva) / static_cast<double>(true_interactions_ccLArFV) << std::endl;


// Calculate total tracks from histogram
double total_tracks_corrected = 0.0;
for (int i = 1; i <= reco_mult_prim_total->GetNbinsX(); ++i) { // Start at 1 (ROOT convention)
    int adjusted_index = i - 1; // Convert ROOT index to zero-based index
    total_tracks_corrected += reco_mult_prim_total->GetBinContent(i) * adjusted_index;
    // total_tracks_corrected += reco_mult_prim_total->GetBinContent(i) * i;
}
std::cout << "Corrected total tracks from histogram: " << total_tracks_corrected << std::endl;

// Calculate total tracks from counters
int total_tracks_from_counters = total_muons + total_pions + total_protons + total_kaons;

std::cout << "Total tracks from counters: " << total_tracks_from_counters << std::endl;

// Print total counts after processing all entries
// std::cout << "Total tracks contributing to reco_mult_prim_total:" << std::endl;
// std::cout << "  Muons: " << total_muons << std::endl;
// std::cout << "  Charged Pions: " << total_pions << std::endl;
// std::cout << "  Protons: " << total_protons << std::endl;
// std::cout << "  Kaons: " << total_kaons << std::endl;
// std::cout << "  Total: " << total_muons + total_pions + total_protons + total_kaons << std::endl;

// Compare results
if (total_tracks_corrected == total_tracks_from_counters) {
    std::cout << "Histogram and counters match!" << std::endl;
} else {
    std::cerr << "Discrepancy detected: histogram = " << total_tracks_corrected
              << ", counters = " << total_tracks_from_counters << std::endl;
}

for (int i = 1; i <= reco_mult_prim_total->GetNbinsX(); ++i) {
    double bin_content = reco_mult_prim_total->GetBinContent(i);
    std::cout << "Bin " << i-1 << ": Content = " << bin_content << std::endl;
}



    return 1;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cout << "\n USAGE: " << argv[0] << " input_caf_file_list output_root_file boolTagIfData \n" << std::endl;
        return 1;
    }

    std::string input_file_list = argv[1];
    std::string output_rootfile = argv[2];
    std::string dataString=argv[3];
  bool dataOnly=true;
  if (dataString=="0") dataOnly=false;

    caf_plotter(input_file_list, output_rootfile,dataOnly);

    return 0;
}

void SaveHistogramAsPNG(TObject* histogram, const std::string& suffix) {
    if (!histogram) return;

    if (histogram->InheritsFrom(TH1::Class())) {
        TH1* hist = static_cast<TH1*>(histogram);

        hist->SetLineColor(kRed);
        hist->SetLineWidth(3);
        hist->SetLineStyle(1);

        TCanvas canvas("canvas", "Canvas", 800, 600);
        gStyle->SetOptStat("emr");

        hist->Draw("HIST E");

        std::string fileName = std::string(hist->GetName()) + "_" + suffix + ".png";
        canvas.SaveAs(fileName.c_str());
    } else if (histogram->InheritsFrom(TH2::Class())) {
        TH2* hist = static_cast<TH2*>(histogram);

        TCanvas canvas("canvas", "Canvas", 800, 600);
        gStyle->SetOptStat("emr");

        hist->Draw("COLZ");

        std::string fileName = std::string(hist->GetName()) + "_" + suffix + ".png";
        canvas.SaveAs(fileName.c_str());
    } else {
        std::cerr << "The provided object is not a histogram." << std::endl;
    }
}


void CustomizeStatBox(TH1D* hist, const std::string& suffix, int total_tracks, double weighted_mean, double weighted_stddev) {
    if (!hist) return;

    TCanvas canvas("canvas", "Canvas", 800, 600);
    gStyle->SetOptStat(0); // Turn off default stats

    hist->SetLineColor(kRed);
    hist->SetLineWidth(3);
    hist->SetLineStyle(1);

    hist->Draw("HIST E");

    // Add a custom stats box
    TPaveText* stats_box = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
    stats_box->SetBorderSize(1);
    stats_box->SetFillColor(0);

    // Set text font, size, and alignment
    stats_box->SetTextFont(42);  // ROOT default font
    stats_box->SetTextSize(0.06); // Increase font size (default is ~0.03)
    stats_box->SetTextAlign(12); // Left-aligned

    // Add custom statistics
    std::stringstream ss;
    ss << "Entries: " << total_tracks;
    stats_box->AddText(ss.str().c_str());

    ss.str("");
    ss << "Mean: " << std::fixed << std::setprecision(2) << weighted_mean;
    stats_box->AddText(ss.str().c_str());

    ss.str("");
    ss << "Std Dev: " << std::fixed << std::setprecision(2) << weighted_stddev;
    stats_box->AddText(ss.str().c_str());

    // Reduce spacing between lines
    stats_box->SetMargin(0.1); // Decrease left/right padding within the box
    stats_box->SetTextSize(0.03);

    stats_box->Draw();

    std::string fileName = std::string(hist->GetName()) + "_" + suffix + ".png";
    canvas.SaveAs(fileName.c_str());
}


std::tuple<double, double, double, double> ProcessHistogram(TH1D* hist, const std::string& suffix) {
    if (!hist) return {0.0, 0.0, 0.0, 0.0};

    double total_tracks_corrected = 0.0;
    double weighted_sum = 0.0;
    double weighted_variance_sum = 0.0;
    double total_entries = 0.0;

    // Loop over bins to calculate total tracks, weighted mean, and variance
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double bin_content = hist->GetBinContent(i);
        double bin_center = hist->GetBinCenter(i);
        total_tracks_corrected += bin_content * (i - 1); // Adjust for zero-based index
        weighted_sum += bin_content * bin_center;
        total_entries += bin_content;
    }

    // Compute weighted mean
    double weighted_mean = total_entries > 0 ? weighted_sum / total_entries : 0.0;

    // Compute weighted variance
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double bin_content = hist->GetBinContent(i);
        double bin_center = hist->GetBinCenter(i);
        weighted_variance_sum += bin_content * std::pow(bin_center - weighted_mean, 2);
    }

    double weighted_stddev = total_entries > 0 ? std::sqrt(weighted_variance_sum / total_entries) : 0.0;

    // Customize the stat box (if needed)
    CustomizeStatBox(hist, suffix, static_cast<int>(total_tracks_corrected), weighted_mean, weighted_stddev);

    // Return calculated statistics
    return {total_tracks_corrected, total_entries, weighted_mean, weighted_stddev};
}


void SaveHistogramToCSV(TH1* histogram, const std::string& suffix) {
    if (!histogram) return;

    std::string fileName = std::string(histogram->GetName()) + "_" + suffix + ".csv";
    std::ofstream csvFile(fileName);

    if (!csvFile.is_open()) {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        return;
    }

    csvFile << "Bin Number,Bin Center,Entries,Error,Lower End,Upper End\n";

    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        int binNumber = i;
        double binCenter = histogram->GetBinCenter(i);
        double binContent = histogram->GetBinContent(i);
        double binError = histogram->GetBinError(i);
        double lowerEnd = binContent - binError;
        double upperEnd = binContent + binError;

        csvFile << binNumber << "," << binCenter << "," << binContent << "," << binError << "," << lowerEnd << "," << upperEnd << "\n";
    }

    csvFile.close();
    std::cout << "Histogram data saved to: " << fileName << std::endl;
}


// Modified function to accept histogram data as a parameter
void SaveStatisticsToCSV(const std::string& suffix, const std::vector<std::pair<TH1*, std::string>>& histograms) {
    std::string fileName = "stats_" + suffix + ".csv";
    std::ofstream csvFile(fileName);

    if (!csvFile.is_open()) {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        return;
    }

    // Write headers
    csvFile << "histogram,entries,mean,std_dev\n";

    // Write histogram statistics
    for (const auto& hist_pair : histograms) {
        TH1* histogram = hist_pair.first;
        const std::string& name = hist_pair.second;

        if (histogram) {
            csvFile << name << ","
                    << std::fixed << std::setprecision(2)
                    << histogram->GetEntries() << ","
                    << histogram->GetMean() << ","
                    << histogram->GetStdDev() << "\n";
        }
    }

    csvFile.close();
    std::cout << "Statistics saved to: " << fileName << std::endl;
}

// Modified function to accept histogram data as a parameter
// void SaveStatisticsToCSV(const std::string& suffix, const std::vector<std::pair<TH1*, std::string>>& histograms) {
    // std::string fileName = "stats_" + suffix + ".csv";
    // std::ofstream csvFile(fileName);

    // if (!csvFile.is_open()) {
        // std::cerr << "Failed to open file: " << fileName << std::endl;
        // return;
    // }

    // Write headers
    // csvFile << "histogram,entries,mean,std_dev\n";

    // for (const auto& hist_pair : histograms) {
        // TH1* histogram = hist_pair.first;
        // const std::string& name = hist_pair.second;

        // if (histogram) {
            // Attempt to compute weighted statistics
            // auto histD = dynamic_cast<TH1D*>(histogram);
            // if (histD) {
                // Replace with weighted statistics if available
                // auto [weighted_total, entries, weighted_mean, weighted_stddev] = ProcessHistogram(histD, suffix);

                // csvFile << name << ","
                        // << std::fixed << std::setprecision(2)
                        // << weighted_total << ","  // Replace Entries
                        // << weighted_mean << ","   // Replace Mean
                        // << weighted_stddev << "\n"; // Replace StdDev
            // } else {
                // Fall back to basic statistics
                // csvFile << name << ","
                        // << std::fixed << std::setprecision(2)
                        // << histogram->GetEntries() << ","
                        // << histogram->GetMean() << ","
                        // << histogram->GetStdDev() << "\n";
            // }
        // }
    // }

    // csvFile.close();
    // std::cout << "Statistics (including weighted) saved to: " << fileName << std::endl;
// }

