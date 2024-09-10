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

void SaveHistogramAsPNG(TObject* histogram);
void Save2DHistogramAsPNG(TH2* histogram, const char* fileName);
void SaveHistogramToCSV(TH1D* histogram, const char* fileName);

int caf_plotter(std::string input_file_list, std::string output_rootfile, bool dataOnly) {
    int start_bins = 50;
    int start_bins_L = -70;
    int start_bins_U = 70;
	
	int mult_bins_genie = 20;
    int mult_edge_L_genie = 0;
    int mult_edge_U_genie = 20;

    int totalInteractions = 0;
	int CCinteractions = 0;
	int CCinteractionsWithinLArFV = 0;
       double mnvOffsetX=0; double mnvOffsetY=0;
       if (dataOnly){ mnvOffsetX=-10; mnvOffsetY=5;}


    TH1D *reco_nu_vtx_x = new TH1D("reco_nu_vtx_x", "Neutrino Vertex-x (sandbox 0.9s_match);Vertex-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_y = new TH1D("reco_nu_vtx_y", "Neutrino Vertex-y (sandbox 0.9s_match);Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_z = new TH1D("reco_nu_vtx_z", "Neutrino Vertex-z (sandbox 0.9s_match);Vertex-z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1D *reco_nu_vtx_x_CC = new TH1D("reco_nu_vtx_x_CC", "Neutrino Vertex-x (CC) (sandbox 0.9s_match);Vertex-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_y_CC = new TH1D("reco_nu_vtx_y_CC", "Neutrino Vertex-y (CC) (sandbox 0.9s_match);Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_z_CC = new TH1D("reco_nu_vtx_z_CC", "Neutrino Vertex-z (CC) (sandbox 0.9s_match);Vertex-z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1D *reco_nu_vtx_x_CC_LArFV = new TH1D("reco_nu_vtx_x_CC_LArFV", "Neutrino Vertex-x (within LArFV) (sandbox 0.9s_match);Vertex-x (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_y_CC_LArFV = new TH1D("reco_nu_vtx_y_CC_LArFV", "Neutrino Vertex-y (within LArFV) (sandbox 0.9s_match);Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U);
    TH1D *reco_nu_vtx_z_CC_LArFV = new TH1D("reco_nu_vtx_z_CC_LArFV", "Neutrino Vertex-z (within LArFV) (sandbox 0.9s_match);Vertex-z (cm);Counts", start_bins, start_bins_L, start_bins_U);

    TH1D *reco_mult_prim_total = new TH1D("reco_mult_prim_total", "Reconstructed Data Multiplicity;Multiplicity;Counts", mult_bins_genie, mult_edge_L_genie, mult_edge_U_genie);

    TH2D *reco_vtx_x_vs_y_CC_LArFV = new TH2D("reco_vtx_x_vs_y_CC_LArFV", "Neutrino Vertex x vs y (within LArFV);Vertex-x (cm);Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U, start_bins, start_bins_L, start_bins_U);
    TH2D *reco_vtx_y_vs_z_CC_LArFV = new TH2D("reco_vtx_y_vs_z_CC_LArFV", "Neutrino Vertex z vs y (within LArFV);Vertex-z (cm);Vertex-y (cm);Counts", start_bins, start_bins_L, start_bins_U, start_bins, start_bins_L, start_bins_U);
    TH2D *reco_vtx_z_vs_x_CC_LArFV = new TH2D("reco_vtx_z_vs_x_CC_LArFV", "Neutrino Vertex z vs x (within LArFV);Vertex-z (cm);Vertex-x (cm);Counts", start_bins, start_bins_L, start_bins_U, start_bins, start_bins_L, start_bins_U);




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


    for (long n = 0; n < Nentries; n++) {
        if (n % 10000 == 0) std::cout << "Processing trigger " << n << " of " << Nentries << std::endl;
        caf_chain->GetEntry(n);

        for (long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++) {
            totalInteractions++; // Increment total number of interactions
			
			
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

            reco_nu_vtx_x->Fill(reco_nu_vtxX);
            reco_nu_vtx_y->Fill(reco_nu_vtxY);
            reco_nu_vtx_z->Fill(reco_nu_vtxZ);

			// Conditions for All Modules			
			if (abs(reco_nu_vtxX) > maxX || abs(reco_nu_vtxX) < minX) continue;
			if (abs(reco_nu_vtxY) > maxY || abs(reco_nu_vtxY) < minY) continue;
			if (abs(reco_nu_vtxZ) > maxZ || abs(reco_nu_vtxZ) < minZ) continue;
			if (reco_nu_vtxX > Xinner_neg && reco_nu_vtxX < Xinner_pos) continue; 
			if (reco_nu_vtxZ > Zinner_neg && reco_nu_vtxZ < Zinner_pos) continue; 

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
				}	// if ((abs(start_pos.z) > maxZ || abs(end_pos.z) > maxZ) ) {

				reco_nu_vtx_x_CC->Fill(reco_nu_vtxX);
				reco_nu_vtx_y_CC->Fill(reco_nu_vtxY);
				reco_nu_vtx_z_CC->Fill(reco_nu_vtxZ);	

				reco_vtx_x_vs_y_CC_LArFV->Fill(reco_nu_vtxX, reco_nu_vtxY);
				reco_vtx_y_vs_z_CC_LArFV->Fill(reco_nu_vtxZ, reco_nu_vtxY);
				reco_vtx_z_vs_x_CC_LArFV->Fill(reco_nu_vtxZ, reco_nu_vtxX);


				if (length > 1) {
					reco_multiplicity++;				
				}	// if (length > 1) {
			}	// for (size_t k = 0; k < sr->common.ixn.dlp[nixn].part.dlp.size(); k++) {

			if (reco_multiplicity > 0 && maxDotProductDS > 0.99) {
				CCinteractionsWithinLArFV++;		
				reco_mult_prim_total->Fill(reco_multiplicity);
			}	// if (reco_multiplicity > 0 && maxDotProductDS > 0.99) {
        }	// for (long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++) {
    }	// for (long n = 0; n < Nentries; n++) {

    TFile *caf_out_file = new TFile(output_rootfile.c_str(), "recreate");

    reco_nu_vtx_x->Write();
    reco_nu_vtx_y->Write();
    reco_nu_vtx_z->Write();

    reco_nu_vtx_x_CC->Write();
    reco_nu_vtx_y_CC->Write();
    reco_nu_vtx_z_CC->Write();

    reco_nu_vtx_x_CC_LArFV->Write();
    reco_nu_vtx_y_CC_LArFV->Write();
    reco_nu_vtx_z_CC_LArFV->Write();
	
	SaveHistogramAsPNG(reco_nu_vtx_x_CC);
	SaveHistogramAsPNG(reco_nu_vtx_y_CC);
	SaveHistogramAsPNG(reco_nu_vtx_z_CC);
	
	reco_vtx_x_vs_y_CC_LArFV->Write();
	reco_vtx_y_vs_z_CC_LArFV->Write();
	reco_vtx_z_vs_x_CC_LArFV->Write();
	Save2DHistogramAsPNG(reco_vtx_x_vs_y_CC_LArFV, "reco_vtx_x_vs_y_CC_LArFV.png");
	Save2DHistogramAsPNG(reco_vtx_y_vs_z_CC_LArFV, "reco_vtx_y_vs_z_CC_LArFV.png");
	Save2DHistogramAsPNG(reco_vtx_z_vs_x_CC_LArFV, "reco_vtx_z_vs_x_CC_LArFV.png");

	
	reco_mult_prim_total->Write();
	SaveHistogramAsPNG(reco_mult_prim_total);
	SaveHistogramToCSV(reco_mult_prim_total, "reco_mult_prim_total.csv");

	
	

    caf_out_file->Close();

    std::cout << "Total number of interactions (All): " << totalInteractions << std::endl;
    std::cout << "Total number of interactions (CC): " << CCinteractions << std::endl;	
    std::cout << "Total number of interactions (CC LArFV): " << CCinteractionsWithinLArFV << std::endl;

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

void SaveHistogramAsPNG(TObject* histogram) {
    if (!histogram) return;

    if (histogram->InheritsFrom(TH1::Class())) {
        TH1* hist = static_cast<TH1*>(histogram);

        hist->SetLineColor(kRed);
        hist->SetLineWidth(3);
        hist->SetLineStyle(1);

        TCanvas canvas("canvas", "Canvas", 800, 600);
        gStyle->SetOptStat("emr");

        hist->Draw();

        std::string fileName = std::string(hist->GetName()) + ".png";
        canvas.SaveAs(fileName.c_str());
    } else if (histogram->InheritsFrom(TH2::Class())) {
        TH2* hist = static_cast<TH2*>(histogram);

        TCanvas canvas("canvas", "Canvas", 800, 600);
        gStyle->SetOptStat("emr");

        hist->Draw("COLZ");

        std::string fileName = std::string(hist->GetName()) + ".png";
        canvas.SaveAs(fileName.c_str());
    } else {
        std::cerr << "The provided object is not a histogram." << std::endl;
    }
}


void Save2DHistogramAsPNG(TH2* histogram, const char* fileName) {
    if (!histogram) return;

    TCanvas canvas("canvas", "Canvas", 800, 600);
    gStyle->SetOptStat(0); // Turn off statistics box

    histogram->Draw("COLZ");

    canvas.SaveAs(fileName);
}


void SaveHistogramToCSV(TH1D* histogram, const char* fileName) {
    // Open a file in write mode
    std::ofstream csvFile(fileName);
    
    // Check if the file is open
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        return;
    }

    // Write the header to the CSV
    csvFile << "Bin Number,Entries\n";

    // Loop over histogram bins and save bin number and bin content
    for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
        int binNumber = i;
        double binContent = histogram->GetBinContent(i);
        csvFile << binNumber << "," << binContent << "\n";
    }

    // Close the file
    csvFile.close();
    std::cout << "Histogram data saved to: " << fileName << std::endl;
}