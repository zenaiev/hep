#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <cassert>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Boost.h>

// Function to generate a two-body decay
std::tuple<ROOT::Math::PxPyPzMVector, ROOT::Math::PxPyPzMVector> two_body_decay(const ROOT::Math::PxPyPzMVector &parent, double m1, double m2, std::default_random_engine &generator) {
    double m0 = parent.M(); // Parent mass

    // Calculate momenta and energies
    double p12 = std::sqrt((m0 * m0 - (m1 + m2) * (m1 + m2)) * (m0 * m0 - (m1 - m2) * (m1 - m2))) / (2 * m0);
    double E1 = (m0 * m0 + m1 * m1 - m2 * m2) / (2 * m0);
    double p = std::sqrt(E1 * E1 - m1 * m1);

    // Generate random angles for isotropic decay
    std::uniform_real_distribution<double> uniform_theta(0, M_PI);
    std::uniform_real_distribution<double> uniform_phi(0, 2 * M_PI);
    double theta = uniform_theta(generator);
    double phi = uniform_phi(generator);

    // Momentum components for daughter 1
    double p1x = p * sin(theta) * cos(phi);
    double p1y = p * sin(theta) * sin(phi);
    double p1z = p * cos(theta);

    // Create 4-momentum vectors for the daughters
    ROOT::Math::PxPyPzMVector daughter1(p1x, p1y, p1z, m1);
    ROOT::Math::PxPyPzMVector daughter2(-p1x, -p1y, -p1z, m2);

    // Boost daughters back to the lab frame if the parent is moving
    ROOT::Math::Boost boost_vector(parent.BoostToCM());
    daughter1 = boost_vector(daughter1);
    daughter2 = boost_vector(daughter2);

    return std::make_tuple(daughter1, daughter2);
}

// Function to smear particle momentum
ROOT::Math::PxPyPzMVector smear_momentum(const ROOT::Math::PxPyPzMVector &particle, std::default_random_engine &generator) {
    std::normal_distribution<double> normal_dist(1.0, 0.05);
    double fluct = normal_dist(generator);
    return ROOT::Math::PxPyPzMVector(particle.X() * fluct, particle.Y() * fluct, particle.Z() * fluct, particle.M());
}

// Function to generate random values for momentum and direction
std::vector<ROOT::Math::PxPyPzMVector> generate_and_decay(double parent_mass, double daughter_mass_1, double daughter_mass_2, std::default_random_engine &generator) {
    std::exponential_distribution<double> exp_dist(1.0);
    std::uniform_real_distribution<double> uniform_theta(0, M_PI);
    std::uniform_real_distribution<double> uniform_phi(0, 2 * M_PI);

    double momentum = exp_dist(generator);
    double theta = uniform_theta(generator);
    double phi = uniform_phi(generator);

    // Create 3-momentum using spherical coordinates
    double px = momentum * sin(theta) * cos(phi);
    double py = momentum * sin(theta) * sin(phi);
    double pz = momentum * cos(theta);

    // Create the parent particle
    ROOT::Math::PxPyPzMVector parent(px, py, pz, parent_mass);

    // Generate daughter particles via two-body decay
    ROOT::Math::PxPyPzMVector daughter1, daughter2;
    std::tie(daughter1, daughter2) = two_body_decay(parent, daughter_mass_1, daughter_mass_2, generator);

    // Smear momentum of daughter tracks
    daughter1 = smear_momentum(daughter1, generator);
    daughter2 = smear_momentum(daughter2, generator);

    return {daughter1, daughter2};
}


// compile using gcc (notice optimisation flag -O2):
// g++ -O2 -o invmass `root-config --libs --cflags` -lGenVector invmass.cpp
// or uncomment int invmass(), comment out int main() and run in ROOT interpretator:
// root -l -q invmass.cpp
// or compile with ROOT interpretator (notice + in the end)
// root -l -q invmass.cpp+
int main() {
//int invmass() {
    // Random seed for reproducibility
    std::default_random_engine generator(42);

    // Particle masses in GeV
    double mass_pi_ch = 0.13957;
    double mass_k_zero = 0.497611;
    double mass_d_zero = 1.86484;

    // Create a histogram to store generated momenta
    TH1F *hInvMass = new TH1F("hInvMass", "Invariant Mass", 300, 0, 3);

    // Create tree (store events)
    //TFile *fileout = new TFile("tracks.root", "recreate");
    //TTree *tree = new TTree("tree", "Tree with tracks");
    //std::vector<ROOT::Math::PxPyPzMVector> tracks_vec;
    //tree->Branch("tracks", &tracks_vec);

    int nParticles = 10000;

    for (int i = 0; i < nParticles; ++i) {
        std::vector<ROOT::Math::PxPyPzMVector> tracks;

        // Generate decays and collect tracks
        auto decay1 = generate_and_decay(mass_k_zero, mass_pi_ch, mass_pi_ch, generator);
        tracks.push_back(decay1[0]);
        tracks.push_back(decay1[1]);

        auto decay2 = generate_and_decay(mass_d_zero, mass_pi_ch, mass_pi_ch, generator);
        tracks.push_back(decay2[0]);
        tracks.push_back(decay2[1]);

        assert(tracks.size() == 4);

        // Fill invariant mass histogram
        for (size_t itr1 = 0; itr1 < tracks.size(); ++itr1) {
            for (size_t itr2 = itr1 + 1; itr2 < tracks.size(); ++itr2) {
                hInvMass->Fill((tracks[itr1] + tracks[itr2]).M());
            }
        }

        // Fill tree with tracks
        //tracks_vec = tracks;
        //tree->Fill();
    }

    // Write TTree to file
    //fileout->cd();
    //tree->Write();
    //fileout->Close();

    // Fit function
    TF1 *fitFunc = new TF1("fitFunc",
                           "[0]/(sqrt(2*TMath::Pi())*[2])*TMath::Exp(-0.5*((x-[1])/[2])^2)"
                           " + [3]/(sqrt(2*TMath::Pi())*[5])*TMath::Exp(-0.5*((x-[4])/[5])^2)"
                           " + [6] + [7]*x + [8]*x^2 + [9]*x^3 + [10]*x^4 + [11]*x^5",
                           -10, 10);

    fitFunc->SetParameters(1, 0.5, 0.01, 1, 1.85, 0.05, 0, 0, 0, 0, 0);
    fitFunc->SetParLimits(1, 0.48, 0.52); // K0 mass
    fitFunc->SetParLimits(2, 0, 0.1);     // K0 width
    fitFunc->SetParLimits(4, 1.8, 1.9);   // D0 mass
    fitFunc->SetParLimits(5, 0, 0.2);     // D0 width

    hInvMass->Fit(fitFunc);

    // Output signal events
    double bin_width = hInvMass->GetBinLowEdge(2) - hInvMass->GetBinLowEdge(1);
    std::cout << "Number of signal events #1 = " << fitFunc->GetParameter(0) / bin_width
              << " +- " << fitFunc->GetParError(0) / bin_width << std::endl;
    std::cout << "Number of signal events #2 = " << fitFunc->GetParameter(3) / bin_width
              << " +- " << fitFunc->GetParError(3) / bin_width << std::endl;

    // Draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Invariant Mass", 600, 600);
    hInvMass->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV]");
    hInvMass->GetYaxis()->SetTitle("Events");
    hInvMass->Draw();

    // Save the canvas as a PDF/PNG
    canvas->SaveAs("invmass.pdf");
    canvas->SaveAs("invmass.png");

    return 0;
}
