// compile:
// g++ unfold.cpp -o unfold `root-config --libs --cflags` -I /home/zenaiev/soft/RooUnfold/build/ -I /home/zenaiev/soft/RooUnfold/src -lRooUnfold -L /home/zenaiev/soft/RooUnfold/build
// export LD_LIBRARY_PATH=/home/zenaiev/soft/RooUnfold/build:$LD_LIBRARY_PATH
// run:
// ./unfold

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldIds.h"
#include <iostream>
#include <random>
#include <vector>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TFile.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>
#endif

using namespace std;

double smear(double momentum) {
    static default_random_engine generator(42); // random seed for reproducibility
    normal_distribution<double> distribution(0.75, 0.25);
    double fluct = distribution(generator);
    return momentum * fluct;
}

void generate(int nevents, double lambda_exp,
              RooUnfoldResponse* response = nullptr,
              TH1F* h_true = nullptr,
              TH1F* h_reco = nullptr) {
    static default_random_engine generator(42); // random seed for reproducibility
    exponential_distribution<double> exp_dist(lambda_exp);

    for (int i = 0; i < nevents; ++i) {
        double momentum_true = exp_dist(generator);
        double momentum_measured = smear(momentum_true);

        if (response) {
            response->Fill(momentum_measured, momentum_true);
        }
        if (h_reco) {
            h_reco->Fill(momentum_measured);
        }
        if (h_true) {
            h_true->Fill(momentum_true);
        }
    }
}

int main() {
//int unfold() {
    // Initialize ROOT application
    //TApplication app("app", nullptr, nullptr);

    // Draw smearing distribution
    TCanvas* fig = new TCanvas("canvas", "Smearing", 800, 600);
    TH1F* h_smear_02 = new TH1F("h_smear_02", "Smear for p = 0.2 GeV", 50, 0, 3.);
    TH1F* h_smear_2 = new TH1F("h_smear_2", "Smear for p = 2 GeV", 50, 0, 3.);

    for (int i = 0; i < 10000; ++i) {
        h_smear_02->Fill(smear(0.2));
        h_smear_2->Fill(smear(2.0));
    }

    h_smear_02->SetLineColor(kBlue);
    h_smear_2->SetLineColor(kRed);
    h_smear_02->Draw();
    h_smear_2->Draw("SAME");

    TLine* line1 = new TLine(0.2, 0, 0.2, h_smear_02->GetMaximum());
    line1->SetLineColor(kBlack);
    line1->Draw();

    TLine* line2 = new TLine(2.0, 0, 2.0, h_smear_2->GetMaximum());
    line2->SetLineColor(kBlack);
    line2->Draw();

    fig->SaveAs("smeared.pdf");
    fig->SaveAs("smeared.png");

    // Generate events
    double lambda_exp = 2.0;
    RooUnfoldResponse response(120, 0., 3., 30, 0., 3.);
    generate(100000, lambda_exp, &response);

    TH1F* h_data = new TH1F("h_data", "Data", 120, 0., 3.);
    TH1F* h_true = new TH1F("h_true", "True", 30, 0., 3.);
    generate(100000, lambda_exp, nullptr, h_true, h_data);

    // Unfolding (unregularised and regularised)
    auto unfold_noreg = new RooUnfoldTUnfold(&response, h_data, 0.);
    auto h_unfold_noreg = unfold_noreg->Hunfold();

    auto unfold_reg = new RooUnfoldTUnfold(&response, h_data);
    //unfold_reg->OptimiseTau();
    auto h_unfold_reg = unfold_reg->Hunfold();

    // Plotting response matrix
    TCanvas* c_2d = new TCanvas("canvas_2d", "Response", 600, 600);
    auto h_response = response.Hresponse();
    h_response->SetStats(0);
    h_response->GetXaxis()->SetTitle("momentum [GeV] (reco)");
    h_response->GetYaxis()->SetTitle("momentum [GeV] (true)");
    h_response->Draw("COLZ");
    c_2d->SaveAs("response.pdf");
    c_2d->SaveAs("response.png");

    // Plot true, reconstructed and unfolded distributions
    TCanvas* c_1d = new TCanvas("canvas_1d", "Unfolding", 300, 600);
    c_1d->Divide(1, 3);

    h_data->SetLineColor(kBlue);
    h_true->SetLineColor(kBlack);
    h_unfold_noreg->SetLineColor(kRed);
    h_unfold_reg->SetLineColor(kMagenta);

    // Rescale data for plotting
    h_data->Scale(h_true->Integral("width") / h_data->Integral("width"));

    auto draw_all = [&]() {
        h_data->Draw("hist");
        h_true->Draw("hist same");
        h_unfold_noreg->Draw("E0 same");
        h_unfold_reg->Draw("E0 same");
    };

    c_1d->cd(1);
    draw_all();

    TLegend* legend = new TLegend(0.5, 0.6, 0.9, 0.9);
    legend->AddEntry(h_true, "True", "l");
    legend->AddEntry(h_data, "Reco", "l");
    legend->AddEntry(h_unfold_noreg, "Reco unfolded (unregularised)", "l");
    legend->AddEntry(h_unfold_reg, "Reco unfolded (regularised)", "l");
    legend->Draw();

    c_1d->cd(2);
    draw_all();
    gPad->SetLogy();

    c_1d->cd(3);
    TH1F* h_unfold_noreg_ratio = static_cast<TH1F*>(h_unfold_noreg->Clone("h_unfold_noreg_ratio"));
    h_unfold_noreg_ratio->Divide(h_true);
    TH1F* h_unfold_reg_ratio = static_cast<TH1F*>(h_unfold_reg->Clone("h_unfold_reg_ratio"));
    h_unfold_reg_ratio->Divide(h_true);

    h_unfold_noreg_ratio->GetXaxis()->SetTitle("momentum [GeV]");
    h_unfold_noreg_ratio->GetYaxis()->SetTitle("Events");
    h_unfold_noreg_ratio->SetStats(0);
    h_unfold_noreg_ratio->Draw("E0");
    h_unfold_reg_ratio->Draw("E0 same");

    c_1d->SaveAs("unfolding.pdf");
    c_1d->SaveAs("unfolding.png");

    //app.Run(); // Run the application

    delete h_smear_02;
    delete h_smear_2;
    delete h_data;
    delete h_true;
    delete h_unfold_noreg;
    delete h_unfold_reg;
    delete c_2d;
    delete c_1d;

    return 0;
}
