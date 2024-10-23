import ROOT
#import RooUnfold
import tqdm
import random
from matplotlib import pyplot as plt

def smear(momentum):
  mean = 0.75
  sigma = 0.25
  fluct = random.normalvariate(mean, sigma)
  return momentum * fluct

def generate(nevents, lambda_exp, response=None, h_true=None, h_reco=None):
  for i in tqdm.tqdm(range(nevents)):
    momentum_true = random.expovariate(lambda_exp)
    momentum_measured = smear(momentum_true)
    #print(momentum_true, momentum_measured)
    if response is not None:
      response.Fill(momentum_measured, momentum_true)
    if h_reco is not None:
      h_reco.Fill(momentum_measured)
    if h_true is not None:
      h_true.Fill(momentum_true)

if __name__ == '__main__':
  # random seed for reproducible results
  random.seed(42)

  # draw smearing distribution
  fig, ax = plt.subplots()
  ax.hist([smear(0.2) for i in range(10000)], bins=50, label='reco p = 0.2 GeV')
  ax.axvline(x=0.2, color='red', label='true p = 0.2 GeV')
  ax.hist([smear(2.0) for i in range(10000)], bins=50, label='reco p = 2 GeV', alpha=0.5)
  ax.axvline(x=2.0, color='orange', label='true p = 2 GeV')
  ax.legend()
  ax.set_xlabel('momentum [GeV]')
  #plt.show() # display immediately
  fig.savefig('smeared.pdf')
  fig.savefig('smeared.png')
    
  # generate events
  lambda_exp = 2.0
  response = ROOT.RooUnfoldResponse(120, 0., 3., 30, 0., 3.)
  generate(100000, lambda_exp, response=response)
  h_data = ROOT.TH1F('h_data', 'Data', 120, 0., 3.)
  h_true = ROOT.TH1F('h_true', 'True', 30, 0., 3.)
  generate(100000, lambda_exp, h_reco=h_data, h_true=h_true)
  
  # unregularised unfolding (TUnfold)
  unfold_noreg = ROOT.RooUnfoldTUnfold (response, h_data)
  unfold_noreg.FixTau(0.)
  h_unfold_noreg = unfold_noreg.Hunfold()

  # regularised unfolding (TUnfold)
  unfold_reg = ROOT.RooUnfoldTUnfold (response, h_data)
  #unfold_noreg.FixTau(0.01)
  unfold_reg.OptimiseTau() # by default it will do L curve scan
  h_unfold_reg = unfold_reg.Hunfold()

  # plot response matrix
  c_2d = ROOT.TCanvas("canvas_2d", "Response", 600, 600)
  c_2d.SetMargin(0.1, 0.2, 0.1, 0.1)
  h_response = response.Hresponse()
  h_response.SetStats(0) # suppress statistics box
  h_response.GetXaxis().SetTitle('momentum [GeV] (reco)')
  h_response.GetYaxis().SetTitle('momentum [GeV] (true)')
  h_response.Draw('colz')
  c_2d.SaveAs('response.pdf')
  c_2d.SaveAs('response.png')

  # plot true, reconstructed and unfolded distributions
  c_1d = ROOT.TCanvas("canvas_1d", "Unfolding", 300, 600)
  c_1d.Divide(1,3)
  h_data.GetXaxis().SetTitle('momentum [GeV]')
  h_data.GetYaxis().SetTitle('Events')
  h_data.SetStats(0)
  h_true.SetLineColor(ROOT.kBlack)
  h_data.SetLineColor(ROOT.kBlue)
  h_unfold_noreg.SetMarkerColor(ROOT.kRed)
  h_unfold_noreg.SetLineColor(ROOT.kRed)
  h_unfold_reg.SetMarkerColor(ROOT.kMagenta)
  h_unfold_reg.SetLineColor(ROOT.kMagenta)
  h_data.Scale(h_true.Integral('width')/h_data.Integral('width')) # rescale data for plot, because the number of bins is different
  #h_data.SetMinimum(-1*h_data.GetMaximum())
  def draw_all():
    h_data.Draw('hist')
    h_true.Draw('hist same')
    h_unfold_noreg.Draw('e0 same')
    h_unfold_reg.Draw('e0 same')
  c_1d.cd(1)
  draw_all()
  legend = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)
  legend.AddEntry(h_true, 'True', 'l')
  legend.AddEntry(h_data, 'Reco', 'l')
  legend.AddEntry(h_unfold_noreg, 'Reco unfolded (unregulaarised)', 'l')
  legend.AddEntry(h_unfold_reg, 'Reco unfolded (regulaarised)', 'l')
  legend.Draw()
  c_1d.cd(2)
  draw_all()
  ROOT.gPad.SetLogy()
  c_1d.cd(3)
  h_unfold_noreg_ratio = h_unfold_noreg.Clone()
  h_unfold_noreg_ratio.Divide(h_true)
  h_unfold_reg_ratio = h_unfold_reg.Clone()
  h_unfold_reg_ratio.Divide(h_true)
  #h_unfold_noreg_ratio.GetYaxis().SetRangeUser(0.5,1.5)
  h_unfold_noreg_ratio.GetXaxis().SetTitle('momentum [GeV]')
  h_unfold_noreg_ratio.GetYaxis().SetTitle('Events')
  h_unfold_noreg_ratio.SetStats(0)
  h_unfold_noreg_ratio.Draw('e0')
  h_unfold_reg_ratio.Draw('e0 same')
  c_1d.SaveAs('unfolding.pdf')
  c_1d.SaveAs('unfolding.png')
