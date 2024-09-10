# hello.py

import ROOT

if __name__ == '__main__':
  print("Hello world")
  hist = ROOT.TH1F("hist", "My histo", 300, 0., 3.)
  hist.Fill(1.2)
  hist.Print()
