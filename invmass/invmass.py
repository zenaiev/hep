import ROOT
import math
import random
import tqdm

def generate_and_decay(parent_mass, daughter_mass_1, daughter_mass_2):
  # Generate random values for momentum and direction
  momentum = random.expovariate(1.0)  # Exponential distribution for momentum
  theta = random.uniform(0, math.pi)       # Uniform in [0, pi]
  phi = random.uniform(0, 2 * math.pi)     # Uniform in [0, 2pi]

  # Create a 3-momentum vector using spherical coordinates
  px = momentum * math.sin(theta) * math.cos(phi)
  py = momentum * math.sin(theta) * math.sin(phi)
  pz = momentum * math.cos(theta)

  # Create parent particle
  parent = ROOT.Math.PxPyPzMVector(px, py, pz, parent_mass) 

  # Decay
  daughter1, daughter2 = two_body_decay(parent, daughter_mass_1, daughter_mass_2)

  # Smear momentum of daughter tracks
  daughter1 = smear_momentum(daughter1)
  daughter2 = smear_momentum(daughter2)

  return [daughter1, daughter2]

# Function to smear particle momentum
def smear_momentum(particle):
  fluct = random.normalvariate(1., 0.05) # Normal distribution with mean = 1 and width = 0.05 (5% resolution)
  return ROOT.Math.PxPyPzMVector(particle.X() * fluct, particle.Y() * fluct, particle.Z() * fluct, particle.M()) 

# Function to generate a two-body decay in the rest frame of the parent particle
def two_body_decay(parent, m1, m2):
  # Calculate the momentum of the daughter particles in the parent's rest frame
  m0 = parent.M()  # Parent mass

  # Energy of the daughters in the rest frame
  # In the rest frame, the daughters have equal momenta (p12)
  p12 = math.sqrt((m0**2 - (m1 + m2)**2) * (m0**2 - (m1 - m2)**2)) / (2 * m0)
  # E0 = E1 + E2
  # E0 = m0
  # E1 = sqrt(m1**2 + p12**2)
  # E2 = sqrt(m2**2 + p12**2)
  # m0 = sqrt(m1**2 + p12**2) + sqrt(m2**2 + p12**2)
  # m0**2 = m1**2 + m2**2 + 2 * p12**2 + 2 * sqrt((m1**2 + p12**2) * (m2**2 + p12**2))

  # p12**2 = (E1**2 + E2**2 - m1**2 - m2**2) / 2
  E1 = (m0**2 + m1**2 - m2**2) / (2 * m0)
  E2 = (m0**2 + m2**2 - m1**2) / (2 * m0)

  # Momentum of the daughters in the rest frame
  p = math.sqrt(E1**2 - m1**2)
  #print(p, p12)

  # Generate random angles for isotropic decay
  theta = random.uniform(0, math.pi)
  phi = random.uniform(0, 2 * math.pi)

  # Momentum components for daughter 1 in the rest frame
  p1x = p * math.sin(theta) * math.cos(phi)
  p1y = p * math.sin(theta) * math.sin(phi)
  p1z = p * math.cos(theta)

  # Create four-momentum vectors for the two daughters in the rest frame
  daughter1 = ROOT.Math.PxPyPzMVector(p1x, p1y, p1z, m1)
  daughter2 = ROOT.Math.PxPyPzMVector(-p1x, -p1y, -p1z, m2)

  # Boost daughters back to the lab frame if the parent is moving
  boost_vector = ROOT.Math.Boost(parent.BoostToCM())
  daughter1 = boost_vector(daughter1)
  daughter2 = boost_vector(daughter2)

  return daughter1, daughter2



if __name__ == '__main__':
  # Set the number of particles to generate
  nParticles = 1000
  #nParticles = 10

  # Create a histogram to store generated momenta
  #hMomentum = ROOT.TH1F("hMomentum", "Particle Momentum", 100, 0, 10)
  hInvMass = ROOT.TH1F("hInvMass", "Invariant Mass", 300, 0, 3)

  # Particle masses
  mass_pi_ch = 0.13957
  mass_k_ch = 0.493677
  mass_k_zero = 0.497611
  mass_d_zero = 1.86484
  #mass_dplus = 1.86966
  mass_bzero = 5.27972

  # seed
  random.seed(42)

  # Loop over particles and generate random kinematics
  #for i in range(nParticles):
  for i in tqdm.tqdm(range(nParticles)):
    tracks = []
    tracks += generate_and_decay(mass_k_zero, mass_pi_ch, mass_pi_ch)
    #tracks += generate_and_decay(mass_d_zero, mass_pi_ch, mass_k_ch)
    tracks += generate_and_decay(mass_d_zero, mass_pi_ch, mass_pi_ch)

    #assert len(tracks) == 2
    for itr1 in range(len(tracks)):
      for itr2 in range(itr1 + 1, len(tracks)):
        hInvMass.Fill((tracks[itr1] + tracks[itr2]).M())
  
  # Fit
  fitFunc = ROOT.TF1("fitFunc", "[0]/(sqrt(2*TMath::Pi())*[2])*TMath::Exp(-0.5*((x-[1])/[2])^2)"   # First Gaussian
                                      " + [3]/(sqrt(2*TMath::Pi())*[5])*TMath::Exp(-0.5*((x-[4])/[5])^2)" # Second Gaussian
#                                      " + [6] + [7]*x + [8]*x^2 + [9]*x^3",               # Polynomial (2nd degree)
                                      " + pol5(6)",               # Polynomial
                          -10, 10)
  fitFunc.SetParameters(1, 0.5, 0.01,  # Amplitude, mean, sigma of first Gaussian
                        1, 1.85, 0.05,  # Amplitude, mean, sigma of second Gaussian
                        0, 0, 0, 0);    # Coefficients of the polynomial
  fitFunc.SetParLimits(1, 0.48, 0.52)
  fitFunc.SetParLimits(2, 0., 0.1)
  fitFunc.SetParLimits(4, 1.8, 1.9)
  fitFunc.SetParLimits(5, 0., 0.2)
  hInvMass.Fit(fitFunc)
  bin_width = (hInvMass.GetBinLowEdge(hInvMass.GetNbinsX() + 1) - hInvMass.GetBinLowEdge(1)) / hInvMass.GetNbinsX()
  print(f'Number of signal events #1 = {fitFunc.GetParameter(0)/bin_width} +- {fitFunc.GetParError(0)/bin_width}')
  print(f'Number of signal events #2 = {fitFunc.GetParameter(3)/bin_width} +- {fitFunc.GetParError(3)/bin_width}')

  # Create a canvas to draw the histogram
  canvas = ROOT.TCanvas("canvas", "Invariant Mass", 800, 600)
  hInvMass.Draw()

  # Save the histogram as an image
  canvas.SaveAs("invmass.pdf")
  canvas.SaveAs("invmass.png")
