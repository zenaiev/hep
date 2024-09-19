import ROOT
import pandas as pd

# Open the ROOT file
file = ROOT.TFile("tracks.root")

# Access the TTree
tree = file["tree"]  # Replace 'myTree' with your actual TTree name

# Create lists to store TLorentzVector components
px_list, py_list, pz_list, m_list = [], [], [], []

# Loop through the events and extract Lorentz vectors
i = 0
for event in tree:
    i += 1
    if(i<=5): print(f'i = {i}')
    px, py, pz, m = [], [], [], []
    for vec in event.tracks:  # vec is a TLorentzVector
        if(i<=5): print(vec)
        px.append(vec.Px())
        py.append(vec.Py())
        pz.append(vec.Pz())
        m.append(vec.M())

    # Store the event-level data
    px_list.append(px)
    py_list.append(py)
    pz_list.append(pz)
    m_list.append(m)

# Create a DataFrame
df = pd.DataFrame({
    'px': px_list,
    'py': py_list,
    'pz': pz_list,
    'M': m_list
})

# Inspect the DataFrame
print(df.head())
