from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, SVG

def show(mol,molSize=(475,175),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    assert mc.GetNumConformers() > 0
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    image = SVG(svg.replace('svg:', ''))
    return svg.replace('svg:', '')