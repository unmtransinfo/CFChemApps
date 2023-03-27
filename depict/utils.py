from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, SVG

from .enums import InputType

def show(mol, molSize=(475, 175), kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    assert mc.GetNumConformers() > 0
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    image = SVG(svg.replace("svg:", ""))
    return svg.replace("svg:", "")

def get_content(type, request):
    content = None
    datas = None

    if type == InputType.INPUT.value:
        print(request)
        content = request.data.get('intxt')
        datas = content.split("\r\n")

    if type == InputType.DEMO.value:
        f = open('web/demo_compounds.txt')
        content = f.read()
        datas = content.split("\n")

    
    print(content, datas)

    return datas