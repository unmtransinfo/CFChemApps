from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, SVG

import os

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
    input_text = None
    datas = None

    if type == InputType.INPUT.value:
        input_text = request.data.get('intxt')
        datas = input_text.split("\r\n")

    if type == InputType.DEMO.value:
        f = open('depict/demo_compounds.txt')
        input_text = f.read()
        datas = input_text.split("\n")
    
    return input_text, datas

def get_content_from_csv(filename):
    f = open("media/{}".format(filename))
    csv_text = f.read()
    datas = csv_text.split("\n")
    print(filename, csv_text, datas)
    return csv_text, datas

def delete_csv(file):
    os.remove(file)

def get_file_type(filename):
    extension = filename.split(".")[1]
    return 1 if extension == "csv" else 0