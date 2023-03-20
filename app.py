from flask import Flask, render_template

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import AllChem
from IPython.display import display, SVG
from rdkit.Chem.Draw import rdMolDraw2D
import csv
from flask import request
import os

template_dir = os.path.abspath('flask/templates')
app = Flask(__name__, template_folder = template_dir)

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

@app.route("/")
def app_start():
    m = Chem.MolFromSmiles('OCCc1ccn2cnccc12')
    s = Chem.MolFromSmiles('c1nccc2n1ccc2')
    AllChem.Compute2DCoords(s)
    svg = show(s)
    # AllChem.GenerateDepictionMatching2DStructure(m, s)
    return render_template("index.html", images = [])

@app.route("/get_mols", methods = ['POST'])
def test():
    data = request.form['intxt']
    datas = data.split("\r\n")
    output = []
    for d in datas:
        smile, name = d.split(" ")
        comp = Chem.MolFromSmiles(smile)
        AllChem.Compute2DCoords(comp)
        svg = show(comp)
        output.append([svg, name])

    print(output)
    return render_template("index.html", images = output)
