from flask import Flask, render_template

from rdkit import Chem
from rdkit.Chem import AllChem
import csv
from flask import request
import os

from flask.utils import show

template_dir = os.path.abspath('templates')
app = Flask(__name__, template_folder = template_dir)

@app.route("/test")
def app_start():
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

    return render_template("index.html", images = output)

if __name__ == "__main__":
    port = int(os.environ.get('PORT', 5000))
    app.run(debug=True, host='0.0.0.0', port=port, TEMPLATES_AUTO_RELOAD = True)