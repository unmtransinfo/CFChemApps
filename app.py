from flask import Flask, render_template

from rdkit import Chem
from rdkit.Chem import AllChem
import csv
from flask import request
import os

from web.utils import show, get_content

template_dir = os.path.abspath("templates")
app = Flask(__name__, template_folder=template_dir)


@app.route("/")
def app_start():
    return render_template("index.html", images=[])

@app.route("/get_mols/<type>", methods = ['POST', "GET"])
def get_mols(type):
    datas = get_content(type, request)
    # data = request.form.get('intxt') or None
    # datas = data.split("\r\n")
    print("size", len(datas))

@app.route("/get_mols", methods=["POST"])
def test():
    data = request.form["intxt"]
    datas = data.split("\r\n")
    output = []
    for d in datas:
        print("da", d)
        smile, name = d.split(" ")
        comp = Chem.MolFromSmiles(smile)
        AllChem.Compute2DCoords(comp)
        svg = show(comp)
        output.append([svg, name])

    return render_template("index.html", images=output)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(debug=True, host="0.0.0.0", port=port, TEMPLATES_AUTO_RELOAD=True)
