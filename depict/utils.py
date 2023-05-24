from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, SVG
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from django.core.files.storage import FileSystemStorage

import os
from PIL import Image

from .enums import InputType, FileType, ImageFormat
from cfchem.Constants import *

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
        input_text = request.data.get(IN_TEXT).strip()
        datas = input_text.split("\r\n")

    if type == InputType.DEMO.value:
        f = open(DEMO_COMPOUNDS_FILE_NAME)
        input_text = f.read()
        datas = input_text.split("\n")

    return input_text, datas

def get_content_from_csv(filename):
    f = open(filename)
    csv_text = f.read()
    datas = csv_text.split("\n")
    return csv_text, datas

def save_file(request):
    myfile = request.FILES[INFILE]

    ensure_directory_exists(MEDIA_FOLDER)

    fs = FileSystemStorage(MEDIA_FOLDER)
    filename = fs.save(myfile.name, myfile)
    _ = fs.url(filename)
    filename = "{}/{}".format(MEDIA_FOLDER, myfile.name)

    return filename
    
def delete_csv(file):
    os.remove(file)

def get_file_type(filename):
    extension = filename.split(".")[1]
    return FileType.CSV if extension == FileType.CSV.value else FileType.MOL

def ensure_directory_exists(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def create_media_filename(filename):
    ensure_directory_exists(MEDIA_FOLDER)
    return "{}/{}".format(MEDIA_FOLDER, filename)

def create_svg(m):
    AllChem.Compute2DCoords(m)
    svg = show(m)

    return svg

def create_png_jpeg_image(m, filename, format):
    extension = ImageFormat.PNG.value if format == ImageFormat.PNG.value else ImageFormat.JPG.value
    # img = Draw.MolToFile(m, filename + ".png")
    if extension == ImageFormat.PNG.value:
        img = Draw.MolToFile(m, filename + ".{}".format(extension))
    else:
        img = Image.open(filename + ".{}".format(extension))
        img.convert("RGB").save(filename + ".{}".format(format))

    name = filename + ".{}".format(format)

    return name


def get_svgs_from_mol_file(filename, format):
    output = []
    row_output = []
    counter = 0
    suppl = Chem.SDMolSupplier(filename)

    if len(suppl) == 1:
        mol = Chem.MolFromMolFile(filename)
        name = create_png_jpeg_image(mol, filename, format)

        output.append([[name, NO_COMPOUND_NAME]])

        return output
    
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("PUBCHEM_COMPOUND_CID") if mol.HasProp("PUBCHEM_COMPOUND_CID") else mol.GetProp("PubChem CID")
        image_name = create_png_jpeg_image(mol, create_media_filename(name), format)
        counter += 1
        row_output.append([image_name, name])
        if counter == 3:
            output.append(row_output)
            row_output = []
            counter = 0
    output.append(row_output)
    return output

def get_svgs_from_data(datas, format):
    output = []
    row_output = []
    counter = 0
    for d in datas:     
        d = d.strip() 
        try:
            smile, name = d.split(" ")
        except Exception as e:
            print(e)
            smile, name = d.rstrip(), NO_COMPOUND_NAME
        comp = Chem.MolFromSmiles(smile)
        image_name = create_png_jpeg_image(comp, create_media_filename(name), format)
        counter += 1
        row_output.append([image_name, name])
        if counter == 3:
            output.append(row_output)
            row_output = []
            counter = 0
    output.append(row_output)

    return output