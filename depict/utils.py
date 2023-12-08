from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, SVG
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from django.core.files.storage import FileSystemStorage

import os
from PIL import Image
import secrets

from .enums import InputType, FileType, ImageFormat, ImageSize
from cfchem.Constants import *

def generate_random_name():
    name = secrets.token_hex(16)
    return name

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

def get_content_from_file(filename):
    f = open(filename)
    text = f.read()

    return text

def get_content_from_smi(filename):
    f = open(filename)
    text = f.read()
    datas = text.split("\n")

def save_file(request):
    myfile = request.FILES[INFILE]

    ensure_directory_exists(MEDIA_FOLDER)

    fs = FileSystemStorage(MEDIA_FOLDER)
    filename = fs.save(myfile.name, myfile)
    _ = fs.url(filename)
    filename = "{}/{}".format(MEDIA_FOLDER, myfile.name)

    return filename
    
def delete_file(file):
    os.remove(file)

def get_file_type(filename):
    extension = filename.split(".")[-1]
    file_types = {
        FileType.CSV.value: FileType.CSV,
        FileType.MOL.value: FileType.MOL,
        FileType.SDF.value: FileType.SDF,
        FileType.SMI.value: FileType.SMI
    }
    

    return file_types[extension]

def get_image_size(size):
    
    image_sizes = {
        ImageSize.xs.name: ImageSize.xs.value,
        ImageSize.s.name: ImageSize.s.value,
        ImageSize.m.name: ImageSize.m.value,
        ImageSize.l.name: ImageSize.l.value,
        ImageSize.xl.name: ImageSize.xl.value 
    }

    return image_sizes[size]


def get_images_per_row(size) -> int:
    """
    Return the number of images to display per row based on the size of each image.
    :param Tuple[int,int] size: size of each image
    :return int: images to display per row on web page
    """
    # TODO: try to account for screen size/width, will be useful for more devices
    row_dict = {
        ImageSize.xs.value : 8,
        ImageSize.s.value : 5,
        ImageSize.m.value: 5,
        ImageSize.l.value: 3,
        ImageSize.xl.value: 2
    }
    return row_dict[size]


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

def create_image(m, filename, format, size, smarts):
    substructure = Chem.MolFromSmarts(smarts)
    all_atom_matches = m.GetSubstructMatches(substructure)
    bond_matches = []
    for atom_matches in all_atom_matches:
        for bond in substructure.GetBonds():
            idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond_matches.append(m.GetBondBetweenAtoms(atom_matches[idx1], atom_matches[idx2]).GetIdx())
    all_atom_matches = sum(all_atom_matches, ()) # this just combines the tuple of tuples into a single tuple
    img_name = "{}.{}".format(filename, format)
    if format in [ImageFormat.JPG.value, ImageFormat.PNG.value]:
        pil_image = Draw.MolToImage(m, size=size, highlightAtoms=all_atom_matches, highlightBonds=bond_matches)
        pil_image.save(img_name)
    elif format == ImageFormat.SVG.value:
        m = Draw.rdMolDraw2D.PrepareMolForDrawing(m)
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        drawer.DrawMolecule(m, highlightAtoms=all_atom_matches, highlightBonds=bond_matches)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        with open(img_name, 'w') as f:
            f.write(svg)
    else:
        raise ValueError("Unsupported image format: {}".format(format))
    return img_name


def get_svgs_from_mol_file(filename, format, size, smarts):
    output = []
    counter = 0
    suppl = Chem.SDMolSupplier(filename)
    for mol in suppl:
        if mol is None:
            continue

        name = mol.GetProp("_Name")
        image_name = create_image(mol, create_media_filename(name), format, size, smarts)
        counter += 1
        output.append([image_name, name])
    return output

def get_svgs_from_data(datas, format, size, smarts):
    output = []
    counter = 0
    for d in datas:     
        d = d.strip() 
        if not d:
            continue

        d_split = d.split(maxsplit = 1)
        print(d_split)
        try:
            if len(d_split) == 1:
                smile, name = d_split[0], NO_COMPOUND_NAME
            else:
                smile, name = d_split
            name = name.split("\t")[0]
        except Exception as e:
            print(e, d)
            smile, name = d.split("\t") #, NO_COMPOUND_NAME
        comp = Chem.MolFromSmiles(smile)

        filename = name if name != NO_COMPOUND_NAME else generate_random_name()
        image_name = create_image(comp, create_media_filename(filename), format, size, smarts) 
        counter += 1
        print(image_name, name)
        output.append([image_name, name])
    return output