from rdkit import Chem, Geometry
from rdkit.Chem import rdDepictor
from rdkit.Chem import Draw

from django.core.files.storage import FileSystemStorage

import os
import csv
import secrets

from .enums import InputType, FileType, ImageFormat, ImageSize
from cfchem.Constants import *


def get_delimiter(
    file_path=None,
    data=None,
    bytes=4096,
):
    # determine the delimiter of an input csv/txt/tsv/smi file
    # from: https://stackoverflow.com/a/69796836
    sniffer = csv.Sniffer()
    delimiters = [",", "\t", " "]  # allowed delimiters
    if data is None:
        data = open(file_path, "r").read(bytes)
    delimiter = sniffer.sniff(data, delimiters=delimiters).delimiter
    return delimiter


def get_mol_supplier(file_type: str, file_path=None, file_data=None):
    if file_path is None and file_data is None:
        # must provide at least one
        raise ValueError("No file_path or file_data provided")

    suppl = None
    if file_type == FileType.SDF.value or file_type == FileType.MOL.value:
        if file_path is not None and os.path.exists(file_path):
            suppl = Chem.SDMolSupplier(file_path, sanitize=True, removeHs=True)
        else:
            suppl = Chem.SDMolSupplier()
            suppl.SetData(file_data)
    elif file_type in [ft.value for ft in FileType]:
        # file is one of tsv, csv, smi, or txt
        delimiter = get_delimiter(file_path, data=file_data)
        if file_path is not None and os.path.exists(file_path):
            suppl = Chem.SmilesMolSupplier(
                file_path,
                delimiter=delimiter,
                smilesColumn=0,
                nameColumn=1,
                titleLine=False,
                sanitize=True,
            )
        else:
            suppl = Chem.SmilesMolSupplierFromText(
                file_data,
                delimiter=delimiter,
                smilesColumn=0,
                nameColumn=1,
                titleLine=False,
                sanitize=True,
            )
    else:
        raise ValueError("Unsupported file type: {}".format(file_type))
    return suppl


def generate_random_name():
    name = secrets.token_hex(16)
    return name


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


def ensure_directory_exists(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def create_media_filename(filename):
    ensure_directory_exists(MEDIA_FOLDER)
    return "{}/{}".format(MEDIA_FOLDER, filename)


def get_atom_bond_matches(m, smarts):
    substructure = Chem.MolFromSmarts(smarts)
    all_atom_matches = m.GetSubstructMatches(substructure)
    bond_matches = []
    for atom_matches in all_atom_matches:
        for bond in substructure.GetBonds():
            idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond_matches.append(m.GetBondBetweenAtoms(atom_matches[idx1], atom_matches[idx2]).GetIdx())
    all_atom_matches = sum(all_atom_matches, ()) # this just combines the tuple of tuples into a single tuple
    return all_atom_matches, bond_matches, substructure


def create_image(m, filename, format, size, smarts, first_match_coords=None, align_smarts: bool = False):
    all_atom_matches, bond_matches, substructure = get_atom_bond_matches(m, smarts)
    if align_smarts and (len(all_atom_matches) > 0 or len(bond_matches) > 0):
        # have match to smarts
        match = m.GetSubstructMatch(substructure)
        if first_match_coords is not None:
            # align coordinates to first match
            coord_dict={}
            for i, coord in enumerate(first_match_coords):
                coord_dict[match[i]] = coord
            rdDepictor.Compute2DCoords(m, coordMap=coord_dict)
        else:
            # is first match, so compute coordinates
            rdDepictor.Compute2DCoords(m)
            coords = [m.GetConformer().GetAtomPosition(x) for x in match]
            first_match_coords = [Geometry.Point2D(pt.x,pt.y) for pt in coords] 
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
    return img_name, first_match_coords


def get_svgs_from_mol_file(filename, format, size, smarts, align_smarts: bool, file_text: str = None):
    output = []
    counter = 0
    suppl = None
    if os.path.exists(filename):
        suppl = Chem.SDMolSupplier(filename)
    elif file_text is not None:
        suppl = Chem.SDMolSupplier()
        suppl.SetData(file_text)
    else:
        # this should be unreachable
        raise ValueError("No file or file text provided")
    first_match_coords = None
    for mol in suppl:
        if mol is None:
            continue

        name = mol.GetProp("_Name")
        image_name, first_match_coords = create_image(mol, create_media_filename(name), format, size, smarts, first_match_coords, align_smarts)
        counter += 1
        output.append([image_name, name])
    return output


def get_svgs_from_mol_supplier(mol_supplier, format, size, smarts, align_smarts: bool):
    output = []
    first_match_coords = None
    for mol in mol_supplier:
        if mol is None:
            continue

        name = mol.GetProp("_Name")
        image_name, first_match_coords = create_image(
            mol,
            create_media_filename(name),
            format,
            size,
            smarts,
            first_match_coords,
            align_smarts,
        )
        output.append([image_name, name])
    return output


def get_svgs_from_data(datas, format, size, smarts, align_smarts: bool):
    output = []
    counter = 0
    # track coordinates of first smarts match to use for subsequent matches
    # no effect if align_smarts == False
    first_match_coords = None
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
        image_name, first_match_coords = create_image(comp, create_media_filename(filename), format, size, smarts, first_match_coords, align_smarts) 
        counter += 1
        print(image_name, name)
        output.append([image_name, name])
    return output
