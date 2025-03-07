import logging
import os
import secrets
import sys
from io import StringIO

from django.core.files.storage import FileSystemStorage
from rdkit import Chem, Geometry, rdBase
from rdkit.Chem import Draw, rdDepictor

from cfchem.Constants import *

from .enums import ImageFormat, ImageSize, InputType, MolsRow

rdBase.WrapLogs()  # log rdkit errors to stderr


def generate_random_name():
    name = secrets.token_hex(16)
    return name


def get_mol_supplier(
    input_format: str,
    file_path=None,
    file_data=None,
    has_header=False,
    smiles_col=0,
    names_col=1,
    sanitize_mols=True,
    delimiter="\t",
):
    if file_path is None and file_data is None:
        # must provide at least one
        raise ValueError("No file_path or file_data provided")

    suppl = None
    if input_format == MOLFILE:
        if file_path is not None and os.path.exists(file_path):
            suppl = Chem.SDMolSupplier(file_path, sanitize=sanitize_mols, removeHs=True)
        else:
            suppl = Chem.SDMolSupplier()
            suppl.SetData(file_data, sanitize=sanitize_mols, removeHs=True)
    elif input_format == SMILES_FILE:
        # file is one of tsv, csv, smi, or txt
        if file_path is not None and os.path.exists(file_path):
            suppl = Chem.SmilesMolSupplier(
                file_path,
                delimiter=delimiter,
                smilesColumn=smiles_col,
                nameColumn=names_col,
                titleLine=has_header,
                sanitize=sanitize_mols,
            )
        else:
            suppl = Chem.SmilesMolSupplierFromText(
                file_data,
                delimiter=delimiter,
                smilesColumn=smiles_col,
                nameColumn=names_col,
                titleLine=has_header,
                sanitize=sanitize_mols,
            )
    else:
        raise ValueError("Unsupported file type: {}".format(input_format))
    return suppl


def get_input_text(type, request):
    input_text = None
    if type == InputType.INPUT.value:
        input_text = request.data.get(IN_TEXT)
    elif type == InputType.DEMO.value:
        f = open(DEMO_COMPOUNDS_FILE_NAME)
        input_text = f.read()

    return input_text


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
    if os.path.exists(file):
        os.remove(file)


def get_image_size(size):
    image_sizes = {
        ImageSize.xs.name: ImageSize.xs.value,
        ImageSize.s.name: ImageSize.s.value,
        ImageSize.m.name: ImageSize.m.value,
    }
    return image_sizes[size]

def get_mols_row(row):
    mols_row = {
        MolsRow.ten.name: MolsRow.ten.value,
        MolsRow.eight.name: MolsRow.eight.value,
        MolsRow.six.name: MolsRow.six.value,
        MolsRow.four.name: MolsRow.four.value,
    }
    return mols_row[row]


def ensure_directory_exists(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


def create_media_filename(filename):
    ensure_directory_exists(MEDIA_FOLDER)
    return "{}/{}".format(MEDIA_FOLDER, filename)


def get_atom_bond_matches(m, smarts):
    substructure = Chem.MolFromSmarts(smarts)
    if substructure is None:
        # bad SMARTS given
        raise ValueError("Invalid SMARTS pattern: {}".format(smarts))
    all_atom_matches = m.GetSubstructMatches(substructure)
    bond_matches = []
    for atom_matches in all_atom_matches:
        for bond in substructure.GetBonds():
            idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond_matches.append(
                m.GetBondBetweenAtoms(atom_matches[idx1], atom_matches[idx2]).GetIdx()
            )
    all_atom_matches = sum(
        all_atom_matches, ()
    )  # this just combines the tuple of tuples into a single tuple
    return all_atom_matches, bond_matches, substructure


def create_image(
    m,
    filename,
    format,
    size,
    smarts,
    first_match_coords=None,
    align_smarts: bool = False,
    kekulize: bool = True,
):
    all_atom_matches, bond_matches, substructure = get_atom_bond_matches(m, smarts)
    if align_smarts and (len(all_atom_matches) > 0 or len(bond_matches) > 0):
        # have match to smarts
        match = m.GetSubstructMatch(substructure)
        if first_match_coords is not None:
            # align coordinates to first match
            coord_dict = {}
            for i, coord in enumerate(first_match_coords):
                coord_dict[match[i]] = coord
            rdDepictor.Compute2DCoords(m, coordMap=coord_dict)
        else:
            # is first match, so compute coordinates
            rdDepictor.Compute2DCoords(m)
            coords = [m.GetConformer().GetAtomPosition(x) for x in match]
            first_match_coords = [Geometry.Point2D(pt.x, pt.y) for pt in coords]
    img_name = "{}.{}".format(filename, format)
    if format in [ImageFormat.JPG.value, ImageFormat.PNG.value]:
        pil_image = Draw.MolToImage(
            m,
            size=size,
            highlightAtoms=all_atom_matches,
            highlightBonds=bond_matches,
            kekulize=kekulize,
        )
        pil_image.save(img_name)
    elif format == ImageFormat.SVG.value:
        m = Draw.rdMolDraw2D.PrepareMolForDrawing(m, kekulize=kekulize)
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        dopts = drawer.drawOptions()
        dopts.prepareMolsBeforeDrawing = kekulize
        drawer.DrawMolecule(
            m,
            highlightAtoms=all_atom_matches,
            highlightBonds=bond_matches,
        )
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        with open(img_name, "w") as f:
            f.write(svg)
    else:
        raise ValueError("Unsupported image format: {}".format(format))
    return img_name, first_match_coords


def validate_smarts(smarts: str, failures: list, sio: StringIO) -> str:
    """
    Ensure that the given smarts string is valid.
    If not, return "" and record error.
    :param str smarts: given SMARTS pattern
    :param list failures: list of errors
    :param StringIO sio: error logger
    :return str: valid SMARTS pattern ("" if smarts not valid)
    """
    substructure = Chem.MolFromSmarts(smarts)
    if substructure is None:
        msg = sio.getvalue().strip()
        failures.append(msg)
        logging.warning(msg)
        sio = sys.stderr = StringIO()  # reset the error logger
        return ""
    return smarts


def get_svgs_from_mol_supplier(
    mol_supplier,
    format,
    size,
    row,
    smarts,
    align_smarts: bool,
    start_idx: int = 0,
    max_mols=None,
    kekulize_mols: bool = True,
):
    output = []
    first_match_coords = None
    failures = []  # mols which rdkit could not interpret
    sio = sys.stderr = StringIO()  # redirect error messages to string
    if mol_supplier is None:
        return output, failures
    if max_mols is None:
        max_mols = len(mol_supplier)
    # avoid OOB error
    end_idx = min(len(mol_supplier), start_idx + max_mols)
    # ensure smarts valid, if not record error and set smarts to ""
    smarts = validate_smarts(smarts, failures, sio)
    for i in range(start_idx, end_idx):
        mol = mol_supplier[i]
        if mol is None:
            msg = sio.getvalue().strip()
            failures.append(msg)
            logging.warning(msg)
            sio = sys.stderr = StringIO()  # reset the error logger
            continue
        name = NO_COMPOUND_NAME
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        # edge case for sdf files where whitespace used in place of name
        n_whitespace = sum(1 for char in name if char.isspace())
        name = NO_COMPOUND_NAME if (len(name) == n_whitespace) else name
        fname = generate_random_name()
        image_name, first_match_coords = create_image(
            mol,
            create_media_filename(fname),
            format,
            size,
            smarts,
            first_match_coords,
            align_smarts,
            kekulize_mols,
        )
        output.append([image_name, fname, name])
    return output, failures