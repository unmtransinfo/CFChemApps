from enum import Enum


class InputType(Enum):
    DEMO = "demo"
    INPUT = "input"
    FILE = "file"


class FileType(Enum):
    CSV = "csv"
    MOL = "mol"
    SDF = "sdf"
    SMI = "smi"
    SMILES = "smiles"
    TSV = "tsv"
    TXT = "txt"


class ImageFormat(Enum):
    PNG = "png"
    JPG = "jpeg"
    SVG = "svg"


class ImageSize(Enum):
    s = (450, 300)
    m = (600, 400)
    l = (900, 600)
    xl = (1200, 800)

class MolsRow(Enum):
    ten = 10
    eight = 8
    six = 6
    four = 4

