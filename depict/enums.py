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
    xs = (800, 500)
    s = (1000, 600)
    m = (1400, 800)

class MolsRow(Enum):
    ten = 10
    eight = 8
    six = 6
    four = 4

