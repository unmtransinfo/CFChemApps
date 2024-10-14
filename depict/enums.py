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
    xs = (120, 80)
    s = (150, 100)
    m = (300, 200)
    l = (450, 300)
    xl = (600, 400)
