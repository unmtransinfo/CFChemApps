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
    TSV = "tsv"
    TXT = "txt"

class ImageFormat(Enum):
    PNG = "png"
    JPG = "jpeg"
    SVG = "svg"

class ImageSize(Enum):
    xs = (96, 96)
    s = (160, 160)
    m = (180, 260)
    l = (280, 380)
    xl = (480, 640)