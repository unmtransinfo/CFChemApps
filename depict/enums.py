from enum import Enum
 
class InputType(Enum):
    DEMO = "demo"
    INPUT = "input"
    FILE = "file"
    # WINTER = 4

class FileType(Enum):
    CSV = "csv"
    MOL = "mol"
    SDF = "sdf"
    SMI = "smi"

class ImageFormat(Enum):
    PNG = "png"
    JPG = "jpeg"

class ImageSize(Enum):
    xs = (96, 96)
    s = (160, 160)
    m = (180, 260)
    l = (280, 380)
    xl = (480, 640)