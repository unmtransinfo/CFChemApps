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

class ImageFormat(Enum):
    PNG = "png"
    JPG = "jpg"