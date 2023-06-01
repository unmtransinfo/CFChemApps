from rdkit.Chem import AllChem  

from django.shortcuts import render

from rest_framework.decorators import api_view, permission_classes
from rest_framework.response import Response
from rest_framework import status

from .utils import *

def index(request):
    return render(request, "depict/index.html")

@api_view(['GET', 'POST'])
def get_mols(request, type):
    format = request.POST.get("imgfmt")
    if request.method == 'POST':
        if INFILE in request.FILES:
            filename = save_file(request)
            type = get_file_type(filename)

            if type == FileType.CSV or type == FileType.SMI:
                input_text, datas = get_content_from_csv(filename)
                delete_file(filename)
            elif type == FileType.MOL or type == FileType.SDF:
                output = get_svgs_from_mol_file(filename, format)
                input_text = ""
                context = {
                    IMAGES: output,
                    INPUT_TEXT: input_text
                }
                delete_file(filename)

                return render(request, "depict/index.html", context = context)
        else:
            input_text, datas = get_content(type, request)

    output = get_svgs_from_data(datas, format)
    context = {
        IMAGES: output,
        INPUT_TEXT: input_text
    }

    return render(request, "depict/index.html", context = context)

def help(request):
    return render(request, "depict/help.html")
