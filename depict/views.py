from rdkit.Chem import AllChem

from django.shortcuts import render

from rest_framework.decorators import api_view, permission_classes
from rest_framework.response import Response
from rest_framework import status

from .utils import *

def index(request):
    return render(request, "depict/test.html")

@api_view(['GET', 'POST'])
def get_mols(request, type):
    input_text, datas = get_content(type, request)

    output = []
    for d in datas:      
        smile, name = d.split(" ")
        comp = Chem.MolFromSmiles(smile)
        AllChem.Compute2DCoords(comp)
        svg = show(comp)
        output.append([svg, name])
    
    context = {
        "images": output,
        "input_text": input_text
    }

    return render(request, "depict/test.html", context = context)