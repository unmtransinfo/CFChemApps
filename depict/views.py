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
    input_text, datas = get_content(type, request)
    output = []
    row_output = []
    counter = 0
    for d in datas:      
        try:
            smile, name = d.split(" ")
        except Exception as e:
            print(e)
            smile, name = d.rstrip(), "---"
        comp = Chem.MolFromSmiles(smile)
        AllChem.Compute2DCoords(comp)
        svg = show(comp)
        counter += 1
        row_output.append([svg, name])
        if counter == 3:
            output.append(row_output)
            row_output = []
            counter = 0
    output.append(row_output)
    
    context = {
        "images": output,
        "input_text": input_text
    }

    return render(request, "depict/index.html", context = context)
