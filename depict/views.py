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
    datas = get_content(type, request)

    print("datas", datas)
    output = []
    for d in datas:
        print("da", d)
        
        smile, name = d.split(" ")
        comp = Chem.MolFromSmiles(smile)
        AllChem.Compute2DCoords(comp)
        svg = show(comp)
        output.append([svg, name])
    
    context = {
        "images": output
    }

    return render(request, "depict/index.html", context = context)