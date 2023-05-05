from rdkit.Chem import AllChem

from django.shortcuts import render

from django.core.files.storage import FileSystemStorage

from rest_framework.decorators import api_view, permission_classes
from rest_framework.response import Response
from rest_framework import status

from .utils import *

def index(request):
    return render(request, "depict/index.html")

@api_view(['GET', 'POST'])
def get_mols(request, type):
    if request.method == 'POST' and request.FILES['infile']:
        myfile = request.FILES['infile']
        fs = FileSystemStorage()
        filename = fs.save(myfile.name, myfile)
        uploaded_file_url = fs.url(filename)

        type = get_file_type(filename)
        if type == 1:
            input_text, datas = get_content_from_csv(filename)
            delete_csv("media/{}".format(filename))
        elif type == 0:
            m = Chem.MolFromMolFile("media/{}".format(filename))
            AllChem.Compute2DCoords(m)
            svg = show(m)
            print(m)
            output = [[[svg, "---"]]]
            input_text = ""
            context = {
                "images": output,
                "input_text": input_text
            }
            print(output)

            return render(request, "depict/index.html", context = context)
    else:
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
    print(output)
    context = {
        "images": output,
        "input_text": input_text
    }

    return render(request, "depict/index.html", context = context)
