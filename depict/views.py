from django.shortcuts import render

from rest_framework.decorators import api_view

from .utils import *


def index(request):
    return render(request, "depict/index.html")


@api_view(["GET", "POST"])
def get_mols(request, request_type):
    format = request.POST.get("imgfmt")
    size = request.POST.get("size")
    smarts = request.POST.get("smarts")
    align_smarts = request.POST.get("alignSmarts", "off") == "on"

    size = get_image_size(size)
    input_text = None
    mol_supplier = None

    if request.method == "POST":
        if INFILE in request.FILES:
            filename = save_file(request)
            file_type = get_file_type(filename).value

            # store current file type in case of using input_text later on
            # (else statement below)
            request.session["file_type"] = file_type
            mol_supplier = get_mol_supplier(file_type, file_path=filename)
            if file_type == FileType.CSV.value or file_type == FileType.SMI.value:
                input_text, _ = get_content_from_csv(filename)
            elif file_type == FileType.MOL.value or file_type == FileType.SDF.value:
                input_text = get_content_from_file(filename)
            delete_file(filename)
        else:
            if request_type == InputType.DEMO.value:
                # in case earlier file was sdf/mol
                request.session["file_type"] = FileType.CSV.value
            file_type = request.session.get("file_type")
            if file_type != FileType.SDF.value and file_type != FileType.MOL.value:
                input_text, _ = get_content(request_type, request)
            else:
                # text is from sdf/mol file
                input_text = request.data.get(IN_TEXT).strip()
            print("input_text", input_text)
            mol_supplier = get_mol_supplier(file_type, file_data=input_text)

    output = get_svgs_from_mol_supplier(
        mol_supplier, format, size, smarts, align_smarts
    )
    context = {IMAGES: output, INPUT_TEXT: input_text}

    return render(request, "depict/index.html", context=context)


def help(request):
    return render(request, "depict/help.html")
