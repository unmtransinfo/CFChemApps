from django.shortcuts import render
from django.contrib import messages

from rest_framework.decorators import api_view

from .utils import *

accepted_filetypes = [f".{filetype.value}" for filetype in FileType]


def index(request):
    return render(
        request,
        "depict/index.html",
        context={
            ACCEPTABLE_FILETYPES: accepted_filetypes,
        },
    )


@api_view(["GET", "POST"])
def get_mols(request, request_type):
    format = request.POST.get("imgfmt")
    size = request.POST.get("size")
    smarts = request.POST.get("smarts")
    align_smarts = request.POST.get("alignSmarts", "off") == "on"
    # options for tsv/csv/smiles file
    smiles_col = int(request.POST.get("smiles_col"))
    names_col = int(request.POST.get("name_col"))
    has_header = request.POST.get("has_header", "off") == "on"
    sanitize_mols = request.POST.get("sanitize_mols", "off") == "on"

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
            try:
                mol_supplier = get_mol_supplier(
                    file_type,
                    file_path=filename,
                    has_header=has_header,
                    smiles_col=smiles_col,
                    names_col=names_col,
                    sanitize_mols=sanitize_mols,
                )
                input_text = get_content_from_file(filename)
            except Exception as e:
                msg = f"Error reading file: {str(e)}"
                messages.error(request, str(e))
            if input_text[0] == "\n":
                # this covers an edge case, middleware removes "\n" otherwise
                # which leads to error in processing certain SDF/MOL files
                input_text = " " + input_text
            delete_file(filename)
        else:
            if request_type == InputType.DEMO.value:
                # in case earlier file was sdf/mol
                request.session["file_type"] = FileType.CSV.value
            file_type = request.session.get("file_type")
            if file_type == FileType.MOL.value or file_type == FileType.SDF.value:
                # text is from sdf/mol file
                input_text = request.data.get(IN_TEXT)
            else:
                input_text, _ = get_content(request_type, request)
            if len(input_text) > 0:
                try:
                    mol_supplier = get_mol_supplier(
                        file_type,
                        file_data=input_text,
                        has_header=has_header,
                        smiles_col=smiles_col,
                        names_col=names_col,
                        sanitize_mols=sanitize_mols,
                    )
                except Exception as e:
                    msg = f"Error reading data: {str(e)}"
                    messages.error(request, msg)
    output, failures = get_svgs_from_mol_supplier(
        mol_supplier, format, size, smarts, align_smarts
    )
    context = {
        IMAGES: output,
        INPUT_TEXT: input_text,
        ACCEPTABLE_FILETYPES: accepted_filetypes,
        FAILURES: failures,
    }

    return render(request, "depict/index.html", context=context)


def help(request):
    return render(request, "depict/help.html")
