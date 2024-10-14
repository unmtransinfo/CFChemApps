from django.contrib import messages
from django.shortcuts import render
from rest_framework.decorators import api_view

from .enums import FileType
from .utils import *

ACCEPTED_FILE_TYPES_LIST = [f".{filetype.value}" for filetype in FileType]


def index(request):
    return render(
        request,
        "depict/index.html",
        context={
            ACCEPTABLE_FILETYPES: ACCEPTED_FILE_TYPES_LIST,
        },
    )


@api_view(["GET", "POST"])
def get_mols(request, request_type):
    input_format = request.POST.get("molfmt")  # either MolFile or SMILES
    # options for tsv/csv/smiles file
    delimiter = request.POST.get("delimiter")
    delimiter = "\t" if delimiter == "\\t" else delimiter
    smiles_col = int(request.POST.get("smiles_col"))
    names_col = int(request.POST.get("name_col"))
    has_header = request.POST.get("has_header", "off") == "on"
    # output options
    sanitize_mols = request.POST.get("sanitize_mols", "off") == "on"
    kekulize_mols = request.POST.get("kekulize_mols", "off") == "on"
    image_format = request.POST.get("imgfmt")
    size = request.POST.get("size")
    start_idx = int(request.POST.get("start_idx"))
    max_mols = int(request.POST.get("max_mols"))
    # SMARTS options
    smarts = request.POST.get("smarts")
    align_smarts = request.POST.get("alignSmarts", "off") == "on"

    size = get_image_size(size)
    input_text = None
    mol_supplier = None

    if request.method == "POST":
        if INFILE in request.FILES:
            filename = save_file(request)

            # store current file type in case of using input_text later on
            # (else statement below)
            try:
                mol_supplier = get_mol_supplier(
                    input_format,
                    file_path=filename,
                    has_header=has_header,
                    smiles_col=smiles_col,
                    names_col=names_col,
                    sanitize_mols=sanitize_mols,
                    delimiter=delimiter,
                )
                input_text = get_content_from_file(filename)
            except Exception as e:
                msg = f"Error reading file: {str(e)}"
                messages.error(request, str(e))
            if input_text and input_text[0] == "\n":
                # this covers an edge case, middleware removes "\n" otherwise
                # which leads to error in processing certain SDF/MOL files
                input_text = " " + input_text
            delete_file(filename)
        else:
            if input_format == MOLFILE:
                # text is from sdf/mol file
                input_text = request.data.get(IN_TEXT)
            else:
                input_text = get_input_text(request_type, request)
            if len(input_text) > 0:
                try:
                    mol_supplier = get_mol_supplier(
                        input_format,
                        file_data=input_text,
                        has_header=has_header,
                        smiles_col=smiles_col,
                        names_col=names_col,
                        sanitize_mols=sanitize_mols,
                        delimiter=delimiter,
                    )
                except Exception as e:
                    msg = f"Error reading data: {str(e)}"
                    messages.error(request, msg)
    output, failures = [], []
    output, failures = get_svgs_from_mol_supplier(
        mol_supplier,
        image_format,
        size,
        smarts,
        align_smarts,
        start_idx,
        max_mols,
        kekulize_mols,
    )
    context = {
        IMAGES: output,
        INPUT_TEXT: input_text,
        ACCEPTABLE_FILETYPES: ACCEPTED_FILE_TYPES_LIST,
        FAILURES: failures,
    }
    prev_image_paths = request.session.get(IMAGE_PATHS, [])
    # o[0] = image path
    request.session[IMAGE_PATHS] = prev_image_paths + [o[0] for o in output]
    return render(request, "depict/index.html", context=context)


def help(request):
    return render(request, "depict/help.html")
