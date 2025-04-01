import math

from django.contrib import messages
from django.http import JsonResponse
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
    molsRow = request.POST.get("mols-row")
    start_idx = int(request.POST.get("start_idx"))
    max_mols = int(request.POST.get("max_mols"))
    # SMARTS options
    smarts = request.POST.get("smarts")
    align_smarts = request.POST.get("alignSmarts", "off") == "on"

    size = get_image_size(size)
    row = get_mols_row(molsRow)
    input_text = None
    mol_supplier = None

    # image width calculation
    screen_width = request.GET.get('screen_width', 1000)  # Default width
    selected_rows = request.GET.get('selected_rows', 6)   # Default rows
    image_width = int(screen_width) // int(selected_rows)  # Calculate width
    output, failures = [], []
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
                    failures=failures,
                )
                input_text = get_content_from_file(filename)
            except Exception as e:
                    msg = f"Error reading file: {str(e)}"
                    failures.append(msg)
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
                        failures=failures
                    )
                except Exception as e:
                    msg = f"Error reading data: {str(e)}"
                    failures.append(msg)


    output, svg_failures = get_svgs_from_mol_supplier(
            mol_supplier,
            image_format,
            size,
            smarts,
            align_smarts,
            start_idx,
            max_mols,
            kekulize_mols,
    )
    failures = failures + svg_failures
    context = {
        IMAGES: output,
        INPUT_TEXT: input_text,
        ACCEPTABLE_FILETYPES: ACCEPTED_FILE_TYPES_LIST,
        FAILURES: failures,
        'image_width': image_width,
    }

    prev_image_paths = request.session.get(IMAGE_PATHS, [])
    # o[0] = image path
    request.session[IMAGE_PATHS] = prev_image_paths + [o[0] for o in output]
    return render(request, "depict/index.html", context=context)

@api_view(["GET", "POST"])
def calculate_image_width(request):
    try:
        # Get the selected row value and screen width from the request
        selected_rows = int(request.GET.get('selected_rows', 1))
        screen_width = int(request.GET.get('screen_width', 800))  # default fallback
        print(selected_rows)
    except (ValueError, TypeError):
        return JsonResponse({'error': 'Invalid parameters'}, status=400)

    # The logic from your JS:
    # imageWidth = Math.floor(((screenSize-90)/(integerValue))-33);
    computed_width = math.floor((screen_width - 90) / selected_rows) - 33

    # If you want to apply the further subtraction as in the JS:
    final_width = computed_width - 10

    return JsonResponse({'image_width': final_width})

def help(request):
    return render(request, "depict/help.html")