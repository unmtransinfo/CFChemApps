from django.shortcuts import render

from rest_framework.decorators import api_view

from .utils import *

def index(request):
    return render(request, "depict/index.html")

@api_view(['GET', 'POST'])
def get_mols(request, request_type):
    format = request.POST.get("imgfmt")
    size = request.POST.get("size")
    smarts = request.POST.get("smarts")
    align_smarts = request.POST.get('alignSmarts', 'off') == 'on'


    size = get_image_size(size)

    if request.method == 'POST':
        if INFILE in request.FILES:
            filename = save_file(request)
            file_type = get_file_type(filename)

            # store current file type in case of using input_text later on
            # (else statement below)
            request.session['file_type'] = file_type.value

            if file_type == FileType.CSV or file_type == FileType.SMI:
                input_text, datas = get_content_from_csv(filename)
                delete_file(filename)
            elif file_type == FileType.MOL or file_type == FileType.SDF:
                output = get_svgs_from_mol_file(filename, format, size, smarts, align_smarts)
                input_text = get_content_from_file(filename)
                context = {
                    IMAGES: output,
                    INPUT_TEXT: input_text
                }
                delete_file(filename)

                return render(request, "depict/index.html", context = context)
        else:
            if request_type == InputType.DEMO.value:
                # in case earlier file was sdf/mol
                request.session['file_type'] = FileType.CSV.value
            file_type = request.session.get('file_type')

            if file_type != FileType.SDF.value and file_type != FileType.MOL.value:
                input_text, datas = get_content(request_type, request)
            else:
                # text is from sdf/mol file
                input_text = request.data.get(IN_TEXT).strip()
                output = get_svgs_from_mol_file("", format, size, smarts, align_smarts, input_text)
                context = {
                    IMAGES: output,
                    INPUT_TEXT: input_text
                }
                return render(request, "depict/index.html", context = context)

    output = get_svgs_from_data(datas, format, size, smarts, align_smarts)
    context = {
        IMAGES: output,
        INPUT_TEXT: input_text
    }

    return render(request, "depict/index.html", context = context)

def help(request):
    return render(request, "depict/help.html")
