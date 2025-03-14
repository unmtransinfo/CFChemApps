from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name = "index"),
    path("mols/<str:request_type>/", views.get_mols, name = "get_mols"),
    path("mols/<str:request_type>/calculate_image_width/", views.calculate_image_width, name="calculate_image_width"),
    path("help/", views.help, name = "help"),
]
