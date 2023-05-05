from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name = "index"),
    path("mols/<str:type>/", views.get_mols, name = "get_mols"),
    path("help/", views.help, name = "help")
]
