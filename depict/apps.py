import logging

from django.apps import AppConfig

from cfchem.Constants import LOGGING_FNAME


class DepictConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "depict"
    def ready(self):
        # startup code, configure logging
        logging.basicConfig(filename=LOGGING_FNAME, level=logging.WARNING) 
