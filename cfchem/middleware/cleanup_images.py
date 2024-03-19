import os

from django.contrib.sessions.models import Session
from django.utils import timezone

from cfchem.Constants import IMAGE_PATHS


class CleanupImagesMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def _cleanup_session(self, session):
        # remove images if session expired
        if session.expire_date < timezone.now():
            data = session.get_decoded()
            image_paths = data.pop(IMAGE_PATHS, [])
            for f in image_paths:
                if os.path.exists(f):  # should always exist, just to be safe
                    os.remove(f)
            # session no longer needed
            session.delete()

    def __call__(self, request):
        # Code to be executed for each request before
        # the view (and later middleware) are called.
        response = self.get_response(request)

        # Code to be executed for each request/response after
        # the view is called.
        sessions = Session.objects.all()
        for session in sessions:
            self._cleanup_session(session)

        return response
