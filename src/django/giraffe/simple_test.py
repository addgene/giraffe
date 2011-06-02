
def setup():
    """Load the generic environment settings to make tests work"""
    from django.core.management import setup_environ
    import giraffe.settings
    setup_environ(giraffe.settings)
