# ===============================================================================
# Exception class, simple error handling
# ===============================================================================
class LocalError(Exception):
    def __init__(self, e, title):
        self.msg = str(e)
        self.title = title

    def __str__(self):
        return self.msg
