class LoggedException(Exception):
    """Exception class to log error messages relevant for the user."""

    message: str

    def __init__(self, message):
        self.message = message
