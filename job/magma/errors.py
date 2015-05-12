class FileFormatError(Exception):
    """Raised when an error is found in an input file"""
    pass


class DataProcessingError(Exception):
    """Raised when an error occurs during processing of the data"""
    pass
