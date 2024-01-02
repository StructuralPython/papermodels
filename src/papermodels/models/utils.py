import csv


def str_to_int(s: str) -> int | str:
    """
    Returns the string 's' converted into an integer if possible.
    Returns the original string, 's', otherwise.
    """
    try:
        return int(s)
    except ValueError:
        return s


def str_to_float(s: str) -> float | str:
    """
    Returns the string 's' converted into a float if possible.
    Returns the original string, 's', otherwise.
    """
    try:
        return float(s)
    except ValueError:
        return s


def read_csv_file(filename: str) -> str:
    """
    Returns data contained in the file, 'filename' as a list of lists.
    It is assumed that the data in the file is "csv-ish", meaning with
    comma-separated values.
    """
    with open(filename, "r") as file:
        csv_data = list(csv.reader(file))
    return csv_data
