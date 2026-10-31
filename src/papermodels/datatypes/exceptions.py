"""Exceptions raised while parsing markups and building geometry."""


class LegendError(Exception):
    """Raised when legend entries are missing, malformed or ambiguous."""

    pass


class GeometryError(Exception):
    """
    Raised when element geometry is invalid or cannot be resolved into a load path.
    """

    pass


class AnnotationError(Exception):
    """Raised when an annotation cannot be parsed into a structural element."""

    pass
