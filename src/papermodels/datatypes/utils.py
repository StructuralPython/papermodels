from typing import Any


def class_representation(object: Any) -> str:
    """
    Returns a generic repr string for a given 'object' where 'object' is a
    "data class"-like object.
    """
    class_name = str(type(object)).split(".")[-1].replace("'", "").replace(">", "")
    attrs = []
    for k, v in object.__dict__.items():
        if k.startswith("_"):
            continue
        if isinstance(v, str):
            v = f"'{v}'"
        attrs.append(f"{k}={v}")
    return f"{class_name}({', '.join(attrs)})"
