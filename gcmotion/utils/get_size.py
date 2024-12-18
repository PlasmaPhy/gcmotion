import sys
from pint import UnitRegistry


def get_size(obj, seen=None):
    """Recursively finds size of objects.

    All of the code is gracefully borrowed from Wissam Jarjoui
    from `<https://goshippo.com/blog/measure-real-size-any-python-object>`_.

    Parameters
    ----------
    object : object
        The object to find the size of.
    seen : set
        Set that keeps track of object ids that have already been measured.
        Should be left to default since the fuction uses it when calling itself
        recursively to handle self-referential objects. Defaults to None.

    Returns
    -------
    int
        The size of the object in bytes.
    """

    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0

    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, "__dict__"):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, "__iter__") and not isinstance(
        obj, (str, bytes, bytearray)
    ):
        try:
            size += sum([get_size(i, seen) for i in obj])
        except TypeError:
            pass
    return size


def get_iter_size(obj):

    sizes = {}
    for key, value in vars(obj).items():
        size = get_size(value)
        sizes[key] = size

    # Sort by size
    sizes = dict(sorted(sizes.items(), key=lambda item: item[1], reverse=True))

    # Quantify
    ureg = UnitRegistry()
    for key, value in sizes.items():
        sizes[key] = value * ureg.bytes

    # Print total size of obj
    total_size = get_size(obj) * ureg.bytes
    print(f"obj size = {total_size:.4g~P#}\n")

    # Print results
    for key, value in sizes.items():
        print(f"{key} size = {value:.4g~P#}")
