"""Validate public NumPy-style docstrings against Python signatures."""

import inspect
from typing import get_type_hints

from numpydoc.validate import validate

import MomentRT

# An extended summary, example, and See Also section are useful when needed,
# but do not add information to every simple property or conversion method.
OPTIONAL_SECTIONS = {"ES01", "EX01", "SA01"}


def main() -> None:
    """Check public class documentation, parameter coverage, and type resolution."""
    errors = []

    checked = 0

    for name in MomentRT.__all__:
        cls = getattr(MomentRT, name)

        if not inspect.isclass(cls) or issubclass(cls, Exception):
            continue

        objects = []

        if not getattr(cls, "_is_protocol", False):
            objects.append((f"{cls.__module__}.{cls.__name__}", cls))

        for member_name, member in vars(cls).items():
            if member_name.startswith("_") and member_name != "__init__":
                continue

            if isinstance(member, property):
                function = member.fget
            elif isinstance(member, classmethod):
                function = member.__func__
            else:
                function = member

            if (
                not inspect.isfunction(function)
                or function.__module__ != cls.__module__
            ):
                continue

            path = f"{cls.__module__}.{cls.__name__}.{member_name}"

            objects.append((path, function))

            signature = inspect.signature(function)

            for parameter in signature.parameters.values():
                if (
                    parameter.name not in {"self", "cls"}
                    and parameter.annotation is inspect.Parameter.empty
                ):
                    errors.append(f"{path}: missing annotation for {parameter.name}")

            if signature.return_annotation is inspect.Signature.empty:
                errors.append(f"{path}: missing return annotation")

            try:
                get_type_hints(function)
            except (NameError, TypeError, AttributeError) as error:
                errors.append(f"{path}: unresolved annotation: {error}")

        for path, _ in objects:
            checked += 1

            for code, message in validate(path)["errors"]:
                if code not in OPTIONAL_SECTIONS:
                    errors.append(f"{path}: {code}: {message}")

    if errors:
        raise SystemExit("\n".join(errors))

    print(
        f"Validated {checked} public class and member docstrings and their type hints."
    )


if __name__ == "__main__":
    main()
