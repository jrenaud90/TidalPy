"""Promoting a configured truncation level that is not tabulated, shared by the eccentricity and obliquity drivers.

A configuration file written before the tabulated levels changed can hold a level that no longer exists. Rather than
failing, the level is promoted to the next tabulated one (or, past the highest, to the exact or general functions) with
a once-per-session warning, so the file keeps working and the accuracy never silently decreases. A level passed
directly to a function is checked strictly instead, by each driver's ``validate_*_truncation``.
"""
import warnings

from TidalPy.configurations import warning_enabled


def promote_truncation(
        level,
        tabulated: tuple,
        names: dict,
        beyond_highest: int,
        family: str,
        warned_levels: set,
        warn=None,
        name_of=str) -> int:
    """Resolve a configured truncation to a tabulated level, a named code, or the code past the highest level.

    Parameters
    ----------
    level : int or str
        The configured level, or one of ``names``.
    tabulated : tuple of int
        The tabulated levels, ascending.
    names : dict
        Lower-case names mapped to their codes (``"exact"``; ``"off"``, ``"gen"``).
    beyond_highest : int
        The code a level past the highest tabulated one is promoted to (the exact or general functions).
    family : str
        ``"Eccentricity"`` or ``"Obliquity"``, for the messages.
    warned_levels : set
        Levels already warned about this session; a promoted level is added.
    warn : bool, optional
        Whether to warn; the ``[warnings] truncation_promotion`` switch of the TidalPy configuration when None.
    name_of : callable, optional
        A code as a configuration writes it, for the warning.

    Returns
    -------
    int
        The tabulated level or code to use.

    Raises
    ------
    ValueError
        For a negative level, a bool, or an unknown name.
    """
    supported = f"Supported levels: {tabulated}, or {', '.join(repr(name) for name in names)}."
    if isinstance(level, str):
        text = level.strip().lower()
        if text in names:
            return names[text]
        try:
            level = int(text)
        except ValueError:
            raise ValueError(f"{family} truncation {level!r} is not supported. {supported}") from None
    if isinstance(level, bool):
        raise ValueError(f"{family} truncation must be a level or a name, not {level!r}. {supported}")
    level = int(level)
    if (level in tabulated) or (level == beyond_highest):
        return level
    if level < 0:
        raise ValueError(f"{family} truncation {level} is not supported. {supported}")
    promoted = next((candidate for candidate in tabulated if candidate > level), beyond_highest)
    if warn is None:
        warn = warning_enabled("truncation_promotion")
    if warn and (level not in warned_levels):
        warned_levels.add(level)
        warnings.warn(
            f"{family} truncation {level} is not tabulated; using {name_of(promoted)!r} instead. {supported}")
    return promoted
