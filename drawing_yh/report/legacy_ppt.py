"""转发到 report_toolkit.legacy_ppt(引擎唯一源)."""
from report_toolkit import legacy_ppt as _src  # noqa: E402

globals().update({k: v for k, v in vars(_src).items() if not k.startswith("__")})
