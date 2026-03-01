"""Compatibility shims for optional dependencies."""

try:
    from tqdm import tqdm
except ImportError:
    # Minimal tqdm-compatible fallback
    class tqdm:
        def __init__(self, iterable=None, *args, desc=None, disable=False, **kwargs):
            self._iter = iterable
            self._disable = disable
            self._desc = desc
        def __iter__(self):
            return iter(self._iter)
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass
        def update(self, n=1):
            pass
        @classmethod
        def write(cls, s):
            pass
