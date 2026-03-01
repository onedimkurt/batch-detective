"""Utility functions for batch-detective."""

import socket
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def find_available_port(start_port: int) -> int:
    """Find next available TCP port starting from start_port.

    Args:
        start_port: Starting port number.

    Returns:
        First available port number.

    Raises:
        RuntimeError: If no port found in range.
    """
    from .exceptions import BatchDetectiveError
    port = start_port
    while port < start_port + 100:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(("localhost", port)) != 0:
                return port
        port += 1
    raise BatchDetectiveError("No available port found in range.")
