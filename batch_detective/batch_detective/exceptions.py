"""Custom exceptions for batch-detective."""


class BatchDetectiveError(Exception):
    """Base exception. All user-facing errors inherit from this."""


class BatchDetectiveValidationError(BatchDetectiveError):
    """Input file validation failures."""


class BatchDetectiveDependencyError(BatchDetectiveError):
    """Missing optional dependency."""


class BatchDetectiveDataError(BatchDetectiveError):
    """Data quality issues that prevent analysis."""
