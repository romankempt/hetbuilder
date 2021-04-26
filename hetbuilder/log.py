import logging
import pretty_errors
from rich.logging import RichHandler


pretty_errors.configure(
    separator_character="*",
    filename_display=pretty_errors.FILENAME_EXTENDED,
    line_number_first=True,
    display_link=True,
    lines_before=2,
    lines_after=1,
    line_color=pretty_errors.RED + "> " + pretty_errors.default_config.line_color,
    code_color="  " + pretty_errors.default_config.line_color,
    truncate_code=True,
    display_locals=True,
)
pretty_errors.blacklist("c:/python")


class DuplicateFilter(logging.Filter):
    def filter(self, record):
        # add other fields if you need more granular comparison, depends on your app
        current_log = (record.module, record.levelno, record.msg)
        if current_log != getattr(self, "last_log", None):
            self.last_log = current_log
            return True
        return False


def setup_custom_logger(name):
    formatter = logging.Formatter(fmt="{message:s}", style="{")
    handler = RichHandler(
        show_time=False, markup=True, rich_tracebacks=True, show_path=False
    )
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    logger.addFilter(DuplicateFilter())
    return logger


logger = setup_custom_logger("root")


def set_verbosity_level(verbosity):
    logger = logging.getLogger("root")
    if verbosity == 1:
        level = "WARNING"
    elif verbosity == 2:
        level = "INFO"
    else:
        level = "DEBUG"
        formatter = logging.Formatter(
            fmt="{levelname:8s} {module:20s} {funcName:20s} |\n {message:s}", style="{"
        )
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(level)

