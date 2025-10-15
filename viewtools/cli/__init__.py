"""CLI module initialization."""

import click
import functools
import sys
from .. import __version__
import logging
from .._logging import get_logger


CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}


@click.version_option(version=__version__)
@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--post-mortem", help="Post mortem debugging", is_flag=True, default=False
)
@click.option(
    "--output-profile",
    help="Profile performance with Python cProfile and dump the statistics "
    "into a binary file",
    type=str,
    default="",
)
@click.option("-v", "--verbose", help="Verbose logging.", count=True)
@click.option(
    "-d",
    "--debug",
    help="On error, drop into the post-mortem debugger shell.",
    is_flag=True,
    default=False,
)
def cli(post_mortem, output_profile, verbose, debug):
    """Flexible tools for Hi-C data processing.

    All pairtools have a few common options, which should be typed _before_
    the command name.

    """
    if post_mortem:
        import traceback

        try:
            import ipdb as pdb
        except ImportError:
            import pdb

        def _excepthook(exc_type, value, tb):
            traceback.print_exception(exc_type, value, tb)
            print()
            pdb.pm()

        sys.excepthook = _excepthook

    if output_profile:
        import cProfile
        import atexit

        pr = cProfile.Profile()
        pr.enable()

        def _atexit_profile_hook():
            pr.disable()
            pr.dump_stats(output_profile)

        atexit.register(_atexit_profile_hook)

    # Initialize logging to stderr
    logging.basicConfig(stream=sys.stderr)
    logging.captureWarnings(True)
    root_logger = get_logger()

    # Set verbosity level
    if verbose > 0:
        root_logger.setLevel(logging.DEBUG)
        if verbose > 1:  # pragma: no cover
            try:
                import psutil
                import atexit

                @atexit.register
                def process_dump_at_exit():
                    process_attrs = [
                        "cmdline",
                        # 'connections',
                        "cpu_affinity",
                        "cpu_num",
                        "cpu_percent",
                        "cpu_times",
                        "create_time",
                        "cwd",
                        # 'environ',
                        "exe",
                        # 'gids',
                        "io_counters",
                        "ionice",
                        "memory_full_info",
                        # 'memory_info',
                        # 'memory_maps',
                        "memory_percent",
                        "name",
                        "nice",
                        "num_ctx_switches",
                        "num_fds",
                        "num_threads",
                        "open_files",
                        "pid",
                        "ppid",
                        "status",
                        "terminal",
                        "threads",
                        # 'uids',
                        "username",
                    ]
                    p = psutil.Process()
                    info_ = p.as_dict(process_attrs, ad_value="")
                    for key in process_attrs:
                        root_logger.debug("PSINFO:'{}': {}".format(key, info_[key]))

            except ImportError:
                root_logger.warning("Install psutil to see process information.")

    else:
        root_logger.setLevel(logging.INFO)

    # Set hook for postmortem debugging
    if debug:  # pragma: no cover
        import traceback

        try:
            import ipdb as pdb
        except ImportError:
            import pdb

        def _excepthook(exc_type, value, tb):
            traceback.print_exception(exc_type, value, tb)
            print()
            pdb.pm()

        sys.excepthook = _excepthook


def common_io_options(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


from . import rearrange_genome

# Register commands
cli.add_command(rearrange_genome.cli, name="rearrange-genome")
