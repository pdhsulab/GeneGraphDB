"""Util methods for IPython notebook."""

from IPython.display import HTML, display


def disp_notebook_full_width():
    """Display notebook at full screen width."""
    # https://stackoverflow.com/a/34058270/3893740
    display(HTML("<style>.container { width:100% !important; }</style>"))
