"""
Basic checks for arg parsing and command dispatch.
"""

from contextlib import redirect_stdout
from io import StringIO

from arcasHLA import main


class TestMain:
    def test_basic(self):
        """
        Basic check, do we see help printed.
        """

        with StringIO() as err_con:
            with redirect_stdout(err_con):
                try:
                    main(["--help"])

                except SystemExit as e:
                    assert e.code == 0

                    err_con.seek(0)
                    assert "usage:" in err_con.readline()
