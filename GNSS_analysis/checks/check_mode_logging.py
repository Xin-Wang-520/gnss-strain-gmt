"""Verify modal rounding against local GMT and live/file log agreement."""
from contextlib import redirect_stdout
import io
from pathlib import Path
import subprocess
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from GNSS_analysis.gnss_processing import (
    locate_gmt, gmt_environment, requested_mode_count, run_gpsgridder, report_used_modes)


def main():
    env = gmt_environment()
    with tempfile.TemporaryDirectory(prefix='test_mode_logging_') as tmp:
        folder = Path(tmp)
        source = folder / 'toy.txt'
        source.write_text('80 30 1 3 1 1\n81 30 2 1 1 1\n80 31 -1 4 1 1\n')
        for percent, expected, verbose in ((25, 2, False), (75, 4, True)):
            # 6*0.25=1.5 rounds to 2; 6*0.75=4.5 rounds to 4.
            # Three points plus -L avoid removal of an exactly fitting plane.
            command = [locate_gmt(None), 'gpsgridder', str(source), '-R79/82/29/32',
                       '-I1', '-fg', '-W', '-r', '-L', f'-Cn{percent}%',
                       f'-G{folder / "toy.nc"}', '-V']
            terminal = io.StringIO()
            log = folder / 'gpsgridder.log'
            with redirect_stdout(terminal):
                run_gpsgridder(command, env, log, verbose=verbose)
            assert terminal.getvalue() == (log.read_text() if verbose else '')
            report = report_used_modes(folder, 3, percent)
            assert report['actual_modes_used'] == expected
            assert report['requested_mode_count'] == expected
        terminal = io.StringIO()
        with redirect_stdout(terminal):
            try:
                run_gpsgridder([sys.executable, '-c', 'print("failed run"); raise SystemExit(7)'],
                              env, folder / 'failed.log', verbose=True)
            except subprocess.CalledProcessError as error:
                assert error.returncode == 7
            else:
                raise AssertionError('Failed command did not propagate its exit code')
        assert terminal.getvalue() == (folder / 'failed.log').read_text()
        assert 'failed run' in terminal.getvalue()
    assert [requested_mode_count(1685, p) for p in (100, 75, 50, 25)] == [3370, 2528, 1685, 842]
    print('PASS: GMT rounding, actual log counts, optional live output, saved log and failure propagation.')


if __name__ == '__main__':
    main()
