import subprocess
import sys
from pathlib import Path

import pandas as pd

from reprint_tool import __version__

DATA = Path(__file__).parent / "data" / "COSMIC_v2_SBS_GRCh37.txt"


def run_cli(*args):
    return subprocess.run(
        [sys.executable, "-m", "reprint_tool.cli", *args],
        capture_output=True,
        check=False,
        text=True,
    )


def test_version():
    result = run_cli("--version")
    assert result.returncode == 0
    assert __version__ in result.stdout


def test_requires_input():
    result = run_cli()
    assert result.returncode != 0
    assert "--input" in result.stderr


def test_save_reprint_and_pdf(tmp_path):
    out_dir = tmp_path / "out"
    reprint_path = tmp_path / "reprint.tsv"
    result = run_cli(
        "--input", str(DATA),
        "--output_dir", str(out_dir),
        "--save_reprint", str(reprint_path),
    )
    assert result.returncode == 0, result.stderr
    assert reprint_path.exists()
    assert (out_dir / "reprint_all_reprints.pdf").exists()

    df_in = pd.read_csv(DATA, sep="\t", index_col=0)
    df_out = pd.read_csv(reprint_path, sep="\t", index_col=0)
    assert list(df_out.columns) == [f"{c}_reprint" for c in df_in.columns]


def test_export_png_creates_missing_output_dir(tmp_path):
    out_dir = tmp_path / "does" / "not" / "exist"
    result = run_cli("--input", str(DATA), "--output_dir", str(out_dir), "--export_png")
    assert result.returncode == 0, result.stderr
    pngs = list(out_dir.glob("*_seamless.png"))
    assert len(pngs) == pd.read_csv(DATA, sep="\t", index_col=0).shape[1]
