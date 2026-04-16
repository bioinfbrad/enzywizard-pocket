from __future__ import annotations

import argparse

from .commands.pocket import add_pocket_parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzywizard-pocket",
        description="EnzyWizard-Pocket: Detect and characterize binding pockets from a cleaned protein structure and generate a detailed JSON report."
    )
    add_pocket_parser(parser)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)