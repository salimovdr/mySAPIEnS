import os
import json
import pandas as pd
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser(description="Collect scOpen metrics.json into a single CSV")
    parser.add_argument(
        "--base_dir",
        required=True,
        help="Path to clustering_scopen_cicero folder containing matrix_* subfolders"
    )
    parser.add_argument(
        "--output",
        required=False,
        default=None,
        help="Path to output CSV file (default: base_dir/metrics.csv)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    base_dir = args.base_dir

    if args.output is None:
        output_csv = os.path.join(base_dir, "metrics.csv")
    else:
        output_csv = args.output

    rows = []

    for name in sorted(os.listdir(base_dir)):
        folder = os.path.join(base_dir, name)
        if not os.path.isdir(folder):
            continue

        metrics_path = os.path.join(folder, "metrics.json")
        if not os.path.exists(metrics_path):
            continue

        with open(metrics_path, "r") as f:
            metrics = json.load(f)

        row = {"matrix_name": name}

        # flatten metrics dict
        def flatten(prefix, d):
            for k, v in d.items():
                if isinstance(v, dict):
                    flatten(f"{prefix}{k}_", v)
                else:
                    row[f"{prefix}{k}"] = v

        flatten("", metrics)
        rows.append(row)

    if not rows:
        print("Не найдено ни одного metrics.json в", base_dir)
        return

    df = pd.DataFrame(rows)
    df = df.sort_values("matrix_name")

    df.to_csv(output_csv, index=False)
    print(f"Сохранено {len(df)} строк в {output_csv}")


if __name__ == "__main__":
    main()

