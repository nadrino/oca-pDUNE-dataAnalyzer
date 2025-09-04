import os
import argparse
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description='Convert .dat files in specified folder')
    parser.add_argument('-i', type=str, required=True,
                        help='Input folder path containing .dat files')
    parser.add_argument("-o", type=str, required=True, help="Output folder path")

    args = parser.parse_args()
    input_folder = Path(args.i)
    output_folder = Path(args.o)

    if not input_folder.exists():
        print(f"Error: Input folder '{input_folder}' does not exist.")
        return

    dat_files = list(input_folder.glob('*.dat'))

    if not dat_files:
        print(f"No .dat files found in '{input_folder}'")
        return

    cal_files = [f for f in dat_files if "_CAL_" in f.name]
    beam_files = [f for f in dat_files if "_BEAM_" in f.name]

    print(f"Found {len(dat_files)} .dat files in '{input_folder}':")
    print(f"Processing CAL files ({len(cal_files)}):")
    for file_path in cal_files:
        print(f"- {file_path.name}")

        if os.path.exists(os.path.join(output_folder, file_path.name.replace('.dat', '.root'))):
            print("skipping")
            continue

        subprocess.run([
            'bmRawToRootConverter',
            '-i', str(file_path),
            '--is-calib', '--skip-event',
            '-o', output_folder
        ])

    # print(f"Processing BEAM files ({len(beam_files)}):")
    # for file_path in beam_files:
    #     print(f"- {file_path.name}")




if __name__ == '__main__':
    main()
