import os
import argparse
from pathlib import Path
import subprocess

# python \
#   /nfs/sw/beam-plug-monitor/repo/oca-pDUNE-dataAnalyzer/python/convertRawFiles.py
#   -i /data3/np02-beam-monitor-data
#   -o /nfs/sw/beam-plug-monitor/data
#   --override

def main():
    parser = argparse.ArgumentParser(description='Convert .dat files in specified folder')
    parser.add_argument('-i', type=str, required=True,
                        help='Input folder path containing .dat files')
    parser.add_argument("-o", type=str, required=True, help="Output folder path")
    parser.add_argument("--override", action="store_true", help="Overwrite BEAM files")

    args = parser.parse_args()
    input_folder = Path(args.i)
    output_folder = Path(args.o)
    overwrite_beam_files = bool(args.override)

    if not input_folder.exists():
        print(f"Error: Input folder '{input_folder}' does not exist.")
        return

    dat_files = list(input_folder.glob('*.dat'))

    if not dat_files:
        print(f"No .dat files found in '{input_folder}'")
        return

    cal_files = sorted([f for f in dat_files if "_CAL_" in f.name], key=lambda x: x.name)
    beam_files = sorted([f for f in dat_files if "_BEAM_" in f.name], key=lambda x: x.name)

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


    def getRunNumber(filePath_):
        # SCD_RUN00283_CAL_20250812_083520
        return int(filePath_.name[7:12])

    print(f"Processing BEAM files ({len(beam_files)}):")
    for file_path in beam_files:
        print(f"- {file_path.name}")

        if not overwrite_beam_files and os.path.exists(os.path.join(output_folder, file_path.name.replace('.dat', '.root'))):
            print("skipping")
            continue

        beamRunNum = getRunNumber(file_path)
        selected_cal_file = None
        minRunOffset = -1
        for cal_file in cal_files:
            if minRunOffset == -1 or abs(beamRunNum - getRunNumber(cal_file)) < minRunOffset:
                minRunOffset = abs(beamRunNum - getRunNumber(cal_file))
                selected_cal_file = cal_file

        print(f"  - Using CAL: {selected_cal_file.name}")
        subprocess.run([
            'bmRawToRootConverter',
            '-i', str(file_path),
            '-c', os.path.join(output_folder, str(selected_cal_file.name).replace(".dat", ".root")),
            '-o', output_folder,
            '-t', str(3),
            '-if' # if beam .txt output
        ])
        
        




if __name__ == '__main__':
    main()
