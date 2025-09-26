"""
This program takes two arguments:
- The first is required and is the location of the outputs directory in which simulation slices are located.
- The second is optional. If any argument, e.g. show, is given here then the figure will be shown instead.
"""

from pytools import read_h5_array
import h5py
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import json

def extract_em_energy(fnames):
    "Returns a dictionary of lists specifying the amount of energy stored in each field component of each simulation slice."
    fields = ['bx', 'by', 'bz', 'ex', 'ey', 'ez']
    # Create empty dict to store results
    field_energies = { field: [] for field in fields}
    # Iteratie over files and field components, measure energy therein
    for fname in fnames:
        f5 = h5py.File(fname, "r")
        for field in fields:
            field_energies[field].append(np.sum(read_h5_array(f5, field)**2))
    # Convert lists of energies to numpy arrays
    field_energies = { f: np.array(e) for f, e in field_energies.items()}
    # Add a total energy statistic
    field_energies['total'] = sum(field_energies.values())
    # Done
    return field_energies

def make_cherenkov_plot(laps, field_energies, outname, total_kinetic_energy=None):
    for field, energies in field_energies.items():
        plt.plot(laps, energies, label=field)
        plt.yscale('log')
    plt.xlabel('Laps $\\rightarrow$')
    plt.ylabel('Energy per field component $\\rightarrow$')
    if total_kinetic_energy:
        plt.axhline(y=total_kinetic_energy, color='r', linestyle='-', label='KE limit')
    plt.legend()
    if outname != 'show':
        plt.savefig(outname)
    else:
        plt.show()

def read_json(fname):
    file = open(fname, 'r')
    return json.load(file)


# Get the output folder argument
outputs_dir = sys.argv[1]
print(f'Analysing simulation output located in "{outputs_dir}".')
# Remove any trailing / from the output folder
outputs_dir = outputs_dir[:-1] if (outputs_dir[-1]=='/') else outputs_dir
# Load the parameters json
#params_dict = read_json(f'{outputs_dir}/params.json')
# Get a list of all the field files in the folder and extract a list of their lap numbers.
prefix = f"{outputs_dir}/flds_"
postfix = f".h5"
field_filenames = glob.glob(f"{prefix}*{postfix}")
laps = [ int(fname[len(prefix):-len(postfix)]) for fname in field_filenames]
# Make sure this list is sorted in order of ascending time
laps, field_filenames = zip(*sorted(zip(laps, field_filenames), key=lambda tuple : tuple[0]))
# Print a statistic:
file_completeness_fraction = (laps[1]-laps[0])*(len(laps) - 1) / (laps[-1]-laps[0])
print(f"Found {len(laps)} files, ranging from lap {laps[0]} to lap {laps[-1]}.")
print(f"{100*file_completeness_fraction:n}% of files present within this range.")

# Next, extract the total EM energy per file:
field_energies = extract_em_energy(field_filenames)

# Do outputs
if len(sys.argv) > 2:
    output_filename = 'show'
else:
    output_filename = outputs_dir+'/field_energies.pdf'
total_kinetic_energy = None
#total_kinetic_energy = params_dict['total_kinetic_energy']
make_cherenkov_plot(laps, field_energies, output_filename, total_kinetic_energy)
print(f"Wrote summary figure to {output_filename}.")