import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
from PIL import Image
import csv
import numpy as np
import sys
import glob
import re
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import Draw 
#import tools

# Call the function when the script runs
if __name__ == "__main__":
    print("Starting the QCxMS2 plotter")


def read_csv(file_path):
    """
    Reads a CSV file and returns two lists: m/z and intensity.
    """
    mz_values = []
    intensities = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row:  # Skip empty rows
                mz_values.append(float(row[0]))
                intensities.append(float(row[1]))
    return mz_values, intensities

def normalize_intensities(intensities):
    """
    Normalizes the intensities to 100%.
    """
    max_intensity = max(intensities)
    return [i / max_intensity * 100 for i in intensities]

def find_high_intensity_peaks(file_path, imp_frags, threshold=2000):
    """
    Searches for lines in allpeaks.dat where the intensity (third column) exceeds the threshold.

    Parameters:
    - file_path: Path to the allpeaks.dat file.
    - threshold: Intensity threshold (default 2000).
    - imp_frags: Array of important fragments (directory, m/z value, intensity).
    """
    class fragment:
        def __init__(self, dir, mz, intensity):
            self.dir = dir
            self.mz = mz
            self.intensity = intensity 
   


    # imp_frags is passed as a parameter, no need to re-initialize it here

    with open(file_path, 'r') as file:
        for line in file:
            # Skip empty lines or improperly formatted lines
            if not line.strip():
                continue

            parts = line.split()
            if len(parts) >= 3:  # Ensure there are at least 3 columns
                try:
                    intensity = float(parts[2])  # Parse the intensity (third column)
                    if intensity > threshold:
                        print(line.strip())  # Print the entire line
                        imp_frag = fragment(parts[0],float(parts[1]),float(parts[2]))
                        imp_frags.append(imp_frag)
                except ValueError:
                    # Skip lines where parsing fails (e.g., header lines or invalid data)
                    continue

def plot_mass_spectrum(peaks_file, exp_file):
    """
    Plots the mass spectrum from peaks.csv and exp.csv files with normalization.
    """
    # Read data from files
    peaks_mz, peaks_intensity = read_csv(peaks_file)
    exp_mz, exp_intensity = read_csv(exp_file)

    # Convert lists to NumPy arrays
    peaks_mz = np.array(peaks_mz)
    peaks_intensity = np.array(peaks_intensity)
    exp_mz = np.array(exp_mz)
    exp_intensity = np.array(exp_intensity)

    # Normalize the intensities to 100%
    peaks_intensity = normalize_intensities(peaks_intensity)
    exp_intensity = normalize_intensities(exp_intensity)

    # Invert the intensities of the experimental data
    inverted_exp_intensity = [-i for i in exp_intensity]

    # Plotting
    plt.figure(figsize=(10, 6))
    
    stems1, _, _ = plt.stem(peaks_mz, peaks_intensity, 
                         linefmt='#07529a', markerfmt='', basefmt=' ', label="QCxMS2")  # No markers, no baseline
    plt.setp(stems1, linewidth=3)  # Thicker stems doesnt work?

    # Plot experiment data with inverted intensities
    stems2, _, _ = plt.stem(exp_mz, inverted_exp_intensity, 
                         linefmt='#909085', markerfmt='', basefmt=' ', label="experiment")  # No markers, no baseline
    plt.setp(stems2, linewidth=3)   # Thicker stems doesnt work?


    # Add a single horizontal baseline from (0,0) to (max_X,0)
    plt.axhline(y=0, xmin=0, xmax=1, color='black', linewidth=2, linestyle='-')
    # Customizing the plot
    #plt.title('Mass Spectrum (Normalized to 100%)')
    plt.title('')
    plt.xlabel('m/z')
    plt.ylabel('Relative Intensity (%)')
    plt.legend(loc='upper left')
    plt.grid(alpha=0.3)
    
    # Find important fragments
    fragment_file = "allpeaks.dat"
    imp_frags = []  
    find_high_intensity_peaks(fragment_file,imp_frags,threshold=imp_frags_threshold)

    # Add important fragments to the plot
    if (depict_fragments):
        add_lewis_formulas(imp_frags,peaks_mz, peaks_intensity)

    # Add reaction barriers o fimportant fragments the the plot
    if (depict_barriers):
        for fragment in imp_frags:
            # check if fragment is fragment or isomer
            if os.path.isfile(fragment.dir+'/fragment.xyz'):
                reactiondir = fragment.dir[:-2] # remove fx from fragment dir
            elif os.path.isfile(fragment.dir+'/isomer.xyz'):
                reactiondir = fragment.dir # isomer is also reaction dir

            barrier = read_barrier_file(reactiondir)
            precursor_dir = find_precursordir(reactiondir)

            fragment_file = "allpeaks.dat"
            precursor_mz = find_fragment_mass('allpeaks.dat',precursor_dir)
            print("Precursor m/z",precursor_mz)
            
            # set for now to 50
            arrow_heigth = get_arrow_height()
            
            successor_mz = fragment.mz
            successor_intensity = fragment.intensity / 100

            print(f"Adding reaction arrow for fragment {fragment.mz} with barrier {barrier}",precursor_mz,arrow_heigth,successor_mz,successor_intensity)
            add_reaction_arrow(barrier,precursor_mz,arrow_heigth,successor_mz,successor_intensity)
            


    # Show the plot
    plt.tight_layout()
    plt.show()

def arrow_height_generator():
    """Generator that yields alternating positive and negative arrow heights."""
    height = 10
    sign = 1
    while True:
        yield height * sign
        sign *= -1
        height += 10


def get_arrow_height():
    """Fetch the next arrow height from the generator."""
    return next(arrow_height_gen)



def find_fragment_mass(filename: str, fragment_dir: str) -> float:
    print(f"Searching for fragment {fragment_dir} in {filename}")
    # if fragment dir is isomer go one level up
    fragment_dir = re.sub(r'p\d+$', '', fragment_dir)
    print(f"Searching for fragment {fragment_dir} in {filename}")
    if (fragment_dir == ''):
        fragment_dir = '.'
    if (fragment_dir == '.'):
        fragment_dir = 'input'
        with open(filename, 'r') as file:
            for line in file:
                parts = line.split()  # Split by whitespace
                if len(parts) >= 4 and parts[0] == fragment_dir:

                    return float(parts[2])  # Convert the third column to float and return it
    else:
                   
        with open(filename, 'r') as file:
            for line in file:
                parts = line.split()  # Split by whitespace
                if len(parts) >= 3 and parts[0] == fragment_dir:
                    return float(parts[1])  # Convert the second column to float and return it

def find_precursordir(dir: str) -> str:

    p_count = dir.count('p')
    if p_count == 1:
        return '.'

    precursordir = re.sub(r'p\d$', '', dir)
    return precursordir
   

def read_barrier_file(reactiondir: str) -> str:
    """
    Searches for a file in the given directory whose name starts with 'barrier_'
    and returns its content as a string. If multiple files match, the first one found is used.

    Args:
        reactiondir (str): The directory to search in.

    Returns:
        str: The content of the file.

    Raises:
        FileNotFoundError: If no matching file is found.
    """
    # Create the pattern for files starting with "barrier_" in the reactiondir
    pattern = os.path.join(reactiondir, "barrier_*")
    # Use glob to find all matching files
    matching_files = glob.glob(pattern)
    
    if not matching_files:
        raise FileNotFoundError(f"No file starting with 'barrier_' was found in {reactiondir}")
    
    # For example, we take the first matching file.
    file_to_read = matching_files[0]
    
    with open(file_to_read, 'r') as file:
        content = float(file.read().strip())
    
    barrier = round(content, 1)
    barrier = str(barrier) +" eV"
    return barrier
   
def add_lewis_formulas(imp_frags,peaks_mz, peaks_intensity):

    for fragment in imp_frags:
        #print(f"Important fragment: {fragment.dir} {fragment.mz} {fragment.intensity}")
    # Generate Lewis structure from XYZ file and save to mol.png

        # check if fragment is fragment or isomer
        if os.path.isfile(fragment.dir+'/fragment.xyz'):
            fname = fragment.dir+'/fragment.xyz'
        else:
            fname = fragment.dir+'/isomer.xyz'

        
        depict_xyz(fname)

    
        # add Lewis structures
        molpic = Image.open("mol.png") # Open the image
        molpic_array = np.array(molpic) 
        # Define the position in data coordinates (m/z = 32, intensity = 100)
        molpic_mz = fragment.mz  # m/z value
        molpic_mz = molpic_mz  # shift a little bit for better visualization
        
        # get intensity of the peak corresponding to the fragment
    
        
        mz_matching = []
        mz_matching = np.where(abs(peaks_mz - fragment.mz) <=  frag_mztol)[0]  # Get the index

        # depict only fragment with highest intensity if multiple fragments for one peak are found
        max_intensity = 0
        if (len(mz_matching) > 0):
            for index in mz_matching:
                if (peaks_intensity[index] > max_intensity):
                    max_intensity = peaks_intensity[index]  # Get the corresponding intensity
        else:
            print(f"Fragment {fragment.mz} not found in peaks")
            continue        

    
        molpic_intensity = max_intensity  # Get the corresponding intensity

        molpic_intensity = molpic_intensity  # shift a little bit for better visualization
        # Create an OffsetImage object
        molpic_box = OffsetImage(molpic_array, zoom=0.2)  # Adjust zoom to control image size

        # Create a text area with the fragment name
        text_area = TextArea(str(molpic_mz)+" m/z", textprops=dict(color='black', size=12))

        molpic_box_named = VPacker(children=[molpic_box,text_area],
                align="center",  # Center the children horizontally
                pad=0,
                sep=5)  
        # Create an AnnotationBbox to place the image at the specified data coordinates
        print(f"Adding Lewis structure at m/z {molpic_mz} with intensity {molpic_intensity}")
        ab = AnnotationBbox(molpic_box_named, (molpic_mz, molpic_intensity), 
                    box_alignment=(0, 0),  # Box alignment relative to the point
                    xycoords='data',       # Coordinates for the arrow and box
                    boxcoords="offset points", # Position the box relative to the anchor point
                    pad=0.5,               # Padding between the box and anchor
                    arrowprops=dict(arrowstyle="->", lw=1.5, color="black"),
                    bboxprops=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"))
        

        # Add the AnnotationBbox to the plot
        plt.gca().add_artist(ab)


def add_reaction_arrow(barrier,precursor_mz,precursor_intensity,successor_mz,successor_intensity):
            
            text_x = (precursor_mz + successor_mz) / 2
            text_y = precursor_intensity
            # Use annotate() to connect two points with an arrow and add text above the arrow.
            plt.annotate(
                barrier,         # The text you want to display
                xy=(successor_mz, successor_intensity),                  # The anchor point where the arrow points (target)
                xytext=(text_x, text_y),              # The position where the text will be placed (source)
                arrowprops=dict(
                    arrowstyle="->",          # Arrow style
                    lw=1.5,                   # Line width of the arrow
                    color="red"             # Color of the arrow
                ),
                horizontalalignment='center',
                verticalalignment='center'
            )
            plt.annotate(
                '',         # The text you want to display
                xy=(text_x+2, text_y),                  # The anchor point where the arrow points (target)
                xytext=(precursor_mz, precursor_intensity),              # The position where the text will be placed (source)
                arrowprops=dict(
                    arrowstyle="->",          # Arrow style
                    lw=1.5,                   # Line width of the arrow
                    color="red"             # Color of the arrow
                ),
                horizontalalignment='center',
                verticalalignment='center'
            )

def depict_xyz(infile):
    """
    Reads the molecule from the XYZ file and saves it as a PNG file.
    """
    # Read the molecule from the XYZ file
    raw_mol = Chem.MolFromXYZFile(infile)

    # Determine bond orders with explicit charge
    mol = Chem.Mol(raw_mol)

    mol.UpdatePropertyCache(strict=False)

    
    #DetermineBonds is rather unstable. For now we try it and 
    # if it fails fall back to always(?) working connectivity
    # TODO switch to obabel for this 
    try: 
        rdkit.Chem.rdDetermineBonds.DetermineBonds(mol, charge=1)
        print("Bonds successfully detected")
    except:
        print("Bonds could not be detected, depiciting only connectivity")
        rdkit.Chem.rdDetermineBonds.DetermineConnectivity(mol, charge=1)

    # Draw the molecule and save it to a file
    Chem.Draw.MolToFile(mol, 'mol.png', size=(300, 300))

if __name__ == "__main__":

    # some global variables

    # Create a generator instance
    arrow_height_gen = arrow_height_generator()

    depict_fragments = True  # Depict important fragments
    depict_barriers = True  # Depict reaction barriers
    imp_frags_threshold = 3000  # Define the intensity threshold for important fragments normalized to 10000
    frag_mztol = 0.3  # m/z tolerance for matching fragments to peaks


    # Input files (update file names if needed)
    peaks_file = 'peaks.csv'
    exp_file = 'exp.csv'

   

    # Plot the mass spectrum
    plot_mass_spectrum(peaks_file, exp_file)
    
    

    



