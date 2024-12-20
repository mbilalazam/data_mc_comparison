import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import csv

# Initialize a list to store scale factors for writing to CSV
scale_factors = []

# Define bin centers for each histogram type
bin_centers_dict = {
    "mult": [0.5 + i for i in range(20)],  # 0.5, 1.5, ..., 19.5
    "length": [1 + i * 2 for i in range(50)],  # 1, 3, 5, ..., 99
    "cosTheta": [-0.99 + i * 0.02 for i in range(100)],  # -0.99, -0.97, ..., 0.99
    "energy": [0.01 + i * 0.02 for i in range(50)],  # 0.01, 0.03, ..., 0.99
    "start_xyz": [round(-68.6 + i * 2.8, 1) for i in range(50)], # -68.6, -65.8, -63, -60.2, ... ,  63, 65.8, 68.6
    "nu": [round(-68.6 + i * 2.8, 1) for i in range(50)] # -68.6, -65.8, -63, -60.2, ... ,  63, 65.8, 68.6
}

# Define particle types and corresponding names for histogram types
particle_titles = {
    "reco_mult_prim_total": "Charged Track", 
    "reco_mult_prim_shower": "Shower", 
    "reco_mult_muon": "Muons", 
    "reco_mult_pion": "Charged Pions", 
    "reco_mult_proton": "Protons",
    "reco_mult_kaon": "Charged Kaons", 
    "reco_length_prim_muon": "Muons", 
    "reco_length_prim_pion": "Charged Pions", 
    "reco_length_prim_proton": "Protons", 
    "reco_length_prim_kaon": "Charged Kaons",
    "reco_energy_prim_muon": "Muons", 
    "reco_energy_prim_pion": "Charged Pions", 
    "reco_energy_prim_proton": "Protons", 
    "reco_energy_prim_kaon": "Charged Kaons",
    "reco_cosTheta_prim_muon": "Muons", 
    "reco_cosTheta_prim_pion": "Charged Pions", 
    "reco_cosTheta_prim_proton": "Protons", 
    "reco_cosTheta_prim_kaon": "Charged Kaons",
    "reco_start_x_muon": "Muons",
    "reco_start_y_muon": "Muons",
    "reco_start_z_muon": "Muons",
    "reco_start_x_pion": "Charged Pions",
    "reco_start_y_pion": "Charged Pions",
    "reco_start_z_pion": "Charged Pions",
    "reco_start_x_proton": "Protons",
    "reco_start_y_proton": "Protons",
    "reco_start_z_proton": "Protons",
    "reco_start_x_kaon": "Charged Kaons",
    "reco_start_y_kaon": "Charged Kaons",
    "reco_start_z_kaon": "Charged Kaons",
    "reco_nu_vtx_x_CC_LArFV": "Neutrino",
    "reco_nu_vtx_y_CC_LArFV": "Neutrino",
    "reco_nu_vtx_z_CC_LArFV": "Neutrino"
}


# Function to get bin centers based on histogram type
def get_bin_centers(histogram_name):
    if "mult" in histogram_name:
        return bin_centers_dict["mult"]
    elif "length" in histogram_name:
        return bin_centers_dict["length"]
    elif "energy" in histogram_name:
        return bin_centers_dict["energy"]
    elif "cosTheta" in histogram_name:
        return bin_centers_dict["cosTheta"]
    elif "start_x" in histogram_name or "start_y" in histogram_name or "start_z" in histogram_name:
        return bin_centers_dict["start_xyz"]
    elif "nu" in histogram_name:
        return bin_centers_dict["nu"]
    else:
        raise ValueError(f"Unknown histogram type for {histogram_name}")

# Function to get x-axis label based on histogram type
def get_xlabel(histogram_name):
    if "mult" in histogram_name:
        return "Multiplicity"
    elif "length" in histogram_name:
        return "Track Length (cm)"
    elif "energy" in histogram_name:
        return "Energy (GeV)"
    elif "cosTheta" in histogram_name:
        return "Cosine Theta"
    elif "start_x" in histogram_name:
        return "Start-x (cm)"
    elif "start_y" in histogram_name:
            return "Start-y (cm)" 
    elif "start_z" in histogram_name:
            return "Start-z (cm)"         
    elif "nu_vtx_x_CC_LArFV" in histogram_name:
            return "Vertex-x (cm)"         
    elif "nu_vtx_y_CC_LArFV" in histogram_name:
            return "Vertex-y (cm)"         
    elif "nu_vtx_z_CC_LArFV" in histogram_name:
            return "Vertex-z (cm)"         
    else:
        return "Bin Center"

# Function to generate the plot title based on histogram type and module
def get_plot_title(histogram_name, module):
    particle = particle_titles.get(histogram_name, "")
    module_name = module.replace("modules_", "Module ").replace("modules_0123", "Modules 0123")
    
    if "mult" in histogram_name:
        return f"{particle} Multiplicity Distribution ({module_name})"
    elif "length" in histogram_name:
        return f"{particle} Track Length Distribution ({module_name})"
    elif "energy" in histogram_name:
        return f"{particle} Energy Distribution ({module_name})"
    elif "cosTheta" in histogram_name:
        return f"{particle} Angular Distribution ({module_name})"
    elif "start_x" in histogram_name:
        return f"{particle} Start-x Distribution ({module_name})"        
    elif "start_y" in histogram_name:
        return f"{particle} Start-y Distribution ({module_name})" 
    elif "start_z" in histogram_name:
        return f"{particle} Start-z Distribution ({module_name})"         
    elif "nu_vtx_x_CC_LArFV" in histogram_name:
        return f"{particle} Vertex-x Distribution ({module_name})"         
    elif "nu_vtx_y_CC_LArFV" in histogram_name:
        return f"{particle} Vertex-y Distribution ({module_name})"         
    elif "nu_vtx_z_CC_LArFV" in histogram_name:
        return f"{particle} Vertex-z Distribution ({module_name})"         
    else:
        return f"{histogram_name} ({module_name}) Comparison Across Datasets"

# Function to plot histogram with error bars
def plot_histogram_with_errorbars(bin_centers, entries, lower_errors, upper_errors, color, label):
    lower_error = [entry - lower for entry, lower in zip(entries, lower_errors)]
    upper_error = [upper - entry for entry, upper in zip(entries, upper_errors)]
    lower_error = [max(0, err) for err in lower_error]
    upper_error = [max(0, err) for err in upper_error]
    asymmetric_error = [lower_error, upper_error]
    plt.step(bin_centers, entries, where='mid', color=color, linewidth=2, label=label)
    plt.errorbar(bin_centers, entries, yerr=asymmetric_error, fmt='none', color=color, capsize=2)

# Function to read CSV data
def read_csv_data(csv_file):
    data = pd.read_csv(csv_file)
    return data['Bin Center'], data['Entries'], data['Lower End'], data['Upper End']
    
def read_stats(stats_file, histogram_name):
    # Attempt to load the CSV file with default settings
    try:
        stats_data = pd.read_csv(stats_file)
        # Ensure 'histogram' is in columns; otherwise, print an error and column names
        if 'histogram' not in stats_data.columns:
            print(f"Error: 'histogram' column not found in {stats_file}. Columns available: {stats_data.columns}")
            return None, None, None
        stats_row = stats_data[stats_data["histogram"] == histogram_name].iloc[0]
        return stats_row["entries"], stats_row["mean"], stats_row["std_dev"]
    except Exception as e:
        print(f"Failed to read {stats_file} with error: {e}")
        return None, None, None
    

# Function to calculate normalization factor based on 'set1' total entries
# def calculate_normalization_factor(set1_entries, comparison_entries):
    # set1_total = sum(set1_entries)
    # comparison_total = sum(comparison_entries)
    # return set1_total / comparison_total if comparison_total > 0 else 1

# Function to calculate normalization factor for set3 (fixed at 0.01)
def calculate_normalization_factor(set1_entries, comparison_entries):
    return 0.0091364902506964

# Function to apply normalization based on set1 to entries and errors
def normalize_to_set1(entries, lower_end, upper_end, normalization_factor):
    normalized_entries = [entry * normalization_factor for entry in entries]
    normalized_lower_end = [entry * normalization_factor for entry in lower_end]
    normalized_upper_end = [entry * normalization_factor for entry in upper_end]
    return normalized_entries, normalized_lower_end, normalized_upper_end

# Function to set custom x-ticks based on histogram type
def set_custom_xticks(histogram_name):
    if "mult" in histogram_name:
        plt.xticks(range(1, 21))  # Ticks at each integer from 1 to 20
    elif "length" in histogram_name:
        plt.xticks(range(0, 101, 10))  # Major ticks at intervals of 10
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(5))  # Minor ticks at intervals of 5
    elif "cosTheta" in histogram_name:
        plt.xticks([round(x, 2) for x in [-0.99, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 0.99]])  # Ticks at specified intervals
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(0.10))  # Minor ticks at intervals of 0.10
    elif "energy" in histogram_name:
        plt.xticks([round(x, 2) for x in [i / 10.0 for i in range(11)]])  # Ticks at 0.1, 0.2, ..., 1.0
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(0.05))  # Minor ticks at intervals of 0.05
    elif "start_x" in histogram_name:
        plt.xticks(range(-70, 71, 10))  # Major ticks at intervals of 10
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(5))  # Minor ticks at intervals of 5
    elif "start_y" in histogram_name:
        plt.xticks(range(-70, 71, 10))  # Major ticks at intervals of 10
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(5))  # Minor ticks at intervals of 5
    elif "start_z" in histogram_name:
        plt.xticks(range(-70, 71, 10))  # Major ticks at intervals of 10
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(5))  # Minor ticks at intervals of 5
    elif "nu_vtx_x_CC_LArFV" in histogram_name:
        plt.xticks(range(-70, 71, 10))  # Major ticks at intervals of 10
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(5))  # Minor ticks at intervals of 5
    elif "nu_vtx_y_CC_LArFV" in histogram_name:
        plt.xticks(range(-70, 71, 10))  # Major ticks at intervals of 10
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(5))  # Minor ticks at intervals of 5
    elif "nu_vtx_z_CC_LArFV" in histogram_name:
        plt.xticks(range(-70, 71, 10))  # Major ticks at intervals of 10
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(5))  # Minor ticks at intervals of 5
    else:
        print(f"Histogram type for {histogram_name} not recognized.")




# Function to plot histograms of the same type across datasets and only for modules_0123
def plot_histogram_for_modules_0123(histogram_name, dataset_files, stats_files, output_file_prefix):
# def plot_histogram_for_modules_0123(histogram_name, dataset_files, stats_files, scaled_files, output_file_prefix):
    # Retrieve the correct bin centers for the histogram type
    bin_centers = get_bin_centers(histogram_name)
    xlabel = get_xlabel(histogram_name)  # Get x-axis label based on histogram type

    # Extract and use 'modules_0123' for normalization factors and plotting
    _, set1_entries, set1_lower, set1_upper = read_csv_data(dataset_files['set1']['modules_0123'])
    _, set3_entries, set3_lower, set3_upper = read_csv_data(dataset_files['set3']['modules_0123'])
    # _, set5_entries, set5_lower, set5_upper = read_csv_data(dataset_files['set5']['modules_0123'])

    # Calculate normalization factors for set3 and scaled relative to set1
    set3_normalization_factor = calculate_normalization_factor(set1_entries, set3_entries)
    # set5_normalization_factor = calculate_normalization_factor(set1_entries, set5_entries)

    plt.figure(figsize=(8, 6))  # Figure size

    # Read and plot set1 data for modules_0123 in red
    set1_stats = read_stats(stats_files['set1']['modules_0123'], histogram_name)
    set1_label = f"Sandbox Data (v4) ({int(set1_stats[0])}, {set1_stats[1]:.2f}, {set1_stats[2]:.2f})"
    plot_histogram_with_errorbars(bin_centers, set1_entries, set1_lower, set1_upper, 'red', set1_label)

    # Read, normalize, and plot set3 data for modules_0123 in green
    normalized_set3_entries, normalized_set3_lower, normalized_set3_upper = normalize_to_set1(set3_entries, set3_lower, set3_upper, set3_normalization_factor)
    set3_stats = read_stats(stats_files['set3']['modules_0123'], histogram_name)
    set3_label = f"Simulation (MR6.1) ({int(set3_stats[0])}, {set3_stats[1]:.2f}, {set3_stats[2]:.2f})\nNormalization Factor: {set3_normalization_factor:.5f}"
    plot_histogram_with_errorbars(bin_centers, normalized_set3_entries, normalized_set3_lower, normalized_set3_upper, 'green', set3_label)

    # Read, normalize, and plot set5 data for modules_0123 in green
    # normalized_set5_entries, normalized_set5_lower, normalized_set5_upper = normalize_to_set1(set5_entries, set5_lower, set5_upper, set5_normalization_factor)
    # set5_stats = read_stats(stats_files['set5']['modules_0123'], histogram_name)
    # set5_label = f"Scaled Simulation (MR6.1) ({int(set5_stats[0])}, {set5_stats[1]:.2f}, {set5_stats[2]:.2f})\nNormalization Factor: {set5_normalization_factor:.2f}"
    # plot_histogram_with_errorbars(bin_centers, normalized_set5_entries, normalized_set5_lower, normalized_set5_upper, 'blue', set5_label)

    # Title, labels, and legend
    plt.title(get_plot_title(histogram_name, 'modules_0123'), fontsize=16)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel('Normalized Counts', fontsize=14)
    plt.legend(loc='upper right')

    # Set y-axis limits based on the maximum value of set1 data for this module
    max_y_value = max(set1_entries) * 1.6  # Add 20% buffer
    plt.ylim(0, max_y_value)
    plt.xlim(min(bin_centers), max(bin_centers))  # Set x-axis range based on bin centers

    # Apply custom x-ticks based on histogram type
    set_custom_xticks(histogram_name)

    # Disable grid
    plt.grid(False)

    # Save the plot as a PNG
    output_file = f"{output_file_prefix}_modules_0123.png"
    plt.savefig(output_file)
    plt.close()

# Main function to iterate over histogram types and plot only for modules_0123
if __name__ == "__main__":
    # Define the histogram types and datasets
    histogram_types = [
        "reco_mult_prim_total", 
        "reco_mult_prim_shower", 
        "reco_mult_muon", 
        "reco_mult_pion", 
        "reco_mult_proton",
        "reco_mult_kaon", 
        "reco_length_prim_muon", 
        "reco_length_prim_pion", 
        "reco_length_prim_proton", 
        "reco_length_prim_kaon",
        "reco_energy_prim_muon", 
        "reco_energy_prim_pion", 
        "reco_energy_prim_proton", 
        "reco_energy_prim_kaon",
        "reco_cosTheta_prim_muon", 
        "reco_cosTheta_prim_pion", 
        "reco_cosTheta_prim_proton", 
        "reco_cosTheta_prim_kaon",
        "reco_start_x_muon",
        "reco_start_y_muon",
        "reco_start_z_muon",
        "reco_start_x_pion",
        "reco_start_y_pion",
        "reco_start_z_pion",
        "reco_start_x_proton",
        "reco_start_y_proton",
        "reco_start_z_proton",
        "reco_start_x_kaon",
        "reco_start_y_kaon",
        "reco_start_z_kaon",
        "reco_nu_vtx_x_CC_LArFV",
        "reco_nu_vtx_y_CC_LArFV",
        "reco_nu_vtx_z_CC_LArFV"
    ]
    
 

    # Define datasets and categories (only using modules_0123)
    # sets = ["sandboxV4", "MR6p1", "MR6p1_scaled"]
    sets = ["sandboxV4", "MR6p1"]
    categories = ["modules_0123"]

    # Set up file paths for each dataset and category
    for histogram_name in histogram_types:
        dataset_files = {
            'set1': {"modules_0123": f"{histogram_name}_sandboxV4_modules_0123.csv"},
            'set3': {"modules_0123": f"{histogram_name}_MR6p1_modules_0123.csv"},
            # 'set5': {"modules_0123": f"{histogram_name}_MR6p1_modules_0123_scaled.csv"}
        }
        stats_files = {
            'set1': {"modules_0123": f"stats_sandboxV4_modules_0123.csv"},
            'set3': {"modules_0123": f"stats_MR6p1_modules_0123.csv"},
            # 'set5': {"modules_0123": f"stats_MR6p1_modules_0123_scaled.csv"},
        }

        # Output file prefix for each histogram type
        output_file_prefix = f"{histogram_name}_comparison"

        # Calculate normalization factors for MR6.1 and scaled MR6.1 relative to Sandbox Data (v4)
        set1_0123_bin_centers, set1_0123_entries, set1_0123_lower, set1_0123_upper = read_csv_data(dataset_files['set1']['modules_0123'])
        set3_0123_bin_centers, set3_0123_entries, set3_0123_lower, set3_0123_upper = read_csv_data(dataset_files['set3']['modules_0123'])
        # set5_0123_bin_centers, set5_0123_entries, set5_0123_lower, set5_0123_upper = read_csv_data(dataset_files['set5']['modules_0123'])

        set3_normalization_factor = calculate_normalization_factor(set1_0123_entries, set3_0123_entries)
        # set5_normalization_factor = calculate_normalization_factor(set1_0123_entries, set5_0123_entries)        

        # Use a fixed normalization factor for set3
        set3_normalization_factor = 0.0091364902506964
        # set5_normalization_factor = 0.77

        # Print scale factors to the console
        print(f"Scale factors for {histogram_name} (modules_0123):")
        print(f"  - MR6.1: {set3_normalization_factor:.5f}")
        # print(f"  - Scaled MR6.1: {set5_normalization_factor:.2f}")

        # Append scale factors to the list for CSV writing
        # scale_factors.append([histogram_name, f"{set3_normalization_factor:.2f}", f"{set5_normalization_factor:.2f}"])
        scale_factors.append([histogram_name, f"{set3_normalization_factor:.5f}"])

        # Plot histogram for modules_0123 for the current histogram type
        plot_histogram_for_modules_0123(histogram_name, dataset_files, stats_files, output_file_prefix)

    # Write scale factors to CSV
    with open("scale_factors.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # writer.writerow(["histogram", "MR6.1", "Scaled MR6.1"])  # Write header
        writer.writerow(["histogram", "MR6.1"])  # Write header
        writer.writerows(scale_factors)  # Write scale factor data
