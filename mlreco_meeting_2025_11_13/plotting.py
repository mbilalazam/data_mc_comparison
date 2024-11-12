import pandas as pd
import matplotlib.pyplot as plt

# Define bin centers for each histogram type
bin_centers_dict = {
    "mult": [0.5 + i for i in range(20)],  # 0.5, 1.5, ..., 19.5
    "length": [1 + i * 2 for i in range(50)],  # 1, 3, 5, ..., 99
    "cosTheta": [-0.99 + i * 0.02 for i in range(100)],  # -0.99, -0.97, ..., 0.99
    "energy": [0.01 + i * 0.02 for i in range(50)]  # 0.01, 0.03, ..., 0.99
}

# Define particle types and corresponding names for histogram types
particle_titles = {
    "reco_length_prim_muon": "Muons",
    "reco_length_prim_pion": "Charged Pions",
    "reco_length_prim_proton": "Protons",
    "reco_cosTheta_prim_muon": "Muons",
    "reco_cosTheta_prim_pion": "Charged Pions",
    "reco_cosTheta_prim_proton": "Protons",
    "reco_energy_prim_muon": "Muons",
    "reco_energy_prim_pion": "Charged Pions",
    "reco_energy_prim_proton": "Protons",
    "reco_mult_prim_total": "Charged Track",
    "reco_mult_prim_shower": "Shower"
}

# Function to get bin centers based on histogram type
def get_bin_centers(histogram_name):
    if "mult" in histogram_name:
        return bin_centers_dict["mult"]
    elif "length" in histogram_name:
        return bin_centers_dict["length"]
    elif "cosTheta" in histogram_name:
        return bin_centers_dict["cosTheta"]
    elif "energy" in histogram_name:
        return bin_centers_dict["energy"]
    else:
        raise ValueError(f"Unknown histogram type for {histogram_name}")

# Function to get x-axis label based on histogram type
def get_xlabel(histogram_name):
    if "mult" in histogram_name:
        return "Multiplicity"
    elif "length" in histogram_name:
        return "Track Length (cm)"
    elif "cosTheta" in histogram_name:
        return "Cosine Theta (in degrees)"
    elif "energy" in histogram_name:
        return "Energy (GeV)"
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
    elif "cosTheta" in histogram_name:
        return f"{particle} Angular Distribution ({module_name})"
    elif "energy" in histogram_name:
        return f"{particle} Energy Distribution ({module_name})"
    else:
        return f"{histogram_name} ({module_name}) Comparison Across Datasets"

# Function to plot histogram with error bars
def plot_histogram_with_errorbars(bin_centers, entries, lower_errors, upper_errors, color, label):
    lower_error = [entry - lower for entry, lower in zip(entries, lower_errors)]
    upper_error = [upper - entry for entry, upper in zip(entries, upper_errors)]
    lower_error = [max(0, err) for err in lower_error]
    upper_error = [max(0, err) for err in upper_error]
    asymmetric_error = [lower_error, upper_error]
    plt.step(bin_centers, entries, where='mid', color=color, linewidth=1, label=label)
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
def calculate_normalization_factor(set1_entries, comparison_entries):
    set1_total = sum(set1_entries)
    comparison_total = sum(comparison_entries)
    return set1_total / comparison_total if comparison_total > 0 else 1

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
    elif "energy" in histogram_name:
        plt.xticks([i * 0.1 for i in range(11)])  # Major ticks at intervals of 0.1 from 0 to 1
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(2))  # Minor ticks at intervals of 2
    elif "cosTheta" in histogram_name:
        plt.xticks([i * 0.1 for i in range(11)])  # Major ticks at intervals of 0.1 from 0 to 1
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(0.02))  # Minor ticks every 0.02 within each major interval

# Function to plot histograms of the same type across datasets and modules
def plot_histogram_type_across_modules(histogram_name, dataset_files, stats_files, output_file_prefix):
    # Retrieve the correct bin centers for the histogram type
    bin_centers = get_bin_centers(histogram_name)
    xlabel = get_xlabel(histogram_name)  # Get x-axis label based on histogram type
    
    # Extract and use 'modules_0123' for normalization factors
    set1_0123_bin_centers, set1_0123_entries, set1_0123_lower, set1_0123_upper = read_csv_data(dataset_files['set1']['modules_0123'])
    set2_0123_bin_centers, set2_0123_entries, set2_0123_lower, set2_0123_upper = read_csv_data(dataset_files['set2']['modules_0123'])
    set3_0123_bin_centers, set3_0123_entries, set3_0123_lower, set3_0123_upper = read_csv_data(dataset_files['set3']['modules_0123'])

    # Calculate normalization factors for set2 and set3 relative to set1
    set2_normalization_factor = calculate_normalization_factor(set1_0123_entries, set2_0123_entries)
    set3_normalization_factor = calculate_normalization_factor(set1_0123_entries, set3_0123_entries)

    # Iterate over each category (module) and plot on separate canvases
    for module in dataset_files['set1'].keys():
        plt.figure(figsize=(8, 6))  # Figure size

        # Read and plot set1 data for the current module in red
        _, set1_entries, set1_lower, set1_upper = read_csv_data(dataset_files['set1'][module])
        set1_stats = read_stats(stats_files['set1'][module], histogram_name)
        set1_label = f"Sandbox Data (v4) ({int(set1_stats[0])}, {set1_stats[1]:.2f}, {set1_stats[2]:.2f})"
        plot_histogram_with_errorbars(bin_centers, set1_entries, set1_lower, set1_upper, 'red', set1_label)

        # Read, normalize, and plot set2 data for the current module in blue
        _, set2_entries, set2_lower, set2_upper = read_csv_data(dataset_files['set2'][module])
        normalized_set2_entries, normalized_set2_lower, normalized_set2_upper = normalize_to_set1(set2_entries, set2_lower, set2_upper, set2_normalization_factor)
        set2_stats = read_stats(stats_files['set2'][module], histogram_name)
        set2_label = f"Simulation (MR6) ({int(set2_stats[0])}, {set2_stats[1]:.2f}, {set2_stats[2]:.2f})"
        plot_histogram_with_errorbars(bin_centers, normalized_set2_entries, normalized_set2_lower, normalized_set2_upper, 'blue', set2_label)

        # Read, normalize, and plot set3 data for the current module in green
        _, set3_entries, set3_lower, set3_upper = read_csv_data(dataset_files['set3'][module])
        normalized_set3_entries, normalized_set3_lower, normalized_set3_upper = normalize_to_set1(set3_entries, set3_lower, set3_upper, set3_normalization_factor)
        set3_stats = read_stats(stats_files['set3'][module], histogram_name)
        set3_label = f"Simulation (MR6.1) ({int(set3_stats[0])}, {set3_stats[1]:.2f}, {set3_stats[2]:.2f})"
        plot_histogram_with_errorbars(bin_centers, normalized_set3_entries, normalized_set3_lower, normalized_set3_upper, 'green', set3_label)

        # Title, labels, and legend
        plt.title(get_plot_title(histogram_name, module), fontsize=16)
        plt.xlabel(xlabel, fontsize=14)
        plt.ylabel('Normalized Counts', fontsize=14)
        plt.legend(loc='upper right')

        # Set y-axis limits based on the maximum value of set1 data for this module
        max_y_value = max(set1_entries) * 1.2  # Add 20% buffer
        plt.ylim(0, max_y_value)
        plt.xlim(min(bin_centers), max(bin_centers))  # Set x-axis range based on bin centers

        # Apply custom x-ticks based on histogram type
        set_custom_xticks(histogram_name)

        # Disable grid
        plt.grid(False)

        # Save the plot as a PNG
        output_file = f"{output_file_prefix}_{module}.png"
        plt.savefig(output_file)
        plt.close()

# Main function to iterate over histogram types and plot
if __name__ == "__main__":
    # Define the histogram types and datasets
    histogram_types = [
        "reco_mult_prim_total", "reco_mult_prim_shower", "reco_length_prim_muon",
        "reco_length_prim_pion", "reco_length_prim_proton", "reco_cosTheta_prim_muon",
        "reco_cosTheta_prim_pion", "reco_cosTheta_prim_proton", "reco_energy_prim_muon",
        "reco_energy_prim_pion", "reco_energy_prim_proton"
    ]
    
    # Define datasets and categories
    sets = ["sandboxV4", "MR6", "MR6p1"]
    categories = ["modules_0123", "modules_0", "modules_1", "modules_2", "modules_3"]

    # Set up file paths for each dataset and category
    for histogram_name in histogram_types:
        dataset_files = {
            'set1': {category: f"{histogram_name}_sandboxV4_{category}.csv" for category in categories},
            'set2': {category: f"{histogram_name}_MR6_{category}.csv" for category in categories},
            'set3': {category: f"{histogram_name}_MR6p1_{category}.csv" for category in categories}
        }
        stats_files = {
            'set1': {category: f"stats_sandboxV4_{category}.csv" for category in categories},
            'set2': {category: f"stats_MR6_{category}.csv" for category in categories},
            'set3': {category: f"stats_MR6p1_{category}.csv" for category in categories}
        }

        # Output file prefix for each histogram type
        output_file_prefix = f"{histogram_name}_comparison"

        # Plot histograms across modules for the current histogram type
        plot_histogram_type_across_modules(histogram_name, dataset_files, stats_files, output_file_prefix)
