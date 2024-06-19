#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Title: Crosslink Distance Analysis in Protein Structures
Author: Zachary Berndsen, ChatGPT
Date: April 2024
Description: 
    This module contains all the functions used in the analysis script
"""

import pandas as pd

#######################################################################################
def read_csv_with_header(csv_file_path):
    """
    Reads a CSV file into a DataFrame using the first row of the CSV as column headers.

    Args:
        csv_file_path (str): The path to the CSV file.

    Returns:
        pd.DataFrame: A DataFrame with the data from the CSV file, using the first row as column names.
    """
    # Read the CSV file, assuming the first row includes the column names
    df = pd.read_csv(csv_file_path)

    return df

#######################################################################################    
def filter_by_spectral_count(df, spectral_count_threshold):
    """
    Filters a DataFrame to include only rows where the spectral count is greater than or equal to a specified threshold.

    Args:
        df (pd.DataFrame): The DataFrame to be filtered, expected to have a 'Spectral Count' column.
        spectral_count_threshold (int): The cutoff value for spectral counts.

    Returns:
        pd.DataFrame: A new filtered DataFrame with rows where 'Spectral Count' is greater than or equal to the threshold.
    """
    if 'Spectral Count' not in df.columns:
        raise ValueError("DataFrame must include a 'Spectral Count' column.")

    # Filter the DataFrame
    filtered_df = df[df['Spectral Count'] >= spectral_count_threshold]
    return filtered_df

#######################################################################################
def domain_crosslink_stats(df, sequence):
    """
    Calculates statistics for predefined domains, including the total number of residues, total lysine residues,
    lysine residues participating in crosslinks, and total unique crosslinks based on a protein sequence.

    Args:
        df (pd.DataFrame): DataFrame containing residue pairs and possibly residue types.
        sequence (str): A string representing the amino acid sequence of the protein.

    Returns:
        pd.DataFrame: A new DataFrame with domain statistics.
    """
    # Predefined domain ranges
    domain_ranges = {
        'NTD': [(1, 1011)],
        'beta-belt': [(1011, 1186), (1276, 1289), (1355, 2016), (2050, 2062), (2757, 3123), (3166, 3179), (3336, 3668), (3700, 3713), (3880, 4058), (4551, 4564)],
        'insert 1': [(1186, 1276)],
        'insert 2': [(1289, 1355)],
        'insert 3': [(2016, 2050)],
        'insert 4': [(2062, 2757)],
        'insert 5': [(3123, 3166)],
        'insert 6': [(3179, 3336)],
        'insert 7': [(3668, 3700)],
        'insert 8': [(3713, 3880)],
        'insert 9': [(4058, 4551)]
    }

    # Find lysine residues by their positions in the sequence
    lysine_residues = {index + 1 for index, residue in enumerate(sequence) if residue == 'K'}

    output_df = pd.DataFrame(columns=['Domain', 'Total Residues', 'Total Lysine Residues',
                                      'Lysine Residues in Crosslinks', 'Total Unique Crosslinks'])

    for domain, ranges in domain_ranges.items():
        domain_residues = set()
        for start, end in ranges:
            domain_residues.update(range(start, end + 1))

        lysines_in_domain = domain_residues.intersection(lysine_residues)

        domain_crosslinks = df[(df['Residue1'].isin(domain_residues)) | (df['Residue2'].isin(domain_residues))]
        
        unique_crosslinks = set()
        for _, row in domain_crosslinks.iterrows():
            pair = tuple(sorted([row['Residue1'], row['Residue2']]))
            unique_crosslinks.add(pair)

        lysines_in_crosslinks = {res for pair in unique_crosslinks for res in pair if res in lysines_in_domain}

        output_df = output_df.append({
            'Domain': domain,
            'Total Residues': len(domain_residues),
            'Total Lysine Residues': len(lysines_in_domain),
            'Lysine Residues in Crosslinks': len(lysines_in_crosslinks),
            'Total Unique Crosslinks': len(unique_crosslinks)
        }, ignore_index=True)

    return output_df

#######################################################################################
def summarize_by_domain_association(df):
    """
    Prepares a DataFrame by adding a 'Canonical Pair' column that uniquely identifies crosslinks and then
    summarizes it by the 'Domain Association' column, aggregating various statistics.
    
    Args:
        df (pd.DataFrame): The DataFrame to be prepared and summarized.
        
    Returns:
        pd.DataFrame: A DataFrame with aggregated statistics for each domain association.
    """
    # Prepare the DataFrame by adding 'Canonical Pair'
    df['Canonical Pair'] = df.apply(lambda x: tuple(sorted([x['Residue1'], x['Residue2']])), axis=1)

    # Group by 'Domain Association' and perform aggregations
    summary_df = df.groupby('Domain Association').agg({
        'Canonical Pair': 'nunique',
        'CA Distance': [
            list, 
            lambda x: (x <= 20).mean() * 100, 
            lambda x: (x <= 50).mean() * 100
        ],
        'Spectral Count': [
            list, 
            lambda x: (x >= 10).mean() * 100, 
            lambda x: (x >= 20).mean() * 100
        ],
        'Sequence Distance': list
    })

    # Flatten the MultiIndex in columns created by aggregations and rename columns to be more descriptive
    summary_df.columns = [
        'Unique Crosslinks', 
        'All CA Distances', 
        'Percentage with CA Distance <= 20', 
        'Percentage with CA Distance <= 50', 
        'All Spectral Counts', 
        'Percentage with Spectral Count >= 10', 
        'Percentage with Spectral Count >= 20', 
        'All Sequence Distances'
    ]

    # Reset index to make 'Domain Association' a column
    summary_df.reset_index(inplace=True)

    return summary_df

#######################################################################################
def save_dataframe_to_csv(df, output_filename):
    """
    Saves a DataFrame to a CSV file.

    Args:
        df (pd.DataFrame): The DataFrame to be saved.
        output_filename (str): The name of the output CSV file including the path.

    Returns:
        None
    """
    try:
        df.to_csv(output_filename, index=False)  # Save DataFrame without the index
        print(f"DataFrame successfully saved to {output_filename}")
    except Exception as e:
        print(f"An error occurred while saving the DataFrame: {e}")
      
#######################################################################################
def filter_by_domain(df, domain_string=None, exclude_intra_domain=True, exclude=False, domains_to_exclude=None):
    """
    Returns a DataFrame containing rows where "Domain Association" column matches specified conditions
    for inter-domain associations or any domain association, optionally excluding specified domains.

    Args:
        df (pd.DataFrame): The DataFrame to be filtered.
        domain_string (str, optional): The domain name or domain pair to match in "Domain Association" column.
                                      If None, all inter-domain associations are considered.
        exclude_intra_domain (bool): If True, excludes intra-domain associations.
        exclude (bool): If True, excludes the domains listed in domains_to_exclude.
        domains_to_exclude (list, optional): List of domain names to be excluded.

    Returns:
        pd.DataFrame: A new DataFrame with rows fitting the filtering criteria.
    """
    column_name = "Domain Association"

    if column_name not in df.columns:
        raise ValueError(f"The DataFrame must include a '{column_name}' column.")

    # Function to filter inter-domain associations
    def filter_inter_domain(x):
        domains = x.split(" to ")
        # Check for intra-domain links if required
        if exclude_intra_domain and len(set(domains)) == 1:
            return False
        # Exclude specified domains if required
        if exclude and domains_to_exclude and any(domain in domains_to_exclude for domain in domains):
            return False
        # Handle both specific domain pairs and single domain filters
        if domain_string:
            # Split the domain_string to handle directional domain pairs
            domain_parts = domain_string.split(" to ")
            if len(domain_parts) == 2:  # Specific domain pair (bi-directional)
                inverse_string = " to ".join(domain_parts[::-1])
                if x != domain_string and x != inverse_string:
                    return False
            elif not any(domain_string == domain.strip() for domain in domains):  # Single domain in any association
                return False
        return True

    # Apply filtering based on the defined criteria
    filtered_df = df[df[column_name].apply(filter_inter_domain)]

    return filtered_df

#######################################################################################        
def format_chimera_dist_selection(df, mod_number):
    """
    Formats each pair of residues into a string format suitable for generating distance commands in ChimeraX.

    Args:
        df (pd.DataFrame): DataFrame containing pairs of residues.
                           Expected columns are 'Residue1' and 'Residue2'.

    Returns:
        str: A single string with each command separated by a newline.
    """
    # Check for required columns
    if not {'Residue1', 'Residue2'}.issubset(df.columns):
        raise ValueError("DataFrame must include 'Residue1' and 'Residue2' columns.")

    # Generate formatted strings for each row
    formatted_strings = [
        f"distance #{mod_number}:{int(row['Residue1'])}@CA #{mod_number}:{int(row['Residue2'])}@CA" for index, row in df.iterrows()
    ]

    return '\n'.join(formatted_strings)

#######################################################################################   
def filter_by_residue(df, residue):
    """
    Filters a DataFrame to include only unique crosslinks that involve a specified residue, returning all original columns.

    Args:
        df (pd.DataFrame): The DataFrame to be filtered, expected to have 'Residue1' and 'Residue2' columns.
        residue (int): The specific residue number to filter for crosslinks.

    Returns:
        pd.DataFrame: A new filtered DataFrame with unique crosslinks involving the specified residue, including all original columns.
    """
    if 'Residue1' not in df.columns or 'Residue2' not in df.columns:
        raise ValueError("DataFrame must include 'Residue1' and 'Residue2' columns.")

    # Use .copy() to ensure the operations are performed on a copy of the data, avoiding SettingWithCopyWarning
    df = df.copy()

    # Create a canonical form of crosslink pairs to identify unique links
    df['Canonical Pair'] = df.apply(lambda x: tuple(sorted([x['Residue1'], x['Residue2']])), axis=1)

    # Remove duplicates based on the canonical pair
    df_deduplicated = df.drop_duplicates(subset='Canonical Pair')

    # Filter for pairs that include the specified residue
    filtered_df = df_deduplicated[
        (df_deduplicated['Residue1'] == residue) | (df_deduplicated['Residue2'] == residue)
    ]

    # Drop the 'Canonical Pair' column to clean up the DataFrame before returning
    filtered_df = filtered_df.drop(columns=['Canonical Pair'])

    return filtered_df

#######################################################################################
def filter_by_ca_distance(df, ca_distance_cutoff, comparison_type):
    """
    Filters a DataFrame to include only rows where the C-alpha distance meets a specified comparison (less than or greater than) to a cutoff.

    Args:
        df (pd.DataFrame): The DataFrame to be filtered, expected to have a 'CA Distance' column.
        ca_distance_cutoff (float): The cutoff value for C-alpha distances.
        comparison_type (str): Type of comparison, 'less' for <= and 'greater' for >=.

    Returns:
        pd.DataFrame: A new filtered DataFrame with rows where 'CA Distance' meets the comparison criteria.
    """
    if 'CA Distance' not in df.columns:
        raise ValueError("DataFrame must include a 'CA Distance' column.")
    
    # Determine the type of filter based on the comparison_type
    if comparison_type == 'less':
        filtered_df = df[df['CA Distance'] <= ca_distance_cutoff]
    elif comparison_type == 'greater':
        filtered_df = df[df['CA Distance'] >= ca_distance_cutoff]
    else:
        raise ValueError("comparison_type must be 'less' or 'greater'.")

    return filtered_df

#######################################################################################      
def calculate_average_and_stdev_of_column(df, column_name):
    """
    Calculates the average and standard deviation of all numeric entries in a specified column of the provided DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame from which to calculate the average and standard deviation.
        column_name (str): The name of the column to analyze.

    Returns:
        tuple: A tuple containing the average and standard deviation of the specified column, or None if the column does not exist or contains non-numeric data.
    """
    if column_name not in df.columns:
        print(f"Error: The column '{column_name}' does not exist in the DataFrame.")
        return None

    # Ensure the column is numeric and calculate average and standard deviation
    try:
        average_value = df[column_name].mean()
        stdev_value = df[column_name].std()
        return average_value, stdev_value
    except TypeError:
        print(f"Error: The column '{column_name}' contains non-numeric data.")
        return None