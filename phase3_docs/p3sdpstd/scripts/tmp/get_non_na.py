import pandas as pd
import sys

def get_non_na_values(csv_file, column_name):
    """
    Reads the given CSV file and returns a list of values from the specified column
    that are not NA.
    """
    # Read CSV file into a DataFrame.
    # df = pd.read_csv(csv_file, delimiter='|')
    df = pd.read_fwf(csv_file, delimiter='|')
    print(df.columns)

    # Check if the column exists.
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in {csv_file}.")
    
    # Drop NA values and return the non-NA values as a list.
    return df[column_name].dropna().tolist()

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <column_name>")
        sys.exit(1)
        
    column_name = sys.argv[1]
    
    primary_file = "table_fixedwidth_primary.csv"
    extension_file = "table_fixedwidth_extention.csv"
    
    print(f"Non-NA values from '{primary_file}' for column '{column_name}':")
    try:
        primary_values = get_non_na_values(primary_file, column_name)
        print(primary_values)
    except Exception as e:
        print(f"Error reading {primary_file}: {e}")
    
    print("\n" + "-"*50 + "\n")
    
    print(f"Non-NA values from '{extension_file}' for column '{column_name}':")
    try:
        extension_values = get_non_na_values(extension_file, column_name)
        print(extension_values)
    except Exception as e:
        print(f"Error reading {extension_file}: {e}")

if __name__ == "__main__":
    main()