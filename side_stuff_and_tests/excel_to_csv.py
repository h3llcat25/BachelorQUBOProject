import pandas as pd

def excel_to_csv(excel_file_path, csv_file_path):
    # Read the Excel file
    df = pd.read_excel(excel_file_path)

    # Save the DataFrame to a CSV file
    df.to_csv(csv_file_path, index=False)  # Set index=False to avoid writing row indices


# Usage
excel_file_path = 'C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\large_dots_stats.xlsx'  # Update this to the path of your Excel file
csv_file_path = 'C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\large_dots_stats.csw'  # Specify the desired path for the output CSV file

# Convert the Excel file to a CSV file
excel_to_csv(excel_file_path, csv_file_path)

print(f"Excel file {excel_file_path} has been converted to CSV file {csv_file_path}.")

