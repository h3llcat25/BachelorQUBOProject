def collect_input_to_set():
    input_set = set()  # Initialize an empty set
    print("Enter your inputs (type 'done' when finished):")

    while True:  # Create an infinite loop to continuously accept input
        user_input = input()  # Take input from the terminal
        if user_input == 'done':  # Check if the input is 'done'
            break  # Exit the loop if the user enters 'done'
        input_set.add(user_input)  # Add the input to the set

    disease_comp_list = ["C00036", "C01165", "C00183", "C00148", "C00074", "C00082", "C05382", "C00169", "C00049",
                         "C00236", "C03287", "C00327", "C00199", "C01005", "C00141", "C00407", "C00233", "C00119"]


    return len(input_set) == len(disease_comp_list)


def extract_substrings(input_string):
    substrings = set()  # Initialize an empty set to store the substrings

    # Iterate through the string to check each substring of length 6
    for i in range(len(input_string) - 5):
        # Extract the current 6-character long substring
        substring = input_string[i:i+6]

        # Check if the substring meets the criteria (6 characters long and starts with "C0")
        if substring.startswith("C0"):
            substrings.add(substring)  # Add the substring to the set

    disease_comp_list = ["C00036", "C01165", "C00183", "C00148", "C00074", "C00082", "C05382", "C00169", "C00049",
                         "C00236", "C03287", "C00327", "C00199", "C01005", "C00141", "C00407", "C00233", "C00119"]

    print(f"length of should be is {len(disease_comp_list)} and is {len(substrings)}")

    return substrings

# Example usage
input_string = ("C01165, C00183 2 urinary C00036, C00148 3 bleeding C00074, C00148, C00082"
"4 digestive C05382, C00148, C00082 5 liver C00036, C00169, C03287 6 vision C00169, C00049, C00236"
"7 hematopoietic C00036, C00169, C05382, C00082 8 immune C00082, C00183, C00327, C01005"
"9 infection C00036, C00148, C03287, C00199 10 muscle C00036, C00169, C00082, C01005"
"11 cardiovascular C00036, C00148, C00074, C00169, C03287, C00199, C00141, C00082, C01005"
"12 inflammation C00036, C00148, C05382, C00407, C00169, C03287, C00233, C00082, C01005"
"13 neural C00036, C00148, C05382, C00074, C00169, C03287, C00199, C00119, C00082, C01005")
resulting_substrings = extract_substrings(input_string)
print("Extracted substrings:", resulting_substrings)

# Function usage example
# user_input_set = collect_input_to_set()
# print("Collected inputs:", user_input_set)