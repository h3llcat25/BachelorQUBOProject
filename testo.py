import random
import randomGraphGenerator

# Define a list of values with associated probabilities
values = [1, 2, 3, 4]
probabilities = [0.2, 0.3, 0.4, 0.1]

# Generate a random value based on the probabilities
random_value = random.choices(values, probabilities)[0]
random_vale = random.choice(values)


def remove_numbers_from_list(input_list, numbers_to_remove):
    # Use list comprehension to filter out numbers not in the set
    result_list = [x for x in input_list if x not in numbers_to_remove]
    return result_list


# Example usage:
my_list = [1, 2, 3, 4, 5]

fyl = my_list.remove(3)

myl = my_list[:]

numbers_set = {2, 4}
# filtered_list = remove_numbers_from_list(my_list, numbers_set)

myl.append(6)

all_reacts_list = [randomGraphGenerator.React(i) for i in range(10)]

dicto = {1: "stringo", "extra": 2.2, "seven": "heyjaw"}

print(dicto)
